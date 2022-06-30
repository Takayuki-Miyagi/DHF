#include <iomanip>
#include <set>
#include <unordered_set>
#include <iostream>
#include <armadillo>
#include "Orbits.hh"
#include "TwoBodySpace.hh"
#include "ModelSpace.hh"


OneBodySpace::~OneBodySpace()
{}

OneBodySpace::OneBodySpace()
{}

OneBodySpace::OneBodySpace(Orbits& orbs)
  : orbits(&orbs)
{
  std::set<int> tmp;
  for (Orbit& o : orbits->orbits) tmp.insert(o.kappa);
  for (auto & it : tmp) kappas.push_back(it);
  std::sort(kappas.begin(), kappas.end(), std::greater<int>());
  for (int channel_idx=0; channel_idx<kappas.size(); channel_idx++){
    std::vector <int> idxs;
    int n_large = 0;
    int n_small = 0;
    for (auto & o : orbits->orbits){
      if (o.kappa != kappas[channel_idx]) continue;
      if (o.ls == 1) n_large += 1;
      if (o.ls ==-1) n_small += 1;
      idxs.push_back(orbits->GetOrbitIndex(o));
      orbit_index_to_channel_index[orbits->GetOrbitIndex(o)] = channel_idx;
    }
    channels.push_back(idxs);
    n_large_channel.push_back(n_large);
    n_small_channel.push_back(n_small);
  }
  number_channels = channels.size();
};

void OneBodySpace::PrintSpace()
{
  for (int idx=0; idx<GetNumberChannels(); idx++)
  {
    std::cout << std::endl;
    std::vector<int> Indices = channels[idx];
    for (int i : Indices){
      Orbit& o = orbits->GetOrbit(i);
      o.Print();
    }
  }
}

ModelSpace::~ModelSpace()
{}

ModelSpace::ModelSpace(int Ne, double Z, double zeta, Orbits orbs)
  : Ne(Ne), Z(Z), zeta(zeta), orbits(orbs)
{
  hole_occ = AssignHoles(Ne);
  one = OneBodySpace(orbits);
  two = TwoBodySpace(orbits);
  //PrintHoleOrbits();
}

ModelSpace::ModelSpace(std::string atom, double zeta, Orbits orbs, std::string valence_space)
  : zeta(zeta), orbits(orbs)
{
  std::vector<std::string> periodic_table = {
    "NA",
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
    "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og" };

  int i = 0;
  while (not isdigit(atom[i]) and atom[i] != '+' and atom[i] != '-') i++;
  std::string element = atom.substr(0,i);
  auto it_elem = find(periodic_table.begin(),periodic_table.end(),element);
  Z = it_elem - periodic_table.begin();
  if(element==atom){
    Ne = Z;
  }
  else if (atom[i]=='+') {
    int n;
    std::stringstream( atom.substr(i+1,atom.size()-i-1)) >> n;
    Ne = Z - n;
  }
  else if (atom[i]=='-') {
    int n;
    std::stringstream( atom.substr(i+1,atom.size()-i-1)) >> n;
    Ne = Z + n;
  }
  else {
    std::cout << "Warning: Unknown format of atom" << std::endl;
  }
  hole_occ = AssignHoles(Ne);
  UpdateOccupation(hole_occ);

  small_components.clear();
  large_components.clear();
  all_orbits.clear();
  holes.clear();
  particles.clear();
  core.clear();
  valence.clear();
  qspace.clear();

  for (int i=0; i<GetNumberOrbits(); i++){
    Orbit & o = GetOrbit(i);
    if(o.ls ==-1) small_components.insert(i);
    if(o.ls == 1) large_components.insert(i);
    all_orbits.insert(i);
    auto it = hole_occ.find(i);
    if(it != hole_occ.end()) holes.insert(i);
    else particles.insert(i);
  }

  if(valence_space == ""){
    for (int i = 0; i<GetNumberOrbits(); i++){
      auto it = hole_occ.find(i);
      if(it == hole_occ.end()) qspace.insert(i);
      else core.insert(i);
    }
  }
  else{
    std::istringstream ss(valence_space);
    std::string orbit_str, core_str;
    getline(ss, core_str, ',');
    auto it = find(periodic_table.begin(),periodic_table.end(),core_str);
    int Ncore = it - periodic_table.begin();
    std::map<int,double> holes = AssignHoles(Ncore);
    for (auto & it : holes) core.insert(it.first);
    while(getline(ss, orbit_str, ',')) valence.insert(GetOrbitIndex(orbit_str));
    for (int i=0; i<GetNumberOrbits(); i++){
      auto it_core = core.find(i);
      auto it_valence = valence.find(i);
      if(it_core == core.end() and it_valence == valence.end()) qspace.insert(i);
    }
  }
  one = OneBodySpace(orbits);
  two = TwoBodySpace(orbits);
  SetKetIndices();
  //PrintHoleOrbits();
  PrintOrbitals();
}

std::map<int,double> ModelSpace::AssignHoles(int N_ele)
{
  //
  // K:  0s
  // L:  1s                                     1p 1p 1p
  // M:  2s                                     2p 2p 2p
  // N:  3s                      2d 2d 2d 2d 2d 3p 3p 3p
  // O:  4s                      3d 3d 3d 3d 3d 4p 4p 4p
  // P:  5s 3f 3f 3f 3f 3f 3f 3f 4d 4d 4d 4d 4d 5p 5p 5p
  // Q:  6s 4f 4f 4f 4f 4f 4f 4f 5d 5d 5d 5d 5d 6p 6p 6p
  //
  std::map<int,double> tmp_holes;
  int N = 0;
  for (int e=0; e<=orbits.GetNmax(); ++e) {
    int d = std::min(N_ele-N, 2); // n=e, l=0
    tmp_holes[orbits.GetOrbitIndex(e,0,1,1)] = d * 0.5;
    N += d;
    if(N==N_ele) return tmp_holes;
    for (int l=std::min(e,orbits.GetLmax()); l>=1; --l){
      int n = e - l - std::max(0,l-1);
      if(n > orbits.GetNmax() or n < 0) continue;
      for (int j2=std::abs(2*l-1); j2<=2*l+1; j2+=2)
      {
        d = std::min(N_ele-N, j2+1);
        tmp_holes[orbits.GetOrbitIndex(n,l,j2,1)] = d / (j2+1.0);
        N += d;
        if(N==N_ele) return tmp_holes;
      }
    }
  }
  std::cout << "Something seems wrong in AssignHoles, N=" << N << std::endl;
  std::cout << "Try again with a larger model space" << std::endl;
  exit(0);
  return tmp_holes;
}


void ModelSpace::PrintModelSpace(bool print_obs, bool print_tbs)
{
  std::cout << " Number of electrons: " << GetElectronNumber() << std::endl;
  std::cout << " Number of protons: " << GetProtonNumber() << std::endl;
  std::cout << " Basis parameter zeta: " << GetZeta() << std::endl;
  if(print_obs) one.PrintSpace();
  if(print_tbs) two.PrintSpace();
}

void ModelSpace::PrintHoleOrbits()
{
  std::cout << std::endl;
  std::cout << " Number of electrons: " << GetElectronNumber() << std::endl;
  std::cout << " List of hole orbits: " << std::endl;
  for (auto & it: hole_occ)
  {
    Orbit& o = orbits.GetOrbit(it.first);
    int width = 4;
    std::cout << " idx =" << std::setw(width) << o.idx <<
      ", n =" << std::setw(width) << o.n << ", l =" << std::setw(width) << o.l << ", j2 =" << std::setw(width) << o.j2 <<
      ", prob = " << std::setw(width) << it.second << ", occ = " << std::setw(width) << it.second*(o.j2+1) << std::endl;
  }
  std::cout << std::endl;
}

void ModelSpace::PrintOrbitals()
{
  std::cout << std::endl;
  for (int i = 0; i<GetNumberOrbits(); i++){
    Orbit & o = orbits.GetOrbit(i);
    std::string ph = "", cvq = "";
    auto it_hole = holes.find(i);
    auto it_particle = particles.find(i);
    auto it_core = core.find(i);
    auto it_valence = valence.find(i);
    auto it_qspace = qspace.find(i);
    if(it_hole != holes.end()) ph = "hole";
    if(it_particle != particles.end()) ph = "particle";
    if(it_core != core.end()) cvq = "core";
    if(it_valence != valence.end()) cvq = "valence";
    if(it_qspace != qspace.end()) cvq = "qspace";
    int w_int = 4, w_str = 10;
    std::cout << " index =" << std::setw(w_int) << o.idx <<
      ", n =" << std::setw(w_int) << o.n << ", l =" << std::setw(w_int) << o.l << ", j2 =" << std::setw(w_int) << o.j2 <<
      ", L/S =" << std::setw(w_int) << o.ls <<
      ", prob = " << std::setw(w_int+2) << o.occ << ", occ = " << std::setw(w_int+2) << o.occ*(o.j2+1) <<
      ", " << std::setw(w_str) << ph;
      if(cvq != "") std::cout << ", " << std::setw(w_str) << cvq;
      std::cout << std::endl;
  }
  std::cout << std::endl;
}

void ModelSpace::UpdateOccupation(std::map<int,double> tmp_holes)
{
  for (auto & o : orbits.orbits) o.occ = 0;
  for (auto & it : tmp_holes){
    Orbit & o = GetOrbit(it.first);
    o.occ = it.second;
  }
}

void ModelSpace::SetKetIndices()
{
  std::vector<int> cvq;
  for (int i = 0; i<GetNumberOrbits(); i++){
    auto it_core = core.find(i);
    auto it_valence = valence.find(i);
    auto it_qspace = qspace.find(i);
    if(it_core != core.end()) cvq.push_back(0);
    if(it_valence != valence.end()) cvq.push_back(1);
    if(it_qspace != qspace.end()) cvq.push_back(2);
  }

  for (TwoBodyChannel & tbc : two.channels){
    tbc.KetIndex_cc.clear();
    tbc.KetIndex_vc.clear();
    tbc.KetIndex_qc.clear();
    tbc.KetIndex_vv.clear();
    tbc.KetIndex_qv.clear();
    tbc.KetIndex_qq.clear();

    for (int i=0; i<tbc.GetNumberStates(); i++){
      int ip = tbc.GetOrbitIndex1(i);
      int iq = tbc.GetOrbitIndex2(i);
      int cvq_p = cvq[ip];
      int cvq_q = cvq[iq];
      if (cvq_p+cvq_q==0)      tbc.KetIndex_cc.push_back(i); // 00
      if (cvq_p+cvq_q==1)      tbc.KetIndex_vc.push_back(i); // 01
      if (std::abs(cvq_p-cvq_q)==2) tbc.KetIndex_qc.push_back(i); // 02
      if (cvq_p*cvq_q==1)      tbc.KetIndex_vv.push_back(i); // 11
      if (cvq_p+cvq_q==3)      tbc.KetIndex_qv.push_back(i); // 12
      if (cvq_p+cvq_q==4)      tbc.KetIndex_qq.push_back(i); // 22
    }
  }
}

//void ModelSpace::UpdateOrbitals(arma::vec SPEs)
//{
//  std::map<int,double> tmp_holes;
//  OneBodySpace& obs = GetOneBodySpace();
//
//  small_components.clear();
//  large_components.clear();
//  all_orbits.clear();
//  holes.clear();
//  particles.clear();
//  for (int i=0; i<GetNumberOrbits(); i++) all_orbits.insert(i);
//  for (int ich=0; ich< obs.GetNumberChannels(); ich++) {
//    std::vector<double> tmp_holes_ch;
//    for (auto & it : hole_occ){
//      Orbit& o_h = GetOrbit(it.first);
//      if (o_h.kappa != obs.kappas[ich]) continue;
//      tmp_holes_ch.push_back(it.second);
//    }
//
//    std::vector<unsigned long long> tmp(obs.channels[ich].begin(), obs.channels[ich].end());
//    arma::uvec sub_idx(std::vector<unsigned long long>(tmp.begin(), tmp.end()));
//    arma::vec SPE_ch = SPEs(sub_idx);
//    arma::uvec sorted_idx = sub_idx(arma::sort_index(SPE_ch));
//    int cnt = 0;
//    for (int i : sorted_idx){
//      cnt += 1;
//      if(cnt <= obs.n_small_channel[ich]){
//        small_components.insert(i);
//        all_orbits.insert(i);
//      }
//      else {
//        large_components.insert(i);
//        all_orbits.insert(i);
//      }
//      if(cnt <= obs.n_small_channel[ich]) continue;
//      if(cnt-obs.n_small_channel[ich] > tmp_holes_ch.size()) break;
//      tmp_holes[i] = tmp_holes_ch[cnt-obs.n_small_channel[ich]-1];
//      holes.insert(i);
//    }
//  }
//  hole_occ = tmp_holes;
//  UpdateOccupation(tmp_holes);
//  for (int i : all_orbits){
//    auto it = holes.find(i);
//    if(it != holes.end()) continue;
//    particles.insert(i);
//  }
//  // TODO: core valence, qspace
//}
