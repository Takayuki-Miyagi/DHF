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
  for (Orbit& o : orbits->orbits){ kappas.insert(o.k); }
  std::vector<int> kappas1(kappas.begin(), kappas.end());
  std::sort(kappas1.begin(), kappas1.end(), std::greater<int>());
  for (int channel_idx=0; channel_idx<kappas1.size(); channel_idx++){
    std::vector <int> idxs;
    for (auto o : orbits->orbits){
      if (o.k != kappas1[channel_idx]) {continue;}
      idxs.push_back(orbits->GetOrbitIndex(o));
      orbit_index_to_channel_index[orbits->GetOrbitIndex(o)] = channel_idx;
    }
    channels.push_back(idxs);
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
  PrintHoleOrbits();
}

ModelSpace::ModelSpace(std::string atom, double zeta, Orbits orbs)
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
  one = OneBodySpace(orbits);
  two = TwoBodySpace(orbits);
  PrintHoleOrbits();
  
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
  // e = n + l
  std::map<int,double> holes;
  int N = 0;
  for (int e=0; e<=200; ++e) {
    int d = std::min(N_ele-N, 2); // n=e, l=0
    holes[orbits.GetOrbitIndex(e,0,1,1)] = d * 0.5;
    N += d;
    if(N==N_ele) return holes;
    for (int l=e; l>=1; --l){
      int n = e - l - std::max(0,l-1);
      if(n < 0) continue;
      for (int j2=std::abs(2*l-1); j2<=2*l+1; j2+=2)
      {
        d = std::min(N_ele-N, j2+1);
        holes[orbits.GetOrbitIndex(n,l,j2,1)] = d / (j2+1.0);
        N += d;
        if(N==N_ele) return holes;
      }
    }
  }
  std::cout << "Something seems wrong in AssignHoles, N=" << N << std::endl;
  return holes;
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
  std::cout << " Number of electrons: " << GetElectronNumber() << std::endl;
  std::cout << " List of hole orbits: " << std::endl;
  for (int idx=0; idx<orbits.GetNumberOrbits(); idx++)
  {
    if(std::abs(hole_occ[idx]) < 1.e-4) continue;
    Orbit& o = orbits.GetOrbit(idx);
    int width = 4;
    std::cout << " idx =" << std::setw(width) << orbits.GetOrbitIndex(o.n, o.l, o.j2, o.e2) <<
      ", n =" << std::setw(width) << o.n << ", l =" << std::setw(width) << o.l << ", j2 =" << std::setw(width) << o.j2 << ", e2 =" << std::setw(width) << o.e2 <<
      ", prob = " << std::setw(width) << hole_occ[idx] << ", occ = " << std::setw(width) << hole_occ[idx]*(o.j2+1) << std::endl;
    
  }
}
