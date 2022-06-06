#include "Orbits.hh"
#include "TwoBodySpace.hh"
#include "ModelSpace.hh"

#include <set>
#include <unordered_set>
//#include <bits/stdc++.h>
#include <armadillo>

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
      o.PrintOrbit();
    }
  }
}

ModelSpace::~ModelSpace()
{}

ModelSpace::ModelSpace(int Ne, double Z, double zeta, Orbits orbs)
  : Ne(Ne), Z(Z), zeta(zeta), orbits(orbs)
{
  one = OneBodySpace(orbits);
  two = TwoBodySpace(orbits);
}

void ModelSpace::PrintModelSpace(bool print_obs, bool print_tbs)
{
  std::cout << " Number of electrons: " << GetElectronNumber() << std::endl;
  std::cout << " Number of protons: " << GetProtonNumber() << std::endl;
  std::cout << " Basis parameter zeta: " << GetZeta() << std::endl;
  if(print_obs) one.PrintSpace();
  if(print_tbs) two.PrintSpace();
}
