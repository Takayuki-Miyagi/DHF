#include <iomanip>
#include <armadillo>
#include "PhysicalConstants.hh"
#include "Orbits.hh"

arma::mat OrthoNormalize(arma::mat H, arma::mat S)
{
  arma::mat L = arma::chol(S,"lower");
  arma::mat H_orth = arma::inv(L) * H * arma::inv(L).t();
  return H_orth;
}

arma::uvec GetElectronStates(arma::vec SPEs, Orbits orbits)
{
  int num_electron_states = 0;
  for (auto o: orbits.orbits) {if(o.e2==1) num_electron_states += 1;}
  arma::uvec indx(num_electron_states);
  int norbs = orbits.GetNumberOrbits();
  arma::vec SPE_diff(norbs-1);
  for (int i = 0; i<norbs-1; i++){
    SPE_diff(i) = std::abs(SPEs(i+1) - SPEs(i));
  }
  int idx_border = SPE_diff.index_max();
  for (int i = 0; i<num_electron_states; i++){
    indx(i) = i+idx_border+1;
  }
  return indx;
}

int main(int argc, char** argv)
{
  std::cout << " Orbit test " << std::endl;
  double Z = 1;
  double zeta = 4;
  int wint = 4;
  int wdouble = 12;
  Orbits orbits = Orbits(4, 0);
  orbits.PrintOrbits();
  int norbs = orbits.GetNumberOrbits();

  arma::mat H = arma::mat(norbs, norbs, arma::fill::zeros);
  arma::mat S = arma::mat(norbs, norbs, arma::fill::zeros);

  for (int i1=0; i1<norbs; i1++){
    for (int i2=i1; i2<norbs; i2++){
      Orbit o1 = orbits.GetOrbit(i1);
      Orbit o2 = orbits.GetOrbit(i2);
      if(o1.j2 != o2.j2) continue;
      //std::cout << std::setw(wint) << i1 << " "
      //  << std::setw(wint) << i2 << " "
      //  << std::setw(wdouble) << MEOverlap(o1, o2, zeta, Z) << " "
      //  << std::setw(wdouble) << MEKinetic(o1, o2, zeta, Z) << " "
      //  << std::setw(wdouble) << MENuclPot(o1, o2, zeta, Z) << " "
      //  << std::setw(wdouble) << NormCheck(o1, o2, zeta, Z) << std::endl;
      double mass_term = 0;
      if(o1.e2 ==-1) {mass_term = -2*PhysConst::c * PhysConst::c*MEOverlap(o1, o2, zeta, Z);}
      H(i1,i2) = MEKinetic(o1, o2, zeta, Z) - MENuclPot(o1, o2, zeta, Z) + mass_term;
      S(i1,i2) = MEOverlap(o1, o2, zeta, Z);
      H(i2,i1) = H(i1,i2);
      S(i2,i1) = S(i1,i2);
    }
  }

  std::cout << S << std::endl;

  arma::mat H_orth = OrthoNormalize(H, S);
  arma::vec eig;
  arma::mat vec_orth;
  arma::eig_sym(eig, vec_orth, H_orth);
  arma::uvec electrons = GetElectronStates(eig, orbits);
  std::cout << eig(electrons(0)) << std::endl;
}

