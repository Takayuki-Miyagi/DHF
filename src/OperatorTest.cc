#include <armadillo>
#include "PhysicalConstants.hh"
#include "Orbits.hh"
#include "ModelSpace.hh"
#include "TwoBodyOperator.hh"
#include "Operator.hh"
int main(int argc, char** argv)
{
  std::cout << " Operator test " << std::endl;
  int ElectronNumber = 2;
  double Z = 2;
  double zeta = 1;
  int wint = 4;
  int wdouble = 12;
  Orbits orbits = Orbits(2, 0);
  ModelSpace ms = ModelSpace(ElectronNumber, Z, zeta, orbits);
  Operator H = Operator(ms);

  H.SetDiracCoulombHamiltonian();
  std::cout << H.OneBody << std::endl;
  H.OrthoNormalize();
  H.TwoBody.Print();
}
