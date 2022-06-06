#include <armadillo>
#include "PhysicalConstants.hh"
#include "Orbits.hh"
#include "ModelSpace.hh"
#include "TwoBodyOperator.hh"
#include "Operator.hh"
#include "HartreeFock.hh"
int main(int argc, char** argv)
{
  std::cout << " HartreeFock test " << std::endl;
  int ElectronNumber = 2;
  double Z = 2;
  double zeta = 0.5;
  int wint = 4;
  int wdouble = 12;
  Orbits orbits = Orbits(2, 0);
  orbits.PrintOrbits();
  ModelSpace ms = ModelSpace(ElectronNumber, Z, zeta, orbits);
  Operator H = Operator(ms);

  H.SetDiracCoulombHamiltonian(true, true);
  H.OrthoNormalize();
  std::map<int, double> holes;
  holes[0] = 1.0;
  HartreeFock HF = HartreeFock(H, holes);
  HF.Solve();
}
