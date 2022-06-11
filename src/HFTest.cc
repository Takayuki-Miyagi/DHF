#include <armadillo>
#include "Parameters.hh"
#include "PhysicalConstants.hh"
#include "Orbits.hh"
#include "ModelSpace.hh"
#include "TwoBodyOperator.hh"
#include "Operator.hh"
#include "HartreeFock.hh"
int main(int argc, char** argv)
{
  Parameters parameters(argc,argv);
  std::cout << " HartreeFock test " << std::endl;

  std::string orbitals = parameters.s("orbitals");
  std::string atom = parameters.s("atom");
  std::string radial_function_type = parameters.s("radial_function_type");
  double zeta_inv = parameters.d("zeta_inv");
  Orbits orbits = Orbits(orbitals, radial_function_type); 
  orbits.Print();
  ModelSpace ms = ModelSpace(atom, 1/zeta_inv, orbits);
  Operator H = Operator(ms);

  if(orbits.relativistic) H.SetDiracCoulombHamiltonian(true, true);
  else H.SetCoulombHamiltonian(true, true);
  H.OrthoNormalize();
  //H.Print();
  HartreeFock HF = HartreeFock(H);
  HF.Solve();
  std::cout << HF.EHF << std::endl;
}
