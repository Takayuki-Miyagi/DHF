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
  std::cout << " HartreeFock test " << std::endl;
  std::cout << " Number of OpenMP threads: " <<  omp_get_max_threads() << std::endl;
  Parameters parameters(argc,argv);

  std::string orbitals = parameters.s("orbitals");
  std::string atom = parameters.s("atom");
  std::string radial_function_type = parameters.s("radial_function_type");
  int NMesh = parameters.i("NMesh");
  double zeta_inv = parameters.d("zeta_inv");
  Orbits orbits = Orbits(orbitals, radial_function_type); 
  orbits.Print();
  ModelSpace ms = ModelSpace(atom, 1/zeta_inv, orbits);
  Operator H = Operator(ms);
  H.TwoBody.SetNMesh(NMesh);

  if(orbits.relativistic) H.SetDiracCoulombHamiltonian(true, true);
  else H.SetCoulombHamiltonian(true, true);
  H.OrthoNormalize();
  //H.Print();
  HartreeFock HF = HartreeFock(H);
  HF.Solve();
  std::cout << " HF energy: " << std::setw(16) << std::setprecision(8) << HF.EHF << std::endl;
}
