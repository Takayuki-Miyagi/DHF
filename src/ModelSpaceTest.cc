#include "PhysicalConstants.hh"
#include "Orbits.hh"
#include "ModelSpace.hh"
int main(int argc, char** argv)
{
  std::cout << " ModelSpace test " << std::endl;
  double Z = 1;
  double zeta = 1;
  int wint = 4;
  int wdouble = 12;
  Orbits orbits = Orbits(1, 2);
  ModelSpace ms = ModelSpace(2, 2, 2, orbits);
  ms.PrintModelSpace(true, true);
}


