#include "PhysicalConstants.hh"
#include "Orbits.hh"
#include "TwoBodySpace.hh"

int main(int argc, char** argv)
{
  std::cout << " TwoBodySpace test " << std::endl;
  double Z = 1;
  double zeta = 1;
  int wint = 4;
  int wdouble = 12;
  Orbits orbits = Orbits(4,0);
  orbits.PrintOrbits();

  TwoBodySpace tbs = TwoBodySpace(orbits);
  std::cout << " Channel numbers: " << tbs.GetNumberChannels() << std::endl;
  tbs.PrintSpace();
}

