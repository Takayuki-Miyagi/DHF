#include "PhysicalConstants.hh"
#include "Orbits.hh"
#include "ModelSpace.hh"
#include "TwoBodyOperator.hh"
int main(int argc, char** argv)
{
  std::cout << " TwoBodyOperator test " << std::endl;
  double Z = 1;
  double zeta = 1;
  int wint = 4;
  int wdouble = 12;
  Orbits orbits = Orbits(0, 1);
  ModelSpace ms = ModelSpace(2, 2, 2, orbits);
  TwoBodyOperator op = TwoBodyOperator(ms);

  op.Set2BME_J(0, 0, 0, 0, 0, 0, 1);
  op.Set2BME_J(0, 0, 1, 1, 1, 1, 2);
  op.Set2BME_J(0, 0, 1, 1, 0, 0, 3);
  op.Set2BME_J(0, 0, 2, 2, 1, 1, 4);
  op.Set2BME_J(0, 0, 2, 2, 2, 2, 5);
  op.Print();
}



