#include <armadillo>
//#include <boost/range/irange.hpp>
//#include <cmath>
//#include <wignerSymbols.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_monte.h>
//#include <gsl/gsl_monte_plain.h>
//#include <gsl/gsl_monte_vegas.h>
//#include <gsl/gsl_monte_miser.h>
//#include <gsl/gsl_sf_coupling.h>
// #include <wignerSymbols>
// #include "wignerSymbols/wignerSymbols-cpp.h"
// #include<gsl>
#include "PhysicalConstants.hh"
#include "TwoBodyOperator.hh"
#include "Operator.hh"

// Wrapper to convert lambda with capture to gsl function for integration
//template< typename F >  class gsl_function_pp : public gsl_function {
//  public:
//    gsl_function_pp(const F& func) : _func(func) {
//      function = &gsl_function_pp::invoke;
//      params = this;
//    }
//  private:
//    const F& _func;
//    static double invoke(double x, void* params) {
//      return static_cast<gsl_function_pp*>(params)->_func(x);
//    }
//};

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Operator::~Operator()
{}

Operator::Operator()
{}

Operator::Operator(ModelSpace& ms, int rankJ, int rankP, int rankTz)
  : modelspace(&ms), rankJ(rankJ), rankP(rankP), rankTz(rankTz)
{
  Orbits& orbits = modelspace->orbits;
  int norbs = orbits.GetNumberOrbits();
  OneBody = arma::mat(norbs, norbs, arma::fill::zeros);
  S = arma::mat(norbs, norbs, arma::fill::eye);
  TwoBody = TwoBodyOperator(*modelspace, rankJ, rankP, rankTz);
  for (int i=0; i<norbs; i++){
    for (int j=i; j<norbs; j++){
      Orbit& oi = orbits.GetOrbit(i);
      Orbit& oj = orbits.GetOrbit(j);
      S(i,j) = MEOverlap(oi, oj, modelspace->GetZeta(), modelspace->GetProtonNumber());
      S(j,i) = S(i,j);
    }
  }
}

void Operator::SetDiracCoulombHamiltonian(bool OneBodyTerm, bool TwoBodyTerm)
{
  if(OneBodyTerm){ SetOneBodyDiracHamiltonian();}
  if(TwoBodyTerm){ SetTwoBodyCoulombInteraction();}
}

void Operator::SetCoulombHamiltonian(bool OneBodyTerm, bool TwoBodyTerm)
{
  if(OneBodyTerm){ SetOneBodyHamiltonian();}
  if(TwoBodyTerm){ SetTwoBodyCoulombInteraction();}
}

void Operator::SetOneBodyHamiltonian()
{
  Orbits orbits = modelspace->orbits;
  int norbs = orbits.GetNumberOrbits();
  double zeta = modelspace->GetZeta();
  double Z = modelspace->GetProtonNumber();
  for (int i1=0; i1<norbs; i1++){
    for (int i2=i1; i2<norbs; i2++){
      Orbit& o1 = orbits.GetOrbit(i1);
      Orbit& o2 = orbits.GetOrbit(i2);
      if(o1.kappa != o2.kappa) continue;
      OneBody(i1,i2) = MEKinetic(o1, o2, zeta, Z) - Z * MENuclPot(o1, o2, zeta, Z);
      OneBody(i2,i1) = OneBody(i1,i2);
    }
  }
}

void Operator::SetOneBodyDiracHamiltonian()
{
  Orbits orbits = modelspace->orbits;
  int norbs = orbits.GetNumberOrbits();
  double zeta = modelspace->GetZeta();
  double Z = modelspace->GetProtonNumber();
  for (int i1=0; i1<norbs; i1++){
    for (int i2=i1; i2<norbs; i2++){
      Orbit& o1 = orbits.GetOrbit(i1);
      Orbit& o2 = orbits.GetOrbit(i2);
      if(o1.kappa != o2.kappa) continue;
      double mass_term = 0;
      if(o1.ls ==-1) mass_term = -2*PhysConst::c * PhysConst::c*MEOverlap(o1, o2, zeta, Z);
      OneBody(i1,i2) = MEKinetic(o1, o2, zeta, Z) - Z * MENuclPot(o1, o2, zeta, Z) + mass_term;
      OneBody(i2,i1) = OneBody(i1,i2);
    }
  }
}

void Operator::SetTwoBodyCoulombInteraction()
{
  TwoBody.SetTwoBodyCoulombTerm();
}

/*
   (ab:J|U1 x U2|cd:J) = [ (a|U|c) (b|U|d)
   - (phase) (a|U|d) (b|U|c)
   - (phase) (a|U|d) (b|U|c)
   + (phase) (a|U|c) (b|U|d) ] / 2 sqrt((1+del_ab) (1+del_cd))
   */
arma::mat Operator::EmbedBasisTrans2(arma::mat T, TwoBodyChannel& tbc)
{
  Orbits& orbits = modelspace->GetOrbits();
  arma::mat U = arma::mat(tbc.GetNumberStates(), tbc.GetNumberStates(), arma::fill::zeros);
  for (int ileft=0; ileft < tbc.GetNumberStates(); ileft++){
    for (int iright=0; iright < tbc.GetNumberStates(); iright++){
      int i = tbc.GetOrbitIndex1(ileft);
      int j = tbc.GetOrbitIndex2(ileft);
      int k = tbc.GetOrbitIndex1(iright);
      int l = tbc.GetOrbitIndex2(iright);

      Orbit& oi = orbits.GetOrbit(i);
      Orbit& oj = orbits.GetOrbit(j);
      Orbit& ok = orbits.GetOrbit(k);
      Orbit& ol = orbits.GetOrbit(l);
      int phase_ij = pow( (-1), (oi.j2+oj.j2)/2-tbc.J );
      int phase_kl = pow( (-1), (ok.j2+ol.j2)/2-tbc.J );
      U(ileft, iright) = ((1+phase_ij*phase_kl) * T(i,k) * T(j,l) - (phase_ij+phase_kl) * T(i,l) * T(j,k)) * 0.5;
      if (i==j) {U(ileft, iright) /= sqrt(2);}
      if (k==l) {U(ileft, iright) /= sqrt(2);}
    }
  }
  return U;
}

void Operator::OrthoNormalize()
{
  if(orthonormalized){
    std::cout << "The operator is already orthonormalized; exiting..." << std::endl;
    return;
  }
  Orbits& orbits = modelspace->GetOrbits();
  TwoBodySpace& tbs = modelspace->GetTwoBodySpace();
  arma::mat L = arma::chol(S, "lower");
  arma::mat T = arma::inv(L);
  OneBody = T * OneBody * T.t();
  for (auto it : TwoBody.Channels){
    int ichbra = it.first[0];
    int ichket = it.first[1];
    TwoBodyOperatorChannel& opch = it.second;
    TwoBodyChannel * chbra = opch.chbra;
    TwoBodyChannel * chket = opch.chket;
    arma::mat Ubra = EmbedBasisTrans2(T, *chbra);
    arma::mat Uket;
    if(ichbra==ichket) {Uket = Ubra;}
    else {Uket = EmbedBasisTrans2(T, *chket);}
    auto& out = TwoBody.Channels[{ichbra,ichket}].MEs;
    out = Ubra * opch.MEs * Uket.t();
  }
  orthonormalized = true;
}

void Operator::Print()
{
  std::cout << "One-body part" << std::endl;
  std::cout << OneBody << std::endl;
  std::cout << "Two-body part" << std::endl;
  TwoBody.Print();
}
