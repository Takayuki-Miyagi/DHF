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
      if(o1.k != o2.k) continue;
      double mass_term = 0;
      if(o1.e2 ==-1) {mass_term = -2*PhysConst::c * PhysConst::c*MEOverlap(o1, o2, zeta, Z);}
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

//
//    void Operator::set_two_body_ham() {
//        ModelSpace ms = modelspace;
//        two_body_space = ms.two;
//        Orbits orbs = ms.orbits;
//
//        for (auto channels : two.Channels) {
//            chbra = two_body_space.get_channel(channels.first[0]);
//            chket = two_body_space.get_channel(channels.first[1]);
//            for (int idxbra=0; idxbra<chbra.get_number_states(); idxbra++) {
//                for (int idxket=0; idxket<chket.get_number_states(); idxket++) {
//                    if (idxbra > idxket) {continue;}
//                    int a = chbra.get_indices(idxbra)[0];
//                    int b = chbra.get_indices(idxbra)[1];
//                    int c = chket.get_indices(idxket)[0];
//                    int d = chket.get_indices(idxket)[1];
//                    Orbit oa = orbs.get_orbit(a);
//                    Orbit ob = orbs.get_orbit(b);
//                    Orbit oc = orbs.get_orbit(c);
//                    Orbit od = orbs.get_orbit(d);
//                    double norm = 1.0;
//                    if (a==b) {norm /= sqrt(2);}
//                    if (c==d) {norm /= sqrt(2);}
//                    double v;
//                    v = _coulomb(oa, ob, oc, od, chket.J);
//                    v += _coulomb(oa, ob, od, oc, chket.J) * pow( (-1), ( floor( (oc.j+od.j)/2) - chket.J + 1) );
//                    v *= norm;
//                    two.set_2bme(channels.first[0], channels.first[1], idxbra, idxket, v);
//                    two.set_2bme(channels.first[0], channels.first[1], idxket, idxbra, v);
//                }
//            }
//        }
//    }
//
//    double Operator::_coulomb(Orbit& oa, Orbit& ob, Orbit& oc, Orbit& od, int J) {
//        if (oa.e != oc.e) {return 0.0;}
//        if (ob.e != od.e) {return 0.0;}
//        double zeta = modelspace.zeta;
//        double Z = modelspace.Z;
//        int Lmin = floor( max( abs(oa.j-oc.j), abs(ob.j-od.j) )/2 );
//        int Lmax = floor( min(    (oa.j+oc.j),    (ob.j+od.j) )/2 );
//        double rmax = 20; // 20
//        int NMesh = 100;  // 100
//        gsl_integration_fixed_workspace *workspace;
//        const gsl_integration_fixed_type *T = gsl_integration_fixed_legendre;
//        workspace = gsl_integration_fixed_alloc(T, NMesh, 0.0, rmax, 0.0, 0.0);
//        double r = 0.0;
//        for (auto L : irange(Lmin, Lmax+1)) {
//            if ( (oa.l+oc.l+L)%2 == 1 ) {continue;}
//            if ( (ob.l+od.l+L)%2 == 1 ) {continue;}
//            if (abs(oa.l-oc.l) > L or oa.l+oc.l < L) continue;
//            if (abs(ob.l-od.l) > L or ob.l+od.l < L) continue;
//            double angular = gsl_sf_coupling_6j(oa.j, ob.j, 2*J, od.j, oc.j, 2*L) * \
//                             gsl_sf_coupling_3j(oa.j, 2*L, oc.j, -1, 0, 1) * \
//                             gsl_sf_coupling_3j(ob.j, 2*L, od.j, -1, 0, 1);
//            double Integral = 0.0;
//            for (int i=0; i<NMesh; i++){
//              for (int j=0; j<NMesh; j++){
//                double x = workspace->x[i];
//                double wx = workspace->weights[i];
//                double y = workspace->x[j];
//                double wy = workspace->weights[j];
//                Integral += wx * wy * oa.eval_radial_function_rspace(x, zeta, Z, oa.e) * \
//                            ob.eval_radial_function_rspace(y, zeta, Z, ob.e) * \
//                            oc.eval_radial_function_rspace(x, zeta, Z, oc.e) * \
//                            od.eval_radial_function_rspace(y, zeta, Z, od.e) * \
//                            pow( min(x,y), L) / pow( max(x,y), (L+1) ) * oa.e * ob.e;
//              }
//            }
//            r += angular * Integral;
//        }
//        r *= sqrt( (oa.j+1) * (ob.j+1) * (oc.j+1) * (od.j+1)) * pow( (-1), floor((oa.j+oc.j)/2)+J );
//        return r;
//    }
//
//    void Operator::orthogonalize(bool scalar) {
//        ModelSpace ms = modelspace;
//        TwoBodySpace two_body_space = ms.two;
//        Orbits orbs = ms.orbits;
//        mat L = arma::chol(S, "lower");
//        Mat<double> T = arma::inv(L);
//        one = T * one * T.t();
//        if (scalar) {
//            /*
//            (ab:J|U1 x U2|cd:J) = [ (a|U|c) (b|U|d)
//            - (phase) (a|U|d) (b|U|c)
//            - (phase) (a|U|d) (b|U|c)
//            + (phase) (a|U|c) (b|U|d) ] / 2 sqrt((1+del_ab) (1+del_cd))
//            */
//            for (auto channels : two.Channels) {
//                TwoBodyChannel tbc = two_body_space.get_channel(channels.first[0]);
//                Mat<double> U;
//                U.zeros(tbc.get_number_states(), tbc.get_number_states());
//                for (int idxbra=0; idxbra<tbc.get_number_states(); idxbra++) {
//                    for (int idxket=0; idxket<tbc.get_number_states(); idxket++) {
//                        int a = tbc.get_indices(idxbra)[0];
//                        int b = tbc.get_indices(idxbra)[1];
//                        int c = tbc.get_indices(idxket)[0];
//                        int d = tbc.get_indices(idxket)[1];
//
//                        Orbit oa = orbs.get_orbit(a);
//                        Orbit ob = orbs.get_orbit(b);
//                        Orbit oc = orbs.get_orbit(c);
//                        Orbit od = orbs.get_orbit(d);
//
//                        int phase_ab = pow( (-1), floor(oa.j+ob.j)/2 - tbc.J );
//                        int phase_cd = pow( (-1), floor(oc.j+od.j)/2 - tbc.J );
//                        U(idxbra, idxket) = ((1+phase_ab*phase_cd) * T(a,c) * T(b,d) - (phase_ab+phase_cd) * T(a,d) * T(b,c)) * 0.5;
//                        if (a==b) {U(idxbra, idxket) /= sqrt(2);}
//                        if (c==d) {U(idxbra, idxket) /= sqrt(2);}
//                    }
//                }
//                two.Channels.at(channels.first).MEs = U * two.Channels.at(channels.first).MEs * U.t();
//           }
//
//        } else {
//            cout << "Not impletmented yet" << endl;
//            S = eye(orbs.get_num_orbits(), orbs.get_num_orbits());
//        }
//    }
//
//
//
//    void Operator::print_operator(){
//        Orbits orbs = modelspace.orbits;
//        OneBodySpace one_body_space = modelspace.one;
//        for (int ichbra=0; ichbra<one_body_space.number_channels; ichbra++) {
//            for (int ichket=0; ichket<one_body_space.number_channels; ichket++) {
//                // code runs really slow here ??
//                for (auto i : one_body_space.channels[ichbra]) {
//                    for (auto j : one_body_space.channels[ichket]) {
//                        // cout << i << " " << j << endl;
//                        Orbit oi = orbs.get_orbit(i);
//                        Orbit oj = orbs.get_orbit(j);
//                        if (!( abs(oi.j-oj.j) <= 2*rankJ && 2*rankJ <= oi.j+oj.j ) ) {continue;}
//                        if ( pow( (-1), (oi.l+oj.l) ) * rankP == -1 ) {continue;}
//                        printf( "   %-3d  %-3d    %- 3.3e      %- 3f \n", i, j, one(i,j), S(i,j) );
//                    }
//                }
//            }
//        }
//        two.print_two_body_operator();
//    }




