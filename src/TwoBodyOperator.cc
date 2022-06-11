#include <iomanip>
#include <armadillo>
#include <gsl/gsl_sf_coupling.h>
#include "Orbits.hh"
#include "TwoBodySpace.hh"
#include "TwoBodyOperator.hh"


TwoBodyOperatorChannel::~TwoBodyOperatorChannel()
{}

TwoBodyOperatorChannel::TwoBodyOperatorChannel()
{}

TwoBodyOperatorChannel::TwoBodyOperatorChannel(TwoBodyChannel& chbra_in, TwoBodyChannel& chket_in)
  : chbra(&chbra_in), chket(&chket_in)
{
  MEs = arma::mat(chbra->GetNumberStates(), chket->GetNumberStates(), arma::fill::zeros);
}

TwoBodyOperator::~TwoBodyOperator()
{}

TwoBodyOperator::TwoBodyOperator()
{}

TwoBodyOperator::TwoBodyOperator(ModelSpace& ms, int rankJ, int rankP, int rankTz)
  : modelspace(&ms), rankJ(rankJ), rankP(rankP), rankTz(rankTz)
{
  TwoBodySpace& tbs = modelspace->two;
  for (int ichbra=0; ichbra<tbs.GetNumberChannels(); ichbra++){
    for (int ichket=ichbra; ichket<tbs.GetNumberChannels(); ichket++){
      TwoBodyChannel& chbra = tbs.GetChannel(ichbra);
      TwoBodyChannel& chket = tbs.GetChannel(ichket);
      if(std::abs(chbra.J - chket.J) > rankJ or chbra.J+chket.J < rankJ) continue;
      if(chbra.Prty * chket.Prty * rankP == -1) continue;
      if(std::abs(chbra.Tz - chket.Tz) != rankTz) continue;
      Channels[{ichbra,ichket}] = TwoBodyOperatorChannel(chbra, chket);
    }
  }
}

double TwoBodyOperator::Get2BME(int ichbra, int ichket, int i, int j, int k, int l)
{
  TwoBodySpace& tbs = modelspace->two;
  TwoBodyChannel& chbra = tbs.GetChannel(ichbra);
  TwoBodyChannel& chket = tbs.GetChannel(ichket);
  int idxbra = chbra.GetIndex(i,j);
  int idxket = chket.GetIndex(k,l);
  if(i==j and chbra.J%2==1) return 0;
  if(k==l and chket.J%2==1) return 0;
  if(chbra.Prty * chket.Prty * rankP == -1) return 0;
  if(std::abs(chbra.Tz - chket.Tz) != rankTz) return 0;
  if(std::abs(chbra.J - chket.J) > rankJ or chbra.J+chket.J < rankJ) return 0;
  double phase = 1;
  phase *= chbra.GetPhaseFactor(i,j);
  phase *= chket.GetPhaseFactor(k,l);
  if(ichbra>ichket){
    return Channels[{ichket,ichbra}].MEs(idxket,idxbra) * pow((-1), chbra.J-chket.J) * phase;
  }
  return Channels[{ichbra,ichket}].MEs(idxbra,idxket) * phase;
}

double TwoBodyOperator::Get2BME(int ichbra, int ichket, Orbit& oi, Orbit& oj, Orbit& ok, Orbit& ol)
{
  Orbits& orbits = modelspace->orbits;
  int i = orbits.GetOrbitIndex(oi);
  int j = orbits.GetOrbitIndex(oj);
  int k = orbits.GetOrbitIndex(ok);
  int l = orbits.GetOrbitIndex(ol);
  return Get2BME(ichbra, ichket, i, j, k, l);
}

double TwoBodyOperator::Get2BME_J(int Jbra, int Jket, Orbit& oi, Orbit& oj, Orbit& ok, Orbit& ol)
{
  TwoBodySpace & tbs = modelspace->two;
  int Prty_bra = pow((-1), oi.l+oj.l);
  int Prty_ket = pow((-1), ok.l+ol.l);
  int Tz_bra = (oi.e2 + oj.e2)/2;
  int Tz_ket = (ok.e2 + ol.e2)/2;
  int ichbra = tbs.GetChannelIndex(Jbra, Prty_bra, Tz_bra);
  int ichket = tbs.GetChannelIndex(Jket, Prty_ket, Tz_ket);
  return Get2BME(ichbra, ichket, oi, oj, ok, ol);
}

double TwoBodyOperator::Get2BME_J(int Jbra, int Jket, int i, int j, int k, int l)
{
  Orbits& orbits = modelspace->GetOrbits();
  Orbit& oi = orbits.GetOrbit(i);
  Orbit& oj = orbits.GetOrbit(j);
  Orbit& ok = orbits.GetOrbit(k);
  Orbit& ol = orbits.GetOrbit(l);
  return Get2BME_J(Jbra, Jket, oi, oj, ok, ol);
}

void TwoBodyOperator::Set2BME(int ichbra, int ichket, int i, int j, int k, int l, double me)
{
  TwoBodySpace& tbs = modelspace->GetTwoBodySpace();
  TwoBodyChannel& chbra = tbs.GetChannel(ichbra);
  TwoBodyChannel& chket = tbs.GetChannel(ichket);
  if(i==j and chbra.J%2==1) return;
  if(k==l and chket.J%2==1) return;
  if(chbra.Prty * chket.Prty * rankP == -1) return;
  if(std::abs(chbra.Tz - chket.Tz) != rankTz) return;
  if(std::abs(chbra.J - chket.J) > rankJ or chbra.J+chket.J < rankJ) return;
  int idxbra = chbra.GetIndex(i,j);
  int idxket = chket.GetIndex(k,l);
  double phase = 1;
  phase *= chbra.GetPhaseFactor(i,j);
  phase *= chket.GetPhaseFactor(k,l);
  if(ichbra>ichket){
    Channels[{ichket,ichbra}].MEs(idxket,idxbra)  = me * pow((-1), chbra.J-chket.J) * phase;
    return;
  }
  Channels[{ichbra,ichket}].MEs(idxbra,idxket) = me * phase;
  if(ichbra != ichket) return;
  Channels[{ichbra,ichket}].MEs(idxket,idxbra) = me * phase;
}

void TwoBodyOperator::Set2BME(int ichbra, int ichket, Orbit& oi, Orbit& oj, Orbit& ok, Orbit& ol, double me)
{
  Orbits& orbits = modelspace->orbits;
  int i = orbits.GetOrbitIndex(oi);
  int j = orbits.GetOrbitIndex(oj);
  int k = orbits.GetOrbitIndex(ok);
  int l = orbits.GetOrbitIndex(ol);
  Set2BME(ichbra, ichket, i, j, k, l, me);
}

void TwoBodyOperator::Set2BME_J(int Jbra, int Jket, Orbit& oi, Orbit& oj, Orbit& ok, Orbit& ol, double me)
{
  int Prty_bra = pow((-1), oi.l+oj.l);
  int Prty_ket = pow((-1), ok.l+ol.l);
  int Tz_bra = (oi.e2 + oj.e2)/2;
  int Tz_ket = (ok.e2 + ol.e2)/2;
  TwoBodySpace & tbs = modelspace->two;
  int ichbra = tbs.GetChannelIndex(Jbra, Prty_bra, Tz_bra);
  int ichket = tbs.GetChannelIndex(Jket, Prty_ket, Tz_ket);
  Set2BME(ichbra, ichket, oi, oj, ok, ol, me);
}

void TwoBodyOperator::Set2BME_J(int Jbra, int Jket, int i, int j, int k, int l, double me)
{
  Orbits& orbits = modelspace->orbits;
  Orbit& oi = orbits.GetOrbit(i);
  Orbit& oj = orbits.GetOrbit(j);
  Orbit& ok = orbits.GetOrbit(k);
  Orbit& ol = orbits.GetOrbit(l);
  Set2BME_J(Jbra, Jket, oi, oj, ok, ol, me);
}

void TwoBodyOperator::Print()
{
  for (auto it : Channels){
    int ichbra = it.first[0];
    int ichket = it.first[1];
    TwoBodyOperatorChannel& op_ch = it.second;
    TwoBodyChannel * chbra = op_ch.chbra;
    TwoBodyChannel * chket = op_ch.chket;
    std::cout
      << "Jbra=" << std::setw(4) << chbra->J
      << ", Prty bra=" << std::setw(4) << chbra->Prty
      << ", Tz bra=" << std::setw(4) << chbra->Tz
      << " | Jket=" << std::setw(4) << chket->J
      << ", Prty ket=" << std::setw(4) << chket->Prty
      << ", Tz ket=" << std::setw(4) << chket->Tz << std::endl;
    arma::mat& m = op_ch.MEs;
    std::cout << m << std::endl;
  }
}

void TwoBodyOperator::SetTwoBodyCoulombTerm()
{
  TwoBodySpace& tbs = modelspace->GetTwoBodySpace();
  Orbits& orbits = modelspace->GetOrbits();
  for (int ich=0; ich<tbs.GetNumberChannels(); ich++){
    TwoBodyChannel& tbc = tbs.GetChannel(ich);
    for(int idxbra=0; idxbra<tbc.GetNumberStates(); idxbra++){
      for(int idxket=idxbra; idxket<tbc.GetNumberStates(); idxket++){
        int i = tbc.GetOrbitIndex1(idxbra);
        int j = tbc.GetOrbitIndex2(idxbra);
        int k = tbc.GetOrbitIndex1(idxket);
        int l = tbc.GetOrbitIndex2(idxket);
        Orbit& oi = orbits.GetOrbit(i);
        Orbit& oj = orbits.GetOrbit(j);
        Orbit& ok = orbits.GetOrbit(k);
        Orbit& ol = orbits.GetOrbit(l);
        double norm = 1.0;
        if (i==j) {norm /= sqrt(2);}
        if (k==l) {norm /= sqrt(2);}
        double v = 0;
        v = MECoulomb(oi, oj, ok, ol, tbc.J);
        v += MECoulomb(oi, oj, ol, ok, tbc.J) * pow((-1), ( (ok.j2+ol.j2)/2 - tbc.J + 1));
        v *= norm;
        Set2BME(ich, ich, i, j, k, l, v);
      }
    }
  }
}

double TwoBodyOperator::MECoulomb(Orbit& o1, Orbit& o2, Orbit& o3, Orbit& o4, int J)
{
  if (o1.e2 != o3.e2) return 0.0;
  if (o2.e2 != o4.e2) return 0.0;
  double zeta = modelspace->GetZeta();
  double Z = modelspace->GetProtonNumber();
  int Lmin = std::max(std::abs(o1.j2-o3.j2), std::abs(o2.j2-o4.j2))/2;
  int Lmax = std::min(        (o1.j2+o3.j2),         (o2.j2+o4.j2))/2;
  double rmax = 50;
  int NMesh = 200;
  gsl_integration_fixed_workspace *workspace;
  //const gsl_integration_fixed_type *T = gsl_integration_fixed_legendre;
  const gsl_integration_fixed_type *T = gsl_integration_fixed_laguerre;
  workspace = gsl_integration_fixed_alloc(T, NMesh, 0.0, 2.0/zeta, 0.0, 0.0);
  double r = 0.0;
  for (int L=Lmin; L<=Lmax; L++) {
    if ( (o1.l+o3.l+L)%2 == 1 ) continue;
    if ( (o2.l+o4.l+L)%2 == 1 ) continue;
    if (abs(o1.l-o3.l) > L or o1.l+o3.l < L) continue;
    if (abs(o2.l-o4.l) > L or o2.l+o4.l < L) continue;
    double angular = gsl_sf_coupling_6j(o1.j2, o2.j2, 2*J, o4.j2, o3.j2, 2*L) *
      gsl_sf_coupling_3j(o1.j2, 2*L, o3.j2, -1, 0, 1) *
      gsl_sf_coupling_3j(o2.j2, 2*L, o4.j2, -1, 0, 1);
    if(std::abs(angular) < 1.e-8) continue;
    double Integral = 0.0;
    for (int i=0; i<NMesh; i++){
      for (int j=0; j<NMesh; j++){
        double x = workspace->x[i];
        double wx = workspace->weights[i];
        double y = workspace->x[j];
        double wy = workspace->weights[j];
        Integral += wx * wy * o1.RadialFunction(x, zeta, Z, o1.e2, true) *
          o2.RadialFunction(y, zeta, Z, o2.e2, true) *
          o3.RadialFunction(x, zeta, Z, o3.e2, true) *
          o4.RadialFunction(y, zeta, Z, o4.e2, true) *
          pow( std::min(x,y), L) / pow( std::max(x,y), (L+1) ) * o1.e2 * o2.e2;
      }
    }
    r += angular * Integral;
  }
  r *= sqrt( (o1.j2+1) * (o2.j2+1) * (o3.j2+1) * (o4.j2+1)) * pow( (-1), (o1.j2+o3.j2)/2+J );
  return r;
}
