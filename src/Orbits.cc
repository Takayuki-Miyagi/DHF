//#include <iostream>
//#include <fstream>
//#include <math.h>
//#include <cmath>
//#include <string>
//#include <regex>
//#include <string_view>
//#include <gsl/gsl_sf.h>
//#include <gsl/gsl_math.h>
//#include <armadillo>

#include <iomanip>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_sf_gamma.h>
#include "Orbits.hh"
#include "PhysicalConstants.hh"
Orbit::~Orbit()
{}

Orbit::Orbit()
  : n(-1), l(-1), j2(-1), k(0), e2(-100), radial_function_type("default")
{}

Orbit::Orbit(int n, int k, int e2)
  : n(n), k(k), e2(e2), radial_function_type("default")
{
  std::array<int,2> lj = k2lj(k);
  l = lj[0];
  j2 = lj[1];
}

Orbit::Orbit(int n, int l, int j2, int e2)
  : n(n), l(l), j2(j2), e2(e2), radial_function_type("default")
{
  k = lj2k(l,j2);
}

std::array<int,2> Orbit::k2lj(int k)
{
  if (k > 0) {
    return { k, 2 * k - 1 };
  }
  else if (k < 0) {
    return { -k - 1, -2 * k - 1 };
  }
  else {return {0,0};}
}

void Orbit::PrintOrbit()
{
  std::cout << " n =" << std::setw(6) << n << ", l =" << std::setw(6) << l << ", j2 =" << std::setw(6) << j2 << ", e2 =" << std::setw(6) << e2 <<
    ", radial wf type: " << radial_function_type << std::endl;
}

double Orbit::RadialFunction(double x, double zeta, double par, int PQ) {
  if (radial_function_type == "default") {
    if (PQ == 1) { return _LSpinor_P_wave_function_rspace(x, zeta, par); }
    if (PQ ==-1) { return _LSpinor_Q_wave_function_rspace(x, zeta, par); }
  }
  return 0;
}

double Orbit::RadialFunctionD(double x, double zeta, double par, int PQ) {
  if (radial_function_type == "default") {
    if (PQ ==  1) { return _LSpinor_P_wave_function_rspace_dr(x, zeta, par); }
    if (PQ == -1) { return _LSpinor_Q_wave_function_rspace_dr(x, zeta, par); }
  }
  return 0;
}

double Orbit::_LSpinor_P_wave_function_rspace(double x, double zeta, double Z) {
  double eta = 2.0 * x / zeta;
  double gam = _get_pars_Lspinor(zeta, Z)[0];
  double N = _get_pars_Lspinor(zeta, Z)[1];
  double Norm = _get_pars_Lspinor(zeta, Z)[2];
  if (Norm == 0.0) { return 0.0; }
  double T = gsl_sf_laguerre_n(n, 2 * gam, eta) * (N - k) / (n + 2 * gam);
  if (n > 0) { T -= gsl_sf_laguerre_n(n - 1, 2 * gam, eta); }
  return T * pow(eta, gam) * exp(-0.5 * eta) * Norm;
}

double Orbit::_LSpinor_Q_wave_function_rspace(double x, double zeta, double Z) {
  double eta = 2.0 * x / zeta;
  double gam = _get_pars_Lspinor(zeta, Z)[0];
  double N = _get_pars_Lspinor(zeta, Z)[1];
  double Norm = _get_pars_Lspinor(zeta, Z)[2];
  if (Norm == 0.0) { return 0.0; }
  double T = -gsl_sf_laguerre_n(n, 2 * gam, eta) * (N - k) / (n + 2 * gam);
  if (n > 0) { T -= gsl_sf_laguerre_n(n - 1, 2 * gam, eta); }
  return T * pow(eta, gam) * exp(-0.5 * eta) * Norm;
}

double Orbit::_LSpinor_P_wave_function_rspace_dr(double x, double zeta, double Z) {
  double eta = 2.0 * x / zeta;
  double gam = _get_pars_Lspinor(zeta, Z)[0];
  double N = _get_pars_Lspinor(zeta, Z)[1];
  double Norm = _get_pars_Lspinor(zeta, Z)[2];
  if (Norm == 0.0) { return 0.0; }
  double T1 = _LSpinor_P_wave_function_rspace(x, zeta, gam) * ((gam / eta) - 0.5);
  double T2;
  if (n > 0) {
    T2 = -gsl_sf_laguerre_n(n - 1, 2 * gam + 1, eta) * (N - k) / (n + 2 * gam);
    if (n > 1) {
      T2 += gsl_sf_laguerre_n(n - 2, 2 * gam + 1, eta);
    }
  } else {
    T2 = 0.0;
  }
  return 2.0 / zeta * (T1 + Norm * T2 * pow(eta, gam) * exp(-0.5 * eta));
}

double Orbit::_LSpinor_Q_wave_function_rspace_dr(double x, double zeta, double Z) {
  double eta = 2.0 * x / zeta;
  double gam = _get_pars_Lspinor(zeta, Z)[0];
  double N = _get_pars_Lspinor(zeta, Z)[1];
  double Norm = _get_pars_Lspinor(zeta, Z)[2];
  if (Norm == 0.0) { return 0.0; }
  double T1 = _LSpinor_Q_wave_function_rspace(x, zeta, gam) * ((gam / eta) - 0.5);
  double T2;
  if (n > 0) {
    T2 = gsl_sf_laguerre_n(n - 1, 2 * gam + 1, eta) * (N - k) / (n + 2 * gam);
    if (n > 1) {
      T2 += gsl_sf_laguerre_n(n - 2, 2 * gam + 1, eta);
    }
  } else {
    T2 = 0.0;
  }
  return 2.0 / zeta * (T1 + Norm * T2 * pow(eta, gam) * exp(-0.5 * eta));
}

std::array<double,3> Orbit::_get_pars_Lspinor(double zeta, double Z) {
  double gam = sqrt(pow(k, 2) - pow(Z, 2) / pow(PhysConst::c, 2));
  double N = sqrt(pow(n, 2) + 2.0 * n * gam + pow(k, 2));
  double Norm;
  if (N - k == 0) {
    Norm = 0.0;
  }
  else {
    Norm = sqrt((tgamma(n + 1) * (2 * gam + n)) / (2 * N * (N - k) * tgamma(2 * gam + n) * zeta));
  }
  return { gam, N, Norm };
}

//
// Orbits class
//
Orbits::~Orbits()
{}

Orbits::Orbits()
  : emax(-1), lmax(-1), norbs(0), radial_function_type("default"), verbose(true),
  labels_orbital_angular_momentum({ 's', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n', 'o', 'q', 'r', 't', 'u', 'v', 'w', 'x', 'y', 'z' })
{}

Orbits::Orbits(int nmax, int lmax, std::string radial_function_type)
  : emax(nmax+lmax), lmax(lmax), radial_function_type(radial_function_type), verbose(true),
  labels_orbital_angular_momentum({ 's', 'p', 'd', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n', 'o', 'q', 'r', 't', 'u', 'v', 'w', 'x', 'y', 'z' })
{
  std::array<int,2> klist = {1,-1};
  for ( int li=0; li < lmax+1; li++ ){
    for ( int n =li; n < nmax+1; n++){
      for (int j : {2*li-1, 2*li+1} ) {
        if ( j<0 ){ continue; }
        for ( int e : klist ) {
          AddOrbit(n, li, j, e);
        }
      }
    }
  }
}

void Orbits::AddOrbit(int n, int l, int j, int e){
  std::array<int,4> nlje = {n,l,j,e};
  for (int i : nlje ){
  //  if (verbose==true) { std::cout << "The orbit " << nlje[0] << nlje[1] << nlje[2] << nlje[3] << " is already there." << std::endl; }
  }
  norbs = orbits.size();
  int idx = norbs;
  nlje_idx[nlje] = idx;
  Orbit orb = Orbit(nlje[0], nlje[1], nlje[2], nlje[3]);
  orb.SetRadialFunctionType(radial_function_type);
  orbits.push_back(orb);
  norbs = orbits.size();
  emax = std::max(emax, nlje[0]+nlje[1]);
  lmax = std::max(lmax, nlje[1]);
}

void Orbits::AddOrbit(std::string value) {

  // string format should be like l0s1 => large-compoennt 0s1/2, s0s1 => small-component 0s1/2

  int z;
  char pn = value.at(0);
  if ( pn == 's' ) { z = -1; }
  else if ( pn == 'l' ) { z = 1; }
  else{
    std::cout << "parse error in add_orbit_from_label: " << value << std::endl;
    return;
  }
  std::string nlj_str = value.erase(0,1);

  char n_str = value.at(0);
  char l_str = value.at(1);
  char j_str = value.at(2);

  n = int(n_str);
  int l = 0;
  for ( char l_label : labels_orbital_angular_momentum ) {
    if ( l_str == l_label ) { break; }
    l += 1;
  }
  int j = int(j_str);
  AddOrbit(n, l, j, z);
}

void Orbits::AddOrbits(std::vector <std::array<int,4>> nlje_list){
  for (auto nlje : nlje_list) { AddOrbit(nlje[0], nlje[1], nlje[2], nlje[3]); }
}

void Orbits::AddOrbits(std::vector <std::string> strings){
  for (std::string label : strings) { AddOrbit(label); }
}

void Orbits::PrintOrbits()
{
  for (auto o: orbits)
  {
    o.PrintOrbit();
  }
}

int Orbits::GetOrbitIndex(int n, int k, int e2)
{
  int l = 0;
  int j2 = 0;
  if (k > 0) {
    l = k;
    j2 = 2*k-1;
  }
  else if (k < 0) {
    l = -k-1;
    j2 = -2*k-1;
  }
  return GetOrbitIndex(n, l, j2, e2);
}



size_t limit = 100;
double epsabs = 1.e-6;
double epsrel = 1.e-10;

//overlap
double MEOverlap(Orbit& o1, Orbit& o2, double zeta, double Z) {
  double result, error;
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);

  if(o1.l != o2.l) return 0;
  if(o1.j2 != o2.j2) return 0;
  if (o1.e2 == 1 && o2.e2 == 1) {
    auto lambda = [&](double x) { return (o1.RadialFunction(x, zeta, Z, o1.e2) * o2.RadialFunction(x, zeta, Z, o2.e2)); };
    gsl_function_pp<decltype(lambda)> Fp(lambda);
    gsl_function* F = static_cast<gsl_function*>(&Fp);
    gsl_integration_qagiu(F, 0, epsabs, epsrel, limit, w, &result, &error);
    gsl_integration_workspace_free(w);
    return result;

  } else if ((o1.e2 == -1 && o2.e2 == -1)){
    auto lambda = [&](double x) {return (o1.RadialFunction(x, zeta, Z, o1.e2) * o2.RadialFunction(x, zeta, Z, o2.e2)); };
    gsl_function_pp<decltype(lambda)> Fp(lambda);
    gsl_function* F = static_cast<gsl_function*>(&Fp);
    gsl_integration_qagiu(F, 0, epsabs, epsrel, limit, w, &result, &error);
    gsl_integration_workspace_free(w);
    return result;

  } else { return 0.0; }
}

// 1 / r
double MENuclPot(Orbit& o1, Orbit& o2, double zeta, double Z) {
  double result, error;
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);

  if(o1.l != o2.l) return 0;
  if(o1.j2 != o2.j2) return 0;
  if (o1.e2 == 1 && o2.e2 == 1) {

    auto lambda = [&](double x) { return (o1.RadialFunction(x, zeta, Z, o1.e2) * o2.RadialFunction(x, zeta, Z, o2.e2) / x); };
    gsl_function_pp<decltype(lambda)> Fp(lambda);
    gsl_function* F = static_cast<gsl_function*>(&Fp);
    gsl_integration_qagiu(F, 0, epsabs, epsrel, limit, w, &result, &error);
    gsl_integration_workspace_free(w);
    return result;
  } else if (o1.e2 == -1 && o2.e2 == -1) {

    auto lambda = [&](double x) { return (o1.RadialFunction(x, zeta, Z, o1.e2) * o2.RadialFunction(x, zeta, Z, o2.e2) / x); };
    gsl_function_pp<decltype(lambda)> Fp(lambda);
    gsl_function* F = static_cast<gsl_function*>(&Fp);
    gsl_integration_qagiu(F, 0, epsabs, epsrel, limit, w, &result, &error);
    gsl_integration_workspace_free(w);
    return result;
  }
  else { return 0.0; }
}


//s dot p
double MEKinetic(Orbit& o1, Orbit& o2, double zeta, double Z) {
  double result, error;
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
  if(o1.l != o2.l) return 0;
  if(o1.j2 != o2.j2) return 0;

  if (o1.e2 == 1 && o2.e2 == -1) {
    auto lambda = [&](double x) { return (-o1.RadialFunction(x, zeta, Z, o1.e2) * (o2.RadialFunctionD(x, zeta, Z, o2.e2) - o2.k * o2.RadialFunction(x, zeta, Z, o2.e2) / x)); };
    gsl_function_pp<decltype(lambda)> Fp(lambda);
    gsl_function* F = static_cast<gsl_function*>(&Fp);
    gsl_integration_qagiu(F, 0, epsabs, epsrel, limit, w, &result, &error);
    gsl_integration_workspace_free(w);
    return result * PhysConst::c;

  } else if (o1.e2 == -1 && o2.e2 == 1) {
    auto lambda = [&](double x) { return (o1.RadialFunction(x, zeta, Z, o1.e2) * (o2.RadialFunctionD(x, zeta, Z, o2.e2) + o2.k * o2.RadialFunction(x, zeta, Z, o2.e2) / x)); };
    gsl_function_pp<decltype(lambda)> Fp(lambda);
    gsl_function* F = static_cast<gsl_function*>(&Fp);
    gsl_integration_qagiu(F, 0, epsabs, epsrel, limit, w, &result, &error);
    gsl_integration_workspace_free(w);
    return result * PhysConst::c;

  } else { return 0.0; }
}

double NormCheck(Orbit& o1, Orbit& o2, double zeta, double par) {
  if (o1.l  != o2.l ) return 0.0;
  if (o1.j2 != o2.j2) return 0.0;
  if (o1.e2 != o2.e2) return 0.0;

  int PQ;
  double result1, error1, result2, error2;
  gsl_integration_workspace* w = gsl_integration_workspace_alloc(limit);
  gsl_integration_workspace* m = gsl_integration_workspace_alloc(limit);
  auto lambda1 = [&](double x) { return (o1.RadialFunction(x, zeta, par, o1.e2) * o2.RadialFunction(x, zeta, par, o2.e2)); };
  gsl_function_pp<decltype(lambda1)> Fp(lambda1);
  gsl_function* F = static_cast<gsl_function*>(&Fp);
  gsl_integration_qagiu(F, 0, epsabs, epsrel, limit, w, &result1, &error1);
  gsl_integration_workspace_free(w);

  auto lambda2 = [&](double x) { return (o1.RadialFunction(x,zeta,par,-o1.e2) * o2.RadialFunction(x,zeta,par,-o2.e2)); };
  gsl_function_pp<decltype(lambda2)> Xp(lambda2);
  gsl_function* X = static_cast<gsl_function*>(&Xp);
  gsl_integration_qagiu(F, 0, epsabs, epsrel, limit, w, &result2, &error2);
  gsl_integration_workspace_free(m);

  return result1+result2;
}

