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

void Orbit::Print()
{
  int width = 6;
  std::cout << " n =" << std::setw(width) << n << ", l =" << std::setw(width) << l << 
    ", j2 =" << std::setw(width) << j2 << 
    ", kappa =" << std::setw(width) << k << 
    ", e2 =" << std::setw(width) << e2 << 
    ", radial wf type: " << radial_function_type << std::endl;
}

double Orbit::RadialFunction(double x, double zeta, double par, int PQ, bool exp_weight) {
  if (radial_function_type == "default") {
    if(PQ == 1) return _LSpinor_P_wave_function_rspace(x, zeta, par, exp_weight); 
    if(PQ ==-1) return _LSpinor_Q_wave_function_rspace(x, zeta, par, exp_weight); 
  }
  else if (radial_function_type=="NonRel_Laguerre"){
    return _Laguerre_wave_function(x, zeta, exp_weight);
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

double Orbit::_Laguerre_wave_function(double x, double zeta, bool exp_weight) {
  //Laguerre function, see[A.E.McCoy and M.A.Caprio, J.Math.Phys. 57, (2016).] for details
  double eta = 2.0 * x / zeta;
  if(exp_weight) return sqrt(2.0 * tgamma(n + 1) / (zeta * tgamma(n + 2 * l + 3))) * 2.0 * pow(eta, l) * gsl_sf_laguerre_n(n, 2.0 * l + 2.0, eta) / zeta * x;
  return sqrt(2.0 * tgamma(n + 1) / (zeta * tgamma(n + 2 * l + 3))) * 2.0 * pow(eta, l) * exp(-0.5 * eta) *  gsl_sf_laguerre_n(n, 2.0 * l + 2.0, eta) / zeta * x;
}

double Orbit::_LSpinor_P_wave_function_rspace(double x, double zeta, double Z, bool exp_weight) {
  double eta = 2.0 * x / zeta;
  double gam = _get_pars_Lspinor(zeta, Z)[0];
  double N = _get_pars_Lspinor(zeta, Z)[1];
  double Norm = _get_pars_Lspinor(zeta, Z)[2];
  if (Norm == 0.0) return 0.0;
  double T = gsl_sf_laguerre_n(n+l, 2 * gam, eta) * (N - k) / (n + l + 2 * gam);
  if (n+l > 0) T -= gsl_sf_laguerre_n(n + l - 1, 2 * gam, eta); 
  if(exp_weight) return T * pow(eta, gam) * Norm;
  return T * pow(eta, gam) * exp(-0.5 * eta) * Norm;
}

double Orbit::_LSpinor_Q_wave_function_rspace(double x, double zeta, double Z, bool exp_weight) {
  double eta = 2.0 * x / zeta;
  double gam = _get_pars_Lspinor(zeta, Z)[0];
  double N = _get_pars_Lspinor(zeta, Z)[1];
  double Norm = _get_pars_Lspinor(zeta, Z)[2];
  if (Norm == 0.0) return 0.0; 
  double T = -gsl_sf_laguerre_n(n+l, 2 * gam, eta) * (N - k) / (n+l + 2 * gam);
  if (n+l > 0) T -= gsl_sf_laguerre_n(n+l - 1, 2 * gam, eta); 
  if(exp_weight) return T * pow(eta, gam) * Norm;
  return T * pow(eta, gam) * exp(-0.5 * eta) * Norm;
}

double Orbit::_LSpinor_P_wave_function_rspace_dr(double x, double zeta, double Z) {
  double eta = 2.0 * x / zeta;
  double gam = _get_pars_Lspinor(zeta, Z)[0];
  double N = _get_pars_Lspinor(zeta, Z)[1];
  double Norm = _get_pars_Lspinor(zeta, Z)[2];
  if (Norm == 0.0) return 0.0; 
  double T1 = _LSpinor_P_wave_function_rspace(x, zeta, gam) * ((gam / eta) - 0.5);
  double T2;
  if (n+l > 0) {
    T2 = -gsl_sf_laguerre_n(n + l - 1, 2 * gam + 1, eta) * (N - k) / (n + l + 2 * gam);
    if (n+l > 1) {
      T2 += gsl_sf_laguerre_n(n + l - 2, 2 * gam + 1, eta);
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
  if (Norm == 0.0) return 0.0;
  double T1 = _LSpinor_Q_wave_function_rspace(x, zeta, gam) * ((gam / eta) - 0.5);
  double T2;
  if (n+l > 0) {
    T2 = gsl_sf_laguerre_n(n + l - 1, 2 * gam + 1, eta) * (N - k) / (n + l + 2 * gam);
    if (n+l > 1) {
      T2 += gsl_sf_laguerre_n(n + l - 2, 2 * gam + 1, eta);
    }
  } else {
    T2 = 0.0;
  }
  return 2.0 / zeta * (T1 + Norm * T2 * pow(eta, gam) * exp(-0.5 * eta));
}

std::array<double,3> Orbit::_get_pars_Lspinor(double zeta, double Z) {
  double gam = sqrt(pow(k, 2) - pow(Z, 2) / pow(PhysConst::c, 2));
  double N = sqrt(pow(n+l, 2) + 2.0 * (n+l) * gam + pow(k, 2));
  double Norm;
  if (N - k == 0) {
    Norm = 0.0;
  }
  else {
    Norm = sqrt((tgamma(n + l + 1) * (2 * gam + n + l)) / (2 * N * (N - k) * tgamma(2 * gam + n + l) * zeta));
  }
  return { gam, N, Norm };
}

//double Orbit::_LSpinor_P_wave_function_rspace(double x, double zeta, double Z) {
//  double eta = 2.0 * x / zeta;
//  double gam = _get_pars_Lspinor(zeta, Z)[0];
//  double N = _get_pars_Lspinor(zeta, Z)[1];
//  double Norm = _get_pars_Lspinor(zeta, Z)[2];
//  if (Norm == 0.0) { return 0.0; }
//  double T = gsl_sf_laguerre_n(n, 2 * gam, eta) * (N - k) / (n + 2 * gam);
//  if (n > 0) { T -= gsl_sf_laguerre_n(n - 1, 2 * gam, eta); }
//  return T * pow(eta, gam) * exp(-0.5 * eta) * Norm;
//}
//
//double Orbit::_LSpinor_Q_wave_function_rspace(double x, double zeta, double Z) {
//  double eta = 2.0 * x / zeta;
//  double gam = _get_pars_Lspinor(zeta, Z)[0];
//  double N = _get_pars_Lspinor(zeta, Z)[1];
//  double Norm = _get_pars_Lspinor(zeta, Z)[2];
//  if (Norm == 0.0) { return 0.0; }
//  double T = -gsl_sf_laguerre_n(n, 2 * gam, eta) * (N - k) / (n + 2 * gam);
//  if (n > 0) { T -= gsl_sf_laguerre_n(n - 1, 2 * gam, eta); }
//  return T * pow(eta, gam) * exp(-0.5 * eta) * Norm;
//}
//
//double Orbit::_LSpinor_P_wave_function_rspace_dr(double x, double zeta, double Z) {
//  double eta = 2.0 * x / zeta;
//  double gam = _get_pars_Lspinor(zeta, Z)[0];
//  double N = _get_pars_Lspinor(zeta, Z)[1];
//  double Norm = _get_pars_Lspinor(zeta, Z)[2];
//  if (Norm == 0.0) { return 0.0; }
//  double T1 = _LSpinor_P_wave_function_rspace(x, zeta, gam) * ((gam / eta) - 0.5);
//  double T2;
//  if (n > 0) {
//    T2 = -gsl_sf_laguerre_n(n - 1, 2 * gam + 1, eta) * (N - k) / (n + 2 * gam);
//    if (n > 1) {
//      T2 += gsl_sf_laguerre_n(n - 2, 2 * gam + 1, eta);
//    }
//  } else {
//    T2 = 0.0;
//  }
//  return 2.0 / zeta * (T1 + Norm * T2 * pow(eta, gam) * exp(-0.5 * eta));
//}
//
//double Orbit::_LSpinor_Q_wave_function_rspace_dr(double x, double zeta, double Z) {
//  double eta = 2.0 * x / zeta;
//  double gam = _get_pars_Lspinor(zeta, Z)[0];
//  double N = _get_pars_Lspinor(zeta, Z)[1];
//  double Norm = _get_pars_Lspinor(zeta, Z)[2];
//  if (Norm == 0.0) { return 0.0; }
//  double T1 = _LSpinor_Q_wave_function_rspace(x, zeta, gam) * ((gam / eta) - 0.5);
//  double T2;
//  if (n > 0) {
//    T2 = gsl_sf_laguerre_n(n - 1, 2 * gam + 1, eta) * (N - k) / (n + 2 * gam);
//    if (n > 1) {
//      T2 += gsl_sf_laguerre_n(n - 2, 2 * gam + 1, eta);
//    }
//  } else {
//    T2 = 0.0;
//  }
//  return 2.0 / zeta * (T1 + Norm * T2 * pow(eta, gam) * exp(-0.5 * eta));
//}
//
//std::array<double,3> Orbit::_get_pars_Lspinor(double zeta, double Z) {
//  double gam = sqrt(pow(k, 2) - pow(Z, 2) / pow(PhysConst::c, 2));
//  double N = sqrt(pow(n, 2) + 2.0 * n * gam + pow(k, 2));
//  double Norm;
//  if (N - k == 0) {
//    Norm = 0.0;
//  }
//  else {
//    Norm = sqrt((tgamma(n + 1) * (2 * gam + n)) / (2 * N * (N - k) * tgamma(2 * gam + n) * zeta));
//  }
//  return { gam, N, Norm };
//}



//
// Orbits class
//
Orbits::~Orbits()
{}

Orbits::Orbits()
  : lmax(-1), radial_function_type("default"), relativistic(true), orbitals("none"), 
  labels_orbital_angular_momentum({ "s", "p", "d", "f", "g", "h", "i", "k", "l", "m", "n", "o", "q", "r", "t", "u", "v", "w", "x", "y", "z" })
{}

Orbits::Orbits(int nmax, int l, std::string radial_function_type, bool relativistic)
  : lmax(-1), radial_function_type(radial_function_type), relativistic(relativistic),
  labels_orbital_angular_momentum({ "s", "p", "d", "f", "g", "h", "i", "k", "l", "m", "n", "o", "q", "r", "t", "u", "v", "w", "x", "y", "z" })
{
  std::vector<int> klist = {1,-1};
  if(not relativistic) klist = {1};
  for ( int n = 0; n < nmax+1; n++){
    for (int j : {2*l-1, 2*l+1} ) {
      if ( j<0 ){ continue; }
      for ( int e : klist ) {
        AddOrbit(n, l, j, e);
      }
    }
  }
}

Orbits::Orbits(std::string orbitals, std::string radial_function_type)
  : lmax(-1), radial_function_type(radial_function_type), orbitals(orbitals), 
  labels_orbital_angular_momentum({ "s", "p", "d", "f", "g", "h", "i", "k", "l", "m", "n", "o", "q", "r", "t", "u", "v", "w", "x", "y", "z" })
{
  // orbitals format : rel-s(num)-p(num)-d(num)-f(num)-g(num)-... or nonrel-s(num)-p(num)-d(num)-...
  
  std::string sep = "-";
  int sep_length = sep.length();
  std::vector<std::string> ss;
  auto offset = std::string::size_type(0);
  while (1) {
    auto pos = orbitals.find(sep, offset);
    if (pos == std::string::npos) {
      ss.push_back(orbitals.substr(offset));
      break;
    }
    ss.push_back(orbitals.substr(offset, pos - offset));
    offset = pos + sep_length;
  }

  auto tmp = ss.begin();
  if(*tmp=="rel" or *tmp=="Rel") relativistic = true;
  if(*tmp=="nonrel" or *tmp=="NonRel") relativistic = false;
  for(std::string str : std::vector<std::string>(std::next(tmp,1), ss.end())){
    auto it = std::find(labels_orbital_angular_momentum.begin(), labels_orbital_angular_momentum.end(), str.substr(0,1));
    int l = it - labels_orbital_angular_momentum.begin();
    int n = stoi(str.substr(1))-1;
    Orbits tmp = Orbits(n, l, radial_function_type, relativistic);
    AddOrbits(tmp);
  }
}

void Orbits::AddOrbit(int n, int l, int j, int e){
  std::array<int,4> nlje = {n,l,j,e};
  //for (int i : nlje ){
  //  if (verbose==true) { std::cout << "The orbit " << nlje[0] << nlje[1] << nlje[2] << nlje[3] << " is already there." << std::endl; }
  //}
  int idx = orbits.size();
  nlje_idx[nlje] = idx;
  Orbit orb = Orbit(nlje[0], nlje[1], nlje[2], nlje[3]);
  orb.SetRadialFunctionType(radial_function_type);
  orbits.push_back(orb);
  lmax = std::max(lmax, nlje[1]);
}

void Orbits::AddOrbit(std::string value) {

  // string format should be like l0s1 => large-compoennt 0s1/2, s0s1 => small-component 0s1/2

  int z;
  std::string pn = value.substr(0,1);
  if ( pn == "s" ) { z = -1; }
  else if ( pn == "l" ) { z = 1; }
  else{
    std::cout << "parse error in add_orbit_from_label: " << value << std::endl;
    return;
  }
  std::string nlj_str = value.erase(0,1);

  std::string n_str = value.substr(0,1);
  std::string l_str = value.substr(1,2);
  std::string j_str = value.substr(2,3);

  n = stoi(n_str);
  int l = 0;
  for ( auto l_label : labels_orbital_angular_momentum ) {
    if ( l_str == l_label ) { break; }
    l += 1;
  }
  int j = stoi(j_str);
  AddOrbit(n, l, j, z);
}

void Orbits::AddOrbits(Orbits orbits){
  for (auto o : orbits.orbits) {AddOrbit(o.n, o.l, o.j2, o.e2);}
}

void Orbits::AddOrbits(std::vector <std::array<int,4>> nlje_list){
  for (auto nlje : nlje_list) { AddOrbit(nlje[0], nlje[1], nlje[2], nlje[3]); }
}

void Orbits::AddOrbits(std::vector <std::string> strings){
  for (std::string label : strings) { AddOrbit(label); }
}

void Orbits::Print()
{
  for (auto o: orbits)
  {
    o.Print();
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

//size_t limit = 1000;
//double epsabs = 1.e-12;
//double epsrel = 1.e-12;

// This should be somewhat same as scipy.integral.quad
size_t limit = 50;
double epsabs = 1.49e-8;
double epsrel = 1.49e-8;

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
    //std::cout << result << " " << error << std::endl;
    return result;

  } else if ((o1.e2 == -1 && o2.e2 == -1)){
    auto lambda = [&](double x) {return (o1.RadialFunction(x, zeta, Z, o1.e2) * o2.RadialFunction(x, zeta, Z, o2.e2)); };
    gsl_function_pp<decltype(lambda)> Fp(lambda);
    gsl_function* F = static_cast<gsl_function*>(&Fp);
    gsl_integration_qagiu(F, 0, epsabs, epsrel, limit, w, &result, &error);
    gsl_integration_workspace_free(w);
    //std::cout << result << " " << error << std::endl;
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
    //std::cout << result << " " << error << std::endl;
    return result;
  } else if (o1.e2 == -1 && o2.e2 == -1) {

    auto lambda = [&](double x) { return (o1.RadialFunction(x, zeta, Z, o1.e2) * o2.RadialFunction(x, zeta, Z, o2.e2) / x); };
    gsl_function_pp<decltype(lambda)> Fp(lambda);
    gsl_function* F = static_cast<gsl_function*>(&Fp);
    gsl_integration_qagiu(F, 0, epsabs, epsrel, limit, w, &result, &error);
    gsl_integration_workspace_free(w);
    //std::cout << result << " " << error << std::endl;
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

  if(o1.radial_function_type=="NonRel_Laguerre"){
    if(o1.n==o2.n) return 0.5*(4*o1.n+2*o1.l+3) / (2*o1.l+3)/zeta/zeta;
    if(o1.n<=o2.n-1) return 0.5*(4*o1.n+4*o1.l+6)/(2*o1.l+3) * 
      sqrt(tgamma(o2.n+1)*tgamma(o1.n+2*o1.l+3) / (tgamma(o1.n+1)*tgamma(o2.n+2*o1.l+3)))/zeta/zeta;
    if(o1.n>=o2.n+1) return 0.5*(4*o2.n+4*o1.l+6)/(2*o1.l+3) * 
      sqrt(tgamma(o1.n+1)*tgamma(o2.n+2*o1.l+3) / (tgamma(o2.n+1)*tgamma(o1.n+2*o1.l+3)))/zeta/zeta;
  }

  if (o1.e2 == 1 && o2.e2 == -1) {
    auto lambda = [&](double x) { return (-o1.RadialFunction(x, zeta, Z, o1.e2) * (o2.RadialFunctionD(x, zeta, Z, o2.e2) - o2.k * o2.RadialFunction(x, zeta, Z, o2.e2) / x)); };
    gsl_function_pp<decltype(lambda)> Fp(lambda);
    gsl_function* F = static_cast<gsl_function*>(&Fp);
    gsl_integration_qagiu(F, 0, epsabs, epsrel, limit, w, &result, &error);
    gsl_integration_workspace_free(w);
    //std::cout << result << " " << error << std::endl;
    return result * PhysConst::c;

  } else if (o1.e2 == -1 && o2.e2 == 1) {
    auto lambda = [&](double x) { return (o1.RadialFunction(x, zeta, Z, o1.e2) * (o2.RadialFunctionD(x, zeta, Z, o2.e2) + o2.k * o2.RadialFunction(x, zeta, Z, o2.e2) / x)); };
    gsl_function_pp<decltype(lambda)> Fp(lambda);
    gsl_function* F = static_cast<gsl_function*>(&Fp);
    gsl_integration_qagiu(F, 0, epsabs, epsrel, limit, w, &result, &error);
    gsl_integration_workspace_free(w);
    //std::cout << result << " " << error << std::endl;
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

