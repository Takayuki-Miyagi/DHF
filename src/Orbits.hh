#ifndef Orbits_h
#define Orbits_h 1

#include <iostream>
#include <stdio.h>
#include <array>
#include <vector>
#include <map>
#include <math.h>
#include <gsl/gsl_integration.h>
class Orbit {
  public:
    ~Orbit();
    Orbit();
    Orbit(int n, int k, int e2);
    Orbit(int n, int l, int j2, int e2);

    // Member variables
    int n, l, j2, k, e2;
    std::string radial_function_type;

    // Class Member Declerations
    int lj2k(int l_in, int j2_in) {return floor(-(4 + j2_in*(j2_in+2) - 4*l_in*(l_in+1)-3)/4);};
    std::array<int,2> k2lj(int k);
    void SetRadialFunctionType(std::string new_radial_function_type) {radial_function_type = new_radial_function_type;};
    void PrintOrbit();

    double RadialFunction(double x, double zeta, double par, int PQ);
    double RadialFunctionD(double x, double zeta, double par, int PQ);
    double _LSpinor_P_wave_function_rspace(double x, double zeta, double Z);
    double _LSpinor_Q_wave_function_rspace(double x, double zeta, double Z);
    double _LSpinor_P_wave_function_rspace_dr(double x, double zeta, double Z);
    double _LSpinor_Q_wave_function_rspace_dr(double x, double zeta, double Z);
    std::array<double,3> _get_pars_Lspinor(double zeta, double Z);
};

class Orbits : public Orbit {
  public:
    // Member Variables
    std::vector <Orbit> orbits;
    int emax, lmax, norbs;
    bool verbose;
    std::string radial_function_type;
    std::map<std::array<int,4>, int> nlje_idx;
    std::vector<char> labels_orbital_angular_momentum;

    // Constructor
    ~Orbits();
    Orbits();
    Orbits(int, int, std::string="default");
    //Orbits(int emax1, int lmax1, std::string radial_function_type1, bool verbose1, bool relativistic1);
    //Orbits(std::string radial_function_type1="default", bool verbose1=false, bool relativistic1=true);

    // Class Function Declerations
    void AddOrbit(int, int, int, int);
    void AddOrbit(std::string value);
    void AddOrbits(std::vector <std::array<int,4>>);
    void AddOrbits(std::vector <std::string> strings);
    void AddOrbits(Orbits);
    Orbit& GetOrbit(int idx) {return orbits[idx];};
    std::string GetOrbitLabel(int);
    int GetOrbitIndex(int, int, int);
    int GetOrbitIndex(Orbit& o) {return nlje_idx[{o.n, o.l, o.j2, o.e2}];};
    int GetOrbitIndex(int n, int l, int j2, int e2) {return nlje_idx[{n,l,j2,e2}];};
    int GetNumberOrbits() {return norbs;};
    void SetOrbits(int nmax, int lmax);
    void PrintOrbits();

};



// Non-Class Function Declerations
template< typename F >  class gsl_function_pp : public gsl_function {
  public:
    gsl_function_pp(const F& func) : _func(func) {
      function = &gsl_function_pp::invoke;
      params = this;
    }
  private:
    const F& _func;
    static double invoke(double x, void* params) {
      return static_cast<gsl_function_pp*>(params)->_func(x);
    }
};
double MEOverlap(Orbit& o1, Orbit& o2, double zeta, double Z);
double MENuclPot(Orbit& o1, Orbit& o2, double zeta, double Z);
double MEKinetic(Orbit& o1, Orbit& o2, double zeta, double Z);
double NormCheck(Orbit& o1, Orbit& o2, double zeta, double par);
#endif
