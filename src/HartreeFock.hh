#ifndef HartreeFock_h
#define HartreeFock_h 1
#include <iostream>
#include <armadillo>
#include <cmath>

#include "Orbits.hh"
#include "Operator.hh"
#include "ModelSpace.hh"
#include "TwoBodySpace.hh"
#include "TwoBodyOperator.hh"

class Monopole
{
  public:
    Operator * H;
    ModelSpace * modelspace;
    std::vector<std::array<int,4>> idx_to_ijkl;
    std::vector<double> vmon2;

    // Constructor
    ~Monopole();
    Monopole();
    Monopole(Operator&);
    void Print();
};

class HartreeFock
{

  public:
    Operator& H;
    std::map<int, double> holes;
    ModelSpace * modelspace;
    Monopole monopole;
    arma::mat C, rho, F, V, S;
    arma::vec SPEs;
    arma::uvec ElectronStates, PositronStates;
    double r;
    double EHF;

    // Constructor
    ~HartreeFock();
    HartreeFock(Operator&, std::map<int, double>);

    void CalcEnergy();
    double UpdateFock(int n_itr = -1);
    void DiagonalizeFock();
    void UpdateDensityMatrix();
    void Solve();
    void PrintStatus(int n_itr, bool detail=true);
};
#endif
