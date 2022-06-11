#ifndef ModelSpace_h
#define ModelSpace_h 1
#include "Orbits.hh"
#include "TwoBodySpace.hh"

#include <set>
#include <unordered_set>

class OneBodySpace
{
  public:
    ~OneBodySpace();
    OneBodySpace();
    OneBodySpace(Orbits&);

    Orbits * orbits;
    std::unordered_set<int> kappas;
    std::vector<std::vector<int>> channels;
    std::map <int, int> orbit_index_to_channel_index;
    int number_channels;
    //
    int GetNumberChannels() {return number_channels;}
    void PrintSpace();
};


class ModelSpace : public OneBodySpace
{
  public:
    ~ModelSpace();
    ModelSpace(int, double, double, Orbits);
    ModelSpace(std::string, double, Orbits);

    Orbits orbits;
    OneBodySpace one;
    TwoBodySpace two;
    std::map<int,double> hole_occ;
    double Z;
    int Ne;
    double zeta;

    int GetElectronNumber() {return Ne;};
    std::map<int,double> AssignHoles(int);
    double GetProtonNumber() {return Z;};
    double GetZeta() {return zeta;};
    void PrintModelSpace(bool, bool);
    void PrintHoleOrbits();
    Orbits& GetOrbits() {return orbits;};
    OneBodySpace& GetOneBodySpace() {return one;};
    TwoBodySpace& GetTwoBodySpace() {return two;};
};
#endif
