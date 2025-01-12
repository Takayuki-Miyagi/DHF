#pragma once
#include "Orbits.h"
#include "TwoBodySpace.h"

using namespace std;

class OneBodySpace : public TwoBodySpace{

    public:

        vector<int> kappas;
        Orbits orbits;
        vector<vector<int>> channels;
        map <int, int> orbit_index_to_channel_index;

        
        OneBodySpace(Orbits orbits=Orbits());
};


class ModelSpace : public OneBodySpace{

    public:
        OneBodySpace one;
        TwoBodySpace two;

        int Z;
        int Ne;
        int zeta;
        double c = 137.035999084;

        ModelSpace(int Ne, int Z, int zeta, Orbits orbs);
        ModelSpace(int Ne=1, int Z=1, int zeta=1);

        // Function Declerations
        void set_model_space_from_orbits(Orbits orbs);
        void set_model_space_from_truncations(int nmax, int lmax);

};
