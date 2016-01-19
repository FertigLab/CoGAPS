#include "GibbsSamplerMap.h"

// ******************** CONSTRUCTOR ********************************************
GibbsSamplerMap::GibbsSamplerMap(unsigned long nEquil, unsigned long nSample, unsigned int nFactor,
                                 double alphaA, double alphaP, double nMaxA, double nMaxP,
                                 unsigned long nIterA, unsigned long nIterP,
                                 double max_gibbsmass_paraA, double max_gibbsmass_paraP,
                                 unsigned long long atomicSize,
                                 char label_A, char label_P, char label_D, char label_S,
                                 vector<vector<double> > &DVector, vector<vector<double> > &SVector,
                                 const string &simulation_id,
                                 vector <vector <double> >  &parameters, char the_fixed_matrix)
    : GibbsSampler(nEquil, nSample, nFactor,
                   alphaA, alphaP, nMaxA, nMaxP, nIterA, nIterP, max_gibbsmass_paraA,
                   max_gibbsmass_paraP, atomicSize,
                   label_A, label_P, label_D, label_S, DVector, SVector, simulation_id) {
    //Constructor body
    for (int i = 0; i < parameters.size(); i++) {
        vector <double> newPattern;
        newPattern.clear();
        double PatSum = 0;

        for (int j = 0; j < parameters[0].size(); j++) {
            PatSum += parameters[i][j];
        }

        for (int j = 0; j < parameters[0].size(); j++) {
            newPattern.push_back(parameters[i][j] / PatSum);
        }

        _MapValues.push_back(newPattern);
    }

    _nFixedMaps = parameters.size();
    _the_fixed_matrix = the_fixed_matrix;
}
