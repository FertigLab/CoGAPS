#include "GibbsSamplerPattern.h"

// ******************** CONSTRUCTOR ********************************************
GibbsSamplerPattern::GibbsSamplerPattern(unsigned long nEquil, unsigned long nSample, unsigned int nFactor,
                                         double alphaA, double alphaP, double nMaxA, double nMaxP,
                                         unsigned long nIterA, unsigned long nIterP,
                                         double max_gibbsmass_paraA, double max_gibbsmass_paraP,
                                         unsigned long long atomicSize,
                                         char label_A, char label_P, char label_D, char label_S,
                                         vector<vector<double> > &DVector, vector<vector<double> > &SVector,
                                         const string &simulation_id,
                                         vector <vector <double> >  &parameters, char the_fixed_matrix, int whichPattern) :
    GibbsSamplerMap(nEquil, nSample, nFactor, alphaA, alphaP, nMaxA, nMaxP, nIterA, nIterP,
                    max_gibbsmass_paraA, max_gibbsmass_paraP, atomicSize, label_A, label_P, label_D, label_S,
                    DVector, SVector, simulation_id, parameters, the_fixed_matrix) {
    _whichPattern = whichPattern;
}
