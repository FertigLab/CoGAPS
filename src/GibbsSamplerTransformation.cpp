#include "GibbsSamplerTransformation.h"

// ******************** CONSTRUCTOR ********************************************
GibbsSamplerTransformation::GibbsSamplerTransformation(unsigned long nEquil, unsigned long nSample, unsigned int nFactor,
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

std::vector<double> GibbsSamplerTransformation::logit(std::vector<double> data) {
    std::vector<double> transformation(data.size());

    for (int i = 0; i < data.size(); ++i) {
        transformation[i] = log(data[i] / (1 - data[i]));
    }

    return transformation;
}

std::vector<double> GibbsSamplerTransformation::identity(std::vector<double> data) {
    return data;
}
