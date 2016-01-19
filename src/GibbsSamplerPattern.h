#ifndef _GIBBSSAMPLERPATTER_H_
#define _GIBBSSAMPLERPATTER_H_

#include "GibbsSamplerMap.h"

class GibbsSamplerPattern : public GibbsSamplerMap {
  protected:
    int _whichPattern;

  public:
    GibbsSamplerPattern(unsigned long nEquil, unsigned long nSample, unsigned int nFactor,
                        double alphaA, double alphaP, double nMaxA, double nMaxP,
                        unsigned long nIterA, unsigned long nIterP,
                        double max_gibbsmass_paraA, double max_gibbsmass_paraP,
                        unsigned long long atomicSize,
                        char label_A, char label_P, char label_D, char label_S,
                        vector<vector<double> > &DVector, vector<vector<double> > &SVector,
                        const string &simulation_id,
                        vector <vector <double> >  &parameters, char the_fixed_matrix, int whichPattern);

    ~GibbsSamplerPattern() {};
    
    void update_pattern();
};

#endif
