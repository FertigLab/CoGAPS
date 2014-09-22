#ifndef _GIBBSSAMPLERMAP_H_
#define _GIBBSSAMPLERMAP_H_

#include<iostream>
#include<fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "GibbsSampler.h"
#include<limits>

using namespace gaps;

class GibbsSamplerMap : public GibbsSampler
{
 protected:
  
  vector <vector <double> > _MapValues; // list of stored doubles (size of elements in map)
  unsigned int _nFixedMaps; // the number of maps
  char _the_fixed_matrix; // which matrix has fixed patterns
 
 public:

  // ******************** CONSTRUCTOR ********************************************
			   
  GibbsSamplerMap(unsigned long nEquil, unsigned long nSample, unsigned int nFactor, 
	       double alphaA, double alphaP, double nMaxA, double nMaxP,
	       unsigned long nIterA, unsigned long nIterP, 
	       double max_gibbsmass_paraA, double max_gibbsmass_paraP, 
	       double lambdaA_scale_factor, double lambdaP_scale_factor, 
               unsigned long atomicSize,
	       char label_A,char label_P,char label_D,char label_S,
	       const string & datafile, const string & variancefile,
               const string & simulation_id, 
			   vector <vector <double> >  &parameters, char the_fixed_matrix);
			   
  GibbsSamplerMap(unsigned long nEquil, unsigned long nSample, unsigned int nFactor, 
	       double alphaA, double alphaP, double nMaxA, double nMaxP,
	       unsigned long nIterA, unsigned long nIterP, 
	       double max_gibbsmass_paraA, double max_gibbsmass_paraP, 
	       double lambdaA_scale_factor, double lambdaP_scale_factor, 
               unsigned long atomicSize,
	       char label_A,char label_P,char label_D,char label_S,
	       vector<vector<double> > &DVector, vector<vector<double> > &SVector,
               const string & simulation_id,
			   vector <vector <double> >  &parameters, char the_fixed_matrix);			   
   
  ~GibbsSamplerMap(){};


  // *************** METHODS FOR INITIALIZATION, DISPLAY, OUTPUT ***********************
	
	// for initializing the correct matrix with the fixed pattern
    void init_Mapped_Matrix();
	
	// for keeping the atomic domain consistent with the initialized matrix
	void initialize_atomic_domain_map();
	
	//for debugging
	void print_A_and_P();

  // **************** METHODS FOR COMPUTING LIKELIHOOD FUNCTIONS *****************

  double computeDeltaLLBDMap(char the_matrix_label,
			double const * const * D,
			double const * const * S,
			double const * const * A,
			double const * const * P,
			vector <double> &newPat, unsigned int chPat);
			
  double computeDeltaLLMEMap(char the_matrix_label,
				    double const * const * D,
				    double const * const * S,
				    double const * const * A,
				    double const * const * P,
				    vector <double> &newPat1, unsigned int chPat1,
					vector <double> &newPat2, unsigned int chPat2);


  // *************** METHODS FOR MAKING PROPOSAL *********************************

  void mapUpdate(char the_matrix_label);
  
  bool Q_fixed(unsigned long long location, char the_matrix_label);
    
  void mappedDeath(char the_matrix_label, 	
				      double const * const * D,
				      double const * const * S,
				      double ** AOrig,
				      double ** POrig);
					  
  void mappedBirth(char the_matrix_label, 	
				      double const * const * D,
				      double const * const * S,
				      double ** AOrig,
				      double ** POrig);
  
  void mappedMove(char the_matrix_label, 	
				      double const * const * D,
				      double const * const * S,
				      double ** AOrig,
				      double ** POrig);

  void mappedExchange(char the_matrix_label,
					double const * const * D,
					double const * const * S,
					double ** AOrig,
					double ** POrig);
					  
 // *************** METHODS FOR UPDATING THE MATRIX ***********************************	
 
  bool calc_new_matrix_Pattern(char the_matrix_label, vector <double> &PatternUpdate, 
										unsigned long long location, double mass);
  
  void update_fixed_pattern(char the_matrix_label, vector <double> &newPat, unsigned int thePat);

};
#endif
