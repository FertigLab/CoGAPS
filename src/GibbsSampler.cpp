#include <iostream>
#include <cmath>
#include <limits>
#include <stdexcept>
#include "GibbsSampler.h"

using namespace std;
using namespace gaps;
using std::vector;

// -----------------------------------------------------------------------------
unsigned long long atomicSize = std::numeric_limits<unsigned long long>::max();  
const double DOUBLE_POSINF = std::numeric_limits<double>::max();
const double DOUBLE_NEGINF = -std::numeric_limits<double>::max();
const double epsilon = 1e-10;
// -----------------------------------------------------------------------------



// ******************** CONSTRUCTOR ********************************************
GibbsSampler:: GibbsSampler(unsigned long nEquil, unsigned long nSample, unsigned int nFactor, 
			    double alphaA, double alphaP, double nMaxA, double nMaxP,
			    unsigned long nIterA, unsigned long nIterP, 
                            double max_gibbsmass_paraA, double max_gibbsmass_paraP,
                            unsigned long long atomicSize,
			    char label_A,char label_P,char label_D,char label_S,
			    const string & datafile, const string & variancefile,
                            const string & simulation_id)
  :_DMatrix(datafile.c_str(),label_D),
   _SMatrix(variancefile.c_str(),label_S){
 
  _nEquil = nEquil;
  _nSample = nSample;
  _nFactor = nFactor;
  _alphaA = alphaA;
  _alphaP = alphaP;
  _nMaxA = nMaxA;
  _nMaxP = nMaxP;
  _nIterA = nIterA;
  _nIterP = nIterP;
  _max_gibbsmass_paraA = max_gibbsmass_paraA;
  _max_gibbsmass_paraP = max_gibbsmass_paraP;
  _lambdaA_scale_factor = 1;
  _lambdaP_scale_factor = 1;
  _simulation_id = simulation_id;
  _atomicSize = atomicSize;
  _label_A = label_A;
  _label_P = label_P;
  _iter = 1;  // tmp use
  _annealingTemperature = 0.4; // tmp use
  _sysChi2 = 0.0; // tmp use
 
}

GibbsSampler:: GibbsSampler(unsigned long nEquil, unsigned long nSample, unsigned int nFactor, 
			    double alphaA, double alphaP, double nMaxA, double nMaxP,
			    unsigned long nIterA, unsigned long nIterP, 
                            double max_gibbsmass_paraA, double max_gibbsmass_paraP,
                            unsigned long long atomicSize,
			    char label_A,char label_P,char label_D,char label_S,
			    vector<vector<double> > &DVector, vector<vector<double> > &SVector,
                            const string & simulation_id)
  :_DMatrix(DVector, label_D),
   _SMatrix(SVector, label_S){
  _nEquil = nEquil;
  _nSample = nSample;
  _nFactor = nFactor;
  _alphaA = alphaA;
  _alphaP = alphaP;
  _nMaxA = nMaxA;
  _nMaxP = nMaxP;
  _nIterA = nIterA;
  _nIterP = nIterP;
  _max_gibbsmass_paraA = max_gibbsmass_paraA;
  _max_gibbsmass_paraP = max_gibbsmass_paraP;
  _lambdaA_scale_factor = 1;
  _lambdaP_scale_factor = 1;
  _simulation_id = simulation_id;
  _atomicSize = atomicSize;
  _label_A = label_A;
  _label_P = label_P;
  _iter = 1;  // tmp use
  _annealingTemperature = 0.4; // tmp use
  _sysChi2 = 0.0; // tmp use

}


// *************** METHODS FOR INITIALIZATION, DISPLAY, OUTPUT ***********************
void GibbsSampler::init_AMatrix_and_PMatrix(){
  // extract information from D as parameters
  _nRow = _DMatrix.get_nRow();
  _nCol = _DMatrix.get_nCol();

  // initialize matrices A and p
  _AMatrix.born_matrix(_nRow,_nFactor,_label_A,_alphaA);
  _PMatrix.born_matrix(_nFactor,_nCol,_label_P,_alphaP);
}

void GibbsSampler::init_AAtomicdomain_and_PAtomicdomain(){
  // extract information from D as parameters
  _nRow = _DMatrix.get_nRow();
  _nCol = _DMatrix.get_nCol();
  double D_mean = _DMatrix.cal_mean();

  // calcuate #Bins and lambda for the atomic spaces
  _nBinsA = _nRow*_nFactor;
  _lambdaA = _alphaA*sqrt(_nFactor/D_mean) * _lambdaA_scale_factor;
  _nBinsP = _nFactor*_nCol;
  _lambdaP = _alphaP*sqrt(_nFactor/D_mean) * _lambdaP_scale_factor;

  // calculate the maximum gibbs mass for A and p
  _max_gibbsmassA = _max_gibbsmass_paraA / _lambdaA;
  _max_gibbsmassP = _max_gibbsmass_paraP / _lambdaP;

  // initialize the atomic spaces
  _AAtomicdomain.initializeAtomic(_nBinsA,atomicSize,_alphaA,_lambdaA,_label_A);
  _PAtomicdomain.initializeAtomic(_nBinsP,atomicSize,_alphaP,_lambdaP,_label_P);

  //cout << "_lambdaA = " << _lambdaA << ", _max_gibbsmassA = " << _max_gibbsmassA << endl;
  //cout << "_lambdaP = " << _lambdaP << ", _max_gibbsmassP = " << _max_gibbsmassP << endl << endl;

}
// For fixing one domain
void GibbsSampler::init_AAtomicdomain_and_PAtomicdomain(char fixeddomain, const char input_file_name[]){
  // extract information from D as parameters
  _nRow = _DMatrix.get_nRow();
  _nCol = _DMatrix.get_nCol();
  double D_mean = _DMatrix.cal_mean();

  // calcuate #Bins and lambda for the atomic spaces
  _nBinsA = _nRow*_nFactor;
  _lambdaA = _alphaA*sqrt(_nFactor/D_mean) * _lambdaA_scale_factor;
  _nBinsP = _nFactor*_nCol;
  _lambdaP = _alphaP*sqrt(_nFactor/D_mean) * _lambdaP_scale_factor;

  // calculate the maximum gibbs mass for A and p
  _max_gibbsmassA = _max_gibbsmass_paraA / _lambdaA;
  _max_gibbsmassP = _max_gibbsmass_paraP / _lambdaP;

  // initialize the atomic spaces (fixed or not)
  if (fixeddomain == 'A'){
  _AAtomicdomain.FixedBins_initializeAtomic(_nBinsA,atomicSize,_alphaA,_lambdaA,_label_A, input_file_name);
  _PAtomicdomain.initializeAtomic(_nBinsP,atomicSize,_alphaP,_lambdaP,_label_P);
  }
  else{
  _AAtomicdomain.initializeAtomic(_nBinsA,atomicSize,_alphaA,_lambdaA,_label_A);
  _PAtomicdomain.FixedBins_initializeAtomic(_nBinsP,atomicSize,_alphaP,_lambdaP,_label_P, input_file_name);
  }

  //cout << "_lambdaA = " << _lambdaA << ", _max_gibbsmassA = " << _max_gibbsmassA << endl;
  //cout << "_lambdaP = " << _lambdaP << ", _max_gibbsmassP = " << _max_gibbsmassP << endl << endl;

}

// For fixing two domains
void GibbsSampler::init_AAtomicdomain_and_PAtomicdomain(const char input_file_nameA[],const char input_file_nameP[]){
  // extract information from D as parameters
  _nRow = _DMatrix.get_nRow();
  _nCol = _DMatrix.get_nCol();
  double D_mean = _DMatrix.cal_mean();

  // calcuate #Bins and lambda for the atomic spaces
  _nBinsA = _nRow*_nFactor;
  _lambdaA = _alphaA*sqrt(_nFactor/D_mean) * _lambdaA_scale_factor;
  _nBinsP = _nFactor*_nCol;
  _lambdaP = _alphaP*sqrt(_nFactor/D_mean) * _lambdaP_scale_factor;

  // calculate the maximum gibbs mass for A and p
  _max_gibbsmassA = _max_gibbsmass_paraA / _lambdaA;
  _max_gibbsmassP = _max_gibbsmass_paraP / _lambdaP;

  // initialize the atomic spaces (BOTH FIXED)
  _AAtomicdomain.FixedBins_initializeAtomic(_nBinsA,atomicSize,_alphaA,_lambdaA,_label_A, input_file_nameA);
  _PAtomicdomain.FixedBins_initializeAtomic(_nBinsP,atomicSize,_alphaP,_lambdaP,_label_P, input_file_nameP);
 

  //cout << "_lambdaA = " << _lambdaA << ", _max_gibbsmassA = " << _max_gibbsmassA << endl;
  //cout << "_lambdaP = " << _lambdaP << ", _max_gibbsmassP = " << _max_gibbsmassP << endl << endl;

}

// For fixing one domain in R
void GibbsSampler::init_AAtomicdomain_and_PAtomicdomain(char fixeddomain, vector<vector<double> > ReadBinProbs){
  // extract information from D as parameters
  _nRow = _DMatrix.get_nRow();
  _nCol = _DMatrix.get_nCol();
  double D_mean = _DMatrix.cal_mean();

  // calcuate #Bins and lambda for the atomic spaces
  _nBinsA = _nRow*_nFactor;
  _lambdaA = _alphaA*sqrt(_nFactor/D_mean) * _lambdaA_scale_factor;
  _nBinsP = _nFactor*_nCol;
  _lambdaP = _alphaP*sqrt(_nFactor/D_mean) * _lambdaP_scale_factor;

  // calculate the maximum gibbs mass for A and p
  _max_gibbsmassA = _max_gibbsmass_paraA / _lambdaA;
  _max_gibbsmassP = _max_gibbsmass_paraP / _lambdaP;

  // initialize the atomic spaces (fixed or not)
  if (fixeddomain == 'A'){
  _AAtomicdomain.FixedBins_initializeAtomic(_nBinsA,atomicSize,_alphaA,_lambdaA,_label_A, ReadBinProbs);
  _PAtomicdomain.initializeAtomic(_nBinsP,atomicSize,_alphaP,_lambdaP,_label_P);
  }
  else{
  _AAtomicdomain.initializeAtomic(_nBinsA,atomicSize,_alphaA,_lambdaA,_label_A);
  _PAtomicdomain.FixedBins_initializeAtomic(_nBinsP,atomicSize,_alphaP,_lambdaP,_label_P, ReadBinProbs);
  }

  //cout << "_lambdaA = " << _lambdaA << ", _max_gibbsmassA = " << _max_gibbsmassA << endl;
  //cout << "_lambdaP = " << _lambdaP << ", _max_gibbsmassP = " << _max_gibbsmassP << endl << endl;

}

// For fixing two domains in R
void GibbsSampler::init_AAtomicdomain_and_PAtomicdomain(vector<vector<double> > ReadBinProbsA,vector<vector<double> > ReadBinProbsP){
  // extract information from D as parameters
  _nRow = _DMatrix.get_nRow();
  _nCol = _DMatrix.get_nCol();
  double D_mean = _DMatrix.cal_mean();

  // calcuate #Bins and lambda for the atomic spaces
  _nBinsA = _nRow*_nFactor;
  _lambdaA = _alphaA*sqrt(_nFactor/D_mean) * _lambdaA_scale_factor;
  _nBinsP = _nFactor*_nCol;
  _lambdaP = _alphaP*sqrt(_nFactor/D_mean) * _lambdaP_scale_factor;

  // calculate the maximum gibbs mass for A and p
  _max_gibbsmassA = _max_gibbsmass_paraA / _lambdaA;
  _max_gibbsmassP = _max_gibbsmass_paraP / _lambdaP;

  // initialize the atomic spaces (BOTH FIXED)
  _AAtomicdomain.FixedBins_initializeAtomic(_nBinsA,atomicSize,_alphaA,_lambdaA,_label_A, ReadBinProbsA);
  _PAtomicdomain.FixedBins_initializeAtomic(_nBinsP,atomicSize,_alphaP,_lambdaP,_label_P, ReadBinProbsP);
 

  //cout << "_lambdaA = " << _lambdaA << ", _max_gibbsmassA = " << _max_gibbsmassA << endl;
  //cout << "_lambdaP = " << _lambdaP << ", _max_gibbsmassP = " << _max_gibbsmassP << endl << endl;

}

// clear all quantities related to the local matrix proposal
void GibbsSampler::clear_Proposal(){
  _Row_changed.clear();
  _Col_changed.clear();
  _mass_changed.clear();
  _atomicProposal.clear();
  _matrixElemChange.clear();

}

// clear all quantities related to the new local matrix proposal
void GibbsSampler::clear_new_Proposal(){
  _new_Row_changed.clear();
  _new_Col_changed.clear();
  _new_mass_changed.clear();
  _new_atomicProposal.clear();
  _new_matrixElemChange.clear();
}


void GibbsSampler::display_matrix(char matrix_label){
  switch(matrix_label){
  case 'D':
    {_DMatrix.display_matrix();break;}
  case 'S':
    {_SMatrix.display_matrix();break;}
  case 'A':
    {_AMatrix.display_matrix();break;}
  case 'P':
    {_PMatrix.display_matrix();break;}
  }
}

	// Just to see the matrices in the terminal
void GibbsSampler:: print_A_and_P(){
	  display_matrix('A');
	  display_matrix('P');
	  }

void GibbsSampler::display_atomicdomain(char atomic_label){
  switch(atomic_label){
  case 'A':
    {_AAtomicdomain.printAtomicInfo(); break;}
  case 'P':
    {_PAtomicdomain.printAtomicInfo(); break;}
  }
}

void GibbsSampler::local_display_matrix(vector<vector<double> > Mat, 
					unsigned int n_row, unsigned int n_col)
{
	/*To be uncommented for debugging 
  cout << endl;
  for(unsigned int m=0;m<n_row;++m)
    {
      for (unsigned int n=0; n<n_col;++n)
	{
	  cout << std::setw(10) << std::right;
	  cout << Mat[m][n] << " ";
	}
      cout << endl;
    }
  cout << endl;
  */
}


void GibbsSampler::local_display_matrix2(double ** Mat_ptr, 
					 unsigned int n_row, unsigned int n_col)
{
	/*To be uncommented for debugging
  cout << endl;
  for(unsigned int m=0;m<n_row;++m)
    {
      for (unsigned int n=0; n<n_col;++n)
	{
	  cout << std::setw(10) << std::right; 
	  cout << Mat_ptr[m][n] << " ";
	}
      cout << endl;
    }
  cout << endl;
  */
}

void GibbsSampler::local_display_matrix2F(ofstream& outputFile, double ** Mat_ptr, 
					  unsigned int n_row, unsigned int n_col){
  for(unsigned int m=0;m<n_row;++m)
    {
      for (unsigned int n=0; n<n_col;++n)
	{
	  outputFile << std::setw(10) << std::right; 
	  outputFile << Mat_ptr[m][n] << " ";
	}
      outputFile << endl;
    }

}

// Added for debugging the fixed bins code
void GibbsSampler::print_totMassinBins(char AtomicDomainLabel){
	/*Uncomment when debugging
 if (AtomicDomainLabel=='A'){
  for (int i=0; i < _nRow*_nFactor; i++){
   if (i%8==0 && i < 50){
   cout << _AAtomicdomain.getTotalMass(i) << endl;
   }
   }
   } else{
  for (int i=0; i < _nFactor * _nCol; i++){
  cout << _PAtomicdomain.getTotalMass(i) << endl;
   }
   }
   */
 }
 
// -----------------------------------------------------------------------------
void GibbsSampler::check_results(){
  //EJF double const * const * D = _DMatrix.get_matrix();
  //EJF double const * const * S = _SMatrix.get_matrix();
  double const * const * A = _AMatrix.get_matrix();
  double const * const * P = _PMatrix.get_matrix();

  vector<vector<double> > AP;
  AP.resize(_nRow,vector<double>(_nCol,0.0));

  for (unsigned int m=0; m < _nRow; ++m){
    for (unsigned int n=0; n < _nCol; ++n){
      for (unsigned int k=0; k < _nFactor; ++k){
	AP[m][n] += A[m][k]*P[k][n];
      }
    }
  }

  //cout << "The product matrix AP = A*P is: " << endl;
  local_display_matrix(AP,_nRow,_nCol);

}

void GibbsSampler::check_resultsF(ofstream& outputFile){
  //EJF double const * const * D = _DMatrix.get_matrix();
  //EJF double const * const * S = _SMatrix.get_matrix();
  double const * const * A = _AMatrix.get_matrix();
  double const * const * P = _PMatrix.get_matrix();

  vector<vector<double> > AP;
  AP.resize(_nRow,vector<double>(_nCol,0.0));

  for (unsigned int m=0; m < _nRow; ++m){
    for (unsigned int n=0; n < _nCol; ++n){
      for (unsigned int k=0; k < _nFactor; ++k){
	AP[m][n] += A[m][k]*P[k][n];
      }
    }
  }

  outputFile << "The product matrix AP = A*P is: " << endl;
  outputFile << endl;
  for(unsigned int m=0; m < _nRow;++m)
    {
      for (unsigned int n=0; n< _nCol;++n)
	{
	  outputFile << setiosflags(ios::right) << setw(10) << AP[m][n] << " ";
	}
      outputFile << endl;
    }
  outputFile << endl;

}

// -----------------------------------------------------------------------------
void GibbsSampler::output_atomicdomain(char atomic_label,unsigned long Samp_cycle){

  char outputFilename[80]; 
  switch(atomic_label){
  case 'A':
    { 
      strcpy(outputFilename,_simulation_id.c_str());
      strcat(outputFilename,"_A_atomicdomain.txt");
      _AAtomicdomain.writeAtomicInfo(outputFilename,Samp_cycle);
      break;
    }
  case 'P':
    {     
      strcpy(outputFilename,_simulation_id.c_str());
      strcat(outputFilename,"_P_atomicdomain.txt");
      _PAtomicdomain.writeAtomicInfo(outputFilename,Samp_cycle);
      break;
    }
  }
}

// -----------------------------------------------------------------------------
void GibbsSampler::output_computing_info(char outputFilename[],
                                         unsigned long Equil_cycle, unsigned long nEquil,
					 unsigned long Samp_cycle, unsigned long nSample){

  ofstream outputFile;
  outputFile.open(outputFilename,ios::out|ios::app);

  outputFile << " *************************************************** " << endl;
  outputFile << " --------------- NEW ROUND ------------------------- " << endl;
  outputFile << " *************************************************** " << endl << endl;
  outputFile << "Equilibration cycle index = " << Equil_cycle << endl;
  outputFile << "Total number of equilibrating cycles to perform = " <<  nEquil << endl;
  outputFile << "Sampling cycle index = " << Samp_cycle << endl;
  outputFile << "Total number of sampling cycles to perform = " <<  nSample << endl;
  outputFile << "System Chi2-value = " << _sysChi2 << endl;

  _AMatrix.display_matrixF(outputFile);
  _PMatrix.display_matrixF(outputFile);
  check_resultsF(outputFile);

  outputFile.close();

}


// ********* METHODS TO GO BETWEEN ATOMIC SPACE AND MATRIX  ********************
unsigned int GibbsSampler::getRow( char matrix_label ,unsigned int iBin){
  switch(matrix_label){ 
  case 'A':  // A - horizontal addressing implicit
    { return (floor(iBin / _nFactor)); break;}
  case 'P':  // P - vertical addressing implicit
    { return iBin % _nFactor; break;}
  }
  
  // EJF dummy return to avoid warnings
  return iBin;
}

unsigned int GibbsSampler::getCol(char matrix_label ,unsigned int iBin){
  switch(matrix_label){
  case 'A':  // A - horizontal addressing implicit
    { return iBin % _nFactor; break;}
  case 'P':  // P - vertical addressing implicit
    { return floor(iBin / _nFactor); break;}
  }
  
  // EJF dummy return to avoid warnings
  return iBin;
}


unsigned int GibbsSampler::getTotNumAtoms(char matrix_label){

  switch(matrix_label){
  case 'A':
    { return _AAtomicdomain.getNAtom();
      break;}
  case 'P':
    { return _PAtomicdomain.getNAtom();
      break;}
  }

  // EJF dummy return to avoid warnings
  return 0; 
}
 
 vector <vector <vector <double> > > GibbsSampler::getNormedMatrices(){
  double ** A = _AMatrix.get_matrix();
  double ** P = _PMatrix.get_matrix();
  vector <vector <double> > AMatrixNormed;
  AMatrixNormed.resize(_nRow,vector<double>(_nFactor,0.0));
  vector <vector <double> > PMatrixNormed;
  PMatrixNormed.resize(_nFactor,vector<double>(_nCol,0.0));
  vector<double> k(_nFactor);  // normalized vector

  // compute the normalization vector
  for (int m=0; m< _nFactor; ++m){
    k[m] = 0.;
    for (int n=0; n< _nCol; ++n){
      k[m] += P[m][n];
    }
    if (k[m] == 0){  // when the whole row of P is zero, then don't do anything
      k[m] = 1.0;
    }
  }
        
    for (int m=0; m < _nRow; ++m){
      for (int n=0; n < _nFactor; ++n){
		AMatrixNormed[m][n] = A[m][n] * k[n];
      }
    }

    for (int m=0; m < _nFactor; ++m){
      for (int n=0; n < _nCol; ++n){
		PMatrixNormed[m][n] = P[m][n] / k[m];
      }
    } 
  vector <vector <vector <double> > > NormedMatrices;
  NormedMatrices.push_back(AMatrixNormed);
  NormedMatrices.push_back(PMatrixNormed);
  return NormedMatrices;
 }

// ======CoGAPS Test methods============

// Get the full atomic domain from the atomic space
 map <unsigned long long, double> GibbsSampler::getAtomicDomain(char matrix_label){
 map <unsigned long long, double> zero;
  if (matrix_label == 'A'){
   return _AAtomicdomain.getDomain();
  }
  else if(matrix_label == 'P'){
   return _PAtomicdomain.getDomain();
  }
  else{
   //cout << "Invalid input matrix" << endl;
   }
   return zero;
  }
  
 // Manually calculate the matrix A from the atomic space passed in.  
  vector <vector <double> > GibbsSampler::createSampleAMat(map <unsigned long long, double> ADomain){
   // Make that a matrix
   vector <vector <double> > SampleAMatrix;
   SampleAMatrix.resize(_nRow);
   for (int i=0; i < _nRow; i++){
    SampleAMatrix[i].resize(_nFactor, 0);
	}
   //Copy the parsed domain 
    map<unsigned long long, double>::const_iterator iter;
	for (iter = ADomain.begin(); iter != ADomain.end(); iter++){
	 unsigned int theBin = _AAtomicdomain.getBin(iter->first);
	 // Put the mass in the bin in the matrix
	 unsigned int theRow = getRow('A', theBin);
	 unsigned int theCol = getCol('A', theBin);
	 SampleAMatrix[theRow][theCol]+=iter->second;
	}
    return SampleAMatrix;	
   }
 
 // Manually calculate the matrix P from the atomic space passed in. 	 
  vector <vector <double> > GibbsSampler::createSamplePMat(map <unsigned long long, double> PDomain){
   // Make that a matrix
   vector <vector <double> > SamplePMatrix;
   SamplePMatrix.resize(_nFactor);
   for (int i=0; i < _nFactor; i++){
    SamplePMatrix[i].resize(_nCol, 0);
	}
   //Copy the parsed domain 
    map<unsigned long long, double>::const_iterator iter;
	for (iter=PDomain.begin(); iter!=PDomain.end(); iter++){
	 unsigned int theBin = _PAtomicdomain.getBin(iter->first);
	 // Put the mass in the bin in the matrix
	 unsigned int theRow = getRow('P', theBin);
	 unsigned int theCol = getCol('P', theBin);
	 SamplePMatrix[theRow][theCol]+=iter->second;
	}
    return SamplePMatrix;	
   }
  
 // Manually calculate the chi squared value based on the 2 matrices passed in  
  double GibbsSampler::ManualCalcChiSqu(vector <vector <double> > SampleAMat, vector <vector <double> > SamplePMat){
   Matrix SampleAMatrix(SampleAMat, 'S');
   //cout << "THe sample A matrix: " << endl;
   //SampleAMatrix.display_matrix();
   Matrix SamplePMatrix(SamplePMat, 'S');
   //cout << "THe sample P matrix: " << endl;
   //SamplePMatrix.display_matrix();
   double ** D = _DMatrix.get_matrix();
   double ** S = _SMatrix.get_matrix();
   double ** A = SampleAMatrix.get_matrix();
   double ** P = SamplePMatrix.get_matrix();
   return GAPSNorm::calChi2(D,S,A,P,_nRow,_nCol,_nFactor);
  }


// ************* METHODS FOR COMPUTING LIKELIHOOD FUNCTIONS ******************
// ---------------------------------------------------------------------------
double GibbsSampler::cal_logLikelihood(){
  double ** D = _DMatrix.get_matrix();
  double ** S = _SMatrix.get_matrix();
  double ** A = _AMatrix.get_matrix();
  double ** P = _PMatrix.get_matrix();
  return GAPSNorm::calChi2(D,S,A,P,_nRow,_nCol,_nFactor)/2.;
}

void GibbsSampler::cal_delloglikelihood_example(){

}

// -----------------------------------------------------------------------------
double GibbsSampler::computeDeltaLL(char the_matrix_label,
				    double const * const * D,
				    double const * const * S,
				    double const * const * A,
				    double const * const * P,
				    unsigned int the_nChange_matrixElemChange,
				    const vector<boost::tuple<unsigned int, unsigned int, double> > the_matrixElemChange){

  double DelLL = 0.0;

 
  switch(the_matrix_label){
  case 'A':
    { 
      if (the_nChange_matrixElemChange == 0){
	DelLL = 0.0;
      } else if (the_nChange_matrixElemChange == 1){
	DelLL = GAPSNorm::calcDeltaLL1E('A',D,S,A,P,the_matrixElemChange,_nRow,
                                        _nCol,_nFactor);
      } else if (the_nChange_matrixElemChange == 2){
	DelLL = GAPSNorm::calcDeltaLL2E('A',D,S,A,P,the_matrixElemChange,_nRow,
                                        _nCol,_nFactor);
      } else {
	DelLL = GAPSNorm::calcDeltaLLGen('A',D,S,A,P,the_matrixElemChange,_nRow,
					 _nCol,_nFactor);
      } // end of if-block according to proposal.size()

      break;} // end of switch-block 'A'
  case 'P':
    { 
      if (the_nChange_matrixElemChange == 0){
	DelLL = 0.0;
      } else if (the_nChange_matrixElemChange == 1){
	DelLL = GAPSNorm::calcDeltaLL1E('P',D,S,A,P,the_matrixElemChange,_nRow,
                                        _nCol,_nFactor);
      } else if (the_nChange_matrixElemChange == 2){
	DelLL = GAPSNorm::calcDeltaLL2E('P',D,S,A,P,the_matrixElemChange,_nRow,
                                        _nCol,_nFactor);
      } else {
	DelLL = GAPSNorm::calcDeltaLLGen('P',D,S,A,P,the_matrixElemChange,_nRow,
					 _nCol,_nFactor);
      } // end of if-block according to proposal.size()

      break;} // end of switch-block 'P' 

  } // end of switch block

  return DelLL;

} // end of computeDeltaLL 


// ---------------- For checking against computeDeltaLL
double GibbsSampler::computeDeltaLL2(char the_matrix_label,
				     double const * const * D,
				     double const * const * S,
				     double const * const * A,
				     double const * const * P,
				     unsigned int the_nChange_matrixElemChange,
				     const vector<boost::tuple<unsigned int, unsigned int, double> > the_matrixElemChange){

  double DelLL = 0.0;
  switch(the_matrix_label){
  case 'A':
    {
      DelLL = GAPSNorm::calcDeltaLLGen('A',D,S,A,P,the_matrixElemChange,_nRow,
				       _nCol,_nFactor);
      break;}
  case 'P':
    {
      DelLL = GAPSNorm::calcDeltaLLGen('P',D,S,A,P,the_matrixElemChange,_nRow,
				       _nCol,_nFactor);
      break;}
  }
  return DelLL;

} // end of computeDeltaLL2 
// -------------------------------------


// *************** METHODS FOR MAKING PROPOSAL ********************************
// -----------------------------------------------------------------------------
void GibbsSampler::update_example(char atomic_domain_label){

} // end of update_example


// -----------------------------------------------------------------------------
vector<vector<double> > GibbsSampler::atomicProposal2Matrix(char atomic_domain_label,
							    double const * const * origMatrix)
{  
  unsigned int bin;
  unsigned int chRow, chCol;
  
  switch(atomic_domain_label){
  case 'A':
    { 
      vector<vector<double> > newMatrix(_nRow,vector<double>(_nFactor,0));
      map<unsigned long long, double> proposal = _AAtomicdomain.getProposedAtoms();
      for (map<unsigned long long, double>::const_iterator
	     iter=proposal.begin(); iter != proposal.end(); ++ iter) {
	bin = _AAtomicdomain.getBin(iter->first);
	chRow = getRow('A',bin);
	chCol = getCol('A',bin);
	newMatrix[chRow][chCol] += iter->second;
      }
      return newMatrix;
      break;} // end of case 'A'

  case 'P':
    {
      vector<vector<double> > newMatrix(_nFactor,vector<double>(_nCol,0));
      map<unsigned long long, double> proposal = _PAtomicdomain.getProposedAtoms();
      for (map<unsigned long long, double>::const_iterator
	     iter=proposal.begin(); iter != proposal.end(); ++ iter) {
	bin = _PAtomicdomain.getBin(iter->first);
	chRow = getRow('P',bin);
	chCol = getCol('P',bin);
	newMatrix[chRow][chCol] += iter->second;
      }
      return newMatrix;
      break;} // end of case 'P'

  } // end of switch 

  // EJF dummy return to avoid warning
  vector<vector<double> > newMatrix(0,vector<double>(0,0));
  return newMatrix;
  
} // end of atomicProposal2Matrix

// -----------------------------------------------------------------------------
vector<vector<double> > GibbsSampler::atomicProposal2FullMatrix(char atomic_domain_label,
								double const * const * origMatrix)
{  
  unsigned int bin;
  unsigned int chRow, chCol;
  
  switch(atomic_domain_label){
  case 'A':
    { 
      vector<vector<double> > FullnewMatrix(_nRow,vector<double>(_nFactor,0));
      
      for (unsigned int iRow=0; iRow < _nRow; ++iRow){
	for (unsigned int iCol=0; iCol < _nFactor; ++iCol){
	  if (origMatrix[iRow][iCol] < epsilon){
	    FullnewMatrix[iRow][iCol]=0.;
	  } else {
	    FullnewMatrix[iRow][iCol]=origMatrix[iRow][iCol];
	  }
        }
      } // end of for block
      
      map<unsigned long long, double> proposal = _AAtomicdomain.getProposedAtoms();
      for (map<unsigned long long, double>::const_iterator
	     iter=proposal.begin(); iter != proposal.end(); ++ iter) {
	bin = _AAtomicdomain.getBin(iter->first);
	chRow = getRow('A',bin);
	chCol = getCol('A',bin);
	FullnewMatrix[chRow][chCol] += iter->second;
      }
      return FullnewMatrix;
      break;} // end of case 'A'

  case 'P':
    {
      vector<vector<double> > FullnewMatrix(_nFactor,vector<double>(_nCol,0));
      
      for (unsigned int iRow=0; iRow < _nFactor; ++iRow){
	for (unsigned int iCol=0; iCol < _nCol; ++iCol){
	  if (origMatrix[iRow][iCol] < epsilon){
	    FullnewMatrix[iRow][iCol]=0.;
	  } else {
	    FullnewMatrix[iRow][iCol]=origMatrix[iRow][iCol];
	  }
        }
      } // end of for block
      
      map<unsigned long long, double> proposal = _PAtomicdomain.getProposedAtoms();
      for (map<unsigned long long, double>::const_iterator
	     iter=proposal.begin(); iter != proposal.end(); ++ iter) {
	bin = _PAtomicdomain.getBin(iter->first);
	chRow = getRow('P',bin);
	chCol = getCol('P',bin);
	FullnewMatrix[chRow][chCol] += iter->second;
      }
      return FullnewMatrix;
      break;} // end of case 'P'

  } // end of switch 

  // EJF dummy return to avoid warning
  vector<vector<double> > newMatrix(0,vector<double>(0,0));
  return newMatrix;  
  
} // end of atomicProposal2FullMatrix

// ----------------------------------------------------------------------------
void GibbsSampler::extract_atomicProposal(char the_matrix_label){
  unsigned int bin, chRow, chCol;
  double chmass;
  map<unsigned long long, double>::const_iterator iter;

  _nChange_matrixElemChange = 0;
  _nChange_atomicProposal = _atomicProposal.size(); 

  if (_nChange_atomicProposal == 0) {   // atomic proposal size = 0
    _nChange_matrixElemChange = 0;
  } 

  else if (_nChange_atomicProposal == 1){  // atomic proposal size = 1 
    iter = _atomicProposal.begin();
    switch(the_matrix_label){
    case 'A':
      {
	bin = _AAtomicdomain.getBin(iter->first);
	chRow = getRow('A',bin);	  
	chCol = getCol('A',bin);
	break;}
    case 'P':
      {
	bin = _PAtomicdomain.getBin(iter->first);
	chRow = getRow('P',bin);
	chCol = getCol('P',bin);
	break;}
    } // end of switch-block for atomic proposal size = 1
    chmass = iter->second;
    _Row_changed.push_back(chRow);
    _Col_changed.push_back(chCol);
    _mass_changed.push_back(chmass);
    _nChange_matrixElemChange = 1;
    _matrixElemChange.push_back(boost::make_tuple(chRow,chCol,chmass));

  } // end of if-block for atomic proposal size = 1

  else {     // atomic proposal size = 2 or above
    unsigned int count = 0;

    for (iter = _atomicProposal.begin(); iter != _atomicProposal.end(); ++ iter) {

      switch (the_matrix_label){
      case 'A':
	{
	  bin = _AAtomicdomain.getBin(iter->first);
	  chRow = getRow('A',bin);
	  chCol = getCol('A',bin);
	  break;}
      case 'P':
	{
	  bin = _PAtomicdomain.getBin(iter->first);
	  chRow = getRow('P',bin);
	  chCol = getCol('P',bin);
	  break;}
      } // end of switch-block
      chmass = iter->second;

      if (count == 0){    // nothing to check for the first count
	_Row_changed.push_back(chRow);
	_Col_changed.push_back(chCol);
	_mass_changed.push_back(chmass);
	count += 1;
	_nChange_matrixElemChange += 1;
      } else {
	for (unsigned int m = 0; m < count; ++m){
	  if (chRow == _Row_changed[m] && chCol == _Col_changed[m]){
	    _mass_changed[m] += chmass;
	    if (_mass_changed[m] == 0) {	 
	      _nChange_matrixElemChange -= 1;
	      _Row_changed.erase(_Row_changed.begin()+m);
	      _Col_changed.erase(_Col_changed.begin()+m);
	      _mass_changed.erase(_mass_changed.begin()+m);
	    }
	  } else {
	    _Row_changed.push_back(chRow);
	    _Col_changed.push_back(chCol);
	    _mass_changed.push_back(chmass);
	    _nChange_matrixElemChange += 1;
	  } // end of if-block when chRow and chRol refers to new matrix elements
	} // end of for-block when looping through existing elements in _mass_changed
	count = _nChange_matrixElemChange;
      } // end of if-block for count != 0

    } // end of for-block with iter looping through the atomic proposal

    // make up _matrixElemChange
    for (unsigned int m = 0; m<_nChange_matrixElemChange; ++m){
      _matrixElemChange.push_back(boost::make_tuple(_Row_changed[m],_Col_changed[m],_mass_changed[m]));
    }
 
  } // end of if-block for proposal size = 2 or above  
} // end of extract_atomicProposal


// -----------------------------------------------------------------------------
void GibbsSampler::extract_new_atomicProposal(char the_matrix_label){
  unsigned int bin, chRow, chCol;
  double chmass;
  map<unsigned long long, double>::const_iterator iter;

  _new_nChange_matrixElemChange = 0;
  _new_nChange_atomicProposal = _new_atomicProposal.size(); 

  if (_new_nChange_atomicProposal == 0) {    // atomic proposal size = 0
    _new_nChange_matrixElemChange = 0;
  } 

  else if (_new_nChange_atomicProposal == 1){  // atomic proposal size = 1 
    iter = _new_atomicProposal.begin();
    switch(the_matrix_label){
    case 'A':
      {
	bin = _AAtomicdomain.getBin(iter->first);
	chRow = getRow('A',bin);	  
	chCol = getCol('A',bin);
	break;}
    case 'P':
      {
	bin = _PAtomicdomain.getBin(iter->first);
	chRow = getRow('P',bin);
	chCol = getCol('P',bin);
	break;}
    } // end of switch-block for atomic proposal size = 1
    chmass = iter->second;
    _new_Row_changed.push_back(chRow);
    _new_Col_changed.push_back(chCol);
    _new_mass_changed.push_back(chmass);
    _new_nChange_matrixElemChange = 1;
    _new_matrixElemChange.push_back(boost::make_tuple(chRow,chCol,chmass));

  } // end of if-block for atomic proposal size = 1

  else {     // atomic proposal size = 2 or above
    unsigned int count = 0;

    for (iter = _new_atomicProposal.begin(); iter != _new_atomicProposal.end(); ++ iter) {

      switch (the_matrix_label){
      case 'A':
	{
	  bin = _AAtomicdomain.getBin(iter->first);
	  chRow = getRow('A',bin);
	  chCol = getCol('A',bin);
	  break;}
      case 'P':
	{
	  bin = _PAtomicdomain.getBin(iter->first);
	  chRow = getRow('P',bin);
	  chCol = getCol('P',bin);
	  break;}
      } // end of switch-block
      chmass = iter->second;

      if (count == 0){    // nothing to check for the first count
	_new_Row_changed.push_back(chRow);
	_new_Col_changed.push_back(chCol);
	_new_mass_changed.push_back(chmass);
	count += 1;
	_new_nChange_matrixElemChange += 1;
      } else {
	for (unsigned int m = 0; m < count; ++m){
	  if (chRow == _new_Row_changed[m] && chCol == _new_Col_changed[m]){
	    _new_mass_changed[m] += chmass;
	    if (_new_mass_changed[m] == 0) {	 
	      _new_nChange_matrixElemChange -= 1;
	      _new_Row_changed.erase(_new_Row_changed.begin()+m);
	      _new_Col_changed.erase(_new_Col_changed.begin()+m);
	      _new_mass_changed.erase(_new_mass_changed.begin()+m);
	    }
	  } else {
	    _new_Row_changed.push_back(chRow);
	    _new_Col_changed.push_back(chCol);
	    _new_mass_changed.push_back(chmass);
	    _new_nChange_matrixElemChange += 1;
	  } // end of if-block when chRow and chRol refers to new matrix elements
	} // end of for-block when looping through existing elements in _mass_changed
	count = _new_nChange_matrixElemChange;
      } // end of if-block for count != 0

    } // end of for-block with iter looping through the atomic proposal

    // make up _new_matrixElemChange
    for (unsigned int m = 0; m<_new_nChange_matrixElemChange; ++m){
      _new_matrixElemChange.push_back(boost::make_tuple(_new_Row_changed[m],
						      _new_Col_changed[m],_new_mass_changed[m]));
    }
 
  } // end of if-block for proposal size = 2 or above

} // end of extract_new_atomicProposal



// -----------------------------------------------------------------------------
void GibbsSampler::update(char the_matrix_label){
  double rng = 0.1; // no use, filling up the list?
  double ** D = _DMatrix.get_matrix();
  double ** S = _SMatrix.get_matrix();
  double ** AOrig = _AMatrix.get_matrix();
  double ** POrig = _PMatrix.get_matrix();
  vector<vector<double> > del_matrix;
  map<unsigned long long, double>::iterator iter;

  bool Q_update;
 
  switch(the_matrix_label){
  case 'A':
    {
      // ----------- making a proposal from atomic space A:
      _AAtomicdomain.makeProposal(rng);
      get_oper_type('A');
      _atomicProposal = _AAtomicdomain.getProposedAtoms();
      extract_atomicProposal('A');
	  
	  /*Uncomment when debugging
      if (_nChange_atomicProposal == 1 && (_oper_type =='E' || _oper_type =='M'))
	{cout << "update inconsistency A1! _nChange_atomicProposal = " << _nChange_atomicProposal <<
	         ", _nChange_matrixElemChange = " << _nChange_matrixElemChange << 
	         ", _oper_type = " << _oper_type << endl;}         
      if (_nChange_atomicProposal == 2 && (_oper_type =='D' || _oper_type =='B'))
	{cout << "update inconsistency A2! _nChange_atomicProposal = " << _nChange_atomicProposal << 
 	         ", _nChange_matrixElemChange = " << _nChange_matrixElemChange << 
	         ", _oper_type = " << _oper_type << endl;}   
	*/ 
	
      // ----------------------------------
      // the proposal is translated into a proposal to matrix A:
	  
      // ----------- modify the proposal in a Gibbs way:      
      //EJF unsigned int iRow, iCol, iFactor;
      if ( _nChange_atomicProposal == 0){}
      if ( _nChange_atomicProposal> 2){
	   throw logic_error("GibbsSampler: can't change more than two atoms!!");
      }
	  
	  // Update based on the _oper_type
	  if (_oper_type == 'D') {
	    Q_update = death('A',D,S,AOrig,POrig);
      } else if (_oper_type == 'B') {
	    Q_update = birth('A',D,S,AOrig,POrig);
	  } else if (_oper_type == 'M') {
	    Q_update = move('A',D,S,AOrig,POrig);
      } else {
	    Q_update = exchange('A',D,S,AOrig,POrig);
	  }
	  
	  // Update the matrix with improved update, if there are further updates
	  if (Q_update == true){
	   _AMatrix.matrix_Elem_update(_new_matrixElemChange,_oper_type,_new_nChange_matrixElemChange);
	}
	
      break;
    } // end of case 'A'
  case 'P':
    {
      // ----------- making a proposal from atomic space P:
      _PAtomicdomain.makeProposal(rng);
      get_oper_type('P');
      _atomicProposal = _PAtomicdomain.getProposedAtoms();
      extract_atomicProposal('P');    

	/*Uncomment when debugging
      if (_nChange_atomicProposal == 1 && (_oper_type =='E' || _oper_type =='M'))
	{cout << "update inconsistency P1! _nChange_atomicProposal = " << _nChange_atomicProposal <<
	         ", _nChange_matrixElemChange = " << _nChange_matrixElemChange << 
	         ", _oper_type = " << _oper_type << endl;}         
      if (_nChange_atomicProposal == 2 && (_oper_type =='D' || _oper_type =='B'))
	{cout << "update inconsistency P2! _nChange_atomicProposal = " << _nChange_atomicProposal << 
 	         ", _nChange_matrixElemChange = " << _nChange_matrixElemChange << 
	         ", _oper_type = " << _oper_type << endl;}   
	*/
	
      // ----------------------------------
      // the proposal is translated into a proposal to matrix P:

      // ----------- modify the proposal in a Gibbs way:     
      //EJF unsigned int iRow, iCol, iFactor;
      if (_nChange_atomicProposal== 0){}
      if (_nChange_atomicProposal > 2){
	   throw logic_error("GibbsSampler: can't chnage more than two atoms!!");
      }
	  
	  // Update based on the _oper_type
	  if (_oper_type == 'D') {
	    Q_update = death('P',D,S,AOrig,POrig);
      } else if (_oper_type == 'B') {
	    Q_update = birth('P',D,S,AOrig,POrig);
	  } else if (_oper_type == 'M') {
	    Q_update = move('P',D,S,AOrig,POrig);
      } else {
	    Q_update = exchange('P',D,S,AOrig,POrig);
	  }
	  
	  if (Q_update == true){
	   _PMatrix.matrix_Elem_update(_new_matrixElemChange,_oper_type,_new_nChange_matrixElemChange);
	}

      break;
    } // end of case 'P'
  } // end of switch block

  // clear Proposal for the next run
  clear_Proposal();
  clear_new_Proposal();

} // end of update()

// ----------------------------------------------------------------------------
void GibbsSampler::init_sysChi2(){
  _sysChi2 = 2.*cal_logLikelihood();

}

void GibbsSampler::update_sysChi2(double delsysChi2){
  _sysChi2 -= 2.*delsysChi2;
}

double GibbsSampler::get_sysChi2(){
  return _sysChi2;
}


// -----------------------------------------------------------------------------
void GibbsSampler::get_oper_type(char the_matrix_label){
  switch(the_matrix_label){
  case 'A':
    {
      _oper_type = _AAtomicdomain.get_oper_type();
      break;
    }
  case 'P':
    {
      _oper_type = _PAtomicdomain.get_oper_type();
      break;
    }
  }

}

bool GibbsSampler::death(char the_matrix_label,
						  double const * const * D,
						  double const * const * S,
						  double ** AOrig,
						  double ** POrig)// in progress
{

  double rng = 0.1; // no use, just fill up the list
  double newMass = 0;
  double attemptMass = 0;

  // read in the original _atomicProposal made from the prior
  unsigned long long location = _atomicProposal.begin()->first;
  double origMass = _atomicProposal.begin()->second;
  //EJF unsigned int bin;
  unsigned int iRow = _Row_changed[0];
  unsigned int iCol = _Col_changed[0];
  double delLL = 0;
  double delLLnew = 0;

    // ------------------------- DEATH -------------------------------------------
    // put in the changes to the atomic space, and compute the corresponding change in 
    // the loglikelihood. 
    switch(the_matrix_label){
    case 'A':
      { 
        delLL = computeDeltaLL('A',D,S,AOrig,POrig,_nChange_matrixElemChange,_matrixElemChange); 
        _AAtomicdomain.acceptProposal(false); // "false" only means not to update _iter!
	    update_sysChi2(delLL); 	    // update system Chi2 
	    _AMatrix.matrix_Elem_update(_matrixElemChange,_oper_type,_nChange_matrixElemChange);
	    break;
		}
    case 'P':
      { 
	    delLL = computeDeltaLL('P',D,S,AOrig,POrig,_nChange_matrixElemChange,_matrixElemChange);  
        _PAtomicdomain.acceptProposal(false); // "false" only means not to update _iter!
	    update_sysChi2(delLL);     // update system Chi2
	    _PMatrix.matrix_Elem_update(_matrixElemChange,_oper_type,_nChange_matrixElemChange);
	    break;
		}
    } // end of switch-block
     
    // an attempt to rebirth
    attemptMass = -origMass;
    switch(the_matrix_label){
    case 'A':
      {
	    // Check other matrix to see if we can use Gibbs
        if (!checkOtherMatrix('A',iRow, iCol, POrig)) { 
	      newMass = attemptMass;
	    } else {
	      newMass = getMass('A',attemptMass,iRow,iCol,POrig,AOrig,D,S,rng); //Gibbs birth
	      // ------- Q: think about it
	      if (newMass <= epsilon) {
	       newMass = attemptMass;
	      }
	    } // end of if-block 

	    _new_atomicProposal.insert(pair<unsigned long long,double>(location,newMass));
	    extract_new_atomicProposal('A');
	    _AAtomicdomain.setProposedAtomMass(_new_atomicProposal,true);  
        delLLnew = computeDeltaLL('A',D,S,AOrig,POrig,_new_nChange_matrixElemChange,_new_matrixElemChange);
	    break;
      } // end of switch-block for A
    case 'P':
      {
	    // Check other matrix to see if we can use Gibbs
	    if (!checkOtherMatrix('P',iRow, iCol, AOrig)) { 
	      newMass = attemptMass;
	     } else {
          newMass = getMass('P',attemptMass,iRow,iCol,AOrig,POrig,D,S,rng); //Gibbs birth
	      // ----- Q: think about it
	      if (newMass <= epsilon) {
	       newMass = attemptMass;
	      }
	    } // end of if-block 

	    _new_atomicProposal.insert(pair<unsigned long long,double>(location,newMass));
	    extract_new_atomicProposal('P');
	    _PAtomicdomain.setProposedAtomMass(_new_atomicProposal,true);
        delLLnew = computeDeltaLL('P',D,S,AOrig,POrig,_new_nChange_matrixElemChange,_new_matrixElemChange);
	    break;
      } // end of switch-block for P

    } // end of switch-block

    // M-H sampling to determine whether or not we can accept Gibbs
    if (delLLnew*_annealingTemperature  < log(randgen('U',0,0))) {

     switch(the_matrix_label){
      case 'A':
	  { 
	    _AAtomicdomain.rejectProposal(false);
 	    return false;
	    break;
      }
      case 'P':
	  { 
	    _PAtomicdomain.rejectProposal(false);
	    return false;
	    break;
	  }
     } // end of switch-block
	  
    } else {
	
      switch(the_matrix_label){
      case 'A':
	  { 
          _AAtomicdomain.acceptProposal(false); 
	      update_sysChi2(delLLnew);  // update system Chi2
	      return true;
          break;
	   }
      case 'P':
	  {   _PAtomicdomain.acceptProposal(false);  
	      update_sysChi2(delLLnew);  // update system Chi2
	      return true;
          break;
	   }
      } // end of switch-block
    } // else of if-block for M-H sampling

    return false;    

}  // end of death method

bool GibbsSampler::birth(char the_matrix_label,
						  double const * const * D,
						  double const * const * S,
						  double ** AOrig,
						  double ** POrig)
{

  double rng = 0.1; 
  double newMass = 0;
  //EJF double attemptMass = 0;

  // read in the original _atomicProposal made from the prior
  unsigned long long location = _atomicProposal.begin()->first;
  double origMass = _atomicProposal.begin()->second;
  //EJF unsigned int bin;
  unsigned int iRow = _Row_changed[0];
  unsigned int iCol = _Col_changed[0];
  double delLL = 0;
  double delLLnew = 0;
  
  // -------------- BIRTH ------------------------------------------------------
  
    switch(the_matrix_label){
    case 'A':
      {
	   // checking conditions for update
	   if (iRow >= _nRow || iCol >= _nFactor) {
	    throw logic_error("Cannot update pattern out of range in A.");
	   }
	// Check other matrix to see if we can use Gibbs
	if ( !checkOtherMatrix('A', iRow, iCol, POrig)) {
	
	  _AAtomicdomain.acceptProposal(false); //accept original proposal 
      delLL = computeDeltaLL('A',D,S,AOrig,POrig,_nChange_matrixElemChange,_matrixElemChange);
	  update_sysChi2(delLL);  // update system Chi2
	  _new_atomicProposal.insert(pair<unsigned long long,double>(location,origMass)); //update _new_atomicProposal with original change
      extract_new_atomicProposal('A');
	  
	  return true;
	  }
     
	 //Otherwise, do a Gibbs birth
     newMass = getMass('A',origMass,iRow,iCol,POrig,AOrig,D,S,rng);
	 _new_atomicProposal.insert(pair<unsigned long long,double>(location,newMass));
  	 extract_new_atomicProposal('A');
	 _AAtomicdomain.setProposedAtomMass(_new_atomicProposal,false);
     delLLnew = computeDeltaLL('A',D,S,AOrig,POrig,_new_nChange_matrixElemChange,_new_matrixElemChange);

	break;
     } // end of case-A block
	 
    case 'P':
      {
	   // checking conditions for update
	   if (iRow >= _nFactor || iCol >= _nCol) {
	    throw logic_error("Cannot update pattern out of range in P.");
	   }
	   // Check other matrix to see if we can use Gibbs
	   if ( !checkOtherMatrix('P', iRow, iCol, AOrig)) {
	    
		_PAtomicdomain.acceptProposal(false); // accept original proposal
	    delLL = computeDeltaLL('P',D,S,AOrig,POrig,_nChange_matrixElemChange,_matrixElemChange);  
  	    update_sysChi2(delLL);  // update system Chi2
	    _new_atomicProposal.insert(pair<unsigned long long,double>(location,origMass));
    	extract_new_atomicProposal('P');
		
	    return true;
	   }

	   	 //Otherwise, do a Gibbs birth
        newMass = getMass('P',origMass,iRow,iCol,AOrig,POrig,D,S,rng);
	    _new_atomicProposal.insert(pair<unsigned long long,double>(location,newMass));
  	    extract_new_atomicProposal('P');
	    _PAtomicdomain.setProposedAtomMass(_new_atomicProposal,false);
        delLLnew = computeDeltaLL('P',D,S,AOrig,POrig,_new_nChange_matrixElemChange,_new_matrixElemChange);

	break;
      } // end of case-P block
    } // end of switch-block

    
    // This incorporates the modified change only, if any.
    switch(the_matrix_label){
    case 'A':
      { 
	   _AAtomicdomain.acceptProposal(false);
	   update_sysChi2(delLLnew);  // update system Chi2
	   return true;
	   break;
	   }
    case 'P':
      { 
	   _PAtomicdomain.acceptProposal(false);
	   update_sysChi2(delLLnew);  // update system Chi2
	   return true;
	   break;
	   }
    }

  return false;

}  // end of method birth

bool GibbsSampler::move(char the_matrix_label,
						    double const * const * D,
						    double const * const * S,
						    double ** AOrig,
						    double ** POrig)// in progress
{
  map<unsigned long long, double>::const_iterator atom;
  double chmass1,chmass2;       
  unsigned long long loc1, loc2;
  unsigned int bin1, bin2;
  double mass1, mass2;
  double newMass1, newMass2;
  atom = _atomicProposal.begin();
  chmass1 = atom->second;
  atom++;
  chmass2 = atom->second;
  if (_atomicProposal.size()==1){
  //cout << "Not doing a move due to update inconsistency."<< endl;
  return false;
  }

  // extract location, bin #, mass and changed mass corresponding to the 
  // atomic proposal such that "1" refers to a positive mass change and 
  // "2" a negative one.  
  switch(the_matrix_label){
  case 'A':
    {
      if (chmass1 > chmass2) {
	atom--;
	loc1 = atom->first;
	bin1 = _AAtomicdomain.getBin(loc1);
	mass1 = _AAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom++;
	loc2 = atom->first;
	bin2 = _AAtomicdomain.getBin(loc2);
	mass2 = _AAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      } else {
	loc1 = atom->first;
	bin1 = _AAtomicdomain.getBin(loc1);
	mass1 = _AAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom--;
	loc2 = atom->first;
	bin2 = _AAtomicdomain.getBin(loc2);
	mass2 = _AAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      }  // end of if-block for comparing chmass1 and chmass2
      break;} // end of case 'A' block
  case 'P':
    {
      if (chmass1 > chmass2) {
	atom--;
	loc1 = atom->first;
	bin1 = _PAtomicdomain.getBin(loc1);
	mass1 = _PAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom++;
	loc2 = atom->first;
	bin2 = _PAtomicdomain.getBin(loc2);
	mass2 = _PAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      } else {
	loc1 = atom->first;
	bin1 = _PAtomicdomain.getBin(loc1);
	mass1 = _PAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom--;
	loc2 = atom->first;
	bin2 = _PAtomicdomain.getBin(loc2);
	mass2 = _PAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      }  // end of if-block for comparing chmass1 and chmass2
      break;}  // end of case 'P' block
  } // end of switch-block for extracting the atomic proposal info

  // return false if bin1 == bin2
  if (bin1 == bin2){
    // cout << "Exchanges in the same bin!! So, no change!" << endl;
    //return nullMatrix;
    return false;
  }
  
  // All code having to do with exchange action was deleted here
  // Because we never use Gibbs in move, we automatically move right to
  // Metropolis-Hastings move action 

  // preparing quantities for possible Gibbs computation later.
     double priorLL = 0.;

  // Metropolis-Hasting move action

  //double pold = 0.; These quantities not necessary in move action 
  //double pnew = 0.;
  double lambda;
  switch(the_matrix_label){
  case 'A':
    {lambda = _lambdaA;
      break;}
  case 'P':
    {lambda = _lambdaP;
      break;}
  }

  /*if (pnew == 0. && pold == 0.) {
    priorLL = 0.0;
  } else if(pnew != 0. && pold == 0.) {
    priorLL = DOUBLE_POSINF;
  } else {
    priorLL = log(pnew / pold);
  } */
	 
  double delLLnew;
  switch(the_matrix_label){
  case 'A':
    {delLLnew = computeDeltaLL('A',D,S,AOrig,POrig,_nChange_matrixElemChange,_matrixElemChange);
      break;}
  case 'P':
    {delLLnew = computeDeltaLL('P',D,S,AOrig,POrig,_nChange_matrixElemChange,_matrixElemChange);
      break;}
  }

  //EJF double totalLL = priorLL + delLLnew * _annealingTemperature;

  _new_nChange_matrixElemChange = 2;
  _new_atomicProposal.insert(pair<unsigned long long, double>(loc1,newMass1-mass1));
  _new_atomicProposal.insert(pair<unsigned long long, double>(loc2,newMass2-mass2));
  switch(the_matrix_label){
  case 'A': {
    extract_new_atomicProposal('A');
    break;
   }
  case 'P': {
    extract_new_atomicProposal('P');
    break;
   }
  }

  double tmp;
  /*if (priorLL == DOUBLE_POSINF){ priorLL is always 0 in move
    //return newMatrix;
    return true;
  } else { */
    tmp = priorLL + delLLnew*_annealingTemperature;
  //}

  double rng = log(randgen('U',0,0));
 
    if (tmp < rng) {
      switch(the_matrix_label){
      case 'A':
	  { 
	  _AAtomicdomain.rejectProposal(false);
	  return false;
	  break;
	  }
      case 'P':
	  {
	  _PAtomicdomain.rejectProposal(false);
	  return false;
	  break;
	  }
     } // end of switch-block
    } else { 
      switch(the_matrix_label){
      case 'A':
	  { 	  
	   _AAtomicdomain.acceptProposal(false); 
	   update_sysChi2(delLLnew);  // update system Chi2
	   return true;
	   break;
	  }
      case 'P':
	  {           
	   _PAtomicdomain.acceptProposal(false); 
	   update_sysChi2(delLLnew);  // update system Chi2
	   return true;
	   break;
	  }
     } // end of switch-block       
    }  
 

  // end of M-H sampling

 return false;
 
} // end of method move

bool GibbsSampler::exchange(char the_matrix_label,
						    double const * const * D,
						    double const * const * S,
						    double ** AOrig,
						    double ** POrig)// in progress
{
  map<unsigned long long, double>::const_iterator atom;
  double chmass1,chmass2;       
  unsigned long long loc1, loc2;
  unsigned int bin1, bin2;
  double mass1, mass2;
  double newMass1, newMass2;
  atom = _atomicProposal.begin();
  chmass1 = atom->second;
  atom++;
  chmass2 = atom->second;

  // extract location, bin #, mass and changed mass corresponding to the 
  // atomic proposal such that "1" refers to a positive mass change and 
  // "2" a negative one.  
  switch(the_matrix_label){
  case 'A':
    {
      if (chmass1 > chmass2) {
	atom--;
	loc1 = atom->first;
	bin1 = _AAtomicdomain.getBin(loc1);
	mass1 = _AAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom++;
	loc2 = atom->first;
	bin2 = _AAtomicdomain.getBin(loc2);
	mass2 = _AAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      } else {
	loc1 = atom->first;
	bin1 = _AAtomicdomain.getBin(loc1);
	mass1 = _AAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom--;
	loc2 = atom->first;
	bin2 = _AAtomicdomain.getBin(loc2);
	mass2 = _AAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      }  // end of if-block for comparing chmass1 and chmass2
      break;} // end of case 'A' block
  case 'P':
    {
      if (chmass1 > chmass2) {
	atom--;
	loc1 = atom->first;
	bin1 = _PAtomicdomain.getBin(loc1);
	mass1 = _PAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom++;
	loc2 = atom->first;
	bin2 = _PAtomicdomain.getBin(loc2);
	mass2 = _PAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      } else {
	loc1 = atom->first;
	bin1 = _PAtomicdomain.getBin(loc1);
	mass1 = _PAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom--;
	loc2 = atom->first;
	bin2 = _PAtomicdomain.getBin(loc2);
	mass2 = _PAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      }  // end of if-block for comparing chmass1 and chmass2
      break;}  // end of case 'P' block
  } // end of switch-block for extracting the atomic proposal info

  // return nullMatrix if bin1 == bin2
  if (bin1 == bin2){
    return false;
  }


  // preparing quantities for possible Gibbs computation later.
  //EJF bool exchange = false;
  double priorLL = 0.;

  unsigned int jGene, jSample;
  bool anyNonzero = false;
  bool useGibbs = true;

  unsigned int iGene1, iPattern1, iGene2, iPattern2, iSample1, iSample2;
  switch(the_matrix_label){
  case 'A':
    {
      iGene1 = getRow('A',bin1);
      iPattern1 = getCol('A',bin1);
      iGene2 = getRow('A',bin2);
      iPattern2 = getCol('A',bin2);
      break;}
  case 'P':
    {
      iPattern1 = getRow('P',bin1);
      iSample1 = getCol('P',bin1);
      iPattern2 = getRow('P',bin2);
      iSample2 = getCol('P',bin2);
      break;}
  }

  // ---------------------------------------------
  switch(the_matrix_label){
  case 'A':
    {
      for (jSample = 0; jSample < _nCol; jSample++) {
	if (POrig[iPattern1][jSample] > epsilon) {
	  anyNonzero = true;
	  break;
	}
	if (POrig[iPattern2][jSample] > epsilon) {
	  anyNonzero = true;
	  break;
	}
      }  // end of for-block to determine the existence of corresponding 
      // non-zero elements in P
      if (!anyNonzero)  {  // cannot update in Gibbs way
	useGibbs = false;
      }
      break;}
  case 'P':
    {
      for (jGene = 0; jGene < _nRow; jGene++) {
	if (AOrig[jGene][iPattern1] > epsilon) {
	  anyNonzero = true;
	  break;
	}
	if (AOrig[jGene][iPattern2] > epsilon) {
	  anyNonzero = true;
	  break;
	}
      }  // end of for-block to determine the existence of corresponding 
      // non-zero elements in P
      if (!anyNonzero)  {  // cannot update in Gibbs way
	useGibbs = false;
      }
      break;}
  }  // end of switch-block

  // -------------------------------------------------------------------------
  // EXCHANGE ACTION when initial useGibbs = true and check if Gibbs is usable.

  double s=0.0, su, mean, sd; // EJF -- MFO check
  pair <double, double> alphaparam;

  if (useGibbs == true){

    switch(the_matrix_label){
	
	// ---------- EXCHANGE ACTION WITH A ----------------------------
    case 'A':
	{
	 if (_AAtomicdomain.inDomain(loc1) && _AAtomicdomain.inDomain(loc2))
	 {
      alphaparam = GAPSNorm::calcAlphaParameters('A', _nRow, _nCol, _nFactor, D, S, AOrig, POrig, iGene1, 
	                                              iPattern1, iGene2, iPattern2, iSample1, iSample2);
	   s = alphaparam.first;
	   su = alphaparam.second;

	  s = s * _annealingTemperature;
	  su = su * _annealingTemperature; 
	  mean = su / s;
	  sd = 1./sqrt(s);
	   // end of compute distribution parameters for A	       

	  } // end of if-block for checking whether the changes are in domain (the exchange block)
	
      break;
	  } // end of switch block for EXCHANGE ACTION with A
      
	  // ---------- EXCHANGE ACTION WITH P ----------------------------
    case 'P': 
	{
      if (_PAtomicdomain.inDomain(loc1) && _PAtomicdomain.inDomain(loc2)) 
      {
      alphaparam = GAPSNorm::calcAlphaParameters('P', _nRow, _nCol, _nFactor, D, S, AOrig, POrig, iGene1, 
	                                              iPattern1, iGene2, iPattern2, iSample1, iSample2);
	   s = alphaparam.first;
	   su = alphaparam.second;

	  s = s * _annealingTemperature;
	  su = su * _annealingTemperature; 
	  mean = su / s;
	  sd = 1./sqrt(s);
	  
	   // end of compute distribution parameters for P      

	  } // end of if-block for checking whether the changes are in domain (the exchange block)
	
      break;
	 } // end of switch block for EXCHANGE ACTION with P
 
    }// end of switch block for EXCHANGE ACTION

  } // end of if-block for operations with possibly Gibbs sampling

  if (s == 0. && su == 0.){
    useGibbs = false;
    //cout << "Parameters aren't updated -> useGibbs = false, do M-H" << endl;
  }

  // -------------------------------------------------------------------------
  if(useGibbs == true){

    // set newMass1	 
    // need to retain exponential prior
    double plower = sub_func::pnorm(-mass1, mean, sd, DOUBLE_NEGINF, 0);
    double pupper = sub_func::pnorm(mass2, mean, sd, DOUBLE_NEGINF, 0);
    double u = plower + randgen('U',0,0)*(pupper - plower);
 
    // must sample from prior if the computed parameters are not good for Gibbs
    if (plower >  0.95 || 
	pupper < 0.05 ||
	s < epsilon || 
	newMass1 == DOUBLE_POSINF ||
	newMass1 == DOUBLE_NEGINF) {
      // do not make a change
      useGibbs = false;

    }
	
	double gibbsMass1, gibbsMass2;

    if (useGibbs == true){

      gibbsMass1 = sub_func::qnorm(u, mean, sd, DOUBLE_NEGINF, 0);
      if (gibbsMass1 < -mass1) gibbsMass1 = -mass1;
      if (gibbsMass1 > mass2) gibbsMass1 = mass2;
      gibbsMass2 = - gibbsMass1;

      // update new masses
      double delLLnew;
      _new_nChange_matrixElemChange = 2;
      _new_atomicProposal.insert(pair<unsigned long long, double>(loc1,gibbsMass1));
      _new_atomicProposal.insert(pair<unsigned long long, double>(loc2,gibbsMass2));
      switch(the_matrix_label){
      case 'A':
	 {
	  extract_new_atomicProposal('A');
	  delLLnew = computeDeltaLL('A',D,S,AOrig,POrig,_new_nChange_matrixElemChange,_new_matrixElemChange);
	  _AAtomicdomain.setProposedAtomMass(_new_atomicProposal, false);
	  _AAtomicdomain.acceptProposal(false); 
	  update_sysChi2(delLLnew);  // update system Chi2
	  
	  break;
	  }
      case 'P':
	 {
	  extract_new_atomicProposal('P');
	  delLLnew = computeDeltaLL('P',D,S,AOrig,POrig,_new_nChange_matrixElemChange,_new_matrixElemChange);
	  _PAtomicdomain.setProposedAtomMass(_new_atomicProposal, false);
	  _PAtomicdomain.acceptProposal(false);
	  update_sysChi2(delLLnew);  // update system Chi2
	  break;
	  }
     }  // end of switch-block       
      
      return true;

    } // end of inner if-block for final updating with Gibbs

  } // end of outer if-block for useGibbs == true && s > epsilon


  // ----------------------------
  // We can't use Gibbs, need
  // Metropolis-Hasting exchange action

  double pold = 0.;
  double pnew = 0.;
  double lambda;
  switch(the_matrix_label){
  case 'A':
    {lambda = _lambdaA;
      break;}
  case 'P':
    {lambda = _lambdaP;
      break;}
  }
  
  // Formerly the if(_oper_type == 'E') block in move_exchange
    if (mass1 > mass2) {
      pnew = sub_func::dgamma(newMass1, 2., 1./lambda, false);
      if (newMass1 > newMass2) {
	pold = sub_func::dgamma(mass1, 2., 1./lambda, false);
      } else {
	pold = sub_func::dgamma(mass2, 2., 1./lambda, false);
      }
    } else {
      pnew = sub_func::dgamma(newMass2, 2., 1./lambda, false);
      if (newMass1 > newMass2) {
	pold = sub_func::dgamma(mass1, 2., 1./lambda,  false);
      } else {
	pold = sub_func::dgamma(mass2, 2., 1./lambda, false);
      }
    }  


  
  if (pnew == 0. && pold == 0.) {
    priorLL = 0.0;
  } else if(pnew != 0. && pold == 0.) {
    priorLL = DOUBLE_POSINF;
  } else {
    priorLL = log(pnew / pold);
  }
	 
  double delLLnew;
  switch(the_matrix_label){
  case 'A':
    {delLLnew = computeDeltaLL('A',D,S,AOrig,POrig,_nChange_matrixElemChange,_matrixElemChange);
      break;}
  case 'P':
    {delLLnew = computeDeltaLL('P',D,S,AOrig,POrig,_nChange_matrixElemChange,_matrixElemChange);
      break;}
  }

  //EJF double totalLL = priorLL + delLLnew * _annealingTemperature;

  _new_nChange_matrixElemChange = 2;
  _new_atomicProposal.insert(pair<unsigned long long, double>(loc1,newMass1-mass1));
  _new_atomicProposal.insert(pair<unsigned long long, double>(loc2,newMass2-mass2));
  switch(the_matrix_label){
  case 'A': {
    extract_new_atomicProposal('A');
    break;
  }
  case 'P': {
    extract_new_atomicProposal('P');
    break;
  }
  }

  double tmp;
  if (priorLL == DOUBLE_POSINF){
    return true;
  } else {
    tmp = priorLL + delLLnew*_annealingTemperature;
  }

  double rng = log(randgen('U',0,0));

    if (tmp  < rng) {
      switch(the_matrix_label){
      case 'A':
	{ _AAtomicdomain.rejectProposal(false);
	  return false;
	  break;}
      case 'P':
	{ _PAtomicdomain.rejectProposal(false);
	  return false;
	  break;}
      } // end of switch-block
    } else {
 
      switch(the_matrix_label){
      case 'A':
	{ 	  
      _AAtomicdomain.acceptProposal(false); 
	  update_sysChi2(delLLnew);  // update system Chi2
	  return true;
          break;
        }
      case 'P':
	{           
        _PAtomicdomain.acceptProposal(false); 
	   update_sysChi2(delLLnew);  // update system Chi2
	   return true;
          break;
        }
      } // end of switch-block       
    }  

  // end of M-H sampling
 
  return false;

} // end of method exchange


// ************ METHODS FOR LOOPING AND CONTROL *****************************
void GibbsSampler::set_iter(unsigned long ext_iter){
  _iter = ext_iter;
} 

double GibbsSampler::get_AnnealingTemperature(){
  return _annealingTemperature;
}


// -----------------------------------------------------------------------------
void GibbsSampler::set_AnnealingTemperature() { 
  double SASteps = _nEquil;
  double SATemp = ( (double) _iter + 1. ) / (SASteps / 2.);

  if (SATemp > 1.) SATemp = 1;
  if (SATemp < 0) {
    throw logic_error("Invalid annealing temperature.");
  }

  _annealingTemperature = SATemp;
}


// -----------------------------------------------------------------------------
// 
void GibbsSampler::check_atomic_matrix_consistency(char the_matrix_label)
{
  double total_atom_mass = 0.0;
  double total_matrix_mass = 0.0;
  
  switch(the_matrix_label){
  case 'A': {
    total_atom_mass = _AAtomicdomain.get_atomicDomain_totalmass();
    total_matrix_mass = _AMatrix.cal_totalsum();
    break;
  } 
  case 'P': {
    total_atom_mass = _PAtomicdomain.get_atomicDomain_totalmass();
    total_matrix_mass = _PMatrix.cal_totalsum();
    break;
  }
  } // end of switch-block

  double diff_total_mass = fabs(total_atom_mass - total_matrix_mass);

  if(diff_total_mass > 1.e-5){
	
	/*
    cout << "Mass inconsistency!! Total mass difference = " << diff_total_mass << endl;
    cout << "total atom mass = " << total_atom_mass << endl;
    cout << "total matrix mass = " << total_matrix_mass << endl;
    cout << "Oper_type = " << _oper_type << endl;
	*/
	
    throw logic_error("Mass inconsistency between atomic domain and matrix!");
  } 
  
}

// -----------------------------------------------------------------------------
void GibbsSampler::compute_statistics_prepare_matrices(unsigned long statindx){

  double ** A = _AMatrix.get_matrix();
  double ** P = _PMatrix.get_matrix();
  vector<double> k(_nFactor);  // normalized vector

  // compute the normalization vector
  for (int m=0; m< _nFactor; ++m){
    k[m] = 0.;
    for (int n=0; n< _nCol; ++n){
      k[m] += P[m][n];
    }
    if (k[m] == 0){  // when the whole row of P is zero, then don't do anything
      k[m] = 1.0;
    }
  }


  // construct the mean and var matrices at statindx = 1
  if (statindx == 1){
        
    _Amean = new double * [_nRow];
    for (int m=0; m < _nRow ; ++m)
      {_Amean[m] = new double [_nFactor];}
    for (int m=0; m < _nRow; ++m){
      for (int n=0; n < _nFactor; ++n){
	_Amean[m][n] = A[m][n] * k[n];
      }
    }

    _Asd = new double * [_nRow];
    for (int m=0; m < _nRow ; ++m)
      {_Asd[m] = new double [_nFactor];}
    for (int m=0; m < _nRow; ++m){
      for (int n=0; n < _nFactor; ++n){
	_Asd[m][n] = pow(A[m][n]*k[n],2);
      }
    } 

    _Pmean = new double * [_nFactor];
    for (int m=0; m < _nFactor ; ++m)
      {_Pmean[m] = new double [_nCol];}
    for (int m=0; m < _nFactor; ++m){
      for (int n=0; n < _nCol; ++n){
	_Pmean[m][n] = P[m][n] / k[m];
      }
    } 

    _Psd = new double * [_nFactor];
    for (int m=0; m < _nFactor ; ++m)
      {_Psd[m] = new double [_nCol];}
    for (int m=0; m < _nFactor; ++m){
      for (int n=0; n < _nCol; ++n){
	_Psd[m][n] = pow(P[m][n] / k[m] , 2);
      }
    } 

  } // end of if-block for matrix construction statindx == 1

  // increment the mean and var matrices at statindx != 1
  if (statindx > 1){

    for (int m=0; m < _nRow; ++m){
      for (int n=0; n < _nFactor; ++n){
	_Amean[m][n] += A[m][n] * k[n];
      }
    }

    for (int m=0; m < _nRow; ++m){
      for (int n=0; n < _nFactor; ++n){
	_Asd[m][n] += pow(A[m][n]*k[n],2);
      }
    }

    for (int m=0; m < _nFactor; ++m){
      for (int n=0; n < _nCol; ++n){
	_Pmean[m][n] += P[m][n] / k[m];
      }
    }

    for (int m=0; m < _nFactor; ++m){
      for (int n=0; n < _nCol; ++n){
	_Psd[m][n] += pow(P[m][n] / k[m],2);
      }
    }

  } // end of if-block for matrix incrementation statindx > 1

} // end of method compute_statistics_prepare_matrices

// -----------------------------------------------------------------------------
void GibbsSampler::compute_statistics(char outputFilename[],
                                      char outputAmean_Filename[],char outputAsd_Filename[],
                                      char outputPmean_Filename[],char outputPsd_Filename[],
				                      char outputAPmean_Filename[],
                                      unsigned int Nstat){

  // compute statistics for A
  for (int m=0; m < _nRow ; ++m){
    for (int n=0; n<_nFactor; ++n){
      _Amean[m][n] = _Amean[m][n] / Nstat;
    }
  }

  for (int m=0; m < _nRow ; ++m){
    for (int n=0; n<_nFactor; ++n){
      _Asd[m][n] = sqrt((_Asd[m][n] - Nstat*pow(_Amean[m][n],2)) / (Nstat-1));
    }
  }

  // compute statistics for P
  for (int m=0; m < _nFactor ; ++m){
    for (int n=0; n<_nCol; ++n){
      _Pmean[m][n] = _Pmean[m][n] / Nstat;
    }
  }

  for (int m=0; m < _nFactor ; ++m){
    for (int n=0; n<_nCol; ++n){
      _Psd[m][n] = sqrt((_Psd[m][n] - Nstat*pow(_Pmean[m][n],2)) / (Nstat-1));
    }
  }


  double ** APmean;
  APmean = new double * [_nRow];
  for (int m=0; m < _nRow ; ++m)
    {APmean[m] = new double [_nCol];}
  for (unsigned int m=0; m < _nRow; ++m){
    for (unsigned int n=0; n < _nCol; ++n){
      APmean[m][n]=0.0;
      for (unsigned int k=0; k < _nFactor; ++k){
	APmean[m][n] += _Amean[m][k]*_Pmean[k][n];
      }
    }
  }



  ofstream outputFile;
  outputFile.open(outputFilename,ios::out|ios::app);

  outputFile << " ************************************************* " << endl;
  outputFile << " ---------- OUTPUT FINAL STATISTICS -------------- " << endl;
  outputFile << " ************************************************* " << endl;

  outputFile << " Number of samples for computing statistics = " << Nstat << endl;
  outputFile << " Amean = " << endl << endl;
  local_display_matrix2F(outputFile,_Amean,_nRow,_nFactor);
  outputFile << endl;
  outputFile << " Asd = " << endl << endl;
  local_display_matrix2F(outputFile,_Asd,_nRow,_nFactor);
  outputFile << endl;
  outputFile << " Pmean = " << endl << endl;
  local_display_matrix2F(outputFile,_Pmean,_nFactor,_nCol);
  outputFile << endl;
  outputFile << " Psd = " << endl << endl;
  local_display_matrix2F(outputFile,_Psd,_nFactor,_nCol);
  outputFile << endl;
  outputFile << "The product Amean*Pmean gives " << endl << endl;
  local_display_matrix2F(outputFile,APmean,_nRow,_nCol);
  outputFile << endl;

  outputFile.close();

  // independent output of Amean, Asd, Pmean, Psd, APmean
  ofstream output_Amean;
  output_Amean.open(outputAmean_Filename,ios::out);
  local_display_matrix2F(output_Amean,_Amean,_nRow,_nFactor);
  output_Amean.close();

  ofstream output_Asd;
  output_Asd.open(outputAsd_Filename,ios::out);
  local_display_matrix2F(output_Asd,_Asd,_nRow,_nFactor);
  output_Asd.close();

  ofstream output_Pmean;
  output_Pmean.open(outputPmean_Filename,ios::out);
  local_display_matrix2F(output_Pmean,_Pmean,_nFactor,_nCol);
  output_Pmean.close();

  ofstream output_Psd;
  output_Psd.open(outputPsd_Filename,ios::out);
  local_display_matrix2F(output_Psd,_Psd,_nFactor,_nCol);
  output_Psd.close();

  ofstream output_APmean;
  output_APmean.open(outputAPmean_Filename,ios::out);
  local_display_matrix2F(output_APmean,APmean,_nRow,_nCol);
  output_APmean.close();



}


void GibbsSampler::compute_statistics(unsigned int Nstat,
									vector<vector <double> > &AMeanVect,
									vector<vector <double> > &AStdVect,
									vector<vector <double> > &PMeanVect,
									vector<vector <double> > &PStdVect
									)
{
	//Conor's Code for resizing the vectors corresponding to each matrix
	AMeanVect.resize(_nRow);
	for (int i = 0; i < _nRow; i++)
	{
		AMeanVect[i].resize(_nFactor);
	}
	AStdVect.resize(_nRow);
	for (int i = 0; i < _nRow; i++)
	{
		AStdVect[i].resize(_nFactor);
	}
	PMeanVect.resize(_nFactor);
	for (int i = 0; i < _nFactor; i++)
	{
		PMeanVect[i].resize(_nCol);
	}
	PStdVect.resize(_nFactor);
	for (int i = 0; i < _nFactor; i++)
	{
		PStdVect[i].resize(_nCol);
	}  //End vector resizing
	
	
  // compute statistics for A
	//Variable for holding temp Calculations
	double tempStat;

	//Changed to compute as vectors to pass into R
  for (int m=0; m < _nRow ; ++m){
    for (int n=0; n<_nFactor; ++n){
		tempStat = _Amean[m][n] / Nstat;
      _Amean[m][n] = tempStat;
	  AMeanVect[m][n] = tempStat;
    }
  }
  for (int m=0; m < _nRow ; ++m){
    for (int n=0; n<_nFactor; ++n){
      tempStat = sqrt((_Asd[m][n] - Nstat*pow(_Amean[m][n],2)) / (Nstat-1));
	  _Asd[m][n] = tempStat;
	  AStdVect[m][n] = tempStat;
    }
  }

  // compute statistics for P
  for (int m=0; m < _nFactor ; ++m){
    for (int n=0; n<_nCol; ++n){
      tempStat = _Pmean[m][n] / Nstat;
	  _Pmean[m][n] = tempStat;
	  PMeanVect[m][n] = tempStat;
    }
  }

  for (int m=0; m < _nFactor ; ++m){
    for (int n=0; n<_nCol; ++n){
      tempStat = sqrt((_Psd[m][n] - Nstat*pow(_Pmean[m][n],2)) / (Nstat-1));
	  _Psd[m][n] = tempStat;
	  PStdVect[m][n] = tempStat;
    }
  }
}

// *****************************************************************************
// Adaptation from the original code:

bool GibbsSampler::checkOtherMatrix(char the_matrix_label, unsigned int iRow, unsigned int iCol,
				     double const * const * otherMatrix){

  unsigned int otherDim;

    // check that there is mass for this pattern in the P Matrix
  if (the_matrix_label=='A'){ 
   for (otherDim = 0; otherDim < _nCol; otherDim++) {
	if (otherMatrix[iCol][otherDim] > epsilon) {
	  return true;
	 }
    }
   return false;
   }

    // check that there is mass for this pattern in the A Matrix
  else{
   for (otherDim = 0; otherDim < _nRow; otherDim++) {
	if (otherMatrix[otherDim][iRow] > epsilon) {
	  return true;
	}
   }
   return false;

   }

  return false; 
}

//-----------------------------------------------------------------

//-----------------------------------------------------------------
double GibbsSampler::getMass(char the_matrix_label, double origMass,
			     unsigned int iRow,
			     unsigned int iCol,
			     double const * const * otherMatrix, 
			     double const * const * currentChainMatrix,
			     double const * const * D, double const * const * S,
			     double rng)
{
  double DOUBLE_POSINF = numeric_limits<double>::max();
  unsigned int iGene = 0;
  unsigned int iPattern = 0;
  unsigned int iSample = 0; 
  unsigned int jPattern = 0;
  double lambda;

  switch(the_matrix_label)
    {
    case 'A': 
      {
	iGene = iRow;
	iPattern = iCol;
	lambda = _AAtomicdomain.getLambda();
	break;
      }
    case 'P': 
      {
	iPattern = iRow;
	iSample = iCol;
	lambda = _PAtomicdomain.getLambda();
	break;
      }
    }

  // determine the parameters for finding the mass
  double s  = 0.;
  double su = 0.;
  double mock;

  switch(the_matrix_label){

  case 'A': 
    {
      //double Aeff;
      for (iSample = 0; iSample < _nCol; iSample++) 
	{
	  mock = D[iGene][iSample];
	  for (jPattern = 0; jPattern <_nFactor; jPattern++) 
	    {
              mock -= currentChainMatrix[iGene][jPattern]* otherMatrix[jPattern][iSample];
	    }
	  s += _annealingTemperature * pow(otherMatrix[iPattern][iSample],2) /
	    ( 2* pow(S[iGene][iSample],2) );
	  su += _annealingTemperature * otherMatrix[iPattern][iSample]*mock /
	    ( 2* pow(S[iGene][iSample],2) );	    
	}
      break;
    }

  case 'P': 
    {
      for (iGene = 0; iGene < _nRow; iGene++) 
	{
	  mock = D[iGene][iSample];
	  for (jPattern = 0; jPattern < _nFactor; jPattern++) 
	    {
	      mock -= otherMatrix[iGene][jPattern] * currentChainMatrix[jPattern][iSample];
	    }
	  s += _annealingTemperature * pow(otherMatrix[iGene][iPattern],2) / 
	    ( 2* pow(S[iGene][iSample],2) );
	  su += _annealingTemperature * otherMatrix[iGene][iPattern]*mock / 
	    ( 2* pow(S[iGene][iSample],2) );
	}
      break;
    }

  }

    
  double mean  = (2*su - lambda)/(2*s);
  double sd = 1./sqrt(2*s);
    
  // note: is bounded below by zero so have to use inverse sampling!
  // based upon algorithm in DistScalarRmath.cc (scalarRandomSample)

  double plower = sub_func::pnorm(0., mean, sd, DOUBLE_NEGINF, 0);
  double pupper = 1.;
  double u = plower + randgen('U',0,0) * (pupper - plower);
  // -------------------------------------------------------------
  // this line seems to be misplaced.
  // double newMass = sub_func::qnorm(u, mean, sd, DOUBLE_NEGINF, 0);
  double newMass = 0;
  // ------------------

  // if the likelihood is flat and nonzero, 
  // force to sample strictly from the prior
  if ( plower == 1 || s < 1.e-5 || newMass == DOUBLE_POSINF || pupper == 0) { 
    if (origMass < 0) {    // death case
      newMass = abs(origMass);
    } else {
      newMass = 0.;  // birth case
    }
  } // end of first comparison
  else if (plower >= 0.99) {
    double tmp1 = sub_func::dnorm(0, mean, sd, false); 
    double tmp2 = sub_func::dnorm(10*lambda, mean, sd, false);
    if ( (tmp1 > epsilon) && (fabs(tmp1-tmp2) < epsilon) )   {
      if (origMass < 0) {   // death case
	return 0.;
      }
      return origMass;   // birth case
    }  
  } // end of second comparison
  else {
    newMass = sub_func::qnorm(u, mean, sd, DOUBLE_NEGINF, 0);  // both death and birth
  }  // end of if-block for the remaining case


  // limit the mass range
    switch(the_matrix_label) {
    case 'A':
      {
	if (newMass > _max_gibbsmassA)
	  newMass = _max_gibbsmassA;
	break;
      }
    case 'P':
      {
	if (newMass > _max_gibbsmassP)
	  newMass = _max_gibbsmassP;
	break;
      }
    }

  if (newMass < 0) newMass = 0;   // due to the requirement that newMass > 0 

  return newMass;

}

// -----------------------------------------------------------------------------
void GibbsSampler::detail_check(char outputchi2_Filename[]){

      double chi2 = 2.*cal_logLikelihood();
	  /*
      cout << "oper_type: " << _oper_type <<
              " ,nA: " << getTotNumAtoms('A') <<
	      " ,nP: " << getTotNumAtoms('P') << 
              " ,_sysChi2 = " << _sysChi2 << endl;
		*/


      ofstream outputchi2_File;
      outputchi2_File.open(outputchi2_Filename,ios::out|ios::app);
      outputchi2_File << chi2 << endl;
      outputchi2_File.close();

}
