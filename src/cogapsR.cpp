// cogapsR.cpp

// =============================================================================
// This is the main code for Cogaps. (7th Sep, 2013)
// =============================================================================

#include <iostream>       // for use with standard I/O
#include <string>		  // for string processing
#include <fstream>        // for output to files
#include <limits>         // for extracting numerical limits of C++


// ------ incorporated to use Cogaps_options ------------
#include <vector>
#include <iomanip>
#include <boost/algorithm/string.hpp>
// ------------------------------------------------------
#include "randgen.h"   // for incorporating a random number generator.
#include "Matrix.h"    // for incorporating a Matrix class
#include "AtomicSupport.h"  // for incorporating an Atomic class
#include "GAPSNorm.h"  // for incorporating calculation of statistics in cogaps.
#include "GibbsSampler.h" // for incorporating the GibbsSampler which
                           // does all the atomic space to matrix conversion
                           // and sampling actions.
#include <Rcpp.h>
// ------------------------------------------------------

//namespace bpo = boost::program_options;
using namespace std;
using namespace gaps;
using std::vector;

boost::mt19937 rng(43);

// [[Rcpp::export]]
Rcpp::List cogaps(Rcpp::DataFrame DFrame, Rcpp::DataFrame SFrame, Rcpp::CharacterVector Config){

  // ===========================================================================
  // Initialization of the random number generator.
  // Different seeding methods:
  // --- fixed seed 
  //std::vector<unsigned long> ve(2);
  //ve[0]=198782;ve[1]=89082;
  //boost::random::seed_seq seq(ve);
  //rng.seed(seq);
  // --- seeded with time
  rng.seed(static_cast<boost::uint32_t>(std::time(0)));
  //---------------------

  // ===========================================================================
  // Part 1) Initialization: 
  // In this section, we read in the system parameters from the paremter file
  // parameter.txt, and matrices D and S from datafile.txt. 
  // Then we initialize A and P in both their atomic domains and 
  // matrix forms.
  // ===========================================================================



//---------------------------------------------
//CONOR'S CODE FOR CHANGING R VARIABLES OF THE CONFIG FILE INTO C VARIABLES
//THIS REPLACES COGAPS_PROGAM_OPTIONS
  string temp;
  temp = Rcpp::as<string>(Config[2]);
  unsigned long nEquil = atol(temp.c_str()); //Cogaps_options.nEquil;    // # outer loop iterations 
   
  temp = Rcpp::as<string>(Config[3]);  // for equilibration
  unsigned long nSample =  atol(temp.c_str());//Cogaps_options.nSample;  // # outer loop iterations 
    
  temp = Rcpp::as<string>(Config[4]);  // for equilibration
  unsigned long nObsR =  atol(temp.c_str());//Cogaps_options.nSample;  // # outer loop iterations 
    
  temp = Rcpp::as<string>(Config[0]);// for sampling
  unsigned int nFactor = atoi(temp.c_str()); //Cogaps_options.nFactor;   // # patterns
  
  temp = Rcpp::as<string>(Config[6]);
  double alphaA = atof(temp.c_str());//Cogaps_options.alphaA;
  
  temp = Rcpp::as<string>(Config[10]);
  double alphaP = atof(temp.c_str());//Cogaps_options.alphaP
  
  temp = Rcpp::as<string>(Config[7]);
  double nMaxA = atof(temp.c_str());//Cogaps_options.nMaxA;             // max. number of atoms in A
  
  temp = Rcpp::as<string>(Config[11]);
  double nMaxP = atof(temp.c_str());//Cogaps_options.nMaxP;             // number of atomic bins for P

  //---------------------------------------------------------------------------------
  //string datafile = Rcpp::as<string>(Config[1]);//Cogaps_options.datafile;        // File for D
  //string variancefile = Rcpp::as<string>(Config[2]); //Cogaps_options.variancefile; // File for S
  //---------------------------------------------------------------------------------
  string simulation_id = Rcpp::as<string>(Config[1]); //Cogaps_options.simulation_id; // simulation id
  
  temp = Rcpp::as<string>(Config[8]);
  double max_gibbsmass_paraA = atof(temp.c_str()); //Cogaps_options.max_gibbsmass_paraA; 
                           // maximum gibbs mass parameter for A 
  temp = Rcpp::as<string>(Config[12]);
  double max_gibbsmass_paraP = atof(temp.c_str());//Cogaps_options.max_gibbsmass_paraP; 
                           // maximum gibbs mass parameter for P 
						   
  temp = Rcpp::as<string>(Config[9]);					
  double lambdaA_scale_factor = atof(temp.c_str()); //Cogaps_options.lambdaA_scale_factor;
                           // scale factor for lambdaA
						
  temp = Rcpp::as<string>(Config[13]);
  double lambdaP_scale_factor = atof(temp.c_str());//Cogaps_options.lambdaP_scale_factor;
                           // scale factor for lambdaP
						   
  temp = Rcpp::as<string>(Config[4]);
  bool Q_output_atomic;
  	if(temp == "TRUE" || temp == "true")
		Q_output_atomic = true;	//Cogaps_options.Q_output_atomic;
                           // whether to output the atomic space info
	else
		Q_output_atomic = false;

	
   //Code to make the D and S matrices read from R into C++ vectors to make into Elana's Matrix Objects in Matrix.cpp
   vector<vector<double> > DVector;
   vector<vector<double> > SVector;
   
   //Code to establish the sizes and initialize the C++ vectors to pass 
    int numC = DFrame.size();
	Rcpp::NumericVector tempFrameCol = DFrame[0];
	int numR = tempFrameCol.size() ;
	double tempFrameElement;
	DVector.resize(numR);
	for (int i = 0; i < numR; i++)
	{
		DVector[i].resize(numC);
	}
	for(int i=0; i<numR; i++)
	{
		for(int j=0; j<numC; j++)
		{
			tempFrameCol = DFrame[j];
			tempFrameElement = tempFrameCol[i];
			DVector[i][j] = tempFrameElement;
			//Rcpp::Rcout << DVector[i][j] << endl;
		}
	}
	
	/* Test output of vector to assure it is the correct D
	for(int i = 0; i < DVector.size(); i++)
	{
		for(int j = 0; j < DVector[0].size(); j++)
		{
			Rcpp::Rcout<< DVector[i][j] << "\t";
		}
		Rcpp::Rcout << endl;
	}
	*/
	
	numC = SFrame.size();
	tempFrameCol = SFrame[0];
	numR = tempFrameCol.size() ;
	SVector.resize(numR);
	for (int i = 0; i < numR; i++)
	{
		SVector[i].resize(numC);
	}
	for(int i=0; i<numR; i++)
	{
		for(int j=0; j<numC; j++)
		{
			tempFrameCol = SFrame[j];
			tempFrameElement = tempFrameCol[i];
			SVector[i][j] = tempFrameElement;
		}
	}
	//--------------------END CREATING D AND S C++ VECTORS


  // Parameters or structures to be calculated or constructed:
  unsigned long nIterA = 10;    // initial inner loop iterations for A
  unsigned long nIterP = 10;    // initial inner loop iterations for P  
  unsigned long atomicSize = 0; // number of atomic points

  char label_A = 'A';  // label for matrix A
  char label_P = 'P';  // label for matrix P
  char label_D = 'D';  // label for matrix D
  char label_S = 'S';// label for matrix S

  // Output parameters and computing info to files:
  
  
  //---------------------------------------------------------------------------
  //ENTIRE SECTION ONLY RELEVANT FOR R VERSION 


  // ---------------------------------------------------------------------------
  // Initialize the GibbsSampler.

  //R Version
  GibbsSampler GibbsSamp(nEquil,nSample,nFactor,   // construct GibbsSampler and 
                         alphaA,alphaP,nMaxA,nMaxP,// Read in D and S matrices
                         nIterA,nIterP,
			 max_gibbsmass_paraA, max_gibbsmass_paraP, 
			 lambdaA_scale_factor, lambdaP_scale_factor,
                         atomicSize,
                         label_A,label_P,label_D,label_S,
			 DVector,SVector,simulation_id);
			 

  // ---------------------------------------------------------------------------
  // Based on the information of D, construct and initialize for A and P both 
  // the matrices and atomic spaces.

  GibbsSamp.init_AMatrix_and_PMatrix(); // initialize A and P matrices
  GibbsSamp.init_AAtomicdomain_and_PAtomicdomain(); // intialize atomic spaces
                                                    // A and P
  GibbsSamp.init_sysChi2(); // initialize the system chi2 value

  // ===========================================================================
  // Part 2) Equilibration:
  // In this section, we let the system eqilibrate with nEquil outer loop 
  // iterations. Within each outer loop iteration, A is iterated nIterA times 
  // and P is iterated nIterP times. After equilibration, we update nIterA and 
  // nIterP according to the expected number of atoms in the atomic spaces 
  // of A and P respectively.
  // ===========================================================================


    double chi2;
    double tempChiSq;
  double tempAtomA;
  double tempAtomP;
	int outCount = 0;
	int numOutputs = nObsR;   //CONOR'S OUTPUT CONTROL
	int totalChiSize = nSample + nEquil;
	Rcpp::NumericVector chiVect(totalChiSize); //INITIALIZE THE VECTOR HOLDING THE CHISQUARE.
	Rcpp::NumericVector nAEquil(nEquil);	   //INITIALIZE THE VECTOR HOLDING THE ATOMS FOR EACH MATRIX FOR EACH EQUIL/SAMP
	Rcpp::NumericVector nASamp(nSample);
	Rcpp::NumericVector nPEquil(nEquil);
	Rcpp::NumericVector nPSamp(nSample);

  for (unsigned long ext_iter=1; ext_iter <= nEquil; ++ext_iter){
    GibbsSamp.set_iter(ext_iter);
    GibbsSamp.set_AnnealingTemperature();


    for (unsigned long iterA=1; iterA <= nIterA; ++iterA){
      GibbsSamp.update('A');
      //GibbsSamp.check_atomic_matrix_consistency('A');
      //GibbsSamp.detail_check(outputchi2_Filename);
    }
    GibbsSamp.check_atomic_matrix_consistency('A');

    for (unsigned long iterP=1; iterP <= nIterP; ++iterP){
      GibbsSamp.update('P');
      //GibbsSamp.check_atomic_matrix_consistency('P');
      //GibbsSamp.detail_check(outputchi2_Filename);
    }
    GibbsSamp.check_atomic_matrix_consistency('P');
	
	//Finds the current ChiSq and places it into the vector to be returned to R (and output on occasion)
	tempChiSq = GibbsSamp.get_sysChi2();
	//Rcpp::Rcout << (ext_iter) - 1 << endl;
	chiVect[(ext_iter)-1] = tempChiSq;
    // ----------- output computing info ---------
	tempAtomA = GibbsSamp.getTotNumAtoms('A');
	tempAtomP = GibbsSamp.getTotNumAtoms('P');
	nAEquil[outCount] = tempAtomA;
	nPEquil[outCount] = tempAtomP;
	outCount++;
    if ( ext_iter % numOutputs == 0){
      //chi2 = 2.*GibbsSamp.cal_logLikelihood();
	  

	  //COMMENTED OUT FOR THE R VERSION
      //GibbsSamp.output_computing_info(outputFilename,ext_iter,nEquil,0,nSample
	  //---------------------------------
      Rcpp::Rcout << "Equil:" << ext_iter << " of " << nEquil << 
              ", Atoms:" << tempAtomA << "("
	      << tempAtomP << ")" <<
	// " ,chi2 = " << chi2 <<
              "  Chi2 = " << tempChiSq << endl;
		
    }

    // -------------------------------------------
    // re-calculate nIterA and nIterP to the expected number of atoms 
    nIterA = (unsigned long) randgen('P',max((double) GibbsSamp.getTotNumAtoms('A'),10.));
    nIterP = (unsigned long) randgen('P',max((double) GibbsSamp.getTotNumAtoms('P'),10.));
    //nIterA = (unsigned long) randgen('P',(double) GibbsSamp.getTotNumAtoms('A')+10.);
    //nIterP = (unsigned long) randgen('P',(double) GibbsSamp.getTotNumAtoms('P')+10.);
    // --------------------------------------------

  }  // end of for-block for equilibration
 


  // ===========================================================================
  // Part 3) Sampling:
  // After the system equilibriates in Part 2, we sample the systems with an 
  // outer loop of nSample iterations. Within each outer loop iteration, A is 
  // iterated nIterA times and P is iterated nIterP times. After sampling, 
  // we update nIterA and nIterP according to the expected number of atoms in 
  // the atomic spaces of A and P respectively.
  // ===========================================================================
  


  unsigned int statindx = 0;
  outCount = 0;
  for (unsigned long i=1; i <= nSample; ++i){
    for (unsigned long iterA=1; iterA <= nIterA; ++iterA){
      GibbsSamp.update('A');
      //GibbsSamp.check_atomic_matrix_consistency('A');
      //GibbsSamp.detail_check(outputchi2_Filename);
    }
    GibbsSamp.check_atomic_matrix_consistency('A');

    for (unsigned long iterP=1; iterP <= nIterP; ++iterP){ 
      GibbsSamp.update('P');
      //GibbsSamp.check_atomic_matrix_consistency('P');
      //GibbsSamp.detail_check(outputchi2_Filename);
    }
    GibbsSamp.check_atomic_matrix_consistency('P');

    if (Q_output_atomic == true){
       GibbsSamp.output_atomicdomain('A',i);
       GibbsSamp.output_atomicdomain('P',i);
    }

     statindx += 1;
     GibbsSamp.compute_statistics_prepare_matrices(statindx);
	
	//Do the same as above.
	tempChiSq = GibbsSamp.get_sysChi2();
	//Rcpp::Rcout << (nEquil + i) - 1 << endl;
	chiVect[(nEquil + i)-1] = tempChiSq;
    // ----------- output computing info ---------
		tempAtomA = GibbsSamp.getTotNumAtoms('A');
		tempAtomP = GibbsSamp.getTotNumAtoms('P');
		nASamp[outCount] = tempAtomA;
		nPSamp[outCount] = tempAtomP;
		outCount++;
    if ( i % numOutputs == 0){
	
      // chi2 = 2.*GibbsSamp.cal_logLikelihood();
      //GibbsSamp.output_atomicdomain('A',(unsigned long) statindx);
      //GibbsSamp.output_atomicdomain('P',(unsigned long) statindx
	  
	  //COMMENTED OUT FOR THE R VERSION
      //GibbsSamp.output_computing_info(outputFilename,nEquil,nEquil,i,nSample);
	  //----------------------------------
        Rcpp::Rcout << "Samp: " << i << " of " << nSample <<
        ", Atoms:" << tempAtomA << "("
        << tempAtomP << ")" <<
        // " ,chi2 = " << chi2 <<
        "  Chi2 = " << tempChiSq << endl;
		
	
		
      if (i == nSample){
         chi2 = 2.*GibbsSamp.cal_logLikelihood();
	 Rcpp::Rcout << " *** Check value of final chi2: " << chi2 << " **** " << endl; 
      }


    }

    // -------------------------------------------
    // re-calculate nIterA and nIterP to the expected number of atoms 
    nIterA = (unsigned long) randgen('P',max((double) GibbsSamp.getTotNumAtoms('A'),10.));
    nIterP = (unsigned long) randgen('P',max((double) GibbsSamp.getTotNumAtoms('P'),10.));
    //nIterA = (unsigned long) randgen('P',(double) GibbsSamp.getTotNumAtoms('A')+10.);
    //nIterP = (unsigned long) randgen('P',(double) GibbsSamp.getTotNumAtoms('P')+10.);
    // --------------------------------------------

  }  // end of for-block for Sampling

 

  // ===========================================================================
  // Part 4) Calculate statistics:
  // In this final section, we calculate all statistics pertaining to the final
  // sample and check the results.
  // ===========================================================================

  	vector<vector <double> > AMeanVector;
	vector<vector <double> > AStdVector;
	vector<vector <double> > PMeanVector;
	vector<vector <double> > PStdVector;


  GibbsSamp.compute_statistics(statindx,
							   AMeanVector, AStdVector, PMeanVector, PStdVector);          // compute statistics like mean and s.d.

//Test To make sure vectors are solid!	
/*
	for(int i = 0; i < AMeanVector.size(); i++)
	{
		for(int j = 0; j < AMeanVector[i].size(); j++)
		{
			Rcpp::Rcout << AMeanVector[i][j] << "  ";
		}
		Rcpp::Rcout << endl;
	}
*/
 
	//CODE FOR CONVERTING ALL THE VECTORS FOR THE DATAFILES INTO NUMERIC VECTORS OR LISTS
	
	int numRow = AMeanVector.size();
	int numCol = AMeanVector[0].size() ;
	Rcpp::NumericMatrix AMeanMatrix( numRow, numCol ) ;
	for(int i=0; i<numRow; i++)
	{
		for(int j=0; j<numCol; j++)
		{
			AMeanMatrix(i,j) = AMeanVector[i][j] ;
		}
	}
	
	numRow = AStdVector.size();
	numCol = AStdVector[0].size() ;
	Rcpp::NumericMatrix AStdMatrix( numRow, numCol ) ;
	for(int i=0; i<numRow; i++)
	{
		for(int j=0; j<numCol; j++)
		{
			AStdMatrix(i,j) = AStdVector[i][j] ;
		}
	}
	
	numRow = PMeanVector.size();
	numCol = PMeanVector[0].size() ;
	Rcpp::NumericMatrix PMeanMatrix( numRow, numCol ) ;
	for(int i=0; i<numRow; i++)
	{
		for(int j=0; j<numCol; j++)
		{
			PMeanMatrix(i,j) = PMeanVector[i][j] ;
		}
	}
	
	numRow = PStdVector.size();
	numCol = PStdVector[0].size() ;
	Rcpp::NumericMatrix PStdMatrix( numRow, numCol ) ;
	for(int i=0; i<numRow; i++)
	{
		for(int j=0; j<numCol; j++)
		{
			PStdMatrix(i,j) = PStdVector[i][j] ;
		}
	}
	
	/* Same Code for APMEAN Matrix, not needed
	numRow = APMeanVector.size();
	numCol = APMeanVector[0].size() ;
	Rcpp::NumericMatrix APMeanMatrix( numRow, numCol ) ;
	for(int i=0; i<numRow; i++)
	{
		for(int j=0; j<numCol; j++)
		{
			APMeanMatrix(i,j) = APMeanVector[i][j] ;
		}
	}
	*/


	
	Rcpp::List fileContainer =  Rcpp::List::create(Rcpp::Named("Amean") = AMeanMatrix, 
								Rcpp::Named("Asd") = AStdMatrix, Rcpp::Named("Pmean") = PMeanMatrix, Rcpp::Named("Psd") = PStdMatrix, 
								Rcpp::Named("atomsAEquil") = nAEquil, Rcpp::Named("atomsASamp") = nASamp, 
								Rcpp::Named("atomsPEquil") = nPEquil, Rcpp::Named("atomsPSamp") = nPSamp, Rcpp::Named("chiSqValues") = chiVect);
	
	//fileContainer[0] = AMeanMatrix;
	//fileContainer[1] = AStdMatrix;
	//fileContainer[2] = PMeanMatrix;
	//fileContainer[3] = PStdMatrix;
	//fileContainer[4] = APMeanMatrix;
	//fileContainer[5] = chiVect;
	 
  return (fileContainer);

}
