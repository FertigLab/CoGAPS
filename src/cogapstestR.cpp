// cogapstestR.cpp
// Outputs and saves the atomic domain, A and P matrices, and manually calculated
// chi squared values at each equilibration and sample. Allows for easier
// debugging.
// This code runs cogaps map and cogaps with fixed bin probabilities.
// Specify which in the config file.

// =============================================================================
// This is the main code for CogapsTest. (1th June, 2167)
// This file also works as an interface to R via Rcpp (14th Feb, 2015)
// =============================================================================

#include <iostream>       // for use with standard I/O
#include <string>         // for string processing
#include <fstream>        // for output to files
#include <limits>         // for extracting numerical limits of C++


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

using namespace std;
using namespace gaps;
using std::vector;

// [[Rcpp::export]]
Rcpp::List cogapsTest(Rcpp::DataFrame DFrame, Rcpp::DataFrame SFrame, Rcpp::DataFrame ABinsFrame,
                      Rcpp::DataFrame PBinsFrame, Rcpp::CharacterVector Config, Rcpp::NumericVector ConfigNums, int seed=-1) {
    // ===========================================================================
    // Initialization of the random number generator.
    // Different seeding methods:
    // --- fixed seed
    //std::vector<unsigned long> ve(2);
    //ve[0]=198782;ve[1]=89082;
    //boost::random::seed_seq seq(ve);
    //rng.seed(seq);
    // --- seeded with time
    if (seed <= 0) {
        rng.seed(static_cast<boost::uint32_t>(std::time(0)));
    } else {
        rng.seed(static_cast<boost::uint32_t>(seed));
    }
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
    double tempNumInput;
    tempNumInput = (ConfigNums[1]);
    unsigned long nEquil = (tempNumInput);
    tempNumInput = (ConfigNums[2]);
    unsigned long nSample = (tempNumInput);
    tempNumInput = (ConfigNums[3]);
    unsigned long nObsR = tempNumInput;
    tempNumInput = (ConfigNums[0]);
    unsigned int nFactor = tempNumInput;
    tempNumInput = (ConfigNums[4]);
    double alphaA = tempNumInput;
    tempNumInput = (ConfigNums[7]);
    double alphaP = tempNumInput;
    tempNumInput = (ConfigNums[5]);
    double nMaxA = tempNumInput;
    tempNumInput = (ConfigNums[8]);
    double nMaxP = tempNumInput;
    //---------------------------------------------------------------------------------
    string simulation_id = Rcpp::as<string>(Config[0]);
    tempNumInput = (ConfigNums[6]);
    double max_gibbsmass_paraA = tempNumInput;
    tempNumInput = (ConfigNums[9]);
    double max_gibbsmass_paraP = tempNumInput;
    temp = Rcpp::as<string>(Config[1]);
    bool Q_output_atomic;

    if (temp == "TRUE" || temp == "true") {
        Q_output_atomic = true;

    } else {
        Q_output_atomic = false;
    }

    temp = Rcpp::as<string>(Config[2]);
    bool fixBinProbs;

    if (temp == "TRUE" || temp == "true") {
        fixBinProbs = true;

    } else {
        fixBinProbs = false;
    }

    temp = Rcpp::as<string>(Config[3]);
    string fixedDomainStr = temp;
    //Code to make the D and S matrices read from R into C++ vectors to make into Elana's Matrix Objects in Matrix.cpp
    //Now also allows for variable bin sizes to be created
    vector<vector<double> > DVector;
    vector<vector<double> > SVector;
    vector<vector<double> > ABinsVector;
    vector<vector<double> > PBinsVector;
    //Code to establish the sizes and initialize the C++ vectors to pass
    int numC = DFrame.size();
    Rcpp::NumericVector tempFrameCol = DFrame[0];
    int numR = tempFrameCol.size() ;
    double tempFrameElement;
    DVector.resize(numR);

    for (int i = 0; i < numR; i++) {
        DVector[i].resize(numC);
    }

    for (int i = 0; i < numR; i++) {
        for (int j = 0; j < numC; j++) {
            tempFrameCol = DFrame[j];
            tempFrameElement = tempFrameCol[i];
            DVector[i][j] = tempFrameElement;
        }
    }

    numC = SFrame.size();
    tempFrameCol = SFrame[0];
    numR = tempFrameCol.size() ;
    SVector.resize(numR);

    for (int i = 0; i < numR; i++) {
        SVector[i].resize(numC);
    }

    for (int i = 0; i < numR; i++) {
        for (int j = 0; j < numC; j++) {
            tempFrameCol = SFrame[j];
            tempFrameElement = tempFrameCol[i];
            SVector[i][j] = tempFrameElement;
        }
    }

    //--------------------END CREATING D and S C++ VECTORS
    // Parameters or structures to be calculated or constructed:
    unsigned long nIterA = 10;    // initial inner loop iterations for A
    unsigned long nIterP = 10;    // initial inner loop iterations for P
    unsigned long atomicSize = 0; // number of atomic points
    // TestGAPS structures
    vector <map <unsigned long long, double> > EquilAAtomicDomains;
    // atomic domains of matrix A during equilibration
    vector <map <unsigned long long, double> > EquilPAtomicDomains;
    // atomic domains of matrix P during equilibration
    vector <map <unsigned long long, double> > SampleAAtomicDomains;
    // atomic domains of matrix A during sampling
    vector <map <unsigned long long, double> > SamplePAtomicDomains;
    // atomic domains of matrix P during sampling
    vector <vector <vector <double> > > AMatsbyEquil;
    // list of A matrices during equilibration (calculated
    // manually from list EquilAAtomicDomains)
    vector <vector <vector <double> > > PMatsbyEquil;
    // list of P matrices during equilibration (calculated
    // manually from list EquilPAtomicDomains)
    vector <vector <vector <double> > > AMatsbySample;
    // list of A matrices during sampling (calculated
    // manually from list SamplingAAtomicDomains)
    vector <vector <vector <double> > > PMatsbySample;
    // list of P matrices during sampling (calculated
    // manually from list SamplingAAtomicDomains)
    vector <double> EquilManChiSqs; // list of Chi squared values during equilibration
    // calculated manually from list of A and P
    // matrices during equilibration
    vector <double> SampleManChiSqs; // list of Chi squared values during sampling
    // calculated manually from list of A and P
    // matrices during sampling
    char label_A = 'A';  // label for matrix A
    char label_P = 'P';  // label for matrix P
    char label_D = 'D';  // label for matrix D
    char label_S = 'S';// label for matrix S
    GibbsSampler GibbsSamp(nEquil, nSample, nFactor, // construct GibbsSampler and
                           alphaA, alphaP, nMaxA, nMaxP, // Read in D and S matrices
                           nIterA, nIterP,
                           max_gibbsmass_paraA, max_gibbsmass_paraP,
                           atomicSize,
                           label_A, label_P, label_D, label_S,
                           DVector, SVector, simulation_id);
    // ---------------------------------------------------------------------------
    // Based on the information of D, construct and initialize for A and P both
    // the matrices and atomic spaces.
    GibbsSamp.init_AMatrix_and_PMatrix(); // initialize A and P matrices
    //This Section now is to handle the many possibilities for Variable Bin Sizes (Priors)
    //A for variable A bins, P for variable P Bins, B for both and N for regular uniform bin sizes
    char fixedDomain = fixedDomainStr[0];

    if (fixBinProbs) {
        Rcpp::Rcout << "Running CoGAPS with fixed bin probabilities with fixed domain ";
        Rcpp::Rcout << fixedDomain << endl;

        if (fixedDomain == 'B') {
            numC = ABinsFrame.size();
            tempFrameCol = ABinsFrame[0];
            numR = tempFrameCol.size();
            ABinsVector.resize(numR);

            for (int i = 0; i < numR; i++) {
                ABinsVector[i].resize(numC);
            }

            for (int i = 0; i < numR; i++) {
                for (int j = 0; j < numC; j++) {
                    tempFrameCol = ABinsFrame[j];
                    tempFrameElement = tempFrameCol[i];
                    ABinsVector[i][j] = tempFrameElement;
                }
            }

            numC = PBinsFrame.size();
            tempFrameCol = PBinsFrame[0];
            numR = tempFrameCol.size();
            PBinsVector.resize(numR);

            for (int i = 0; i < numR; i++) {
                PBinsVector[i].resize(numC);
            }

            for (int i = 0; i < numR; i++) {
                for (int j = 0; j < numC; j++) {
                    tempFrameCol = PBinsFrame[j];
                    tempFrameElement = tempFrameCol[i];
                    PBinsVector[i][j] = tempFrameElement;
                }
            }

            GibbsSamp.init_AAtomicdomain_and_PAtomicdomain(ABinsVector, PBinsVector);

        } else if (fixedDomain == 'A') {
            numC = ABinsFrame.size();
            tempFrameCol = ABinsFrame[0];
            numR = tempFrameCol.size() ;
            ABinsVector.resize(numR);

            for (int i = 0; i < numR; i++) {
                ABinsVector[i].resize(numC);
            }

            for (int i = 0; i < numR; i++) {
                for (int j = 0; j < numC; j++) {
                    tempFrameCol = ABinsFrame[j];
                    tempFrameElement = tempFrameCol[i];
                    ABinsVector[i][j] = tempFrameElement;
                }
            }

            GibbsSamp.init_AAtomicdomain_and_PAtomicdomain(fixedDomain, ABinsVector);

        } else if (fixedDomain == 'P') {
            numC = PBinsFrame.size();
            tempFrameCol = PBinsFrame[0];
            numR = tempFrameCol.size() ;
            PBinsVector.resize(numR);

            for (int i = 0; i < numR; i++) {
                PBinsVector[i].resize(numC);
            }

            for (int i = 0; i < numR; i++) {
                for (int j = 0; j < numC; j++) {
                    tempFrameCol = PBinsFrame[j];
                    tempFrameElement = tempFrameCol[i];
                    PBinsVector[i][j] = tempFrameElement;
                }
            }

            GibbsSamp.init_AAtomicdomain_and_PAtomicdomain(fixedDomain, PBinsVector);

        } else {
            throw logic_error("GibbsSampler: Invalid Specification of Fixed Domain(s)!");
        }

    } else {
        GibbsSamp.init_AAtomicdomain_and_PAtomicdomain();   // initialize atomic spaces
        // A and P
    }

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
    int numOutputs = nObsR;                  //CONOR'S OUTPUT CONTROL
    int totalChiSize = nSample + nEquil;
    Rcpp::NumericVector chiVect(totalChiSize); //INITIALIZE THE VECTOR HOLDING THE CHISQUARE.
    Rcpp::NumericVector nAEquil(nEquil);       //INITIALIZE THE VECTOR HOLDING THE ATOMS FOR EACH MATRIX FOR EACH EQUIL/SAMP
    Rcpp::NumericVector nASamp(nSample);
    Rcpp::NumericVector nPEquil(nEquil);
    Rcpp::NumericVector nPSamp(nSample);

    for (unsigned long ext_iter = 1; ext_iter <= nEquil; ++ext_iter) {
        GibbsSamp.set_iter(ext_iter);
        GibbsSamp.set_AnnealingTemperature();

        for (unsigned long iterA = 1; iterA <= nIterA; ++iterA) {
            GibbsSamp.update('A');
        }

        GibbsSamp.check_atomic_matrix_consistency('A');

        for (unsigned long iterP = 1; iterP <= nIterP; ++iterP) {
            GibbsSamp.update('P');
        }

        GibbsSamp.check_atomic_matrix_consistency('P');
        //Finds the current ChiSq and places it into the vector to be returned to R (and output on occasion)
        tempChiSq = GibbsSamp.get_sysChi2();
        chiVect[(ext_iter) - 1] = tempChiSq;
        // ----------- output computing info ---------
        tempAtomA = GibbsSamp.getTotNumAtoms('A');
        tempAtomP = GibbsSamp.getTotNumAtoms('P');
        nAEquil[outCount] = tempAtomA;
        nPEquil[outCount] = tempAtomP;
        outCount++;
        // Add the atomic domains to the list
        EquilAAtomicDomains.push_back(GibbsSamp.getAtomicDomain('A'));
        EquilPAtomicDomains.push_back(GibbsSamp.getAtomicDomain('P'));

        if (ext_iter % numOutputs == 0) {
            Rcpp::Rcout << "Equil:" << ext_iter << " of " << nEquil <<
                        ", Atoms:" << tempAtomA << "("
                        << tempAtomP << ")" <<
                        "  Chi2 = " << tempChiSq << endl;
        }

        // -------------------------------------------
        // re-calculate nIterA and nIterP to the expected number of atoms
        nIterA = (unsigned long) randgen('P', max((double) GibbsSamp.getTotNumAtoms('A'), 10.));
        nIterP = (unsigned long) randgen('P', max((double) GibbsSamp.getTotNumAtoms('P'), 10.));
        // --------------------------------------------
    }  // end of for-block for equilibration

// ===========================================================================
    // Part 2.5) Equilibration Settling:
    // Allow the Equilibration to settle for 10% of the set Equilibrations
    // ===========================================================================
    int nEquilCool = floor(.1 * nEquil);

    for (unsigned long ext_iter = 1; ext_iter <= nEquilCool; ++ext_iter) {
        GibbsSamp.set_iter(ext_iter);

        for (unsigned long iterA = 1; iterA <= nIterA; ++iterA) {
            GibbsSamp.update('A');
        }

        GibbsSamp.check_atomic_matrix_consistency('A');

        for (unsigned long iterP = 1; iterP <= nIterP; ++iterP) {
            GibbsSamp.update('P');
        }

        GibbsSamp.check_atomic_matrix_consistency('P');
    }

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

    for (unsigned long i = 1; i <= nSample; ++i) {
        for (unsigned long iterA = 1; iterA <= nIterA; ++iterA) {
            GibbsSamp.update('A');
        }

        GibbsSamp.check_atomic_matrix_consistency('A');

        for (unsigned long iterP = 1; iterP <= nIterP; ++iterP) {
            GibbsSamp.update('P');
        }

        GibbsSamp.check_atomic_matrix_consistency('P');

        if (Q_output_atomic == true) {
            GibbsSamp.output_atomicdomain('A', i);
            GibbsSamp.output_atomicdomain('P', i);
        }

        statindx += 1;
        GibbsSamp.compute_statistics_prepare_matrices(statindx);
        //Do the same as above.
        tempChiSq = GibbsSamp.get_sysChi2();
        chiVect[(nEquil + i) - 1] = tempChiSq;
        // ----------- output computing info ---------
        tempAtomA = GibbsSamp.getTotNumAtoms('A');
        tempAtomP = GibbsSamp.getTotNumAtoms('P');
        nASamp[outCount] = tempAtomA;
        nPSamp[outCount] = tempAtomP;
        outCount++;
        // Add the atomic domains to the list
        SampleAAtomicDomains.push_back(GibbsSamp.getAtomicDomain('A'));
        SamplePAtomicDomains.push_back(GibbsSamp.getAtomicDomain('P'));

        if (i % numOutputs == 0) {
            Rcpp::Rcout << "Samp: " << i << " of " << nSample <<
                        ", Atoms:" << tempAtomA << "("
                        << tempAtomP << ")" <<
                        "  Chi2 = " << tempChiSq << endl;

            if (i == nSample) {
                chi2 = 2.*GibbsSamp.cal_logLikelihood();
                Rcpp::Rcout << " *** Check value of final chi2: " << chi2 << " **** " << endl;
            }
        }

        // -------------------------------------------
        // re-calculate nIterA and nIterP to the expected number of atoms
        nIterA = (unsigned long) randgen('P', max((double) GibbsSamp.getTotNumAtoms('A'), 10.));
        nIterP = (unsigned long) randgen('P', max((double) GibbsSamp.getTotNumAtoms('P'), 10.));
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

    //Code to fill the chi-Squared vector
    for (int iEquil = 0; iEquil < nEquil; iEquil++) {
        // Create the matrices and add them to the list
        vector <vector <double> > EquilAMat = GibbsSamp.createSampleAMat(EquilAAtomicDomains.at(iEquil));
        vector <vector <double> > EquilPMat = GibbsSamp.createSamplePMat(EquilPAtomicDomains.at(iEquil));
        AMatsbyEquil.push_back(EquilAMat);
        PMatsbyEquil.push_back(EquilPMat);
        // Calculate the corresponding chi square value and add to the list
        EquilManChiSqs.push_back(GibbsSamp.ManualCalcChiSqu(EquilAMat, EquilPMat));
    }

    for (int iSample = 0; iSample < nSample; iSample++) {
        // Create the matrices and add them to the list
        vector <vector <double> > SampleAMat = GibbsSamp.createSampleAMat(SampleAAtomicDomains.at(iSample));
        vector <vector <double> > SamplePMat = GibbsSamp.createSamplePMat(SamplePAtomicDomains.at(iSample));
        AMatsbySample.push_back(SampleAMat);
        PMatsbySample.push_back(SamplePMat);
        // Calculate the corresponding chi squared value and add to the list
        SampleManChiSqs.push_back(GibbsSamp.ManualCalcChiSqu(SampleAMat, SamplePMat));
    }

    //CODE FOR CONVERTING ALL THE VECTORS FOR THE DATAFILES INTO NUMERIC VECTORS OR LISTS
    //Every Chi-Square for Equil and Sample
    Rcpp::NumericVector EquilManualChi(EquilManChiSqs.size());
    Rcpp::NumericVector SampleManualChi(SampleManChiSqs.size());
    EquilManualChi = EquilManChiSqs;
    SampleManualChi = SampleManChiSqs;
    int numRow = AMeanVector.size();
    int numCol = AMeanVector[0].size() ;
    Rcpp::NumericMatrix AMeanMatrix(numRow, numCol) ;

    for (int i = 0; i < numRow; i++) {
        for (int j = 0; j < numCol; j++) {
            AMeanMatrix(i, j) = AMeanVector[i][j] ;
        }
    }

    numRow = AStdVector.size();
    numCol = AStdVector[0].size() ;
    Rcpp::NumericMatrix AStdMatrix(numRow, numCol) ;

    for (int i = 0; i < numRow; i++) {
        for (int j = 0; j < numCol; j++) {
            AStdMatrix(i, j) = AStdVector[i][j] ;
        }
    }

    numRow = PMeanVector.size();
    numCol = PMeanVector[0].size() ;
    Rcpp::NumericMatrix PMeanMatrix(numRow, numCol) ;

    for (int i = 0; i < numRow; i++) {
        for (int j = 0; j < numCol; j++) {
            PMeanMatrix(i, j) = PMeanVector[i][j] ;
        }
    }

    numRow = PStdVector.size();
    numCol = PStdVector[0].size() ;
    Rcpp::NumericMatrix PStdMatrix(numRow, numCol) ;

    for (int i = 0; i < numRow; i++) {
        for (int j = 0; j < numCol; j++) {
            PStdMatrix(i, j) = PStdVector[i][j] ;
        }
    }

    //Code for returning all the Atomic Domains
    //First Create the lists for each equil or sample snapshot of the AtomicDomains
    Rcpp::List EquilAAtomicR(nEquil);
    Rcpp::List EquilPAtomicR(nEquil);
    Rcpp::List SampAAtomicR(nSample);
    Rcpp::List SampPAtomicR(nSample);
    //Now fill the lists with the atomic domain for each iteration
    //To do this make a temp matrix each time of the size of the atomic domain by 2
    int tempRowCount = 0;
    int tempAtomicSize = 0;
    map<unsigned long long, double>::const_iterator iter;

    for (int k = 0; k < nEquil; k++) {
        tempAtomicSize = EquilAAtomicDomains[k].size();
        Rcpp::NumericMatrix tempRDomain(tempAtomicSize, 2);
        tempRowCount = 0;

        for (iter = EquilAAtomicDomains[k].begin(); iter != EquilAAtomicDomains[k].end(); iter++) {
            tempRDomain(tempRowCount, 0) = (iter->first);
            tempRDomain(tempRowCount, 1) = (iter->second);
            tempRowCount++;
        }

        EquilAAtomicR[k] = tempRDomain;
    }

    for (int k = 0; k < nSample; k++) {
        tempAtomicSize = SampleAAtomicDomains[k].size();
        Rcpp::NumericMatrix tempRDomain(tempAtomicSize, 2);
        tempRowCount = 0;

        for (iter = SampleAAtomicDomains[k].begin(); iter != SampleAAtomicDomains[k].end(); iter++) {
            tempRDomain(tempRowCount, 0) = (iter->first);
            tempRDomain(tempRowCount, 1) = (iter->second);
            tempRowCount++;
        }

        SampAAtomicR[k] = tempRDomain;
    }

    for (int k = 0; k < nEquil; k++) {
        tempAtomicSize = EquilPAtomicDomains[k].size();
        Rcpp::NumericMatrix tempRDomain(tempAtomicSize, 2);
        tempRowCount = 0;

        for (iter = EquilPAtomicDomains[k].begin(); iter != EquilPAtomicDomains[k].end(); iter++) {
            tempRDomain(tempRowCount, 0) = (iter->first);
            tempRDomain(tempRowCount, 1) = (iter->second);
            tempRowCount++;
        }

        EquilPAtomicR[k] = tempRDomain;
    }

    for (int k = 0; k < nSample; k++) {
        tempAtomicSize = SamplePAtomicDomains[k].size();
        Rcpp::NumericMatrix tempRDomain(tempAtomicSize, 2);
        tempRowCount = 0;

        for (iter = SamplePAtomicDomains[k].begin(); iter != SamplePAtomicDomains[k].end(); iter++) {
            tempRDomain(tempRowCount, 0) = (iter->first);
            tempRDomain(tempRowCount, 1) = (iter->second);
            tempRowCount++;
        }

        SampPAtomicR[k] = tempRDomain;
    }

    //Code for returning all the Matrices
    //First Create the lists for each equil or sample snapshot of the matrices
    Rcpp::List AMatEquilR(nEquil);
    Rcpp::List PMatEquilR(nEquil);
    Rcpp::List AMatSampR(nSample);
    Rcpp::List PMatSampR(nSample);
    //Now fill the lists with the matrices for each iteration
    //To do this make a temp matrix every time of the dimensions of A and P, filled element by element
    numRow = AMeanVector.size();
    numCol = AMeanVector[0].size() ;
    Rcpp::NumericMatrix tempAEquilRMatrix(numRow, numCol);

    for (int k = 0; k < nEquil; k++) {
        for (int i = 0; i < numRow; i++) {
            for (int j = 0; j < numCol; j++) {
                tempAEquilRMatrix(i, j) = AMatsbyEquil[k][i][j] ;
            }
        }

        AMatEquilR[k] = tempAEquilRMatrix;
    }

    numRow = AMeanVector.size();
    numCol = AMeanVector[0].size();
    Rcpp::NumericMatrix tempASampRMatrix(numRow, numCol);

    for (int k = 0; k < nSample; k++) {
        for (int i = 0; i < numRow; i++) {
            for (int j = 0; j < numCol; j++) {
                tempASampRMatrix(i, j) = AMatsbySample[k][i][j] ;
            }
        }

        AMatSampR[k] = tempASampRMatrix;
    }

    numRow = PMeanVector.size();
    numCol = PMeanVector[0].size() ;
    Rcpp::NumericMatrix tempPEquilRMatrix(numRow, numCol);

    for (int k = 0; k < nEquil; k++) {
        for (int i = 0; i < numRow; i++) {
            for (int j = 0; j < numCol; j++) {
                tempPEquilRMatrix(i, j) = PMatsbyEquil[k][i][j] ;
            }
        }

        PMatEquilR[k] = tempPEquilRMatrix;
    }

    numRow = PMeanVector.size();
    numCol = PMeanVector[0].size() ;
    Rcpp::NumericMatrix tempPSampRMatrix(numRow, numCol);

    for (int k = 0; k < nSample; k++) {
        for (int i = 0; i < numRow; i++) {
            for (int j = 0; j < numCol; j++) {
                tempPSampRMatrix(i, j) = PMatsbySample[k][i][j] ;
            }
        }

        PMatSampR[k] = tempPSampRMatrix;
    }

    Rcpp::List fileContainer =  Rcpp::List::create(Rcpp::Named("Amean") = AMeanMatrix,
                                Rcpp::Named("Asd") = AStdMatrix, Rcpp::Named("Pmean") = PMeanMatrix, Rcpp::Named("Psd") = PStdMatrix,
                                Rcpp::Named("atomsAEquil") = nAEquil, Rcpp::Named("atomsASamp") = nASamp,
                                Rcpp::Named("atomsPEquil") = nPEquil, Rcpp::Named("atomsPSamp") = nPSamp, Rcpp::Named("chiSqValues") = chiVect,
                                Rcpp::Named("matricesAEquil") = AMatEquilR, Rcpp::Named("matricesASamp") = AMatSampR,
                                Rcpp::Named("matricesPEquil") = PMatEquilR, Rcpp::Named("matricesPSamp") = PMatSampR,
                                Rcpp::Named("domainAEquil") = EquilAAtomicR, Rcpp::Named("domainASamp") = SampAAtomicR,
                                Rcpp::Named("domainPEquil") = EquilPAtomicR, Rcpp::Named("domainPSamp") = SampPAtomicR,
                                Rcpp::Named("chiSqEquil") = EquilManualChi, Rcpp::Named("chiSqSamp") = SampleManualChi);
    return (fileContainer);
}
