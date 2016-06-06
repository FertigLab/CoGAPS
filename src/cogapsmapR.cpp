// cogapsmapR.cpp

// =============================================================================
// This is the main code for CogapsMap. (7th August, 2014)
// This file also works as an interface to R via Rcpp (1st June, 2015)
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
#include "GibbsSampler.h" // for incorporating the GibbsSampMapler which
// does all the atomic space to matrix conversion
// and sampling actions.
#include "GibbsSamplerMap.h"
#include <Rcpp.h>
// ------------------------------------------------------

using namespace std;
using namespace gaps;
using std::vector;

// [[Rcpp::export]]
Rcpp::List cogapsMap(Rcpp::DataFrame DFrame, Rcpp::DataFrame SFrame, Rcpp::DataFrame FixedPatt,
                     Rcpp::DataFrame ABinsFrame, Rcpp::DataFrame PBinsFrame, Rcpp::CharacterVector Config, Rcpp::NumericVector ConfigNums, int seed=-1) {
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
    temp = Rcpp::as<string>(Config[4]); //The matrix letter that has fixed mass distributing (map)
    string fixedMatrixStr = temp;
    temp = Rcpp::as<string>(Config[5]);
    bool SampleSnapshots;

    if (temp == "TRUE" || temp == "true") {
        SampleSnapshots = true;

    } else {
        SampleSnapshots = false;
    }

    tempNumInput = (ConfigNums[10]);
    int numSnapshots = tempNumInput;
    //Code to make the D and S matrices read from R into C++ vectors to make into Elana's Matrix Objects in Matrix.cpp
    //Now also allows for variable bin sizes to be created
    vector<vector<double> > DVector;
    vector<vector<double> > SVector;
    vector<vector<double> > FPVector;
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

    // Convert FixedPatt (DataFrame) to vector of vectors (double)
    numC = FixedPatt.size();
    tempFrameCol = FixedPatt[0];
    numR = tempFrameCol.size() ;
    FPVector.resize(numR);

    for (int i = 0; i < numR; i++) {
        FPVector[i].resize(numC);
    }

    for (int i = 0; i < numR; i++) {
        for (int j = 0; j < numC; j++) {
            tempFrameCol = FixedPatt[j];
            tempFrameElement = tempFrameCol[i];
            FPVector[i][j] = tempFrameElement;
        }
    }

    //--------------------END CREATING D and S C++ VECTORS
    // Parameters or structures to be calculated or constructed:
    unsigned long nIterA = 10;    // initial inner loop iterations for A
    unsigned long nIterP = 10;    // initial inner loop iterations for P
    unsigned long atomicSize = 0; // number of atomic points
    char label_A = 'A';  // label for matrix A
    char label_P = 'P';  // label for matrix P
    char label_D = 'D';  // label for matrix D
    char label_S = 'S';// label for matrix S
    char label_FP = fixedMatrixStr[0]; //Label for the Fixed Matrix
    // ---------------------------------------------------------------------------
    // Initialize the GibbsSampMap.
    //R Version
    //Now with Variable Bin capability for priors (Fixed Bins)
    GibbsSamplerMap GibbsSampMap(nEquil, nSample, nFactor, // construct GibbsSampMapler and
                                 alphaA, alphaP, nMaxA, nMaxP, // Read in D and S matrices
                                 nIterA, nIterP,
                                 max_gibbsmass_paraA, max_gibbsmass_paraP,
                                 atomicSize,
                                 label_A, label_P, label_D, label_S,
                                 DVector, SVector, simulation_id, FPVector, label_FP);
    // ---------------------------------------------------------------------------
    // Based on the information of D, construct and initialize for A and P both
    // the matrices and atomic spaces.
    GibbsSampMap.init_AMatrix_and_PMatrix(); // initialize A and P matrices
    //This Section now is to handle the many possibilities for Variable Bin Sizes (Priors)
    //A for variable A bins, P for variable P Bins, B for both and N for regular uniform bin sizes
    char fixedDomain = fixedDomainStr[0];

    //First check for a logic error, cannot fix the same domain as the matrix being fixed
    if (fixedDomain == label_FP) {
        throw logic_error("GibbsSampler: Cannot fix the same domain as the matrix being fixed!");
    }

    if (fixBinProbs) {
        Rcpp::Rcout << "Running CoGAPS Map with fixed bin probabilities with fixed domain ";
        Rcpp::Rcout << fixedDomain << endl;

        if (fixedDomain == 'A') {
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

            //Calls from child class GibbsSampler
            GibbsSampMap.init_AAtomicdomain_and_PAtomicdomain(fixedDomain, ABinsVector);

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

            GibbsSampMap.init_AAtomicdomain_and_PAtomicdomain(fixedDomain, PBinsVector);

        } else {
            throw logic_error("GibbsSampler: Invalid Specification of Fixed Domain!");
        }

    } else {
        GibbsSampMap.init_AAtomicdomain_and_PAtomicdomain();    // initialize atomic spaces
        // A and P
    }

    GibbsSampMap.init_Mapped_Matrix();
    GibbsSampMap.initialize_atomic_domain_map();
    GibbsSampMap.init_sysChi2(); // initialize the system chi2 value
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
    int numOutputs = nObsR;
    int totalChiSize = nSample + nEquil;
    Rcpp::NumericVector chiVect(totalChiSize); //INITIALIZE THE VECTOR HOLDING THE CHISQUARE.
    Rcpp::NumericVector nAEquil(nEquil);       //INITIALIZE THE VECTOR HOLDING THE ATOMS FOR EACH MATRIX FOR EACH EQUIL/SAMP
    Rcpp::NumericVector nASamp(nSample);
    Rcpp::NumericVector nPEquil(nEquil);
    Rcpp::NumericVector nPSamp(nSample);

    for (unsigned long ext_iter = 1; ext_iter <= nEquil; ++ext_iter) {
        GibbsSampMap.set_iter(ext_iter);
        GibbsSampMap.set_AnnealingTemperature();

        for (unsigned long iterA = 1; iterA <= nIterA; ++iterA) {
            GibbsSampMap.mapUpdate('A');
        }

        GibbsSampMap.check_atomic_matrix_consistency('A');

        for (unsigned long iterP = 1; iterP <= nIterP; ++iterP) {
            GibbsSampMap.mapUpdate('P');
        }

        GibbsSampMap.check_atomic_matrix_consistency('P');
        //Finds the current ChiSq and places it into the vector to be returned to R (and output on occasion)
        tempChiSq = GibbsSampMap.get_sysChi2();
        chiVect[(ext_iter) - 1] = tempChiSq;
        // ----------- output computing info ---------
        tempAtomA = GibbsSampMap.getTotNumAtoms('A');
        tempAtomP = GibbsSampMap.getTotNumAtoms('P');
        nAEquil[outCount] = tempAtomA;
        nPEquil[outCount] = tempAtomP;
        outCount++;

        if (ext_iter % numOutputs == 0) {
            Rcpp::Rcout << "Equil:" << ext_iter << " of " << nEquil <<
                        ", Atoms:" << tempAtomA << "("
                        << tempAtomP << ")" <<
                        "  Chi2 = " << tempChiSq << endl;
        }

        // -------------------------------------------
        // re-calculate nIterA and nIterP to the expected number of atoms
        nIterA = (unsigned long) randgen('P', max((double) GibbsSampMap.getTotNumAtoms('A'), 10.));
        nIterP = (unsigned long) randgen('P', max((double) GibbsSampMap.getTotNumAtoms('P'), 10.));
        // --------------------------------------------
    }  // end of for-block for equilibration

// ===========================================================================
    // Part 2.5) Equilibration Settling:
    // Allow the Equilibration to settle for 10% of the set Equilibrations
    // ===========================================================================
    int nEquilCool = floor(.1 * nEquil);

    for (unsigned long ext_iter = 1; ext_iter <= nEquilCool; ++ext_iter) {
        GibbsSampMap.set_iter(ext_iter);

        for (unsigned long iterA = 1; iterA <= nIterA; ++iterA) {
            GibbsSampMap.mapUpdate('A');
        }

        GibbsSampMap.check_atomic_matrix_consistency('A');

        for (unsigned long iterP = 1; iterP <= nIterP; ++iterP) {
            GibbsSampMap.mapUpdate('P');
        }

        GibbsSampMap.check_atomic_matrix_consistency('P');
    }

    // ===========================================================================
    // Part 3) Sampling:
    // After the system equilibriates in Part 2, we sample the systems with an
    // outer loop of nSample iterations. Within each outer loop iteration, A is
    // iterated nIterA times and P is iterated nIterP times. After sampling,
    // we update nIterA and nIterP according to the expected number of atoms in
    // the atomic spaces of A and P respectively.
    // ===========================================================================
    //Initialize Snapshots of A and P
    vector <vector <vector <double> > > ASnap;
    vector <vector <vector <double> > > PSnap;
    unsigned int statindx = 0;
    outCount = 0;

    for (unsigned long i = 1; i <= nSample; ++i) {
        for (unsigned long iterA = 1; iterA <= nIterA; ++iterA) {
            GibbsSampMap.mapUpdate('A');
        }

        GibbsSampMap.check_atomic_matrix_consistency('A');

        for (unsigned long iterP = 1; iterP <= nIterP; ++iterP) {
            GibbsSampMap.mapUpdate('P');
        }

        GibbsSampMap.check_atomic_matrix_consistency('P');

        if (Q_output_atomic == true) {
            GibbsSampMap.output_atomicdomain('A', i);
            GibbsSampMap.output_atomicdomain('P', i);
        }

        statindx += 1;
        GibbsSampMap.compute_statistics_prepare_matrices(statindx);
        //Do the same as above.
        tempChiSq = GibbsSampMap.get_sysChi2();
        chiVect[(nEquil + i) - 1] = tempChiSq;
        // ----------- output computing info ---------
        tempAtomA = GibbsSampMap.getTotNumAtoms('A');
        tempAtomP = GibbsSampMap.getTotNumAtoms('P');
        nASamp[outCount] = tempAtomA;
        nPSamp[outCount] = tempAtomP;
        outCount++;

        if (i % numOutputs == 0) {
            //----------------------------------
            Rcpp::Rcout << "Samp: " << i << " of " << nSample <<
                        ", Atoms:" << tempAtomA << "("
                        << tempAtomP << ")" <<
                        // " ,chi2 = " << chi2 <<
                        "  Chi2 = " << tempChiSq << endl;

            if (i == nSample) {
                chi2 = 2.*GibbsSampMap.cal_logLikelihood();
                Rcpp::Rcout << " *** Check value of final chi2: " << chi2 << " **** " << endl;
            }
        }

        if (SampleSnapshots && (i % (nSample / numSnapshots) == 0)) {
            vector <vector <vector <double> > > NormedMats = GibbsSampMap.getNormedMatrices();
            ASnap.push_back(NormedMats[0]);
            PSnap.push_back(NormedMats[1]);
        }

        // -------------------------------------------
        // re-calculate nIterA and nIterP to the expected number of atoms
        nIterA = (unsigned long) randgen('P', max((double) GibbsSampMap.getTotNumAtoms('A'), 10.));
        nIterP = (unsigned long) randgen('P', max((double) GibbsSampMap.getTotNumAtoms('P'), 10.));
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
    GibbsSampMap.compute_statistics(statindx,
                                    AMeanVector, AStdVector, PMeanVector, PStdVector);          // compute statistics like mean and s.d.
    //CODE FOR CONVERTING ALL THE VECTORS FOR THE DATAFILES INTO NUMERIC VECTORS OR LISTS
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

    //Code for transferring Snapshots in R
    int numSnaps = numSnapshots; //Arbitrary to keep convention

    if (SampleSnapshots == true) {
        Rcpp::List ASnapR(numSnaps);
        Rcpp::List PSnapR(numSnaps);
        numRow = AMeanVector.size();
        numCol = AMeanVector[0].size() ;
        Rcpp::NumericMatrix tempASnapMatrix(numRow, numCol);

        for (int k = 0; k < numSnaps; k++) {
            for (int i = 0; i < numRow; i++) {
                for (int j = 0; j < numCol; j++) {
                    tempASnapMatrix(i, j) = ASnap[k][i][j] ;
                }
            }

            ASnapR[k] = (tempASnapMatrix);
        }

        numRow = PMeanVector.size();
        numCol = PMeanVector[0].size() ;
        Rcpp::NumericMatrix tempPSnapMatrix(numRow, numCol);

        for (int k = 0; k < numSnaps; k++) {
            for (int i = 0; i < numRow; i++) {
                for (int j = 0; j < numCol; j++) {
                    tempPSnapMatrix(i, j) = PSnap[k][i][j] ;
                }
            }

            PSnapR[k] = (tempPSnapMatrix);
        }

        Rcpp::List fileContainer =  Rcpp::List::create(Rcpp::Named("Amean") = AMeanMatrix,
                                    Rcpp::Named("Asd") = AStdMatrix, Rcpp::Named("Pmean") = PMeanMatrix, Rcpp::Named("Psd") = PStdMatrix,
                                    Rcpp::Named("ASnapshots") = ASnapR, Rcpp::Named("PSnapshots") = PSnapR,
                                    Rcpp::Named("atomsAEquil") = nAEquil, Rcpp::Named("atomsASamp") = nASamp,
                                    Rcpp::Named("atomsPEquil") = nPEquil, Rcpp::Named("atomsPSamp") = nPSamp, Rcpp::Named("chiSqValues") = chiVect);
        return (fileContainer);

    } else {
        //Just leave the snapshots as empty lists
        Rcpp::List ASnapR = Rcpp::List::create();
        Rcpp::List PSnapR = Rcpp::List::create();
        Rcpp::List fileContainer =  Rcpp::List::create(Rcpp::Named("Amean") = AMeanMatrix,
                                    Rcpp::Named("Asd") = AStdMatrix, Rcpp::Named("Pmean") = PMeanMatrix, Rcpp::Named("Psd") = PStdMatrix,
                                    Rcpp::Named("ASnapshots") = ASnapR, Rcpp::Named("PSnapshots") = PSnapR,
                                    Rcpp::Named("atomsAEquil") = nAEquil, Rcpp::Named("atomsASamp") = nASamp,
                                    Rcpp::Named("atomsPEquil") = nPEquil, Rcpp::Named("atomsPSamp") = nPSamp, Rcpp::Named("chiSqValues") = chiVect);
        return (fileContainer);
    }
}
