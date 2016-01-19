#ifndef _GIBBSSAMPLER_H_
#define _GIBBSSAMPLER_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include "GAPSNorm.h"
#include "randgen.h"
#include "sub_func.h"
#include "Matrix.h"
#include "AtomicSupport.h"
#include <limits>
#include <cmath>
#include <limits>
#include <stdexcept>

using namespace gaps;
using namespace std;
using std::vector;

// -----------------------------------------------------------------------------
const double DOUBLE_POSINF = std::numeric_limits<double>::max();
const double DOUBLE_NEGINF = -std::numeric_limits<double>::max();
const double epsilon = 1e-10;
// -----------------------------------------------------------------------------

class GibbsSampler {
  protected:
    // Parameters or data that are read in:
    unsigned long _nEquil;    // # outer loop iterations for equilibration
    unsigned long _nSample;   // # outer loop iterations for sampling
    unsigned int _nFactor;    // # patterns
    double _alphaA;
    double _alphaP;
    double _nMaxA;             // max. number of atoms in A
    double _nMaxP;             // number of atomic bins for P
    unsigned long _nIterA;    // initial # of inner-loop iterations for A
    unsigned long _nIterP ;    // initial # of inner-loop iterations for P
    string _simulation_id;   // simulation id for the run
    double _max_gibbsmass_paraA; // max gibbs mass parameter for A
    double _max_gibbsmass_paraP; // max gibbs mass parameter for P
    double _lambdaA_scale_factor; // factor to rescale _lambdaA
    double _lambdaP_scale_factor; // factor to rescale _lambdaP

    // Parameters or structures to be calculated or constructed:
    unsigned int _nRow;       // number of items in observation (= # of genes)
    unsigned int _nCol;       // number of observation (= # of arrays)
    unsigned int _nBinsA;     // number of atomic bins for A
    unsigned int _nBinsP;     // number of atomic bins for P
    double _lambdaA;
    double _lambdaP;
    double _max_gibbsmassA;  // max gibbs mass for A
    double _max_gibbsmassP;  // max gibbs mass for P
    unsigned long _atomicSize; // number of atomic points

    char _label_A;  // label for matrix A
    char _label_P;  // label for matrix P
    char _label_D;  // label for matrix D
    char _label_S;  // label for matrix S

    unsigned long _iter;
    double _annealingTemperature;

    AtomicSupport _AAtomicdomain, _PAtomicdomain;
    Matrix _AMatrix, _PMatrix, _DMatrix, _SMatrix;

    map<unsigned long long, double> _atomicProposal;
    unsigned int _nChange_atomicProposal;
    char _oper_type; // label for the operation type
    double _sysChi2; // the chi2 value for the system

    unsigned int _nChange_matrixElemChange;
    vector<unsigned int> _Row_changed;
    vector<unsigned int> _Col_changed;
    vector<double> _mass_changed;
    vector<boost::tuple<unsigned int, unsigned int, double> > _matrixElemChange;

    map<unsigned long long, double> _new_atomicProposal;
    unsigned int _new_nChange_atomicProposal;

    unsigned int _new_nChange_matrixElemChange;
    vector<unsigned int> _new_Row_changed;
    vector<unsigned int> _new_Col_changed;
    vector<double> _new_mass_changed;
    vector<boost::tuple<unsigned int, unsigned int, double> > _new_matrixElemChange;



    // for computing statistics with matrices A and P
    // unsigned long _statindx_A, _statindx_P;  // counter
    double **_Amean;
    double **_Asd;
    double **_Pmean;
    double **_Psd;

  public:

    // ******************** CONSTRUCTOR ********************************************
    GibbsSampler() {};

    GibbsSampler(unsigned long nEquil, unsigned long nSample, unsigned int nFactor,
                 double alphaA, double alphaP, double nMaxA, double nMaxP,
                 unsigned long nIterA, unsigned long nIterP,
                 double max_gibbsmass_paraA, double max_gibbsmass_paraP,
                 unsigned long long atomicSize,
                 char label_A, char label_P, char label_D, char label_S,
                 vector<vector<double> > &DVector, vector<vector<double> > &SVector,
                 const string &simulation_id);

    ~GibbsSampler() {};


    // *************** METHOS FOR INITIALIZATION, DISPALY, OUTPUT ***********************
    void init_AMatrix_and_PMatrix();

    void init_AAtomicdomain_and_PAtomicdomain();

    void init_AAtomicdomain_and_PAtomicdomain(char fixeddomain, vector<vector<double> > ReadBinProbs);

    void init_AAtomicdomain_and_PAtomicdomain(vector<vector<double> > ReadBinProbsA, vector<vector<double> > ReadBinProbsP);

    void clear_Proposal();

    void clear_new_Proposal();

    void local_display_matrix2F(ofstream &outputFile, double **Mat_ptr,
                                unsigned int n_row, unsigned int n_col);

    void output_atomicdomain(char atomic_label, unsigned long Samp_cycle);

    // ********* METHODS TO GO BETWEEN ATOMIC SPACE AND MATRIX *****************

    unsigned int getRow(char matrix_label , unsigned int iBin);

    unsigned int getCol(char matrix_label , unsigned int iBin);

    unsigned int getTotNumAtoms(char matrix_label);


    // Snapshot methods
    vector <vector <vector <double> > > getNormedMatrices();

    // COGAPS TEST Methods
    map <unsigned long long, double> getAtomicDomain(char matrix_label);

    vector <vector <double> > createSampleAMat(map <unsigned long long, double> ADomain);

    vector <vector <double> > createSamplePMat(map <unsigned long long, double> PDomain);

    double ManualCalcChiSqu(vector <vector <double> > SampleAMat, vector <vector <double> > SamplePMat);

    // **************** METHODS FOR COMPUTING LIKELIHOOD FUNCTIONS *****************
    double cal_logLikelihood();

    /**
      * @short Compute the change in likelihood, DeltaLL, by first getting the proposal from
      *  the atomic space, then it invokes the corresponding methods in GAPSNorm
      *  according to the proposal size.
    * @return the change
    */
    double computeDeltaLL(char the_matrix_label,
                          double const *const *D,
                          double const *const *S,
                          double const *const *A,
                          double const *const *P,
                          unsigned int the_nChange_matrixElemChange,
                          const vector<boost::tuple<unsigned int, unsigned int, double> > the_matrixElemChange);
    /**
      * @short Extract information of the proposal made in the atomic space.
      *  Assuming an _atomicProposal, this method instantiates the
      *  corresponding member variables for _matrixElemChange. In the current version,
      *  we store the new proposal (in matrix space) in two different class variables
      *  of GibbsSampler:
      *  1. vector<unsigned int> _Row_changed; vector<unsigned int> _Col_changed;
      *     vector<double> _mass_changed;
      *  2. _matrixElemChange (this is a boost::tuple)
      */
    void extract_atomicProposal(char the_matrix_label);

    /**
      * @short Extract from _new_atomicProposal
      *  The code is exactly the same as extract_atomicProposal, only using the
      *  _new-quantities for the new proposal.
      */
    void extract_new_atomicProposal(char the_matrix_label);

    /** @short make proposal and update for the matrices
      *  It output "newMatrix", which is the final product of all the steps like
      *  birth, death, move, or exchange. One of these four methods is chosen
      *  based on the oper_type extracted from AtomicSupport using get_oper_type.
      */
    void update(char the_matrix_label);

    void init_sysChi2();

    void update_sysChi2(double delsysChi2);

    double get_sysChi2();

    void get_oper_type(char the_matrix_label);

    // Formerly together as birth_death and move_exchange,
    // these methods have been separated in this version.
    // All of these methods follow a similar pattern:
    /** @short instantiate the matrix element change to
      *  be made, if it can be.
    * @return true if the matrix needs changing, false
    *  if not
    */

    bool death(char the_matrix_label,
               double const *const *D,
               double const *const *S,
               double **AOrig,
               double **POrig);

    bool birth(char the_matrix_label,
               double const *const *D,
               double const *const *S,
               double **AOrig,
               double **POrig);

    bool move(char the_matrix_label,
              double const *const *D,
              double const *const *S,
              double **AOrig,
              double **POrig);

    bool exchange(char the_matrix_label,
                  double const *const *D,
                  double const *const *S,
                  double **AOrig,
                  double **POrig);


    // ************ METHODS FOR LOOPING AND CONTROL ********************************
    void set_iter(unsigned long ext_iter);

    /**
      * @short Calculate the _new annealingTemperature as a function of iteration step _iter.
      * (Note: the annealingTemperature here is really the inverted temperature!)
    */
    void set_AnnealingTemperature();

    void check_atomic_matrix_consistency(char the_matrix_label);

    /**
      * @short Form the matrices _Amean, _Asd, _Pmean, _Psd. They are temporary
      *  matrices to be modified in compute_statistics() to calculate the means and the variances
      *  of individual matrix entries. Note that here we have cumulated to each matrix element all its
      *  realizations in time as:
      *  _Amean_{ij} = \sum_{t} A_{ij}*k_{j}
      *  _Asd_{ij} = \sum_{t}  (A_{ij}*k_{j} )^2
      *  _Pmean_{ij} = \sum_{t} P_{ij}/k_{i}
      *  _Psd_{ij} = \sum_{t}  (P_{ij}/k_{i} )^2
    */
    void compute_statistics_prepare_matrices(unsigned long statindx);

    /** @short Use the matrices prepared in compute_statistics_prepare_matrices()
      * to compute the means and variances of individual elements in A and P.
    */

    void compute_statistics(unsigned int Nstat,
                            vector< vector <double> > &AMeanVect,
                            vector< vector <double> > &AStdVect,
                            vector< vector <double> > &PMeanVect,
                            vector< vector <double> > &PStdVect
                           );

    // ---------------------------------------------------------------------------
    /**
     * @short check whether or not we can use Gibbs Sampling by looking at the other
        matrix. History: formerly two methods, performUpdate and performUpdateKill,
      now consolidated to one method called checkOtherMatrix.
     */
    bool checkOtherMatrix(char the_matrix_label, unsigned int iRow, unsigned int iCol,
                          double const *const *otherMatrix);

    /**
     * @short Give the atom a mass based on Gibbs Sampling, if it can be used
     * @return the new mass of the atom
     */
    double getMass(char the_matrix_label, double origMass,
                   unsigned int iRow,
                   unsigned int iCol,
                   double const *const *otherMatrix,
                   double const *const *currentChainMatrix,
                   double const *const *D, double const *const *S,
                   double rng);
};
#endif
