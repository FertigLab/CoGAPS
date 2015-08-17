#ifndef ATOMIC_SUPPORT_H_
#define ATOMIC_SUPPORT_H_

#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "sub_func.h"
#include "randgen.h"



using std::map;
using namespace std;

namespace gaps {

  class AtomicSupport {
    
  public:
    AtomicSupport();
    ~AtomicSupport();
    
    /**
	 * @short Initialize atoms into the atomic domain -- 
	    adding ONE atom to the location [loc] with mass (mass)
	 */
	void initializeAtomic(unsigned int nBin, unsigned long long NatomLength, 
			  double alpha, double lambda, char atomic_domain_label);
     	
    void FixedBins_initializeAtomic(unsigned int nBin, unsigned long long NatomLength, 
			  double alpha, double lambda, char atomic_domain_label,
			  const char input_file_name[]);
	
	void FixedBins_initializeAtomic(unsigned int nBin, unsigned long long NatomLength, 
			  double alpha, double lambda, char atomic_domain_label,
			  std::vector<std::vector<double> > ReadBinProbs);

    /**
     * @short Find the bin to which the given location refers.
     * @return Bin number.
     */
    unsigned int getBin(unsigned long long location);

    unsigned long long getMidLocation(unsigned int iBin);
    unsigned long long getStartLocation(unsigned int iBin);
    unsigned long long getEndLocation(unsigned int iBin);

    unsigned int getExpectedNAtom() {
      return _alpha * _nBin;
    }

    double getNAtomPriorProb(int delAtom, bool log);

    unsigned long long binToLocation(unsigned int bin);

    double getTotalMass(unsigned int bin);
	
    unsigned int getNumAtoms(unsigned int bin);
	
	map <unsigned long long, double> getDomain();
	
    double get_atomicDomain_totalmass();

    void setInitialAtoms(const map<unsigned long long, double> initAtoms);

    double getLambda();

    /**
     * @short Propose a change to the atomic domain
	   The first method to be called in the Gibbs Sampler. 
       Decide which propose method to call based on value of 
	   double updatestep, which is based on # of atoms and 
	   distribution on number of atoms. Propose methods are
	   below makeProposal. 
    */
    void makeProposal(double rng);
	
	void ProposeBirth();
	void ProposeDeath();
	void ProposeMove();
	void ProposeExchange();


    map<unsigned long long, double> getProposedAtoms() {
      return _proposedAtoms;
    }

    map<unsigned long long, double>::const_iterator getAtomsBegin() {
      return _AtomicDomain.begin();
    }

    map<unsigned long long, double>::const_iterator getAtomsEnd() {
      return _AtomicDomain.end();
    }

    unsigned int getNAtom() {
      return _nAtom;
    }

    bool inDomain(unsigned long long location);
    double getMass(unsigned long long location);
	
	unsigned long long birthAtomLocation();

    void setProposedAtomMass(const map<unsigned long long, double> newProposal,
			     bool isNewProposal);

    void acceptProposal(bool updateIter);

    void rejectProposal(bool updateIter);

    void writeAtomicDiagnostics();

    void writeAtomicDiagnosticsBinary(bool isByCol, bool outByCol,
				      unsigned int nRow, unsigned int nCol);

    void initializeAtomicBinary(char diagnosticFileName[]);

    void writeAtomicHeader(char diagnosticFileName[],
			   int thinDiag);

    void resetAtomicOutputThin(int thinDiag);

    unsigned int getNBin() { 
      return _nBin;
    }

    void updateAtomicBins(double binProbabilities[], unsigned int length,
			  bool onlyUpdateRelativeWidth);

    void printAtomicInfo(); // original
    void printAtomicInfoF(ofstream& outputFile);

    void writeAtomicInfo(char outputFilename[],unsigned long Samp_cycle);

    char get_atomic_domain_label();

    char get_oper_type();


  private:
    // storage of the atomic domain
    map<unsigned long long, double> _AtomicDomain;
    unsigned long long _nAtom;
    int _iter;

    // boundaries of the atomic domain
    map<unsigned int, unsigned long long> _lBoundariesByBin;
    map<unsigned long long, unsigned int> _lBoundaries;

    // proposed changes to the atomic domain
    map<unsigned long long, double> _proposedAtoms;

    // deletion functions that properly clean up the memory
    void cleanDeleteProposal(map<unsigned long long, double>::const_iterator iter);
    void cleanDeleteProposalLocation(unsigned long long location);
    void cleanDeleteAtomic(map<unsigned long long, double>::const_iterator iter);
    void cleanDeleteAtomicLocation(unsigned long long location);
    void cleanClearProposal();
    void cleanClearAtomic();

    // parameters of the distribution
    unsigned int _nBin; // number of bins into which the distribution is divided
    unsigned long long _NatomLength;   // maximum number of slots for atoms
    double _alpha;      // average number of atoms per bin
    double _lambda;     // expected magnitude of each atom
    char _atomic_domain_label;  // label of the atomic domain
    char _oper_type; // the type of operation in makeProposal
    double _epsilon; // small number for setting things to zero

    // output file with diagnostic atomic information
    ofstream atomicDiagFile;
    ofstream atomicDiagFileBinary;
    bool     outputAtomicDiag;
    int      thinAtomicDiag;
    int      _initIterOutput;

  };
}

#endif /* ATOMIC_SUPPORT_H_ */
