#include "AtomicSupport.h"

#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <string.h>
#include <time.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <limits>

using std::vector;
using std::copy;
using std::min;
using std::string;
using std::logic_error;
using std::pair;
using std::map;
using std::endl;

const double DOUB_POSINF = std::numeric_limits<double>::max();
const double DOUB_NEGINF = -std::numeric_limits<double>::max();
double lepsilon = 1.e-10;

namespace gaps {
    
    AtomicSupport::AtomicSupport() {
    }
    
    AtomicSupport::~AtomicSupport() {
    }
    
    
    void AtomicSupport::initializeAtomic(unsigned int nBin,
                                         unsigned long long NatomLength,
                                         double alpha, double lambda,
                                         char atomic_domain_label) {
        // parameter values
        _nBin      = nBin;
        _NatomLength = NatomLength;
        _alpha     = alpha;
        _lambda    = lambda;
        _atomic_domain_label = atomic_domain_label;
        
        // set the small number
        _epsilon = lepsilon;
        
        // check the validity of the parameters
        if (NatomLength <= 0)
            throw logic_error("Invalid value for maximum number of atoms in AtomicSupport.");
        
        if (lambda <= 0)
            throw logic_error("Invalid value for lambda in AtomicSupport.");
        
        if (nBin <= 0)
            throw logic_error("Invalid number of bins in AtomicSupport.");
        
        // initialize the internal atomic domain
        _nAtom = 0;
        cleanClearAtomic();
        cleanClearProposal();
        
        // initialize the diagnostic file
        _iter = 0;
        
        // by default, set boundaries to be uniform -- TODO, check this
        double binProb[nBin];
        for (unsigned int iBin = 0; iBin < nBin; iBin++) {
            binProb[iBin] = 1.;
        }
        
        updateAtomicBins(binProb, nBin, true);
    }
    
    void AtomicSupport::printAtomicInfo() {
        /*cout << endl;
         cout << "nBin: "<<_nBin<<endl;
         cout << "NatomLength: "<< _NatomLength<<endl;
         cout << "lambda: "<<_lambda<<endl;
         cout << "alpha: "<<_alpha<<endl;
         cout << "atomic_domain_label: " << _atomic_domain_label << endl;
         cout << "atomDomain = [";
         
         unsigned int iBin;
         
         for (iBin = 0; iBin < _nBin; iBin++) {
         cout << "[" << iBin << " " << _lBoundariesByBin[iBin]<<" "<<getNumAtoms(iBin)<<" "<<getTotalMass(iBin) << "];";
         }
         cout << "];" << endl;
         cout << endl;*/
    }
    
    
    void AtomicSupport::printAtomicInfoF(ofstream& outputFile) {
        outputFile << endl;
        outputFile << "nBin: "<<_nBin<<endl;
        outputFile << "NatomLength: "<< _NatomLength<<endl;
        outputFile << "lambda: "<<_lambda<<endl;
        outputFile << "alpha: "<<_alpha<<endl;
        outputFile << "atomic_domain_label: " << _atomic_domain_label << endl;
        
        
        outputFile << "atomDomain = [";
        
        unsigned int iBin;
        
        for (iBin = 0; iBin < _nBin; iBin++) {
            outputFile << "[" << iBin << " " << _lBoundariesByBin[iBin]<<" "<<getNumAtoms(iBin)<<" "<<getTotalMass(iBin) << "];";
            
        }
        outputFile << "];" << endl;
        outputFile << endl;
    }
    
    void AtomicSupport::writeAtomicInfo(char outputFilename[],
                                        unsigned long Samp_cycle) {
        
        ofstream outputFile;
        if (Samp_cycle == 1){
            outputFile.open(outputFilename,ios::out);
        } else {
            outputFile.open(outputFilename,ios::out|ios::app);
        }
        
        map<unsigned long long, double>::const_iterator iter;
        for (iter = _AtomicDomain.begin(); iter != _AtomicDomain.end();
             iter++) {
            outputFile << setiosflags(ios::right) << setw(25) << iter->first << " "
            << setw(15) << iter->second << endl;
        }
        
        outputFile.close();
    }
    
    void AtomicSupport::setInitialAtoms(const map<unsigned long long, double>
                                        initAtoms) {
        for (map<unsigned long long, double>::const_iterator iter = initAtoms.begin();
             iter != initAtoms.end(); iter++) {
            
            if (iter->second > _epsilon) {
                _nAtom++;
                _AtomicDomain.insert(pair<unsigned long long, double>(iter->first,
                                                                      iter->second));
            }
        }
    }
    
    
    
    double AtomicSupport::getTotalMass(unsigned int iBin) {
        double massInBin = 0;
        
        if (_AtomicDomain.size() == 0) return 0;
        
        if (iBin >= _nBin) {
            throw logic_error("Cannot determine mass for more bins than in atomic domain in AtomicSupport::getTotalMass");
        }
        
        map<unsigned int, unsigned long long>::const_iterator boundIter =
        _lBoundariesByBin.find(iBin);
        
        unsigned long long lBound = boundIter->second;
        unsigned long long rBound;
        if (boundIter != _lBoundariesByBin.end()) {
            if (boundIter->first != _lBoundariesByBin.rbegin()->first) {
                boundIter++;
                rBound = boundIter->second;
            } else {
                rBound = _NatomLength;
            }
        } else {
            // bin is of size zero -> mass is zero
            return 0.;
        }
        
        map<unsigned long long, double>::const_iterator iter;
        
        
        for (iter = _AtomicDomain.lower_bound(lBound); iter != _AtomicDomain.end();
             iter++) {
            
            if (iter->first < rBound) {
                
                massInBin += iter->second;
                
            } else {
                break;
            }
        }
        
        return massInBin;
    }
    
    unsigned int AtomicSupport::getNumAtoms(unsigned int iBin) {
        unsigned int nAtomInBin = 0;
        
        if (_AtomicDomain.size() == 0) return 0;
        
        if (iBin >= _nBin) {
            throw logic_error("Cannot determine mass for more bins than in atomic domain in AtomicSupport::getTotalMass");
        }
        
        map<unsigned int, unsigned long long>::const_iterator boundIter =
        _lBoundariesByBin.find(iBin);
        
        unsigned long long lBound = boundIter->second;
        unsigned long long rBound;
        if (boundIter != _lBoundariesByBin.end()) {
            if (boundIter->first != _lBoundariesByBin.rbegin()->first) {
                boundIter++;
                rBound = boundIter->second;
            } else {
                rBound = _NatomLength;
            }
        } else {
            // bin is of size zero -> mass is zero
            return 0;
        }
        
        map<unsigned long long, double>::const_iterator iter;
        
        
        for (iter = _AtomicDomain.lower_bound(lBound); iter != _AtomicDomain.end();
             iter++) {
            
            
            if (iter->first < rBound) {
                
                nAtomInBin++;
                
            } else {
                break;
            }
        }
        
        return nAtomInBin;
    }
    
    
    
    double AtomicSupport::get_atomicDomain_totalmass(){
        double atomicDomain_totalmass = 0.;
        map<unsigned long long, double>::const_iterator iter;
        for(iter = _AtomicDomain.begin(); iter != _AtomicDomain.end(); iter++){
            atomicDomain_totalmass += iter->second;
        }
        return atomicDomain_totalmass;
    }
    
    
    bool AtomicSupport::inDomain(unsigned long long location) {
        return _AtomicDomain.find(location) != _AtomicDomain.end();
    }
    
    double AtomicSupport::getMass(unsigned long long location) {
        if (!inDomain(location)) {
            return 0;
        }
        return _AtomicDomain.find(location) -> second;
    }
    
    
    unsigned int AtomicSupport::getBin(unsigned long long location) {
        
        if (location < _lBoundaries.begin()->first || location >= _NatomLength) {
            throw logic_error("Searching for bin in invalid location in AtomicSupport::getBin.");
        }
        
        map<unsigned long long, unsigned int>::iterator binBound =
        _lBoundaries.upper_bound(location);
        
        if (binBound == _lBoundaries.begin()) {
            throw logic_error("Incorrectly finding upper bound in the first bin in AtomicSupport::getBin.");
        }
        
        binBound--;
        
        return binBound->second;
        
    }
    
    unsigned long long AtomicSupport::getMidLocation(unsigned int iBin) {
        unsigned long long lLocation, rLocation;
        lLocation = _lBoundariesByBin[iBin];
        if (iBin < _nBin-1) {
            rLocation = _lBoundariesByBin[iBin+1];
        } else {
            rLocation = _NatomLength;
        }
        
        return (lLocation + rLocation)/2;
    }
    
    unsigned long long AtomicSupport::getStartLocation(unsigned int iBin) {
        return _lBoundariesByBin[iBin];
    }
    
    unsigned long long AtomicSupport::getEndLocation(unsigned int iBin) {
        if (iBin < _nBin-1) {
            return _lBoundariesByBin[iBin+1];
        } else {
            return _NatomLength;
        }
    }
    
    unsigned long long AtomicSupport::binToLocation(unsigned int bin) {
        
        map<unsigned int, unsigned long long>::iterator boundIter =
        _lBoundariesByBin.find(bin);
        if (boundIter != _lBoundariesByBin.end()) {
            return _lBoundariesByBin[bin];
        }
        
        return _NatomLength + 1;
    }
    
    
    double AtomicSupport::getNAtomPriorProb(int delAtom, bool log) {
        
        // for now, only consider the poisson prior
        if (_alpha <= 0) {
            if (log) {
                return 0.;
            }
            return 1.;
        }
        
        if (_nAtom + delAtom > _NatomLength) {
            if (log) {
                return DOUB_NEGINF;
            }
            return 0.;
        }
        
        
        return randgen('P',(double) (_nAtom + delAtom),(double) getExpectedNAtom());
        
    }
    
    void AtomicSupport::makeProposal(double rng){
        
        // initialize the update
        cleanClearProposal();
        
        double minStep = 0.;
        double maxStep = 1.;
        
        if (_nAtom < 2) {
            maxStep = 0.75;
        }
        
        double updateStep = sub_func::runif(minStep, maxStep);
        
        // can only birth if there are zero atoms
        if (_nAtom == 0) {
            updateStep = 0.4;
        }
        
        // cannot birth an atom if there are the maximum number
        if (_nAtom >= _NatomLength) {
            if (updateStep < 0.5 && updateStep >= 0.25) {
                if (updateStep < 1/3) {
                    updateStep = 0.;
                } else {
                    updateStep = 2*(updateStep-0.25) + 0.5;
                }
            }
        }
        
        // select birth / death based on number of atoms
        double pDelete = 1.;
        if (updateStep < 0.5) {
            
            // following the notation of Skilling,
            // _alpha > 0 has a poisson prior
            // _alpha < 0 geometric
            // _alpha = 0 uniform
            if (_alpha > 0 ) {
                
                // prior on number of atoms is a Poisson
                // (note: need to modify better for maximum number of atoms)
                
                if (_nAtom >= _NatomLength) {
                    pDelete = 1.;
                } else if (_nAtom == 0) {
                    pDelete = 0.;
                } else {
                    double maxTerm = ((double) (_NatomLength - _nAtom)) /
                    ((double) _NatomLength);
                    pDelete = (double) _nAtom / ( (double)_nAtom + maxTerm*_alpha*_nBin);
                }
            } else if (_alpha < 0) {
                
                // prior on number of atoms is geometric
                
                double c = -_alpha*_nBin / (-_alpha*_nBin + 1.);
                
                
                if (_nAtom >= _NatomLength) {
                    pDelete = 1.;
                } else if (_nAtom == 0) {
                    pDelete = 0.;
                } else {
                    pDelete = (double) _nAtom / (((double)_nAtom+1.)*c + (double)_nAtom);
                }
            } else {
                
                // prior on number of atoms is uniform
                
                if (_nAtom >= _NatomLength) {
                    pDelete = 1.;
                } else if (_nAtom == 0) {
                    pDelete = 0.;
                } else {
                    pDelete = 0.5;
                }
            }
            
            if (sub_func::runif(0., 1.) < pDelete) {
                updateStep = 0.;
            } else  {
                updateStep = 0.4;
            }
        }
        
        //Determine operation based on updateStep
        if (updateStep < 0.25) {
            ProposeDeath();
        } else if (updateStep < 0.5) {
            ProposeBirth();
        } else if (updateStep < 0.75) {
            ProposeMove();
        } else {
            ProposeExchange();
        }
        
    }
    
    void AtomicSupport::ProposeDeath(){
        
        map<unsigned long long, double>::const_iterator iter;
        unsigned long long proposedLocation = _NatomLength;
        double proposedMass = 0.;
        
        // DEATH STEP
        _oper_type = 'D';
        
        // delete an atom
        
        // select an atom at random
        unsigned int deletionAtom = floor(sub_func::runif(0., 1.)*_nAtom);
        // find the proposed atom
        iter = _AtomicDomain.begin();
        if (deletionAtom > 0) {
            for (unsigned int iAtom = 0; iAtom < deletionAtom; iAtom++) {
                iter++;
                if (iter == _AtomicDomain.end())
                    throw logic_error("Attempting to delete a non-existant atom in AtomicSupport::makeProposal.");
            }
        }
        
        // find the proposed location and mass change
        proposedLocation = iter->first;
        
        // determine new mass based on deletion or resizing
        proposedMass = -iter->second;
        
        // mark the change in the list of atomic changes
        _proposedAtoms.insert(pair<unsigned long long, double>
                              (proposedLocation, proposedMass));
        
	}
	
    void AtomicSupport::ProposeBirth(){
        
        map<unsigned long long, double>::const_iterator iter;
        unsigned long long proposedLocation = _NatomLength;
        double proposedMass = 0.;
        
        // BIRTH STEP
        _oper_type = 'B';
        
        // add an atom
        
        // TODO: loop is inefficient for large number of atoms, REWRITE THIS LOOP
        // find a location not on the list of atoms
        do {
            // note select location between 0 and N-1, for a total of N locations
            proposedLocation = floor(sub_func::runif(0., 1.)*(_NatomLength-1));
            iter = _AtomicDomain.find(proposedLocation);
        } while (iter != _AtomicDomain.end());
        
        //generate a random mass - got rid of normAtomic statement here
        proposedMass = randgen('E', _lambda, 0);
        
        // mark the change in the list of atomic changes
        _proposedAtoms.insert(pair<unsigned long long, double>
                              (proposedLocation, proposedMass));
        
    }
    
    void AtomicSupport::ProposeMove(){
        
        map<unsigned long long, double>::const_iterator iter;
        
        // MOVE STEP
        _oper_type = 'M';
        
        // move an atom
        
        // select the atom to move
        unsigned int moveAtom = floor(sub_func::runif(0., 1.)*_nAtom);
        
        unsigned long long lbound, rbound;
        if (moveAtom == 0) {
            lbound = 0;
        }
        
        // find the proposed atom
        iter = _AtomicDomain.begin();
        if (moveAtom > 0) {
            for (unsigned int iAtom = 0; iAtom < moveAtom; iAtom++) {
                lbound = iter->first;
                if (iter == _AtomicDomain.end())
                    throw logic_error("Attempting to move a non-existent atom in AtomicSupport::makeProposal.");
                iter++;
            }
        }
        unsigned long long moveLocation = iter->first;
        double moveMass = iter->second;
        
        if (moveAtom < _nAtom-1) {
            iter++;
            rbound = iter->first;
        } else {
            rbound = _NatomLength-1;
        }
        
        double newLocation = sub_func::runif((double) lbound, (double) rbound);
        
        // allow the domain to wrap around
        /* if (moveAtom == 0) {
         
         unsigned long long otherLBound;
         for (iter = _AtomicDomain.begin(); iter != _AtomicDomain.end();
	     iter++) {
         otherLBound = iter->first;
         }
         
         double probRBin = (double) (_NatomLength-1 - otherLBound) /
         (double) (rbound + (_NatomLength-1 - otherLBound));
         if (probRBin < sub_func::runif(0., 1.)) {
         newLocation = sub_func::runif((double) otherLBound,
         (double) (_NatomLength-1));
         }
         
         } else if (moveAtom >= _nAtom -1) {
         unsigned long long otherRBound = _AtomicDomain.begin()->first;
         double probLBin = (double) (otherRBound) /
         (double) (otherRBound + (rbound - lbound));
         if (probLBin < sub_func::runif(0., 1.)) {
         newLocation = sub_func::runif(0., (double) otherRBound);
         }
         } */
        
        // add to proposal
        _proposedAtoms.insert(pair<unsigned long long, double>
                              (moveLocation, -moveMass));
        _proposedAtoms.insert(pair<unsigned long long, double>
                              (newLocation, moveMass));
	}
	
    void AtomicSupport::ProposeExchange(){
        
        map<unsigned long long, double>::const_iterator iter;
        
        // EXCHANGE STEP
        _oper_type = 'E';
        
        // exchange atoms (wlog, can always select left most atoms to exchange with
        // right neighbor, where last atom maps to first)
        
        unsigned int exchangeAtom = floor(sub_func::runif(0., 1.)*(_nAtom));
        
        unsigned long long atomLoc1, atomLoc2;
        double atomMass1, atomMass2;
        
        // find the proposed atom
        iter = _AtomicDomain.begin();
        if (exchangeAtom > 0) {
            for (unsigned int iAtom = 0; iAtom < exchangeAtom; iAtom++) {
                
                if (iter == _AtomicDomain.end())
                    throw logic_error("Attempting to move a non-existant atom in AtomicSupport::makeProposal.");
                iter++;
                
            }
        }
        
        atomLoc1 = iter->first;
        atomMass1 = iter->second;
        
        if (exchangeAtom < _nAtom - 1) {
            iter++;
            atomLoc2 = iter->first;
            atomMass2 = iter->second;
        } else {
            atomLoc2 = _AtomicDomain.begin()->first;
            atomMass2 = _AtomicDomain.begin()->second;
        }
        
        double newMass1, newMass2;
        
        if (atomMass1 + atomMass2 == 0) {
            newMass1 = 0;
            newMass2 = 0;
        } else {
            
            
            double pupper = sub_func::pgamma(atomMass1+atomMass2,2., 1./_lambda, DOUB_NEGINF, false);
            double newValue = sub_func::qgamma(sub_func::runif(0., pupper), 2., 1./_lambda, DOUB_NEGINF, false);
            
            if (atomMass1 > atomMass2) {
                //swap with bias to left
                newMass1 = newValue;
                newMass2 = atomMass1 + atomMass2 - newMass1;
            } else {
                //swap with bias to right
                newMass2 = newValue;
                newMass1 = atomMass1 + atomMass2 - newMass2;
            }
            
        }
        
        _proposedAtoms.insert(pair<unsigned long long, double>
                              (atomLoc1, newMass1-atomMass1));
        _proposedAtoms.insert(pair<unsigned long long, double>
                              (atomLoc2, newMass2-atomMass2));
        
    }
    
    void AtomicSupport::acceptProposal(bool updateIterCount) {
        
        // define local variables
        unsigned long long proposedLocation;
        double proposedMass;
        
        // define map iterators
        map<unsigned long long, double>::const_iterator proposeIter;
        map<unsigned long long, double>::iterator updateIter;
        
        // ensures that only one atom is updated at a time
        if (_proposedAtoms.size() > 2 && _iter > 0) {
            throw logic_error("Cannot update more than 2 atoms simultaneously.");
        }
        
        // loop over each new atom that was proposed
        for (map<unsigned long long, double>::iterator
             iter = _proposedAtoms.begin(); iter != _proposedAtoms.end();
             ++iter) {
            
            // find the update
            proposedLocation = iter->first;
            proposedMass     = iter->second;
            
            // adding the atom
            if (proposedMass > _epsilon) {
                
                // allow for addition of mass at an existing location
                proposeIter = _AtomicDomain.find(proposedLocation);
                if (proposeIter != _AtomicDomain.end()) {
                    
                    // find the true total mass at that location
                    proposedMass += proposeIter->second;
                    
                    // remove the atom because one will be added by default
                    _nAtom--;
                    cleanDeleteAtomicLocation(proposedLocation);
                }
                
                // add the atom
                _nAtom++;
                _AtomicDomain.insert(pair<unsigned long long, double>
                                     (proposedLocation, proposedMass));
            } // endif adding atoms
            
            // deleting the atom
            else if (proposedMass < 0) {
                
                // find the atom to remove
                proposeIter = _AtomicDomain.find(proposedLocation);
                if (proposeIter != _AtomicDomain.end()) {
                    
                    // allow for reduction of mass at an existing location
                    if (proposeIter->second + proposedMass > _epsilon) {
                        proposedMass += proposeIter->second;
                        cleanDeleteAtomicLocation(proposedLocation);
                        _AtomicDomain.insert(pair<unsigned long long, double>
                                             (proposedLocation, proposedMass));
                        
                    }
                    
                    // delete the atom
                    else {
                        _nAtom--;
                        cleanDeleteAtomicLocation(proposedLocation);
                    }
                }
                
            } // endif deleting atoms
        } // enddo iterator loop
        
        // erase atoms as all changes have been made
        cleanClearProposal();
        
        if (updateIterCount) {
            _iter++;
        }
    }
    
    void AtomicSupport::rejectProposal(bool updateIter) {
        cleanClearProposal();
        
        if (updateIter) {
            _iter++;
        }
        
    }
    
    void AtomicSupport::cleanClearProposal() {
        
        _proposedAtoms.clear();
    }
    
    void AtomicSupport::cleanClearAtomic() {
        
        _AtomicDomain.clear();
    }
    
    void AtomicSupport::cleanDeleteProposal(map<unsigned long long, double>::
                                            const_iterator iter) {
        cleanDeleteProposalLocation(iter->first);
    }
    
    void AtomicSupport::cleanDeleteAtomic(map<unsigned long long, double>::
                                          const_iterator iter) {
        cleanDeleteAtomicLocation(iter->first);
    }
    
    void AtomicSupport::cleanDeleteProposalLocation(unsigned long long location) {
        _proposedAtoms.erase(location);
    }
    
    void AtomicSupport::cleanDeleteAtomicLocation(unsigned long long location) {
        _AtomicDomain.erase(location);
    }
    
    void AtomicSupport::resetAtomicOutputThin(int thinDiag) {
        thinAtomicDiag = thinDiag;
        if (thinDiag > 0) {
            outputAtomicDiag = true;
            _initIterOutput = _iter;
        }
        atomicDiagFile << "Reset output to start at iteration " << _initIterOutput << " at every "
        << thinDiag << " iterations." << endl;
        //cout << "Reset output to start at iteration " << _initIterOutput << " at every "
        //	 << thinDiag << " iterations." << endl;
    }
    
    void AtomicSupport::writeAtomicHeader(char diagnosticFileName[],
                                          int thinDiag) {
        
        thinAtomicDiag = thinDiag;
        outputAtomicDiag = thinDiag >= 1;
        
        _initIterOutput = _iter;
        
        atomicDiagFile.open(diagnosticFileName, ios::out);
        
        // get the current time for output
        time_t rawtime;
        struct tm * timeinfo;
        time(&rawtime);
        timeinfo = localtime(&rawtime);
        atomicDiagFile << "Atomic diagnostic file created at\t" << asctime(timeinfo) << endl;
        
        // output the model parameters used
        atomicDiagFile << "number of bins\t" << _nBin << endl;
        atomicDiagFile << "maximum number of atoms\t"<< _NatomLength << endl;
        atomicDiagFile << "expected number of atoms per bin (alpha)\t" << _alpha << endl;
        atomicDiagFile << "expected mass of an atom (lambda)\t" << _lambda << endl;
        atomicDiagFile << endl;
    }
    
    void AtomicSupport::writeAtomicDiagnostics() {
        
        if (!outputAtomicDiag) {
            return;
        }
        
        // output results to a file
        if ( ((_iter - _initIterOutput) % thinAtomicDiag) !=  0) {
            return;
        }
        
        //cout << "outputting at: "<<_iter<< "; init: " << _initIterOutput << endl;
        
        unsigned int iBin;
        map<unsigned long long, double>::const_iterator iter;
        unsigned int nAtomPerBin[_nBin];
        double massPerBin[_nBin];
        for (iBin = 0; iBin < _nBin; iBin++) {
            nAtomPerBin[iBin] = 0;
            massPerBin[iBin] = 0.;
        }
        
        // output the location of atoms
        for (iter = _AtomicDomain.begin(); iter != _AtomicDomain.end(); iter++) {
            //atomicDiagFile << "\t" << iter->first;
            iBin = getBin(iter->first);
            nAtomPerBin[iBin]++;
            massPerBin[iBin] += iter->second;
        }
        atomicDiagFile << endl;
        
        for (iBin = 0; iBin < _nBin; iBin++) {
            atomicDiagFile << "\t" << massPerBin[iBin];
        }
        atomicDiagFile << endl;
        
    }
    
    void AtomicSupport::initializeAtomicBinary(char diagnosticFileName[]) {
        atomicDiagFileBinary.open(diagnosticFileName, ios::out);
    }
    
    void AtomicSupport::writeAtomicDiagnosticsBinary(bool isByCol, bool outByCol,
                                                     unsigned int nRow,
                                                     unsigned int nCol) {
        float outputValue;
        unsigned int iBin;
        double massPerBin[_nBin];
        for (iBin = 0; iBin < _nBin; iBin++) {
            massPerBin[iBin] = 0.;
        }
        
        map<unsigned long long, double>::const_iterator iter;
        for (iter = _AtomicDomain.begin(); iter != _AtomicDomain.end(); iter++) {
            iBin = getBin(iter->first);
            massPerBin[iBin] += iter->second;
        }
        
        
        // if both input and output are by row or column, no need to transpose
        if (isByCol == outByCol) {
            
            for (iBin = 0; iBin < _nBin; iBin++) {
                
                outputValue = (float) massPerBin[iBin];
                if (iBin > 0) {
                    atomicDiagFileBinary << "\t";
                }
                atomicDiagFileBinary << outputValue;
            }
            
            
        }
        // atomic data is stored by row by output is by column
        else if (!isByCol && outByCol) {
            
            // outputFirstBin = false;
            for (unsigned int iCol = 0; iCol < nCol; iCol++) {
                for (iBin = 0; iBin < nRow; iBin++) {
                    outputValue = (float) massPerBin[iCol + iBin*nCol];
                    
                    
                    if (iBin + iCol > 0) {
                        atomicDiagFileBinary << "\t";
                    }
                    atomicDiagFileBinary << outputValue;
                }
            }
            
        }
        // atomic data is stored by column and output is by row
        else {
            
            for (unsigned int iRow = 0; iRow < nRow; iRow++) {
                for (iBin = 0; iBin < nCol; iBin++) {
                    
                    
                    outputValue = (float) massPerBin[iRow + iBin*nRow];
                    
                    if (iBin + iRow > 0) {
                        atomicDiagFileBinary << "\t";
                    }
                    atomicDiagFileBinary << outputValue;
                    
                }
            }
            
        }
        atomicDiagFileBinary << endl;
        
        atomicDiagFileBinary.flush();
        
    }
    
    double AtomicSupport::getLambda() {
        return _lambda;
    }
    
    void AtomicSupport::setProposedAtomMass(const map<unsigned long long, double> 
                                            newProposal, bool isNewProposal) {
        
        map<unsigned long long, double>::const_iterator iter_orig;
        map<unsigned long long, double>::const_iterator iter_new;
        
        // check the values
        if (isNewProposal) {
            if (_proposedAtoms.size() > 0) {
                throw logic_error("Cannot create a new proposal when current one is not empty.");
            }
        } else {
            iter_orig = _proposedAtoms.begin();
            for (iter_new = newProposal.begin(); iter_new != newProposal.end();
                 iter_new++) {
                
                if (iter_orig == _proposedAtoms.end()) {
                    throw logic_error("Cannot change more atoms than originally proposed.");
                }
                
            }
            
            
            iter_orig++;
        }
        
        // set the proposal to have the new values
        cleanClearProposal();
        for (iter_new = newProposal.begin(); iter_new != newProposal.end();
             iter_new++) {
            _proposedAtoms.insert(pair<unsigned long long, double>
                                  (iter_new->first, iter_new->second));
        }
        
    }
    
    void AtomicSupport::updateAtomicBins(double binProbabilities[], 
                                         unsigned int length,
                                         bool onlyUpdateRelativeWidth) {
        
        // TODO -- debug this when also updating atomic locations
        
        
        if (length != _nBin) {
            throw logic_error("Must specify same number of bin probabilities as bins to update atomic bins in AtomicSupport::updateAtomicBins.");
        }
        
        unsigned long long origBoundaries[_nBin], newBoundaries[_nBin];
        
        double totalWeight = 0.;
        unsigned int iBin;
        unsigned long long lBoundValue = 0;
        
        // figure out the total value of the probabilities in which to divide bins
        for (iBin = 0; iBin < length; iBin++) {
            totalWeight += binProbabilities[iBin];
        }
        
        // compute threshold probability input to ensure there will be mass in each bin
        // TODO - right now ensures that the minimum is of unit one, when adapting probabilities 
        // as we go should revise to ensure that the threshold is sufficient to account for the number 
        // of atoms which are already in that bin
        vector<double> binProbVector (binProbabilities, binProbabilities+length);
        sort(binProbVector.begin(), binProbVector.end());
        
        
        double minMaxValue = 0, thresholdValue = 0;
        double newTotal = totalWeight, remainingTotal = totalWeight;
        int numRemoved = 0;
        vector<double>::iterator binProbIt, binProbIt2;
        for (binProbIt = binProbVector.begin(); 
             binProbIt != binProbVector.end(); binProbIt++) {
            
            // no need to correct further, all the remaining bins are wide enough
            if (((_NatomLength-1.)*(*binProbIt)/newTotal) >= 1) {
                break;
            }
            
            remainingTotal -= *binProbIt;
            numRemoved++;
            
            // only work last of a given probability value
            if (numRemoved < length) {
                binProbIt2 = binProbIt + 1;
                if (*binProbIt == *binProbIt2) {
                    continue;
                }
            } else {
                throw logic_error("Too many bins to reach threshold in AtomicSupport::updateAtomicBins");
            }
            
            minMaxValue = *binProbIt;
            thresholdValue = remainingTotal / (_NatomLength - 1. - (double) numRemoved);
            newTotal = remainingTotal + ((double) numRemoved) * thresholdValue;
            
        }
        
        // determine the new boundaries, accounting for threshold values
        for (iBin = 0; iBin < length; iBin++) {
            
            // can only obtain original values if already stored
            if (_AtomicDomain.size() > 0) {
                origBoundaries[iBin] = _lBoundariesByBin[iBin];
            }
            
            newBoundaries[iBin] = lBoundValue;
            
            // double check that the bin widths are non zero
            if (iBin > 0 && (newBoundaries[iBin] == newBoundaries[iBin-1])) {
                throw logic_error("AtomicSupport::updateAtomicBins: Attempting to create an atomic bin with zero width.");
            }
            
            if (binProbabilities[iBin] > minMaxValue) {
                
                lBoundValue += (unsigned long long) floor( ( (binProbabilities[iBin]/
                                                              newTotal)*
                                                            (_NatomLength - 1) ) );
            } else {
                
                lBoundValue += (unsigned long long) floor( ( (thresholdValue/newTotal)*
                                                            (_NatomLength-1) ) );
            }
            
        }
        
        // move around atoms to new locations
        // TODO - debug this -- will likely need to be careful to avoid putting two atoms 
        // in the same location in some cases, particularly if the bins are too narrow
        if (_AtomicDomain.size() > 0) {
            
            //cout << "moving existing atoms" << endl;
            
            unsigned long long newLocation[_AtomicDomain.size()];
            double atomSize[_AtomicDomain.size()];
            unsigned int nAtom = 0, iAtom;
            
            for (map<unsigned long long, double>::const_iterator iter = 
                 _AtomicDomain.begin(); iter != _AtomicDomain.end(); iter++) {
                iBin = getBin(iter->first);
                
                if (iBin < _nBin - 1) {
                    newLocation[nAtom] = newBoundaries[iBin] + 
                    (unsigned long long) (
                                          ((double)(iter->first - origBoundaries[iBin])) * 
                                          ((double)(newBoundaries[iBin+1] - newBoundaries[iBin])) / 
                                          ((double)(origBoundaries[iBin+1] - origBoundaries[iBin])));
                } else {
                    newLocation[nAtom] = newBoundaries[iBin] + 
                    (unsigned long long) (
                                          ((double)(iter->first - origBoundaries[iBin])) * 
                                          ((double)(_NatomLength - 1 - newBoundaries[iBin])) / 
                                          ((double)(_NatomLength - 1 - origBoundaries[iBin])));
                }
                atomSize[nAtom] = iter->second;
                nAtom++;
            }
            
            _AtomicDomain.clear();
            for (iAtom = 0; iAtom < nAtom; iAtom++) {
                _AtomicDomain.insert(pair<unsigned long long, double>(newLocation[iAtom],
                                                                      atomSize[iAtom]));
            }
            
            
        }
        
        _lBoundariesByBin.clear();
        _lBoundaries.clear();
        
        for (iBin = 0; iBin < _nBin; iBin++) {
            
            _lBoundariesByBin.insert(pair<unsigned int, unsigned long long>(iBin, newBoundaries[iBin]));
            _lBoundaries.insert(pair<unsigned long long, unsigned int>(newBoundaries[iBin], iBin));
            
            
        }
        
        if (!onlyUpdateRelativeWidth) {
            // TODO - consider changing other parameters to reflect true probabilities
        }
        
    }
    
    // --------- to extract _atomic_domain_label
    char AtomicSupport::get_atomic_domain_label(){
        return _atomic_domain_label;
    }
    
    char AtomicSupport::get_oper_type(){
        return _oper_type;
    }
    
    
}

