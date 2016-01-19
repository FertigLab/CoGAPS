//#include <config.h>

#include "AtomicSupport.h"

#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <string.h>
#include <time.h>
#include <algorithm>
#include <iostream>
#include<vector>
#include<limits>

// ------------ added for reading input file -----------------------------------
#include <sstream>
#include <string>
#include <istream>
#include <iterator>
#include <fstream>

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
// For equal bin sizes
void AtomicSupport::initializeAtomic(unsigned int nBin, unsigned long long NatomLength,
                                     double alpha, double lambda, char atomic_domain_label) {
    // parameter values
    _nBin      = nBin;
    _NatomLength = NatomLength;
    _alpha     = alpha;
    _lambda    = lambda;
    _atomic_domain_label = atomic_domain_label;
    // set the small number
    _epsilon = lepsilon;

    // check the validity of the parameters
    if (lambda <= 0) {
        throw logic_error("Invalid value for lambda in AtomicSupport.");
    }

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

// Use this "constructor" when the bins have different probabilities (priors).
void AtomicSupport:: FixedBins_initializeAtomic(unsigned int nBin, unsigned long long NatomLength,
        double alpha, double lambda, char atomic_domain_label,
        std::vector<std::vector<double> > ReadBinProbs) {
    // parameter values
    _nBin      = nBin;
    _NatomLength = NatomLength;
    _alpha     = alpha;
    _lambda    = lambda;
    _atomic_domain_label = atomic_domain_label;
    // set the small number
    _epsilon = lepsilon;

    // check the validity of the parameters
    if (lambda <= 0) {
        throw logic_error("Invalid value for lambda in AtomicSupport.");
    }

    // initialize the internal atomic domain
    _nAtom = 0;
    cleanClearAtomic();
    cleanClearProposal();
    // initialize the diagnostic file
    _iter = 0;
    // by default, set boundaries to be uniform -- TODO, check this
    double binProb[nBin];
    // Reading in the bin weights from R
    // Populate the list of bin probabilities based on the bin numbers for A and P.
    // For A, bins increment as you go along a row. For P, they increment as you
    // go down a column.
    int n_row = ReadBinProbs.size();
    int n_col = ReadBinProbs[1].size();
    int length = n_row * n_col;

    if (length != nBin) {
        throw logic_error("Number of bins don't line up - you probably passed in bins for the wrong matrix.");
    }

    int iBin = 0;

    if (_atomic_domain_label == 'A') {
        for (int iRow = 0; iRow < n_row; iRow++) {
            for (int iCol = 0; iCol < n_col; iCol++) {
                binProb[iBin] = ReadBinProbs[iRow][iCol];
                iBin++;
            }
        }

    } else { //if (theFixedDomain == 'P')
        for (int iCol = 0; iCol < n_col; iCol++) {
            for (int iRow = 0; iRow < n_row; iRow++) {
                binProb[iBin] = ReadBinProbs[iRow][iCol];
                iBin++;
            }
        }
    }

    updateAtomicBins(binProb, nBin, true); // true means to update only the bin width
}

void AtomicSupport::writeAtomicInfo(char outputFilename[], unsigned long Samp_cycle) {
    ofstream outputFile;

    if (Samp_cycle == 1) {
        outputFile.open(outputFilename, ios::out);

    } else {
        outputFile.open(outputFilename, ios::out | ios::app);
    }

    map<unsigned long long, double>::const_iterator iter;

    for (iter = _AtomicDomain.begin(); iter != _AtomicDomain.end();
            iter++) {
        outputFile << setiosflags(ios::right) << setw(25) << iter->first << " "
                   << setw(15) << iter->second << endl;
    }

    outputFile.close();
}

// Generate a random location for the atom to be birthed in.
// June 2015 - this generates a random 64 bit integer. Different
// than the prior method which used rand()*_NAtomLength
unsigned long long AtomicSupport::birthAtomLocation() {
    boost::mt19937 randGen(rng);
    boost::uniform_int<uint64_t> uInt64Dist(0, std::numeric_limits<uint64_t>::max());
    boost::variate_generator<boost::mt19937 &, boost::uniform_int<uint64_t> > getRand(randGen, uInt64Dist);
    uint64_t location = getRand();
    return (location);
}


double AtomicSupport::get_atomicDomain_totalmass() {
    double atomicDomain_totalmass = 0.;
    map<unsigned long long, double>::const_iterator iter;

    for (iter = _AtomicDomain.begin(); iter != _AtomicDomain.end(); iter++) {
        atomicDomain_totalmass += iter->second;
    }

    return atomicDomain_totalmass;
}

map <unsigned long long, double> AtomicSupport::getDomain() {
    return _AtomicDomain;
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

    if (iBin < _nBin - 1) {
        rLocation = _lBoundariesByBin[iBin + 1];

    } else {
        rLocation = _NatomLength;
    }

    return (lLocation + rLocation) / 2;
}

void AtomicSupport::makeProposal() {
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
            if (updateStep < 1 / 3) {
                updateStep = 0.;

            } else {
                updateStep = 2 * (updateStep - 0.25) + 0.5;
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
        if (_alpha > 0) {
            // prior on number of atoms is a Poisson
            // (note: need to modify better for maximum number of atoms)
            if (_nAtom >= _NatomLength) {
                pDelete = 1.;

            } else if (_nAtom == 0) {
                pDelete = 0.;

            } else {
                double maxTerm = ((double)(_NatomLength - _nAtom)) /
                                 ((double) _NatomLength);
                pDelete = (double) _nAtom / ((double)_nAtom + maxTerm * _alpha * _nBin);
            }

        } else if (_alpha < 0) {
            // prior on number of atoms is geometric
            double c = -_alpha * _nBin / (-_alpha * _nBin + 1.);

            if (_nAtom >= _NatomLength) {
                pDelete = 1.;

            } else if (_nAtom == 0) {
                pDelete = 0.;

            } else {
                pDelete = (double) _nAtom / (((double)_nAtom + 1.) * c + (double)_nAtom);
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

void AtomicSupport::ProposeDeath() {
    map<unsigned long long, double>::const_iterator iter;
    unsigned long long proposedLocation = _NatomLength;
    double proposedMass = 0.;
    // DEATH STEP
    _oper_type = 'D';
    // delete an atom
    // select an atom at random
    unsigned int deletionAtom = floor(sub_func::runif(0., 1.) * _nAtom);
    // find the proposed atom
    iter = _AtomicDomain.begin();

    if (deletionAtom > 0) {
        for (unsigned int iAtom = 0; iAtom < deletionAtom; iAtom++) {
            iter++;

            if (iter == _AtomicDomain.end()) {
                throw logic_error("Attempting to delete a non-existant atom in AtomicSupport::makeProposal.");
            }
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

void AtomicSupport::ProposeBirth() {
    // BIRTH STEP
    _oper_type = 'B';
    map<unsigned long long, double>::const_iterator iter = _AtomicDomain.begin();
    unsigned long long proposedLocation = 0;
    double proposedMass = 0.;

    // Need to make sure the location doesn't already have an atom in it.
    // We use the find function for this. If the location isn't in the list
    // the find function returns the end of the atomic domain and this loop
    // stops.
    while (iter != _AtomicDomain.end()) {
        proposedLocation = birthAtomLocation();
        iter = _AtomicDomain.find(proposedLocation);
    }

    //generate a random mass - got rid of normAtomic statement here
    proposedMass = randgen('E', _lambda, 0);
    // mark the change in the list of atomic changes
    _proposedAtoms.insert(pair<unsigned long long, double>
                          (proposedLocation, proposedMass));
}

void AtomicSupport::ProposeMove() {
    map<unsigned long long, double>::const_iterator iter;
    // MOVE STEP
    _oper_type = 'M';
    // move an atom
    // select the atom to move
    unsigned int moveAtom = floor(sub_func::runif(0., 1.) * _nAtom);
    unsigned long long lbound, rbound;

    if (moveAtom == 0) {
        lbound = 0;
    }

    // find the proposed atom
    iter = _AtomicDomain.begin();

    if (moveAtom > 0) {
        for (unsigned int iAtom = 0; iAtom < moveAtom; iAtom++) {
            lbound = iter->first;

            if (iter == _AtomicDomain.end()) {
                throw logic_error("Attempting to move a non-existent atom in AtomicSupport::makeProposal.");
            }

            iter++;
        }
    }

    unsigned long long moveLocation = iter->first;
    double moveMass = iter->second;

    if (moveAtom < _nAtom - 1) {
        iter++;
        rbound = iter->first;

    } else {
        rbound = _NatomLength - 1;
    }

    double newLocation = sub_func::runif((double) lbound, (double) rbound);
    // add to proposal
    _proposedAtoms.insert(pair<unsigned long long, double>
                          (moveLocation, -moveMass));
    _proposedAtoms.insert(pair<unsigned long long, double>
                          (newLocation, moveMass));
}

void AtomicSupport::ProposeExchange() {
    map<unsigned long long, double>::const_iterator iter;
    // EXCHANGE STEP
    _oper_type = 'E';
    // exchange atoms (wlog, can always select left most atoms to exchange with
    // right neighbor, where last atom maps to first)
    unsigned int exchangeAtom = floor(sub_func::runif(0., 1.) * (_nAtom));
    unsigned long long atomLoc1, atomLoc2;
    double atomMass1, atomMass2;
    // find the proposed atom
    iter = _AtomicDomain.begin();

    if (exchangeAtom > 0) {
        for (unsigned int iAtom = 0; iAtom < exchangeAtom; iAtom++) {
            if (iter == _AtomicDomain.end()) {
                throw logic_error("Attempting to move a non-existant atom in AtomicSupport::makeProposal.");
            }

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
        double pupper = sub_func::pgamma(atomMass1 + atomMass2, 2., 1. / _lambda, DOUB_NEGINF, false);
        double newValue = sub_func::qgamma(sub_func::runif(0., pupper), 2., 1. / _lambda, DOUB_NEGINF, false);

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
                          (atomLoc1, newMass1 - atomMass1));
    _proposedAtoms.insert(pair<unsigned long long, double>
                          (atomLoc2, newMass2 - atomMass2));
}

void AtomicSupport::acceptProposal(bool updateIterCount) {
    // define local variables
    unsigned long long proposedLocation;
    double proposedMass;
    // define map iterators
    map<unsigned long long, double>::const_iterator proposeIter;

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

    unsigned long long newBoundaries[_nBin];
    double totalWeight = 0.;
    unsigned int iBin;
    unsigned long long lBoundValue = 0;

    // figure out the total value of the probabilities in which to divide bins
    for (iBin = 0; iBin < length; iBin++) {
        totalWeight += binProbabilities[iBin];
    }

    // Code deleted here that computed threshold values to make sure
    // there is mass in each bin

    // determine the new boundaries
    for (iBin = 0; iBin < length; iBin++) {
        newBoundaries[iBin] = lBoundValue;

        // double check that the bin widths are non zero
        if (iBin > 0 && (newBoundaries[iBin] == newBoundaries[iBin - 1])) {
            throw logic_error("AtomicSupport::updateAtomicBins: Attempting to create an atomic bin with zero width.");
        }

        lBoundValue += (unsigned long long) floor(((binProbabilities[iBin] / totalWeight) *
                       (_NatomLength - 1)));
    }

    // Code here move around atoms to new locations
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

char AtomicSupport::get_oper_type() {
    return _oper_type;
}


}

