// CoGAPS C++ Version
//
// GibbsSamplerMap:
// Creation of a GibbsSampler that deals with fixed patterns
// This class inherits from GibbsSampler
//
// History: v 1.0  August 7, 2014
//

#include <iostream>
#include <cmath>
#include <limits>
#include <stdexcept>
#include "GibbsSamplerMap.h"


using namespace std;
using namespace gaps;
using std::vector;

// ******************** CONSTRUCTOR ********************************************
GibbsSamplerMap::GibbsSamplerMap(unsigned long nEquil, unsigned long nSample, unsigned int nFactor,
                                 double alphaA, double alphaP, double nMaxA, double nMaxP,
                                 unsigned long nIterA, unsigned long nIterP,
                                 double max_gibbsmass_paraA, double max_gibbsmass_paraP,
                                 unsigned long long atomicSize,
                                 char label_A, char label_P, char label_D, char label_S,
                                 vector<vector<double> > &DVector, vector<vector<double> > &SVector,
                                 const string &simulation_id,
                                 vector <vector <double> >  &parameters, char the_fixed_matrix)
    : GibbsSampler(nEquil, nSample, nFactor,
                   alphaA, alphaP, nMaxA, nMaxP, nIterA, nIterP, max_gibbsmass_paraA,
                   max_gibbsmass_paraP, atomicSize,
                   label_A, label_P, label_D, label_S, DVector, SVector, simulation_id) {
    //Constructor body
    for (int i = 0; i < parameters.size(); i++) {
        vector <double> newPattern;
        newPattern.clear();
        double PatSum = 0;

        for (int j = 0; j < parameters[0].size(); j++) {
            PatSum += parameters[i][j];
        }

        for (int j = 0; j < parameters[0].size(); j++) {
            newPattern.push_back(parameters[i][j] / PatSum);
        }

        _MapValues.push_back(newPattern);
    }

    _nFixedMaps = parameters.size();
    _the_fixed_matrix = the_fixed_matrix;
}


// *************** METHODS FOR INITIALIZATION, DISPLAY, OUTPUT ***********************

// Per the JAGS version of the code, initializing the matrices with
// the given fixed values as follows: read in the patterns from the
// vector _MapValues. For matrix A, start with the first column and assign
// the patterns from left to right. For matrix P,
// start with row nFactor (number of rows in P) - nFixedMaps
// and go down to the end of the matrix
// The problem here now becomes that the matrix gets out of sync with the atomic
// domain. As a result, we call (in main) the method initialize_atomic_domain_map
// to correct this.
void GibbsSamplerMap::init_Mapped_Matrix() {
    vector <double> thefixedPat;

    switch (_the_fixed_matrix) {
        case 'A': {
            for (int iCol = 0; iCol < _nFixedMaps; iCol++) {
                thefixedPat = _MapValues.at(iCol);
                _AMatrix.setCol(thefixedPat, iCol);
            }

            break;
        }

        case 'P': {
            for (int iRow = 0; iRow < _nFixedMaps; iRow++) {
                thefixedPat = _MapValues.at(iRow);
                _PMatrix.setRow(thefixedPat, (_nFactor - _nFixedMaps + iRow));
            }

            break;
        }
    } // end switch block
} //end method

// Keep the atomic domain consistent with the matrix to avoid
// a mass inconsistency. Get the row or column of the fixed pattern
// and birth new atoms across it in the middle of the correct bin.
void GibbsSamplerMap::initialize_atomic_domain_map() {
    map <unsigned long long, double > getNSync;
    unsigned long long updateloc;
    unsigned int updatebin;
    double updatemass;
    vector <double> fixPat;

    switch (_the_fixed_matrix) {
        case 'A': {
            // For A, we start in the 1st column at bin 0. Bins count horizontally
            // so we increment the bin by _nFactor to go down a column.
            for (int iCol = 0; iCol < _nFixedMaps; iCol++) {
                fixPat = _MapValues.at(iCol);
                //reset updatebin
                updatebin = iCol;

                for (int iAtom = 0; iAtom < _nRow; iAtom++) {
                    // Place the atom in the middle of the bin
                    updateloc = _AAtomicdomain.getMidLocation(updatebin);
                    updatemass = fixPat.at(iAtom);
                    getNSync.insert(pair<unsigned long long, double>(updateloc, updatemass));
                    updatebin += _nFactor;
                }
            }

            //Proposal is now finished. Accept it to update the atomic domain.
            _AAtomicdomain.setProposedAtomMass(getNSync, true);
            _AAtomicdomain.acceptProposal(false);
            break;
        } // end case A

        case 'P': {
            // For P, we start with row nFactor (number of rows in P) - nFixedMaps
            // and go down to the end of the matrix. We increment the bin numbers
            // which here count vertically, by nFactor.
            for (int iRow = 0; iRow < _nFixedMaps; iRow++) {
                fixPat = _MapValues.at(iRow);
                //reset updatebin
                updatebin = _nFactor - _nFixedMaps + iRow;

                for (int iAtom = 0; iAtom < _nCol; iAtom++) {
                    // Place the atom in the middle of the bin
                    updateloc = _PAtomicdomain.getMidLocation(updatebin);
                    updatemass = fixPat.at(iAtom);
                    getNSync.insert(pair<unsigned long long, double>(updateloc, updatemass));
                    updatebin += _nFactor;
                }
            }

            //Proposal is now finished. Accept it to update the atomic domain.
            _PAtomicdomain.setProposedAtomMass(getNSync, true);
            _PAtomicdomain.acceptProposal(false);
            break;
        } // end case P
    } // end switch block
} // end initialize_atomic_domain_map


// ************* METHODS FOR COMPUTING LIKELIHOOD FUNCTIONS ******************

// For one pattern change (in birth and death)
double GibbsSamplerMap::computeDeltaLLBDMap(char the_matrix_label,
        double const *const *D,
        double const *const *S,
        double const *const *A,
        double const *const *P,
        vector <double> &newPat, unsigned int chPat) {
    double DelLL;

    switch (the_matrix_label) {
        case 'A': {
            DelLL = GAPSNorm::calcDeltaLLMap('A', D, S, A, P, newPat, chPat, _nRow,
                                             _nCol, _nFactor);
            break;
        } // end of switch-block 'A'

        case 'P': {
            DelLL = GAPSNorm::calcDeltaLLMap('P', D, S, A, P, newPat, chPat, _nRow,
                                             _nCol, _nFactor);
            break;
        } // end of switch-block 'P'
    } // end of switch block

    return DelLL;
} // end of computeDeltaLLBDMap

// For two pattern changes (in move and exchange)
double GibbsSamplerMap::computeDeltaLLMEMap(char the_matrix_label,
        double const *const *D,
        double const *const *S,
        double const *const *A,
        double const *const *P,
        vector <double> &newPat1, unsigned int chPat1,
        vector <double> &newPat2, unsigned int chPat2) {
    double DelLL;

    switch (the_matrix_label) {
        case 'A': {
            DelLL = GAPSNorm::calcDeltaLL2Map('A', D, S, A, P, newPat1, chPat1,
                                              newPat2, chPat2, _nRow, _nCol, _nFactor);
            break;
        } // end of switch-block 'A'

        case 'P': {
            DelLL = GAPSNorm::calcDeltaLL2Map('P', D, S, A, P, newPat1, chPat1,
                                              newPat2, chPat2, _nRow, _nCol, _nFactor);
            break;
        } // end of switch-block 'P'
    } // end of switch block

    return DelLL;
} // end of computeDeltaLLMap


// ******************* METHODS FOR THE UPDATE/PROPOSAL *************************
// -----------------------------------------------------------------------------
// For the modified update method, need to check the location of the atomic
// proposal. If the proposal falls in a fixed area, distribute the mass of the
// proposal across the row such that the percentage of total mass in each fixed
// bin remains the same. Otherwise continue normally.
void GibbsSamplerMap::mapUpdate(char the_matrix_label) {
    // Send directly to the regular Gibbs and break from this method
    // if the matrix differs from the fixed one
    if (_the_fixed_matrix != the_matrix_label) {
        update(the_matrix_label); // (GibbsSampler update method)
        return;
    }

    double **D = _DMatrix.get_matrix();
    double **S = _SMatrix.get_matrix();
    double **AOrig = _AMatrix.get_matrix();
    double **POrig = _PMatrix.get_matrix();
    bool Q_update;

    switch (the_matrix_label) {
        case 'A': {
            // ----------- making a proposal from atomic space A:
            _AAtomicdomain.makeProposal();
            get_oper_type('A');
            _atomicProposal = _AAtomicdomain.getProposedAtoms();
            extract_atomicProposal('A');

            // Check to make sure there are no problems with the proposal
            if (_nChange_atomicProposal == 0) {}

            if (_nChange_atomicProposal > 2) {
                throw logic_error("GibbsSampler: can't change more than two atoms!!");
            }

            // Initialize necessary quantities
            map<unsigned long long, double>::const_iterator atom;
            unsigned long long loc1, loc2;
            bool fixed1, fixed2;
            atom = _atomicProposal.begin();
            loc1 = atom->first;
            fixed1 = Q_fixed(loc1, 'A');

            //------------------------------------------------------
            // Break into 2 cases: birth/death and move/exchange based
            // on the size of the atomic proposal. Within each case,
            // determine if the change maps to a fixed pattern. If yes,
            // send the change to seperate method. Otherwise, send to
            // regular Gibbs methods.
            //-------------------------------------------------------
            // Birth/Death Cases
            if (_nChange_atomicProposal == 1) {
                if (fixed1) {
                    // Do nothing-fixed patterns locked
                }//end if (fixed pattern) block
                else { // not a fixed pattern, call the methods from Gibbs
                    if (_oper_type == 'D') {
                        Q_update = death('A', D, S, AOrig, POrig);

                    } else {
                        Q_update = birth('A', D, S, AOrig, POrig);
                    }

                    // Modify the proposal in normal Gibbs manner
                    // Update the matrix with improved update, if there are further updates
                    if (Q_update == true) {
                        _AMatrix.matrix_Elem_update(_new_matrixElemChange, _oper_type, _new_nChange_matrixElemChange);
                    }
                } // end if block for fixed or not fixed patterns
            } //end of birth/death case

            // -----------------------------------------------------------------
            // Move/Exchange Cases
            // If both of the locations is sent to a fixed pattern,
            // send the change to a seperate method which WILL change the matrix
            // for the fixed locations if allowed by MCMC. If one, but not both, of
            // the locations is in a fixed pattern, do nothing.
            // Otherwise, send to regular Gibbs methods.
            // -----------------------------------------------------------------
            else {
                atom++;
                loc2 = atom->first;
                fixed2 = Q_fixed(loc2, 'A');

                if (fixed1 && fixed2) {
                    // Do nothing-fixed patterns locked
                } // end if one or both locations in a fixed pattern
                else if (!fixed1 && !fixed2) { //both locations not in fixed pattern, normal Gibbs */
                    if (_oper_type == 'M') {
                        Q_update = move('A', D, S, AOrig, POrig);

                    } else {
                        Q_update = exchange('A', D, S, AOrig, POrig);
                    }

                    // Modify the proposal in normal Gibbs manner
                    // Update the matrix with improved update, if there are further updates
                    if (Q_update == true) {
                        _AMatrix.matrix_Elem_update(_new_matrixElemChange, _oper_type, _new_nChange_matrixElemChange);
                    }
                }
            } //end Move/Exchange cases

            break;
        } // end of case 'A'

        case 'P': {
            // ----------- making a proposal from atomic space P:
            _PAtomicdomain.makeProposal();
            get_oper_type('P');
            _atomicProposal = _PAtomicdomain.getProposedAtoms();
            extract_atomicProposal('P');

            // Check to make sure there are no problems with the proposal
            if (_nChange_atomicProposal == 0) {}

            if (_nChange_atomicProposal > 2) {
                throw logic_error("GibbsSampler: can't change more than two atoms!!");
            }

            // Initialize necessary quantities
            map<unsigned long long, double>::const_iterator atom;
            unsigned long long loc1, loc2;
            bool fixed1, fixed2;
            atom = _atomicProposal.begin();
            loc1 = atom->first;
            fixed1 = Q_fixed(loc1, 'P');

            //------------------------------------------------------
            // Break into 2 cases: birth/death and move/exchange based
            // on the size of the atomic proposal. Within each case,
            // determine if the change maps to a fixed pattern. If yes,
            // send the change to seperate method. Otherwise, send to
            // regular Gibbs methods.
            //-------------------------------------------------------
            // Birth/Death Cases
            if (_nChange_atomicProposal == 1) {
                if (fixed1) {
                    // Do nothing-fixed patterns locked
                }//end if (fixed pattern) block
                else { // not a fixed pattern, call the methods from Gibbs
                    if (_oper_type == 'D') {
                        Q_update = death('P', D, S, AOrig, POrig);

                    } else {
                        Q_update = birth('P', D, S, AOrig, POrig);
                    }

                    // Modify the proposal in normal Gibbs manner
                    // Update the matrix with improved update, if there are further updates
                    if (Q_update == true) {
                        _PMatrix.matrix_Elem_update(_new_matrixElemChange, _oper_type, _new_nChange_matrixElemChange);
                    }
                }
            } //end if birth death

            // -----------------------------------------------------------------
            // Move/Exchange Cases
            // If both of the locations is sent to a fixed pattern,
            // send the change to a seperate method which WILL change the matrix
            // for the fixed location(s) ONLY if allowed by MCMC.
            // Otherwise, send to regular Gibbs methods.
            // If only one location is in a fixed pattern, do nothing.
            // -----------------------------------------------------------------
            else {
                atom++;
                loc2 = atom->first;
                fixed2 = Q_fixed(loc2, 'P');

                if (fixed1 && fixed2) {
                    // Do nothing-fixed patterns locked
                } // end if one or both locations in a fixed pattern
                else if (!fixed1 && !fixed2) { //both locations not in fixed pattern, normal Gibbs
                    if (_oper_type == 'M') {
                        Q_update = move('P', D, S, AOrig, POrig);

                    } else {
                        Q_update = exchange('P', D, S, AOrig, POrig);
                    }

                    // Modify the proposal in normal Gibbs manner
                    // Update the matrix with improved update, if there are further updates
                    if (Q_update == true) {
                        _PMatrix.matrix_Elem_update(_new_matrixElemChange, _oper_type, _new_nChange_matrixElemChange);
                    }
                }
            } //end Move/Exchange cases

            break;
        } // end of case 'P'
    } // end of switch block

    // clear Proposal for the next run
    clear_Proposal();
    clear_new_Proposal();
} // end of update()

// ----------------------------------------------------------------------------

// Determine whether or not a location falls in a fixed pattern
bool GibbsSamplerMap::Q_fixed(unsigned long long location,
                              char the_matrix_label) {
    unsigned int theBin, theRow, theCol;

    // If the matrix being tested isn't the fixed matrix, stop
    if (the_matrix_label != _the_fixed_matrix) {
        return false;
    }

    switch (the_matrix_label) {
        case 'A': {
            theBin = _AAtomicdomain.getBin(location);
            theCol = getCol('A', theBin);

            // if the change is in rows 0 to nFixedMaps of A,
            // the change is in a fixed row - return true
            if (theCol < _nFixedMaps) {
                return true;
            }

            break;
        } //end case A

        case 'P': {
            theBin = _PAtomicdomain.getBin(location);
            theRow = getRow('P', theBin);

            // if the change is in rows _nFactor - _nFixedMaps to
            // _nFactor of P,  the change is in a fixed row - return true
            if (theRow >= (_nFactor - _nFixedMaps)) {
                return true;
            }

            break;
        } //end case P
    } //end switch block

    return false;
}


// ==============METHODS FOR TEST GAPS================

// Manually calculate the matrix A from the atomic space passed in.
vector <vector <double> > GibbsSamplerMap::createSampleAMatMap(map <unsigned long long, double> ADomain) {
    // Make that a matrix
    vector <vector <double> > SampleAMatrix;
    SampleAMatrix.resize(_nRow);

    for (int i = 0; i < _nRow; i++) {
        SampleAMatrix[i].resize(_nFactor, 0);
    }

    //Copy the parsed domain
    map<unsigned long long, double>::const_iterator iter;

    for (iter = ADomain.begin(); iter != ADomain.end(); iter++) {
        // Check if the location falls in a fixed pattern
        // If it does, add the mass in the mapped way
        unsigned int theBin = _AAtomicdomain.getBin(iter->first);
        // Put the mass in the bin in the matrix
        unsigned int theRow = getRow('A', theBin);
        unsigned int theCol = getCol('A', theBin);
        bool fixedLoc = Q_fixed(iter->first, 'A');

        if (fixedLoc) {
            for (int iRow = 0; iRow < _nRow; iRow++) {
                SampleAMatrix[iRow][theCol] += (_MapValues[theCol][iRow]) * (iter->second);
            }

        } else {
            // Otherwise add in a regular way
            SampleAMatrix[theRow][theCol] += iter->second;
        }
    }

    return SampleAMatrix;
}

// Manually calculate the matrix P from the atomic space passed in.
vector <vector <double> > GibbsSamplerMap::createSamplePMatMap(map <unsigned long long, double> PDomain) {
    vector <vector <double> > SamplePMatrix;
    SamplePMatrix.resize(_nFactor);

    for (int i = 0; i < _nFactor; i++) {
        SamplePMatrix[i].resize(_nCol, 0);
    }

    //Copy the parsed domain
    map<unsigned long long, double>::const_iterator iter;

    for (iter = PDomain.begin(); iter != PDomain.end(); iter++) {
        // Check if the location falls in a fixed pattern
        // If it does, add the mass in the mapped way
        unsigned int theBin = _PAtomicdomain.getBin(iter->first);
        // Put the mass in the bin in the matrix
        unsigned int theRow = getRow('P', theBin);
        unsigned int theCol = getCol('P', theBin);
        bool fixedLoc = Q_fixed(iter->first, 'P');

        if (fixedLoc) {
            int fixedPatNum = theRow - (_nFactor - _nFixedMaps);

            for (int iCol = 0; iCol < _nCol; iCol++) {
                SamplePMatrix[theRow][iCol] += (_MapValues[fixedPatNum][iCol]) * (iter->second);
            }

        } else {
            SamplePMatrix[theRow][theCol] += iter->second;
        }
    }

    return SamplePMatrix;
}

