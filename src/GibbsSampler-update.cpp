#include "GibbsSampler.h"

#include "MatAlgo.h"

// ************* METHODS FOR COMPUTING LIKELIHOOD FUNCTIONS ******************
// ---------------------------------------------------------------------------
double GibbsSampler::cal_logLikelihood() {
    return GAPSNorm::calChi2(_DMatrix, _SMatrix, _AMatrix, _PMatrix,
        _nFactor) / 2.0;
}

// -----------------------------------------------------------------------------
double GibbsSampler::computeDeltaLL(char the_matrix_label,
                                    const Matrix &D,
                                    const Matrix &S,
                                    const Matrix &A,
                                    const Matrix &P,
                                    unsigned int the_nChange_matrixElemChange,
                                    const vector<ElementChange> the_matrixElemChange) {
    double DelLL = 0.0;

    switch (the_matrix_label) {
        case 'A': {
            if (the_nChange_matrixElemChange == 0) {
                DelLL = 0.0;

            } else if (the_nChange_matrixElemChange == 1) {
                DelLL = GAPSNorm::calcDeltaLL1E('A', D, S, A, P,
                    the_matrixElemChange, _nFactor);

            } else if (the_nChange_matrixElemChange == 2) {
                DelLL = GAPSNorm::calcDeltaLL2E('A', D, S, A, P,
                    the_matrixElemChange, _nFactor);

            } else {
                DelLL = GAPSNorm::calcDeltaLLGen('A', D, S, A, P,
                    the_matrixElemChange, _nFactor);
            } // end of if-block according to proposal.size()

            break;
        } // end of switch-block 'A'

        case 'P': {
            if (the_nChange_matrixElemChange == 0) {
                DelLL = 0.0;

            } else if (the_nChange_matrixElemChange == 1) {
                DelLL = GAPSNorm::calcDeltaLL1E('P', D, S, A, P,
                    the_matrixElemChange, _nFactor);

            } else if (the_nChange_matrixElemChange == 2) {
                DelLL = GAPSNorm::calcDeltaLL2E('P', D, S, A, P,
                    the_matrixElemChange, _nFactor);

            } else {
                DelLL = GAPSNorm::calcDeltaLLGen('P', D, S, A, P,
                    the_matrixElemChange, _nFactor);
            } // end of if-block according to proposal.size()

            break;
        } // end of switch-block 'P'
    } // end of switch block

    return DelLL;
} // end of computeDeltaLL

// ----------------------------------------------------------------------------
void GibbsSampler::extract_atomicProposal(char the_matrix_label) {
    unsigned int bin, chRow, chCol;
    double chmass;
    map<unsigned long long, double>::const_iterator iter;
    _nChange_matrixElemChange = 0;
    _nChange_atomicProposal = _atomicProposal.size();

    if (_nChange_atomicProposal == 0) {   // atomic proposal size = 0
        _nChange_matrixElemChange = 0;

    } else if (_nChange_atomicProposal == 1) { // atomic proposal size = 1
        iter = _atomicProposal.begin();

        switch (the_matrix_label) {
            case 'A': {
                bin = _AAtomicdomain.getBin(iter->first);
                chRow = getRow('A', bin);
                chCol = getCol('A', bin);
                break;
            }

            case 'P': {
                bin = _PAtomicdomain.getBin(iter->first);
                chRow = getRow('P', bin);
                chCol = getCol('P', bin);
                break;
            }
        } // end of switch-block for atomic proposal size = 1

        chmass = iter->second;
        _Row_changed.push_back(chRow);
        _Col_changed.push_back(chCol);
        _mass_changed.push_back(chmass);
        _nChange_matrixElemChange = 1;
        _matrixElemChange.push_back(ElementChange(chRow, chCol, chmass));
    } // end of if-block for atomic proposal size = 1

    else {     // atomic proposal size = 2 or above
        unsigned int count = 0;

        for (iter = _atomicProposal.begin(); iter != _atomicProposal.end(); ++ iter) {
            switch (the_matrix_label) {
                case 'A': {
                    bin = _AAtomicdomain.getBin(iter->first);
                    chRow = getRow('A', bin);
                    chCol = getCol('A', bin);
                    break;
                }

                case 'P': {
                    bin = _PAtomicdomain.getBin(iter->first);
                    chRow = getRow('P', bin);
                    chCol = getCol('P', bin);
                    break;
                }
            } // end of switch-block

            chmass = iter->second;

            if (count == 0) {   // nothing to check for the first count
                _Row_changed.push_back(chRow);
                _Col_changed.push_back(chCol);
                _mass_changed.push_back(chmass);
                count += 1;
                _nChange_matrixElemChange += 1;

            } else {
                for (unsigned int m = 0; m < count; ++m) {
                    if (chRow == _Row_changed[m] && chCol == _Col_changed[m]) {
                        _mass_changed[m] += chmass;

                        if (_mass_changed[m] == 0) {
                            _nChange_matrixElemChange -= 1;
                            _Row_changed.erase(_Row_changed.begin() + m);
                            _Col_changed.erase(_Col_changed.begin() + m);
                            _mass_changed.erase(_mass_changed.begin() + m);
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
        for (unsigned int m = 0; m < _nChange_matrixElemChange; ++m) {
            _matrixElemChange.push_back(ElementChange(_Row_changed[m], _Col_changed[m], _mass_changed[m]));
        }
    } // end of if-block for proposal size = 2 or above
} // end of extract_atomicProposal


// -----------------------------------------------------------------------------
void GibbsSampler::extract_new_atomicProposal(char the_matrix_label) {
    unsigned int bin, chRow, chCol;
    double chmass;
    map<unsigned long long, double>::const_iterator iter;
    _new_nChange_matrixElemChange = 0;
    _new_nChange_atomicProposal = _new_atomicProposal.size();

    if (_new_nChange_atomicProposal == 0) {    // atomic proposal size = 0
        _new_nChange_matrixElemChange = 0;

    } else if (_new_nChange_atomicProposal == 1) { // atomic proposal size = 1
        iter = _new_atomicProposal.begin();

        switch (the_matrix_label) {
            case 'A': {
                bin = _AAtomicdomain.getBin(iter->first);
                chRow = getRow('A', bin);
                chCol = getCol('A', bin);
                break;
            }

            case 'P': {
                bin = _PAtomicdomain.getBin(iter->first);
                chRow = getRow('P', bin);
                chCol = getCol('P', bin);
                break;
            }
        } // end of switch-block for atomic proposal size = 1

        chmass = iter->second;
        _new_Row_changed.push_back(chRow);
        _new_Col_changed.push_back(chCol);
        _new_mass_changed.push_back(chmass);
        _new_nChange_matrixElemChange = 1;
        _new_matrixElemChange.push_back(ElementChange(chRow, chCol, chmass));
    } // end of if-block for atomic proposal size = 1

    else {     // atomic proposal size = 2 or above
        unsigned int count = 0;

        for (iter = _new_atomicProposal.begin(); iter != _new_atomicProposal.end(); ++ iter) {
            switch (the_matrix_label) {
                case 'A': {
                    bin = _AAtomicdomain.getBin(iter->first);
                    chRow = getRow('A', bin);
                    chCol = getCol('A', bin);
                    break;
                }

                case 'P': {
                    bin = _PAtomicdomain.getBin(iter->first);
                    chRow = getRow('P', bin);
                    chCol = getCol('P', bin);
                    break;
                }
            } // end of switch-block

            chmass = iter->second;

            if (count == 0) {   // nothing to check for the first count
                _new_Row_changed.push_back(chRow);
                _new_Col_changed.push_back(chCol);
                _new_mass_changed.push_back(chmass);
                count += 1;
                _new_nChange_matrixElemChange += 1;

            } else {
                for (unsigned int m = 0; m < count; ++m) {
                    if (chRow == _new_Row_changed[m] && chCol == _new_Col_changed[m]) {
                        _new_mass_changed[m] += chmass;

                        if (_new_mass_changed[m] == 0) {
                            _new_nChange_matrixElemChange -= 1;
                            _new_Row_changed.erase(_new_Row_changed.begin() + m);
                            _new_Col_changed.erase(_new_Col_changed.begin() + m);
                            _new_mass_changed.erase(_new_mass_changed.begin() + m);
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
        for (unsigned int m = 0; m < _new_nChange_matrixElemChange; ++m) {
            _new_matrixElemChange.push_back(ElementChange(_new_Row_changed[m],
                                            _new_Col_changed[m], _new_mass_changed[m]));
        }
    } // end of if-block for proposal size = 2 or above
} // end of extract_new_atomicProposal



// -----------------------------------------------------------------------------
void GibbsSampler::update(char the_matrix_label)
{
    Matrix AOrig (_AMatrix);
    Matrix POrig (_PMatrix);
    bool Q_update;

    switch (the_matrix_label)
    {
        case 'A':
        {
            // ----------- making a proposal from atomic space A:
            _AAtomicdomain.makeProposal();
            get_oper_type('A');
            _atomicProposal = _AAtomicdomain.getProposedAtoms();
            extract_atomicProposal('A');

            // ----------------------------------
            // the proposal is translated into a proposal to matrix A:

            // ----------- modify the proposal in a Gibbs way:
            if (_nChange_atomicProposal == 0) {}

            if (_nChange_atomicProposal > 2) {
                throw logic_error("GibbsSampler: can't change more than two atoms!!");
            }

            // Update based on the _oper_type
            if (_oper_type == 'D') {
                Q_update = death('A', _DMatrix, _SMatrix, AOrig, POrig);

            } else if (_oper_type == 'B') {
                Q_update = birth('A', _DMatrix, _SMatrix, AOrig, POrig);

            } else if (_oper_type == 'M') {
                Q_update = move('A', _DMatrix, _SMatrix, AOrig, POrig);

            } else {
                Q_update = exchange('A', _DMatrix, _SMatrix, AOrig, POrig);
            }

            // Update the matrix with improved update, if there are further updates
            if (Q_update == true) {
                _AMatrix.elemUpdate(_new_matrixElemChange);
            }

            break;
        } // end of case 'A'

        case 'P': {
            // ----------- making a proposal from atomic space P:
            _PAtomicdomain.makeProposal();
            get_oper_type('P');
            _atomicProposal = _PAtomicdomain.getProposedAtoms();
            extract_atomicProposal('P');

            // ----------------------------------
            // the proposal is translated into a proposal to matrix P:

            // ----------- modify the proposal in a Gibbs way:
            if (_nChange_atomicProposal == 0) {}

            if (_nChange_atomicProposal > 2) {
                throw logic_error("GibbsSampler: can't chnage more than two atoms!!");
            }

            // Update based on the _oper_type
            if (_oper_type == 'D') {
                Q_update = death('P', _DMatrix, _SMatrix, AOrig, POrig);

            } else if (_oper_type == 'B') {
                Q_update = birth('P', _DMatrix, _SMatrix, AOrig, POrig);

            } else if (_oper_type == 'M') {
                Q_update = move('P', _DMatrix, _SMatrix, AOrig, POrig);

            } else {
                Q_update = exchange('P', _DMatrix, _SMatrix, AOrig, POrig);
            }

            if (Q_update == true) {
                _PMatrix.elemUpdate(_new_matrixElemChange);
            }

            break;
        } // end of case 'P'
    } // end of switch block

    // clear Proposal for the next run
    clear_Proposal();
    clear_new_Proposal();
} // end of update()

// ----------------------------------------------------------------------------
void GibbsSampler::init_sysChi2() {
    _sysChi2 = 2.*cal_logLikelihood();
}

void GibbsSampler::update_sysChi2(double delsysChi2) {
    _sysChi2 -= 2.*delsysChi2;
}

double GibbsSampler::get_sysChi2() {
    return _sysChi2;
}


// -----------------------------------------------------------------------------
void GibbsSampler::get_oper_type(char the_matrix_label) {
    switch (the_matrix_label) {
        case 'A': {
            _oper_type = _AAtomicdomain.get_oper_type();
            break;
        }

        case 'P': {
            _oper_type = _PAtomicdomain.get_oper_type();
            break;
        }
    }
}

bool GibbsSampler::death(char the_matrix_label, const Matrix &D, const Matrix &S,
const Matrix &AOrig, const Matrix &POrig)
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
    switch (the_matrix_label) {
        case 'A': {
            delLL = computeDeltaLL('A', D, S, AOrig, POrig, _nChange_matrixElemChange, _matrixElemChange);
            _AAtomicdomain.acceptProposal(false); // "false" only means not to update _iter!
            update_sysChi2(delLL);      // update system Chi2
            _AMatrix.elemUpdate(_matrixElemChange);
            break;
        }

        case 'P': {
            delLL = computeDeltaLL('P', D, S, AOrig, POrig, _nChange_matrixElemChange, _matrixElemChange);
            _PAtomicdomain.acceptProposal(false); // "false" only means not to update _iter!
            update_sysChi2(delLL);     // update system Chi2
            _PMatrix.elemUpdate(_matrixElemChange);
            break;
        }
    } // end of switch-block

    // an attempt to rebirth
    attemptMass = -origMass;

    switch (the_matrix_label) {
        case 'A': {
            // Check other matrix to see if we can use Gibbs
            if (!checkOtherMatrix('A', iRow, iCol, POrig)) {
                newMass = attemptMass;

            } else {
                newMass = getMass('A', attemptMass, iRow, iCol, POrig, AOrig, D, S, rng); //Gibbs birth

                // ------- Q: think about it
                if (newMass <= epsilon) {
                    newMass = attemptMass;
                }
            } // end of if-block

            _new_atomicProposal.insert(pair<unsigned long long, double>(location, newMass));
            extract_new_atomicProposal('A');
            _AAtomicdomain.setProposedAtomMass(_new_atomicProposal, true);
            delLLnew = computeDeltaLL('A', D, S, AOrig, POrig, _new_nChange_matrixElemChange, _new_matrixElemChange);
            break;
        } // end of switch-block for A

        case 'P': {
            // Check other matrix to see if we can use Gibbs
            if (!checkOtherMatrix('P', iRow, iCol, AOrig)) {
                newMass = attemptMass;

            } else {
                newMass = getMass('P', attemptMass, iRow, iCol, AOrig, POrig, D, S, rng); //Gibbs birth

                // ----- Q: think about it
                if (newMass <= epsilon) {
                    newMass = attemptMass;
                }
            } // end of if-block

            _new_atomicProposal.insert(pair<unsigned long long, double>(location, newMass));
            extract_new_atomicProposal('P');
            _PAtomicdomain.setProposedAtomMass(_new_atomicProposal, true);
            delLLnew = computeDeltaLL('P', D, S, AOrig, POrig, _new_nChange_matrixElemChange, _new_matrixElemChange);
            break;
        } // end of switch-block for P
    } // end of switch-block

    // M-H sampling to determine whether or not we can accept Gibbs
    if (delLLnew * _annealingTemperature  < log(Random::uniform())) {
        switch (the_matrix_label) {
            case 'A': {
                _AAtomicdomain.rejectProposal(false);
                return false;
            }

            case 'P': {
                _PAtomicdomain.rejectProposal(false);
                return false;
            }
        } // end of switch-block

    } else {
        switch (the_matrix_label) {
            case 'A': {
                _AAtomicdomain.acceptProposal(false);
                update_sysChi2(delLLnew);  // update system Chi2
                return true;
            }

            case 'P': {
                _PAtomicdomain.acceptProposal(false);
                update_sysChi2(delLLnew);  // update system Chi2
                return true;
            }
        } // end of switch-block
    } // else of if-block for M-H sampling

    return false;
}  // end of death method

bool GibbsSampler::birth(char the_matrix_label, const Matrix &D, const Matrix &S,
const Matrix &AOrig, const Matrix &POrig)
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

    switch (the_matrix_label) {
        case 'A': {
            // checking conditions for update
            if (iRow >= _nRow || iCol >= _nFactor) {
                throw logic_error("Cannot update pattern out of range in A.");
            }

            // Check other matrix to see if we can use Gibbs
            if (!checkOtherMatrix('A', iRow, iCol, POrig)) {
                _AAtomicdomain.acceptProposal(false); //accept original proposal
                delLL = computeDeltaLL('A', D, S, AOrig, POrig, _nChange_matrixElemChange, _matrixElemChange);
                update_sysChi2(delLL);  // update system Chi2
                _new_atomicProposal.insert(pair<unsigned long long, double>(location,
                                           origMass)); //update _new_atomicProposal with original change
                extract_new_atomicProposal('A');
                return true;
            }

            //Otherwise, do a Gibbs birth
            newMass = getMass('A', origMass, iRow, iCol, POrig, AOrig, D, S, rng);
            _new_atomicProposal.insert(pair<unsigned long long, double>(location, newMass));
            extract_new_atomicProposal('A');
            _AAtomicdomain.setProposedAtomMass(_new_atomicProposal, false);
            delLLnew = computeDeltaLL('A', D, S, AOrig, POrig, _new_nChange_matrixElemChange, _new_matrixElemChange);
            break;
        } // end of case-A block

        case 'P': {
            // checking conditions for update
            if (iRow >= _nFactor || iCol >= _nCol) {
                throw logic_error("Cannot update pattern out of range in P.");
            }

            // Check other matrix to see if we can use Gibbs
            if (!checkOtherMatrix('P', iRow, iCol, AOrig)) {
                _PAtomicdomain.acceptProposal(false); // accept original proposal
                delLL = computeDeltaLL('P', D, S, AOrig, POrig, _nChange_matrixElemChange, _matrixElemChange);
                update_sysChi2(delLL);  // update system Chi2
                _new_atomicProposal.insert(pair<unsigned long long, double>(location, origMass));
                extract_new_atomicProposal('P');
                return true;
            }

            //Otherwise, do a Gibbs birth
            newMass = getMass('P', origMass, iRow, iCol, AOrig, POrig, D, S, rng);
            _new_atomicProposal.insert(pair<unsigned long long, double>(location, newMass));
            extract_new_atomicProposal('P');
            _PAtomicdomain.setProposedAtomMass(_new_atomicProposal, false);
            delLLnew = computeDeltaLL('P', D, S, AOrig, POrig, _new_nChange_matrixElemChange, _new_matrixElemChange);
            break;
        } // end of case-P block
    } // end of switch-block

    // This incorporates the modified change only, if any.
    switch (the_matrix_label) {
        case 'A': {
            _AAtomicdomain.acceptProposal(false);
            update_sysChi2(delLLnew);  // update system Chi2
            return true;
        }

        case 'P': {
            _PAtomicdomain.acceptProposal(false);
            update_sysChi2(delLLnew);  // update system Chi2
            return true;
        }
    }

    return false; // should never be reached?
}  // end of method birth

bool GibbsSampler::move(char the_matrix_label, const Matrix &D, const Matrix &S,
const Matrix &AOrig, const Matrix &POrig)
{
    map<unsigned long long, double>::const_iterator atom;
    double chmass1, chmass2;
    unsigned long long loc1, loc2;
    unsigned int bin1, bin2;
    double mass1, mass2;
    double newMass1, newMass2;
    atom = _atomicProposal.begin();
    chmass1 = atom->second;
    atom++;
    chmass2 = atom->second;

    if (_atomicProposal.size() == 1) {
        return false;
    }

    // extract location, bin #, mass and changed mass corresponding to the
    // atomic proposal such that "1" refers to a positive mass change and
    // "2" a negative one.
    switch (the_matrix_label) {
        case 'A': {
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

            break;
        } // end of case 'A' block

        case 'P': {
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

            break;
        }  // end of case 'P' block
    } // end of switch-block for extracting the atomic proposal info

    // return false if bin1 == bin2
    if (bin1 == bin2) {
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

    switch (the_matrix_label) {
        case 'A': {
            lambda = _lambdaA;
            break;
        }

        case 'P': {
            lambda = _lambdaP;
            break;
        }
    }

    double delLLnew;

    switch (the_matrix_label) {
        case 'A': {
            delLLnew = computeDeltaLL('A', D, S, AOrig, POrig, _nChange_matrixElemChange, _matrixElemChange);
            break;
        }

        case 'P': {
            delLLnew = computeDeltaLL('P', D, S, AOrig, POrig, _nChange_matrixElemChange, _matrixElemChange);
            break;
        }
    }

    //EJF double totalLL = priorLL + delLLnew * _annealingTemperature;
    _new_nChange_matrixElemChange = 2;
    _new_atomicProposal.insert(pair<unsigned long long, double>(loc1, newMass1 - mass1));
    _new_atomicProposal.insert(pair<unsigned long long, double>(loc2, newMass2 - mass2));

    switch (the_matrix_label) {
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
    tmp = priorLL + delLLnew * _annealingTemperature;
    double rng = log(Random::uniform());

    if (tmp < rng) {
        switch (the_matrix_label) {
            case 'A': {
                _AAtomicdomain.rejectProposal(false);
                return false;
            }

            case 'P': {
                _PAtomicdomain.rejectProposal(false);
                return false;
            }
        } // end of switch-block

    } else {
        switch (the_matrix_label) {
            case 'A': {
                _AAtomicdomain.acceptProposal(false);
                update_sysChi2(delLLnew);  // update system Chi2
                return true;
            }

            case 'P': {
                _PAtomicdomain.acceptProposal(false);
                update_sysChi2(delLLnew);  // update system Chi2
                return true;
            }
        } // end of switch-block
    }

    // end of M-H sampling
    return false;
} // end of method move

bool GibbsSampler::exchange(char the_matrix_label, const Matrix &D,
const Matrix &S, const Matrix &AOrig, const Matrix &POrig)
{
    map<unsigned long long, double>::const_iterator atom;
    double chmass1, chmass2;
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
    switch (the_matrix_label) {
        case 'A': {
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

            break;
        } // end of case 'A' block

        case 'P': {
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

            break;
        }  // end of case 'P' block
    } // end of switch-block for extracting the atomic proposal info

    // return nullMatrix if bin1 == bin2
    if (bin1 == bin2) {
        return false;
    }

    // preparing quantities for possible Gibbs computation later.
    //EJF bool exchange = false;
    double priorLL = 0.;
    unsigned int jGene, jSample;
    bool anyNonzero = false;
    bool useGibbs = true;
    unsigned int iGene1, iPattern1, iGene2, iPattern2, iSample1, iSample2;

    switch (the_matrix_label) {
        case 'A': {
            iGene1 = getRow('A', bin1);
            iPattern1 = getCol('A', bin1);
            iGene2 = getRow('A', bin2);
            iPattern2 = getCol('A', bin2);
            break;
        }

        case 'P': {
            iPattern1 = getRow('P', bin1);
            iSample1 = getCol('P', bin1);
            iPattern2 = getRow('P', bin2);
            iSample2 = getCol('P', bin2);
            break;
        }
    }

    // ---------------------------------------------
    switch (the_matrix_label) {
        case 'A': {
            for (jSample = 0; jSample < _nCol; jSample++) {
                if (POrig(iPattern1,jSample) > epsilon) {
                    anyNonzero = true;
                    break;
                }

                if (POrig(iPattern2,jSample) > epsilon) {
                    anyNonzero = true;
                    break;
                }
            }  // end of for-block to determine the existence of corresponding

            // non-zero elements in P
            if (!anyNonzero)  {  // cannot update in Gibbs way
                useGibbs = false;
            }

            break;
        }

        case 'P': {
            for (jGene = 0; jGene < _nRow; jGene++) {
                if (AOrig(jGene,iPattern1) > epsilon) {
                    anyNonzero = true;
                    break;
                }

                if (AOrig(jGene,iPattern2) > epsilon) {
                    anyNonzero = true;
                    break;
                }
            }  // end of for-block to determine the existence of corresponding

            // non-zero elements in P
            if (!anyNonzero)  {  // cannot update in Gibbs way
                useGibbs = false;
            }

            break;
        }
    }  // end of switch-block

    // -------------------------------------------------------------------------
    // EXCHANGE ACTION when initial useGibbs = true and check if Gibbs is usable.
    double s = 0.0, su, mean, sd; // EJF -- MFO check
    pair <double, double> alphaparam;

    if (useGibbs == true) {
        switch (the_matrix_label) {
            // ---------- EXCHANGE ACTION WITH A ----------------------------
            case 'A': {
                if (_AAtomicdomain.inDomain(loc1) && _AAtomicdomain.inDomain(loc2)) {
                    alphaparam = GAPSNorm::calcAlphaParameters('A', _nFactor, D, S, AOrig, POrig, iGene1,
                                 iPattern1, iGene2, iPattern2, iSample1, iSample2);
                    s = alphaparam.first;
                    su = alphaparam.second;
                    s = s * _annealingTemperature;
                    su = su * _annealingTemperature;
                    mean = su / s;
                    sd = 1. / sqrt(s);
                    // end of compute distribution parameters for A
                } // end of if-block for checking whether the changes are in domain (the exchange block)

                break;
            } // end of switch block for EXCHANGE ACTION with A

            // ---------- EXCHANGE ACTION WITH P ----------------------------
            case 'P': {
                if (_PAtomicdomain.inDomain(loc1) && _PAtomicdomain.inDomain(loc2)) {
                    alphaparam = GAPSNorm::calcAlphaParameters('P', _nFactor, D, S, AOrig, POrig, iGene1,
                                 iPattern1, iGene2, iPattern2, iSample1, iSample2);
                    s = alphaparam.first;
                    su = alphaparam.second;
                    s = s * _annealingTemperature;
                    su = su * _annealingTemperature;
                    mean = su / s;
                    sd = 1. / sqrt(s);
                    // end of compute distribution parameters for P
                } // end of if-block for checking whether the changes are in domain (the exchange block)

                break;
            } // end of switch block for EXCHANGE ACTION with P
        }// end of switch block for EXCHANGE ACTION
    } // end of if-block for operations with possibly Gibbs sampling

    if (s == 0. && su == 0.) {
        useGibbs = false;
        //cout << "Parameters aren't updated -> useGibbs = false, do M-H" << endl;
    }

    // -------------------------------------------------------------------------
    if (useGibbs == true) {
        // set newMass1
        // need to retain exponential prior
        double plower = Random::p_norm(-mass1, mean, sd);
        double pupper = Random::p_norm(mass2, mean, sd);
        double u = plower + Random::uniform() * (pupper - plower);

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

        if (useGibbs == true) {
            gibbsMass1 = Random::q_norm(u, mean, sd);

            if (gibbsMass1 < -mass1) {
                gibbsMass1 = -mass1;
            }

            if (gibbsMass1 > mass2) {
                gibbsMass1 = mass2;
            }

            gibbsMass2 = - gibbsMass1;
            // update new masses
            double delLLnew;
            _new_nChange_matrixElemChange = 2;
            _new_atomicProposal.insert(pair<unsigned long long, double>(loc1, gibbsMass1));
            _new_atomicProposal.insert(pair<unsigned long long, double>(loc2, gibbsMass2));

            switch (the_matrix_label) {
                case 'A': {
                    extract_new_atomicProposal('A');
                    delLLnew = computeDeltaLL('A', D, S, AOrig, POrig, _new_nChange_matrixElemChange, _new_matrixElemChange);
                    _AAtomicdomain.setProposedAtomMass(_new_atomicProposal, false);
                    _AAtomicdomain.acceptProposal(false);
                    update_sysChi2(delLLnew);  // update system Chi2
                    break;
                }

                case 'P': {
                    extract_new_atomicProposal('P');
                    delLLnew = computeDeltaLL('P', D, S, AOrig, POrig, _new_nChange_matrixElemChange, _new_matrixElemChange);
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

    switch (the_matrix_label) {
        case 'A': {
            lambda = _lambdaA;
            break;
        }

        case 'P': {
            lambda = _lambdaP;
            break;
        }
    }

    // Formerly the if(_oper_type == 'E') block in move_exchange
    if (mass1 > mass2) {
        pnew = Random::d_gamma(newMass1, 2., 1. / lambda);

        if (newMass1 > newMass2) {
            pold = Random::d_gamma(mass1, 2., 1. / lambda);

        } else {
            pold = Random::d_gamma(mass2, 2., 1. / lambda);
        }

    } else {
        pnew = Random::d_gamma(newMass2, 2., 1. / lambda);

        if (newMass1 > newMass2) {
            pold = Random::d_gamma(mass1, 2., 1. / lambda);

        } else {
            pold = Random::d_gamma(mass2, 2., 1. / lambda);
        }
    }

    if (pnew == 0. && pold == 0.) {
        priorLL = 0.0;

    } else if (pnew != 0. && pold == 0.) {
        priorLL = DOUBLE_POSINF;

    } else {
        priorLL = log(pnew / pold);
    }

    double delLLnew;

    switch (the_matrix_label) {
        case 'A': {
            delLLnew = computeDeltaLL('A', D, S, AOrig, POrig, _nChange_matrixElemChange, _matrixElemChange);
            break;
        }

        case 'P': {
            delLLnew = computeDeltaLL('P', D, S, AOrig, POrig, _nChange_matrixElemChange, _matrixElemChange);
            break;
        }
    }

    //EJF double totalLL = priorLL + delLLnew * _annealingTemperature;
    _new_nChange_matrixElemChange = 2;
    _new_atomicProposal.insert(pair<unsigned long long, double>(loc1, newMass1 - mass1));
    _new_atomicProposal.insert(pair<unsigned long long, double>(loc2, newMass2 - mass2));

    switch (the_matrix_label) {
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

    if (priorLL == DOUBLE_POSINF) {
        return true;

    } else {
        tmp = priorLL + delLLnew * _annealingTemperature;
    }

    double rng = log(Random::uniform());

    if (tmp  < rng) {
        switch (the_matrix_label) {
            case 'A': {
                _AAtomicdomain.rejectProposal(false);
                return false;
            }

            case 'P': {
                _PAtomicdomain.rejectProposal(false);
                return false;
            }
        } // end of switch-block

    } else {
        switch (the_matrix_label) {
            case 'A': {
                _AAtomicdomain.acceptProposal(false);
                update_sysChi2(delLLnew);  // update system Chi2
                return true;
            }

            case 'P': {
                _PAtomicdomain.acceptProposal(false);
                update_sysChi2(delLLnew);  // update system Chi2
                return true;
            }
        } // end of switch-block
    }

    // end of M-H sampling
    return false;
} // end of method exchange


// ************ METHODS FOR LOOPING AND CONTROL *****************************
void GibbsSampler::set_iter(unsigned long ext_iter) {
    _iter = ext_iter;
}

// -----------------------------------------------------------------------------
void GibbsSampler::set_AnnealingTemperature() {
    double SASteps = _nEquil;
    double SATemp = ((double) _iter + 1.) / (SASteps / 2.);

    if (SATemp > 1.) {
        SATemp = 1;
    }

    if (SATemp < 0) {
        throw logic_error("Invalid annealing temperature.");
    }

    _annealingTemperature = SATemp;
}


// -----------------------------------------------------------------------------
//
void GibbsSampler::check_atomic_matrix_consistency(char the_matrix_label) {
    double total_atom_mass = 0.0;
    double total_matrix_mass = 0.0;

    switch (the_matrix_label) {
        case 'A': {
            total_atom_mass = _AAtomicdomain.get_atomicDomain_totalmass();
            total_matrix_mass = MatAlgo::sum(_AMatrix);
            break;
        }

        case 'P': {
            total_atom_mass = _PAtomicdomain.get_atomicDomain_totalmass();
            total_matrix_mass = MatAlgo::sum(_PMatrix);
            break;
        }
    } // end of switch-block

    double diff_total_mass = fabs(total_atom_mass - total_matrix_mass);

    if (diff_total_mass > 0.00001) { // inconsistent since atomic domain uses doubles
        throw logic_error("Mass inconsistency between atomic domain and matrix!");
    }
}

