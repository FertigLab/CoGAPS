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

void AtomicSupport::printAtomicInfoF(ofstream &outputFile) {
    outputFile << endl;
    outputFile << "nBin: " << _nBin << endl;
    outputFile << "NatomLength: " << _NatomLength << endl;
    outputFile << "lambda: " << _lambda << endl;
    outputFile << "alpha: " << _alpha << endl;
    outputFile << "atomic_domain_label: " << _atomic_domain_label << endl;
    outputFile << "atomDomain = [";
    unsigned int iBin;

    for (iBin = 0; iBin < _nBin; iBin++) {
        outputFile << "[" << iBin << " " << _lBoundariesByBin[iBin] << " " << getNumAtoms(iBin) << " " << getTotalMass(
                       iBin) << "];";
    }

    outputFile << "];" << endl;
    outputFile << endl;
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

    if (_AtomicDomain.size() == 0) {
        return 0;
    }

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

    if (_AtomicDomain.size() == 0) {
        return 0;
    }

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

unsigned long long AtomicSupport::getStartLocation(unsigned int iBin) {
    return _lBoundariesByBin[iBin];
}

unsigned long long AtomicSupport::getEndLocation(unsigned int iBin) {
    if (iBin < _nBin - 1) {
        return _lBoundariesByBin[iBin + 1];

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

    return randgen('P', (double)(_nAtom + delAtom), (double) getExpectedNAtom());
}

void AtomicSupport::resetAtomicOutputThin(int thinDiag) {
    thinAtomicDiag = thinDiag;

    if (thinDiag > 0) {
        outputAtomicDiag = true;
        _initIterOutput = _iter;
    }

    atomicDiagFile << "Reset output to start at iteration " << _initIterOutput << " at every "
                   << thinDiag << " iterations." << endl;
}

void AtomicSupport::writeAtomicHeader(char diagnosticFileName[],
                                      int thinDiag) {
    thinAtomicDiag = thinDiag;
    outputAtomicDiag = thinDiag >= 1;
    _initIterOutput = _iter;
    atomicDiagFile.open(diagnosticFileName, ios::out);
    // get the current time for output
    time_t rawtime;
    struct tm *timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    atomicDiagFile << "Atomic diagnostic file created at\t" << asctime(timeinfo) << endl;
    // output the model parameters used
    atomicDiagFile << "number of bins\t" << _nBin << endl;
    atomicDiagFile << "maximum number of atoms\t" << _NatomLength << endl;
    atomicDiagFile << "expected number of atoms per bin (alpha)\t" << _alpha << endl;
    atomicDiagFile << "expected mass of an atom (lambda)\t" << _lambda << endl;
    atomicDiagFile << endl;
}

void AtomicSupport::writeAtomicDiagnostics() {
    if (!outputAtomicDiag) {
        return;
    }

    // output results to a file
    if (((_iter - _initIterOutput) % thinAtomicDiag) !=  0) {
        return;
    }

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
                outputValue = (float) massPerBin[iCol + iBin * nCol];

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
                outputValue = (float) massPerBin[iRow + iBin * nRow];

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

// --------- to extract _atomic_domain_label
char AtomicSupport::get_atomic_domain_label() {
    return _atomic_domain_label;
}

}
