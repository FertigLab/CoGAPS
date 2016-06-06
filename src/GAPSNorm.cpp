
// CoGAPS C++ Version
//
// Functions to compute matrix multiplication and
// Likelihood function
//
// History: v 1.0  Jan 16, 2014
//          Updated to include mapped methods August 7, 2014


#include <iostream>
#include <cmath>
#include <stdexcept>
#include "GAPSNorm.h"
// -- checking
#include <iomanip>
// -----

using std::vector;
using std::logic_error;
using namespace std;

namespace gaps {

// ---------------------------------------------------------------------------
// Calculation of M = A*P (matrix multiplication)
void GAPSNorm::computeMock(double **M, double const *const *A,
                           double const *const *P, unsigned int nRow,
                           unsigned int nCol, unsigned int nFactor) {
    for (unsigned int iRow = 0; iRow < nRow; iRow++) {
        for (unsigned int iCol = 0; iCol < nCol; iCol++) {
            M[iRow][iCol] = 0.;
        }
    }

    for (unsigned int iRow = 0; iRow < nRow; iRow++) {
        for (unsigned int iCol = 0; iCol < nCol; iCol++) {
            for (unsigned int iFactor = 0; iFactor < nFactor; iFactor++) {
                M[iRow][iCol] += A[iRow][iFactor] * P[iFactor][iCol] ;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Calculation of Chi2 = \sum_{i,j} (D_{ij}-(A*P)_{ij})/(S_{ij}^2)
double GAPSNorm::calChi2(double const *const *D,
                         double const *const *S,
                         double const *const *A,
                         double const *const *P,
                         unsigned int nRow, unsigned int nCol,
                         unsigned int nFactor) {
    double **Mock;
    Mock = new double * [nRow];

    for (int m = 0; m < nRow; ++m) {
        Mock[m] = new double [nCol];
    }

    GAPSNorm::computeMock(Mock, A, P, nRow, nCol, nFactor);
    double Chi2 = 0;

    for (unsigned int iRow = 0; iRow < nRow; iRow++) {
        for (unsigned int iCol = 0; iCol < nCol; iCol++) {
            Chi2 += (D[iRow][iCol] - Mock[iRow][iCol]) *
                    (D[iRow][iCol] - Mock[iRow][iCol]) /
                    (S[iRow][iCol] * S[iRow][iCol]);
        }
    }

    for (int m = 0; m < nRow; ++m) {
        delete[] Mock[m];
    }

    delete[] Mock;
    return Chi2;
}


// ---------------------------------------------------------------------------
// Calcuation of the change in the log-likelihood when matrix A (or P) is
// augmented with matrix ElemChange (~delA or delP) when there is only ONE
// change. Mathematically, it calculates changes to the following:
//
// Let the only change be delA_{mn}, define M = D-AP;
// delloglikelihood = \sum_j [2*M_{mj}*(delA_{mn}*P_{nj})-(delA_{mn}*P_{nj})^2]/2
//                          / S_{mj}^2

double GAPSNorm::calcDeltaLL1E
(char matrix_label,
 double const *const *D, double const *const *S,
 double const *const *A, double const *const *P,
 const vector<boost::tuple<unsigned int, unsigned int,
 double> >ElemChange, unsigned int nRow,
 unsigned int nCol, unsigned int nFactor) {
    // ------- Calculate changes in log-likelihood --------
    // Extract where and how large the change is
    unsigned int chRow, chCol;
    double delelem;
    double sqTerm = 0;
    double mockTerm = 0.;
    double delloglikelihood = 0.;
    chRow = ElemChange[0].get<0>();    // chRow = iRow that carries a change
    chCol = ElemChange[0].get<1>();    // chCol = iCol that carries a change
    delelem = ElemChange[0].get<2>();  // delelem = change in the element

    switch (matrix_label) {
        case 'A': {
            // ---- Form M = D - A*P, in particular, we need only M[chRow][]
            double M[1][nCol];

            for (unsigned int iCol = 0; iCol < nCol; ++ iCol) {
                M[0][iCol] = D[chRow][iCol];

                for (unsigned int iPattern = 0; iPattern < nFactor; ++iPattern) {
                    M[0][iCol] -= A[chRow][iPattern] * P[iPattern][iCol];
                }
            }

            for (unsigned int iCol = 0; iCol < nCol; ++iCol) {
                mockTerm = 2.*M[0][iCol] * delelem * P[chCol][iCol];
                sqTerm = pow(delelem * P[chCol][iCol], 2);
                delloglikelihood += (mockTerm - sqTerm) / 2. / pow(S[chRow][iCol], 2);
            }

            break;
        }

        case 'P': {
            // ---- Form M = D - A*P, in particular, we need only M[][chCol]
            double M[nRow][1];

            for (unsigned int iRow = 0; iRow < nRow; ++ iRow) {
                M[iRow][0] = D[iRow][chCol];

                for (unsigned int iPattern = 0; iPattern < nFactor; ++iPattern) {
                    M[iRow][0] -= A[iRow][iPattern] * P[iPattern][chCol];
                }
            }

            for (unsigned int iRow = 0; iRow < nRow; ++iRow) {
                mockTerm = 2.*M[iRow][0] * delelem * A[iRow][chRow];
                sqTerm = pow(delelem * A[iRow][chRow], 2);
                delloglikelihood += (mockTerm - sqTerm) / 2. / pow(S[iRow][chCol], 2);
            }

            break;
        }
    }

    return delloglikelihood;
}


// ---------------------------------------------------------------------------
// Calcuation of the change in the log-likelihood when matrix A (or P) is
// augmented with matrix ElemChange (~delA or delP) when there are TWO changes.
// Mathematically, it calculates changes to the following:
//
// Let the only change be delA_{mn} and delA_{rs}, define M = D-AP;
// Two cases:
// I) m = r
// let Del_j = delA_{mn}*P_{nj} + delA_{ms}*P_{sj}
// delloglikelihood = \sum_j [2*M_{mj}*Del_j-Del_j^2] /2/S_{mj}^2
// II) m != r
// let dellog1 = \sum_j [2*M_{mj}*(delA_{mn}*P_{nj})-(delA_{mn}*P_{nj})^2] /2
//          /S_{mj}^2
//     dellog2 = \sum_j [2*M_{rj}*(delA_{rs}*P_{sj})-(delA_{rs}*P_{sj})^2] /2
//          /S_{rj}^2
//     delloglikelihood = dellog1 + dellog2

double GAPSNorm::calcDeltaLL2E(char matrix_label,
                               double const *const *D, double const *const *S,
                               double const *const *A, double const *const *P,
                               const vector<boost::tuple<unsigned int, unsigned int, double> > ElemChange,
                               unsigned int nRow, unsigned int nCol, unsigned int nFactor) {
    // ------- Calculate changes in log-likelihood --------
    // Extract where and how much the change is.
    unsigned int chRow[2], chCol[2];
    double delelem[2];
    chRow[0] = ElemChange[0].get<0>();  // chRow[0] = iRow for the first change
    chCol[0] = ElemChange[0].get<1>();  // chCol[0] = iCol for the first change
    delelem[0] = ElemChange[0].get<2>(); // delelem[0] = first change
    chRow[1] = ElemChange[1].get<0>();  // chRow[1] = iRow for the second change
    chCol[1] = ElemChange[1].get<1>();  // chCol[1] = iCol for the second change
    delelem[1] = ElemChange[1].get<2>(); // delelem[1] = second change
    // Calculate delloglikelihood.
    double mockTerm0 = 0.; // temp variables
    double sqTerm0 = 0;    // temp variables
    double mockTerm1 = 0;  // temp variables
    double sqTerm1 = 0;    // temp variables
    unsigned int iRow, iCol, iPattern;  // loop counters
    double delloglikelihood = 0.;    // target quantity to compute

    switch (matrix_label) {
        case 'A': {
            // ---- Form M = D - A*P, in particular, we need only M[chRow[0]][]
            // and M[chRow[1]][].
            double M[2][nCol];

            for (iCol = 0; iCol < nCol; ++ iCol) {
                M[0][iCol] = D[chRow[0]][iCol];
                M[1][iCol] = D[chRow[1]][iCol];

                for (iPattern = 0; iPattern < nFactor; ++iPattern) {
                    M[0][iCol] -= A[chRow[0]][iPattern] * P[iPattern][iCol];
                    M[1][iCol] -= A[chRow[1]][iPattern] * P[iPattern][iCol];
                }
            }

            // ---- Two conditions to calculate delloglikelihood.
            if (chRow[0] == chRow[1]) {
                for (iCol = 0; iCol < nCol; ++iCol) {
                    mockTerm0 = 2.*M[0][iCol] * (delelem[0] * P[chCol[0]][iCol] +
                                                 delelem[1] * P[chCol[1]][iCol]);
                    sqTerm0 = pow((delelem[0] * P[chCol[0]][iCol] + delelem[1] * P[chCol[1]][iCol]), 2);
                    delloglikelihood += (mockTerm0 - sqTerm0) / 2. / pow(S[chRow[0]][iCol], 2);
                }

            } else {
                for (iCol = 0; iCol < nCol; ++iCol) {
                    mockTerm0 = 2.*M[0][iCol] * delelem[0] * P[chCol[0]][iCol];
                    sqTerm0 = pow(delelem[0] * P[chCol[0]][iCol], 2);
                    mockTerm1 = 2.*M[1][iCol] * delelem[1] * P[chCol[1]][iCol];
                    sqTerm1 = pow(delelem[1] * P[chCol[1]][iCol], 2);
                    delloglikelihood += (mockTerm0 - sqTerm0) / 2. / pow(S[chRow[0]][iCol], 2) +
                                        (mockTerm1 - sqTerm1) / 2. / pow(S[chRow[1]][iCol], 2);
                }
            }

            break;
        }

        case 'P': {
            // ---- Form M = D - A*P, in particular, we need only M[][chCol[0]]
            // and M[][chCol[1]].
            double M[nRow][2];

            for (iRow = 0; iRow < nRow; ++ iRow) {
                M[iRow][0] = D[iRow][chCol[0]];
                M[iRow][1] = D[iRow][chCol[1]];

                for (iPattern = 0; iPattern < nFactor; ++iPattern) {
                    M[iRow][0] -= A[iRow][iPattern] * P[iPattern][chCol[0]];
                    M[iRow][1] -= A[iRow][iPattern] * P[iPattern][chCol[1]];
                }
            }

            // ------ Two conditions to calculate the delloglikelihood.
            if (chCol[0] == chCol[1]) {
                for (iRow = 0; iRow < nRow; ++iRow) {
                    mockTerm0 = 2.*M[iRow][0] * (A[iRow][chRow[0]] * delelem[0] +
                                                 A[iRow][chRow[1]] * delelem[1]);
                    sqTerm0 = pow(A[iRow][chRow[0]] * delelem[0] + A[iRow][chRow[1]] * delelem[1], 2);
                    delloglikelihood += (mockTerm0 - sqTerm0) / 2. / pow(S[iRow][chCol[0]], 2);
                }

            } else {
                for (iRow = 0; iRow < nRow; ++iRow) {
                    mockTerm0 = 2.*M[iRow][0] * A[iRow][chRow[0]] * delelem[0];
                    sqTerm0 = pow(A[iRow][chRow[0]] * delelem[0], 2);
                    mockTerm1 = 2.*M[iRow][1] * A[iRow][chRow[1]] * delelem[1];
                    sqTerm1 = pow(A[iRow][chRow[1]] * delelem[1], 2);
                    delloglikelihood += (mockTerm0 - sqTerm0) / 2. / pow(S[iRow][chCol[0]], 2) +
                                        (mockTerm1 - sqTerm1) / 2. / pow(S[iRow][chCol[1]], 2);
                }
            }

            break;
        }
    }

    return delloglikelihood;
}

// ---------------------------------------------------------------------------
// Calcuation of the change in the log-likelihood when matrix A (or P) is augmented
// with matrix ElemChange (~delA or delP) generally.
// Mathematically, it calculates changes to the following:
//
// Let the only change be delA, define M = D-AP;
// delloglikelihood = \sum_{ij} [2*M_{ij}*(delA*P)_{ij}-(delA*P)_{ij}^2] /2/S_{ij}^2

double GAPSNorm::calcDeltaLLGen(char matrix_label,
                                double const *const *D, double const *const *S,
                                double const *const *A, double const *const *P,
                                const vector<boost::tuple
                                <unsigned int, unsigned int, double> > ElemChange,
                                unsigned int nRow, unsigned int nCol, unsigned int nFactor) {
    double delloglikelihood = 0;
    // ---- Form M = D - A*P -----
    double M[nRow][nCol];

    for (unsigned int iRow = 0; iRow < nRow; ++iRow) {
        for (unsigned int iCol = 0; iCol < nCol; ++ iCol) {
            M[iRow][iCol] = D[iRow][iCol];

            for (unsigned int iPattern = 0; iPattern < nFactor; ++iPattern) {
                M[iRow][iCol] -= A[iRow][iPattern] * P[iPattern][iCol];
            }
        }
    }

    switch (matrix_label) {
        case 'A': {
            unsigned int chRow, chCol;
            double delelem;
            unsigned int sizeChange = ElemChange.size();
            // ---- Construct delA -----
            double delA[nRow][nFactor];

            for (unsigned int m = 0; m < nRow; m++) {
                for (unsigned int n = 0; n < nFactor ; n++) {
                    delA[m][n] = 0. ;
                }
            }

            for (unsigned int m = 0; m < sizeChange; ++ m) {
                chRow = ElemChange[m].get<0>();
                chCol = ElemChange[m].get<1>();
                delelem = ElemChange[m].get<2>();
                delA[chRow][chCol] += delelem;
            }

            // ----- Compute delA*P ----
            double delAP[nRow][nCol];

            for (unsigned int iRow = 0; iRow < nRow; ++iRow) {
                for (unsigned int iCol = 0; iCol < nCol; ++iCol) {
                    delAP[iRow][iCol] = 0.;

                    for (unsigned int iPattern = 0; iPattern < nFactor; ++iPattern) {
                        delAP[iRow][iCol] += delA[iRow][iPattern] * P[iPattern][iCol];
                    }
                }
            }

            // ------ Compute delloglikelihood -------
            for (unsigned int iRow = 0; iRow < nRow; ++iRow) {
                for (unsigned int iCol = 0; iCol < nCol; ++iCol) {
                    delloglikelihood += (2.*M[iRow][iCol] * delAP[iRow][iCol] - pow(delAP[iRow][iCol], 2))
                                        / 2. / pow(S[iRow][iCol], 2);
                }
            }

            break;
        } // end of the A sub-block

        case 'P': {
            unsigned int chRow, chCol;
            double delelem;
            unsigned int sizeChange = ElemChange.size();
            // ---- Construct delP -----
            double delP[nFactor][nCol];

            for (unsigned int m = 0; m < nFactor; m++) {
                for (unsigned int n = 0; n < nCol ; n++) {
                    delP[m][n] = 0. ;
                }
            }

            for (unsigned int m = 0; m < sizeChange; ++ m) {
                chRow = ElemChange[m].get<0>();
                chCol = ElemChange[m].get<1>();
                delelem = ElemChange[m].get<2>();
                delP[chRow][chCol] += delelem;
            }

            // ----- Compute A*delP ----
            double AdelP[nRow][nCol];

            for (unsigned int iRow = 0; iRow < nRow; ++iRow) {
                for (unsigned int iCol = 0; iCol < nCol; ++iCol) {
                    AdelP[iRow][iCol] = 0.;

                    for (unsigned int iPattern = 0; iPattern < nFactor; ++iPattern) {
                        AdelP[iRow][iCol] += A[iRow][iPattern] * delP[iPattern][iCol];
                    }
                }
            }

            // ------ Compute delloglikelihood -------
            for (unsigned int iRow = 0; iRow < nRow; ++iRow) {
                for (unsigned int iCol = 0; iCol < nCol; ++iCol) {
                    delloglikelihood += (2.*M[iRow][iCol] * AdelP[iRow][iCol] - pow(AdelP[iRow][iCol], 2))
                                        / 2. / pow(S[iRow][iCol], 2);
                }
            }

            break;
        } // end of the P sub-block
    } // end of switch

    return delloglikelihood;
} // end of calcDeltaGen

// ---------------------------------------------------------------------------
// Calcuation of the change in the log-likelihood when matrix A (or P) is augmented
// with a change that occurs in a fixed pattern. This means that an entire
// column/row of A or P is changed.
// Because changing an entire column of A or row of P completely changes the matrix
// A*P, we recalculate A*P completely

double GAPSNorm::calcDeltaLLMap(char matrix_label,
                                double const *const *D, double const *const *S,
                                double const *const *A, double const *const *P,
                                vector <double> &newPat, unsigned int chPat, unsigned int nRow,
                                unsigned int nCol, unsigned int nFactor) {
    double delloglikelihood = 0;
    // ---- Form M = D - A*P -----
    double M[nRow][nCol];

    for (unsigned int iRow = 0; iRow < nRow; ++iRow) {
        for (unsigned int iCol = 0; iCol < nCol; ++ iCol) {
            M[iRow][iCol] = D[iRow][iCol];

            for (unsigned int iPattern = 0; iPattern < nFactor; ++iPattern) {
                M[iRow][iCol] -= A[iRow][iPattern] * P[iPattern][iCol];
            }
        }
    }

    switch (matrix_label) {
        case 'A': {
            // ---- Construct delA -----
            double delA[nRow][nFactor];

            for (unsigned int m = 0; m < nRow; m++) {
                for (unsigned int n = 0; n < nFactor ; n++) {
                    delA[m][n] = 0. ;
                }
            }

            // The change here is along one column, chPat:
            // Remember that the passed in vector represents
            // the actual new row, not the change, hence the
            // adjustment.
            for (unsigned int iRow = 0; iRow < nRow; ++iRow) {
                delA[iRow][chPat] += (newPat.at(iRow) - A[iRow][chPat]);
            }

            // ----- Compute delA*P ----
            double delAP[nRow][nCol];

            for (unsigned int iRow = 0; iRow < nRow; ++iRow) {
                for (unsigned int iCol = 0; iCol < nCol; ++iCol) {
                    delAP[iRow][iCol] = 0.;

                    for (unsigned int iPattern = 0; iPattern < nFactor; ++iPattern) {
                        delAP[iRow][iCol] += delA[iRow][iPattern] * P[iPattern][iCol];
                    }
                }
            }

            // ------ Compute delloglikelihood -------
            for (unsigned int iRow = 0; iRow < nRow; ++iRow) {
                for (unsigned int iCol = 0; iCol < nCol; ++iCol) {
                    delloglikelihood += (2.*M[iRow][iCol] * delAP[iRow][iCol] - pow(delAP[iRow][iCol], 2))
                                        / 2. / pow(S[iRow][iCol], 2);
                }
            }

            break;
        } // end of the A sub-block

        case 'P': {
            // ---- Construct delP -----
            double delP[nFactor][nCol];

            for (unsigned int m = 0; m < nFactor; m++) {
                for (unsigned int n = 0; n < nCol ; n++) {
                    delP[m][n] = 0. ;
                }
            }

            // The change here is along one row, chPat:
            // Remember that the passed in vector represents
            // the actual new row, not the change, hence the
            // adjustment.
            for (unsigned int iCol = 0; iCol < nCol; ++iCol) {
                delP[chPat][iCol] += (newPat.at(iCol) - P[chPat][iCol]);
            }

            // ----- Compute A*delP ----
            double AdelP[nRow][nCol];

            for (unsigned int iRow = 0; iRow < nRow; ++iRow) {
                for (unsigned int iCol = 0; iCol < nCol; ++iCol) {
                    AdelP[iRow][iCol] = 0.;

                    for (unsigned int iPattern = 0; iPattern < nFactor; ++iPattern) {
                        AdelP[iRow][iCol] += A[iRow][iPattern] * delP[iPattern][iCol];
                    }
                }
            }

            // ------ Compute delloglikelihood -------
            for (unsigned int iRow = 0; iRow < nRow; ++iRow) {
                for (unsigned int iCol = 0; iCol < nCol; ++iCol) {
                    delloglikelihood += (2.*M[iRow][iCol] * AdelP[iRow][iCol] - pow(AdelP[iRow][iCol], 2))
                                        / 2. / pow(S[iRow][iCol], 2);
                }
            }

            break;
        } // end of the P sub-block
    } // end of switch

    return delloglikelihood;
} // end of calcDeltaLLGenMap

// For move exchange / multiple pattern changes
double GAPSNorm::calcDeltaLL2Map(char matrix_label,
                                 double const *const *D, double const *const *S,
                                 double const *const *A, double const *const *P,
                                 vector <double> &newPat1, unsigned int chPat1,
                                 vector <double> &newPat2, unsigned int chPat2, unsigned int nRow,
                                 unsigned int nCol, unsigned int nFactor) {
    double delloglikelihood = 0;
    // ---- Form M = D - A*P -----
    double M[nRow][nCol];

    for (unsigned int iRow = 0; iRow < nRow; ++iRow) {
        for (unsigned int iCol = 0; iCol < nCol; ++ iCol) {
            M[iRow][iCol] = D[iRow][iCol];

            for (unsigned int iPattern = 0; iPattern < nFactor; ++iPattern) {
                M[iRow][iCol] -= A[iRow][iPattern] * P[iPattern][iCol];
            }
        }
    }

    switch (matrix_label) {
        case 'A': {
            // ---- Construct delA -----
            double delA[nRow][nFactor];

            for (unsigned int m = 0; m < nRow; m++) {
                for (unsigned int n = 0; n < nFactor ; n++) {
                    delA[m][n] = 0. ;
                }
            }

            // The change here is along 2 columns, chPat1, and chPat2:
            // Remember that the passed in vector represents
            // the actual new row, not the change, hence the
            // adjustment.
            for (unsigned int iRow = 0; iRow < nRow; ++iRow) {
                delA[iRow][chPat1] += (newPat1.at(iRow) - A[iRow][chPat1]);
                delA[iRow][chPat2] += (newPat2.at(iRow) - A[iRow][chPat2]);
            }

            // ----- Compute delA*P ----
            double delAP[nRow][nCol];

            for (unsigned int iRow = 0; iRow < nRow; ++iRow) {
                for (unsigned int iCol = 0; iCol < nCol; ++iCol) {
                    delAP[iRow][iCol] = 0.;

                    for (unsigned int iPattern = 0; iPattern < nFactor; ++iPattern) {
                        delAP[iRow][iCol] += delA[iRow][iPattern] * P[iPattern][iCol];
                    }
                }
            }

            // ------ Compute delloglikelihood -------
            for (unsigned int iRow = 0; iRow < nRow; ++iRow) {
                for (unsigned int iCol = 0; iCol < nCol; ++iCol) {
                    delloglikelihood += (2.*M[iRow][iCol] * delAP[iRow][iCol] - pow(delAP[iRow][iCol], 2))
                                        / 2. / pow(S[iRow][iCol], 2);
                }
            }

            break;
        } // end of the A sub-block

        case 'P': {
            // ---- Construct delP -----
            double delP[nFactor][nCol];

            for (unsigned int m = 0; m < nFactor; m++) {
                for (unsigned int n = 0; n < nCol ; n++) {
                    delP[m][n] = 0. ;
                }
            }

            // The change here is along one row, chPat:
            // Remember that the passed in vector represents
            // the actual new row, not the change, hence the
            // adjustment.
            for (unsigned int iCol = 0; iCol < nCol; ++iCol) {
                delP[chPat1][iCol] += (newPat1.at(iCol) - P[chPat1][iCol]);
                delP[chPat2][iCol] += (newPat2.at(iCol) - P[chPat2][iCol]);
            }

            // ----- Compute A*delP ----
            double AdelP[nRow][nCol];

            for (unsigned int iRow = 0; iRow < nRow; ++iRow) {
                for (unsigned int iCol = 0; iCol < nCol; ++iCol) {
                    AdelP[iRow][iCol] = 0.;

                    for (unsigned int iPattern = 0; iPattern < nFactor; ++iPattern) {
                        AdelP[iRow][iCol] += A[iRow][iPattern] * delP[iPattern][iCol];
                    }
                }
            }

            // ------ Compute delloglikelihood -------
            for (unsigned int iRow = 0; iRow < nRow; ++iRow) {
                for (unsigned int iCol = 0; iCol < nCol; ++iCol) {
                    delloglikelihood += (2.*M[iRow][iCol] * AdelP[iRow][iCol] - pow(AdelP[iRow][iCol], 2))
                                        / 2. / pow(S[iRow][iCol], 2);
                }
            }

            break;
        } // end of the P sub-block
    } // end of switch

    return delloglikelihood;
} // end of calcDeltaLLGenMap


// ---------------------------------------------------------------------------
// Calcuation of the parameters s and su for exchange action
// will be used GibbsSampler to determine whether or not Gibbs exchange
// action is possible. If it is, s and su are adjusted for annealingtemp
// and used to determine distribution on alpha. For full proof, see Fertig (2009)

pair<double, double> GAPSNorm:: calcAlphaParameters(char the_matrix_label, unsigned int nRow, unsigned int nCol,
        unsigned int nFactor,
        double const *const *D, double const *const *S, double **AOrig,
        double **POrig, unsigned int iGene1, unsigned int iPattern1,
        unsigned int iGene2, unsigned int iPattern2, unsigned int iSample1,
        unsigned int iSample2) {
    double s = 0.0;
    double su = 0.0;
    double mock = 0.0;
    double mock1 = 0.0;
    double mock2 = 0.0;

    switch (the_matrix_label) {
        case 'A': {
            // ---------- EXCHANGE ACTION WITH A ----------------------------
            double Aeff;

            // compute the distribution parameters
            if (iGene1 == iGene2) {
                for (int jSample = 0; jSample < nCol; jSample++) {
                    // Calculate the mock term
                    mock = D[iGene1][jSample];

                    for (int jPattern = 0; jPattern < nFactor; jPattern++) {
                        Aeff = AOrig[iGene1][jPattern];
                        mock -= Aeff * POrig[jPattern][jSample];
                    } // end of for-block that make changes to elements in A

                    s += ((POrig[iPattern1][jSample] - POrig[iPattern2][jSample]) *
                          (POrig[iPattern1][jSample] - POrig[iPattern2][jSample])) /
                         pow(S[iGene1][jSample], 2);
                    su += mock * (POrig[iPattern1][jSample] - POrig[iPattern2][jSample]) /
                          pow(S[iGene1][jSample], 2);
                }
            }  // end of case iGene1 = iGene2

            else  { // iGene1 != iGene2
                for (int jSample = 0; jSample < nCol; jSample++) {
                    mock1 = D[iGene1][jSample];
                    mock2 = D[iGene2][jSample];

                    for (int jPattern = 0; jPattern < nFactor; jPattern++) {
                        Aeff = AOrig[iGene1][jPattern];
                        mock1 -= Aeff * POrig[jPattern][jSample];
                        Aeff = AOrig[iGene2][jPattern];
                        mock2 -= Aeff * POrig[jPattern][jSample];
                    } // end of for-block for adding changes to A

                    s  += pow(POrig[iPattern1][jSample] / S[iGene1][jSample], 2) +
                          pow(POrig[iPattern2][jSample] / S[iGene2][jSample], 2);
                    su += mock1 * POrig[iPattern1][jSample] / pow(S[iGene1][jSample], 2) -
                          mock2 * POrig[iPattern2][jSample] / pow(S[iGene2][jSample], 2);
                }
            } // end of if-block for calculating s and su for case 'A'

            break;
        } // end of switch block for EXCHANGE ACTION with A

        case 'P': {
            // ---------- EXCHANGE ACTION WITH P ----------------------------
            double Peff;

            // compute the distribution parameters
            if (iSample1 == iSample2) {
                for (int jGene = 0; jGene < nRow; jGene++) {
                    // Calculate the mock term
                    mock = D[jGene][iSample1];

                    for (int jPattern = 0; jPattern < nFactor; jPattern++) {
                        Peff = POrig[jPattern][iSample1];
                        mock -=  AOrig[jGene][jPattern] * Peff;
                    } // end of for-block that make changes to elements in P

                    s += pow(((AOrig[jGene][iPattern1] - AOrig[jGene][iPattern2]) /
                              S[jGene][iSample1]), 2);
                    su += mock * (AOrig[jGene][iPattern1] - AOrig[jGene][iPattern2]) /
                          pow(S[jGene][iSample1], 2);
                }
            }  // end of case iSample1 = iSample2

            else  { // iSample1 != iSample2
                for (int jGene = 0; jGene < nRow; jGene++) {
                    mock1 = D[jGene][iSample1];
                    mock2 = D[jGene][iSample2];

                    for (int jPattern = 0; jPattern < nFactor; jPattern++) {
                        Peff = POrig[jPattern][iSample1];
                        mock1 -= AOrig[jGene][jPattern] * Peff;
                        Peff = POrig[jPattern][iSample2];
                        mock2 -= AOrig[jGene][jPattern] * Peff;
                    } // end of for-block for adding changes to P

                    s  += pow(AOrig[jGene][iPattern1] / S[jGene][iSample1], 2) +
                          pow(AOrig[jGene][iPattern2] / S[jGene][iSample2], 2);
                    su += mock1 * AOrig[jGene][iPattern1] / pow(S[jGene][iSample1], 2) -
                          mock2 * AOrig[jGene][iPattern2] / pow(S[jGene][iSample2], 2);
                }
            } // end of if-block for calculating s and su for case 'P'

            // end of compute distribution parameters for P
            break;
        } // end of switch block for EXCHANGE ACTION with P
    } // end of switch block for EXCHANGE ACTION

    return make_pair(s, su);
} // end of calcalphaparameters method


} // correspond to namespace gaps
