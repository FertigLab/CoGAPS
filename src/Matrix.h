// Matrix.h
//
// Ref: For reading in files, see:
// http://www.codeproject.com/Articles/21909/Introduction-to-dynamic-two-dimensional-arrays-in

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <boost/tuple/tuple.hpp>
#include <algorithm>

using namespace std;

class Matrix {
  protected:
    unsigned int _n_row;
    unsigned int _n_col;
    unsigned int _length;
    char _label;
    double _alpha;
    double _lambda;
    // ---- added to bridge to Atomic Class ----------
    double **_Matrix;

  public:
    Matrix() {};

    Matrix(unsigned int row_size, unsigned int col_size, char the_matrix_label,
           double the_matrix_alpha);

    Matrix(std::vector<std::vector<double> > &vv, char the_matrix_label);

    Matrix(const char input_file_name[], char the_matrix_label);

    ~Matrix();

// *************** METHODS *******************************************

    void born_matrix(unsigned int row_size, unsigned int col_size,
                     char the_matrix_label, double the_matrix_alpha);

    void matrix_init();


    void setRow(vector <double> theRow, int RowNum);

    void setCol(vector <double> theCol, int ColNum);

    double **get_matrix() const;

    vector <double> get_Row(int rowNum);

    vector <double> get_Col(int colNum);

    unsigned int get_nRow() const;

    unsigned int get_nCol() const;

    unsigned int get_length() const;

    double get_max_given_row(unsigned int row_indx);

    double get_min_given_row(unsigned int row_indx);

    double get_max_given_col(unsigned int col_indx);

    double get_min_given_col(unsigned int col_indx);

    // ********************* OTHER METHODS *****************************************
    // cal_mean calculates the mean of a matrix over all its entries.
    double cal_mean() const;

    double cal_totalsum() const;

    void matrix_update(vector<vector<double> > delMatrix);

    void matrix_Elem_update(vector<boost::tuple<unsigned int, unsigned int, double> > the_matrixElemChange,
                            char oper_type, unsigned int nChange);

    // ********************* DISPLAY METHODS *****************************************

    void display_matrixF(ofstream &outputFile);

};

#endif
