#ifndef __COGAPS_MATRIX_H__
#define __COGAPS_MATRIX_H__

#include <iostream>
#include <vector>
#include <cmath>
#include <boost/tuple/tuple.hpp>
#include <algorithm>

class Matrix
{
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

    Matrix(const std::vector< std::vector<double> > &vv, char the_matrix_label);

    Matrix(const char input_file_name[], char the_matrix_label);

    ~Matrix();

// *************** METHODS *******************************************

    void born_matrix(unsigned int row_size, unsigned int col_size,
                     char the_matrix_label, double the_matrix_alpha);

    void matrix_init();

    double at(int r, int c);

    void setRow(const std::vector<double> &theRow, int RowNum);
    void setCol(const std::vector<double> &theCol, int ColNum);

    std::vector<double> get_Row(int rowNum);

    std::vector<double> get_Col(int colNum);

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

    void matrix_update(const std::vector< std::vector<double> > &delMatrix);

    void matrix_Elem_update(const std::vector< boost::tuple<unsigned int, unsigned int, double> > &the_matrixElemChange,
                            char oper_type, unsigned int nChange);

    // ********************* DISPLAY METHODS *****************************************

    void display_matrixF(std::ofstream &outputFile);

};

#endif
