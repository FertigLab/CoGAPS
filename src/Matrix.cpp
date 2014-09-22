
// CoGAPS C++ Version
// 
// Creation of Matrix class and functions related to it
//
// History: v 1.0  Jan 18, 2014
//          Updated to include mapped methods August 7, 2014

#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>
#include "Matrix.h"

// ------------ added for reading input file -----------------------------------
#include <sstream>
#include <string>
#include <istream>
#include <iterator>
#include <fstream>
// -----------------------------------------------------------------------------

// constants not already in main code
double epsilon = 1.0e-10;

using namespace std;
using std::vector;

// *************** CONSTRUCTORS ************************************************

Matrix::Matrix(unsigned int row_size,unsigned int col_size,
	       char the_matrix_label, double the_matrix_alpha)
{
  _n_row = row_size;
  _n_col = col_size;
  _length = row_size * col_size;
  _label = the_matrix_label;
  _alpha = the_matrix_alpha;
  // allocate memory for the constructed matrix
  _Matrix = new double * [_n_row];
  for (int m=0; m < _n_row ; ++m){
    _Matrix[m] = new double [_n_col];
  }
  // initializtion of the constructed matrix
  matrix_init();

}

Matrix::Matrix(std::vector<std::vector<double> > &vv, char the_matrix_label)
{	
    _label = the_matrix_label;
	_n_row = vv.size();
    _n_col = vv[1].size();
    _length = _n_row * _n_col;
    _Matrix = new double * [_n_row];
    for (int m=0; m < _n_row ; ++m) {
        _Matrix[m] = new double [_n_col];
    }
    for (int m=0; m < _n_row ; ++m) {
        for (int n=0; n < _n_col; ++n) {
            _Matrix[m][n] = vv[m][n];
        }
    }
}

// initialize a matrix from an input file
Matrix::Matrix(const char input_file_name[], char the_matrix_label) {
    _label = the_matrix_label;
    
    // --- read input file
    std::ifstream file;
    file.open(input_file_name);
    std::vector<std::vector<double> > vv;
    std::string line;
    while(std::getline(file, line)) {
        std::stringstream ss(line);
        std::istream_iterator<double> begin(ss), end;
        std::vector<double> v(begin, end);
        vv.push_back(v);
    }
    file.close();
    
    _n_row = vv.size();
    _n_col = vv[1].size();
    _length = _n_row * _n_col;
    _Matrix = new double * [_n_row];
    for (int m=0; m < _n_row ; ++m) {
        _Matrix[m] = new double [_n_col];
    }
    
    for (int m=0; m < _n_row ; ++m) {
        for (int n=0; n < _n_col; ++n) {
            _Matrix[m][n] = vv[m][n];
        }
    }
    
}

Matrix::~Matrix() {
    for (int m=0 ; m < _n_row ; ++m)
        delete [] _Matrix[m];
    delete[] _Matrix;
    
}

// *************** METHODS *******************************************

// create a new matrix and initiatialize it to zero
// really a repeat of a contructor method
void Matrix::born_matrix(unsigned int row_size,unsigned int col_size,
                         char the_matrix_label, double the_matrix_alpha) {
    _n_row = row_size;
    _n_col = col_size;
    _length = row_size * col_size;
    _label = the_matrix_label;
    _alpha = the_matrix_alpha;
    // allocate memory for the constructed matrix
    _Matrix = new double * [_n_row];
    for (int m=0; m < _n_row ; ++m) {
        _Matrix[m] = new double [_n_col];
    }
    
    // initializtion of the constructed matrix
    matrix_init();
    
}


// set all matrix elements to zero
void Matrix::matrix_init() {
    for(int m=0;m<_n_row;++m) {
        for (int n=0; n<_n_col;++n) {        
	  _Matrix[m][n] = 0.;  // for normal operation of Cogaps
        }
    }
}

// set a row of the matrix to passed vector 
// implemented to fix rows at a time with maps
void Matrix::setRow(vector <double> &newRow, int RowNum){
        for (int n=0; n<_n_col;++n) {        
	  _Matrix[RowNum][n] = newRow.at(n); 
        }
    }
	
// set a column of the matrix to passed vector
// implemented to fix columns at a time with maps
void Matrix::setCol(vector <double> &newCol, int ColNum){
        for (int n=0; n<_n_row;++n) {        
	  _Matrix[n][ColNum] = newCol.at(n); 
        }
    } 
	
// set a row of the matrix to passed vector
// implemented to fix rows at a time with maps
void Matrix::setRow(vector <double> const &theRow, int RowNum){
        for (int n=0; n<_n_col;++n) {        
	  _Matrix[RowNum][n] = theRow.at(n); 
        }
    }
	
// set a column of the matrix to passed vector 
// implemented to fix columns at a time with maps
void Matrix:: setCol(vector <double> const &theCol, int ColNum){
        for (int n=0; n<_n_row;++n) {        
	  _Matrix[n][ColNum] = theCol.at(n); 
        }
    } 



// these methods return the matrix or parameters
double ** Matrix::get_matrix() const {
  return _Matrix;
}

void Matrix::get_Row(int rowNum, vector <double> &theRow) const {
  for (int iCol = 0; iCol < _n_col; iCol++){
   theRow.push_back(_Matrix[rowNum][iCol]);
   }
}

void Matrix::get_Col(int colNum, vector <double> &theCol) const {
  for (int iRow=0; iRow < _n_row; iRow++){
   theCol.push_back(_Matrix[iRow][colNum]);
   }
}


unsigned int Matrix::get_nRow() const {
  return _n_row;
}

unsigned int Matrix::get_nCol() const {
  return _n_col;
}

unsigned int Matrix::get_length() const {
  return _length;
}

double Matrix::get_max_given_row(unsigned int row_indx){
  return *std::max_element(_Matrix[row_indx],_Matrix[row_indx]+_n_col);
}

double Matrix::get_min_given_row(unsigned int row_indx){
  return *std::min_element(_Matrix[row_indx],_Matrix[row_indx]+_n_col);
}

double Matrix::get_max_given_col(unsigned int col_indx){
  vector<double> col_vec(_n_row,0.0);
  for (unsigned int m=0; m<_n_row; ++m){
    col_vec[m] = _Matrix[m][col_indx];
  }
  return *std::max_element(col_vec.begin(),col_vec.end());
}

double Matrix::get_min_given_col(unsigned int col_indx){
  vector<double> col_vec(_n_row,0.0);
  for (unsigned int m=0; m<_n_row; ++m){
    col_vec[m] = _Matrix[m][col_indx];
  }
  return *std::min_element(col_vec.begin(),col_vec.end());
}


  // ********************* OTHER METHODS *****************************************
// get the mean value of all matrix elements
double Matrix::cal_mean() const {
    double mean = 0;
    for(int m=0;m < _n_row;++m) {
        for (int n=0; n < _n_col;++n) {
            mean += _Matrix[m][n];
        }
    }
    mean = mean / (_n_row*_n_col);
    return mean;
}

double Matrix::cal_totalsum() const {
    double totalsum = 0.;
    for(int m=0;m < _n_row;++m) {
        for (int n=0; n < _n_col;++n) {
            totalsum += _Matrix[m][n];
        }
    }

    return totalsum;
}


// adds a matrix to the existing matrix
void Matrix::matrix_update(vector<vector<double> > delMatrix){ 
  for (unsigned int m=0; m < _n_row; ++m){
    for (unsigned int n=0; n < _n_col; ++n){
      _Matrix[m][n] = _Matrix[m][n] + delMatrix[m][n];
      if (fabs(_Matrix[m][n]) < epsilon){
	_Matrix[m][n] = 0.0;
      }
    }
  }
}

// update the present matrix with the Element Changes (boost object)
void Matrix::matrix_Elem_update(vector<boost::tuple<unsigned int, unsigned int, double> > the_matrixElemChange,
				char oper_type, unsigned int nChange){
  unsigned int chRow, chCol;
  double delelem;
 
  for (unsigned int m=0; m < nChange; ++m) {
    chRow = the_matrixElemChange[m].get<0>();
    chCol = the_matrixElemChange[m].get<1>();
    delelem = the_matrixElemChange[m].get<2>();
    //cout << "chRow = " << chRow << endl;
    //cout << "chCol = " << chCol << endl;
    //cout << "delelem = " << delelem << endl;
     _Matrix[chRow][chCol] += delelem;
      if (fabs(_Matrix[chRow][chCol]) < epsilon){
	_Matrix[chRow][chCol] = 0.0;
      }
  } // end of for-loop for looping over ElemChange
 
} // end of method matrix_Elem_update 

// ********************* DISPLAY METHODS *****************************************
// print out a matrix on the output
// USED FOR DEBUGGING
void Matrix::display_matrix() {
    /*cout << endl;
    cout << "The matrix " << _label << " is: " << endl;
    for(int m=0;m<_n_row;++m) {
        for (int n=0; n<_n_col;++n) {
	    cout << std::setw(10) << std::right;
            cout << _Matrix[m][n] << " ";
        }
        cout << endl;
    }
    cout << endl;
    */
}

void Matrix::display_matrixF(ofstream& outputFile) {
    outputFile << endl;
    outputFile << "The matrix " << _label << " is: " << endl;
    for(int m=0;m<_n_row;++m) {
        for (int n=0; n<_n_col;++n) {
	  outputFile << std::setw(10) << std::right << _Matrix[m][n] << " ";
        }
        outputFile << endl;
    }
    outputFile << endl;
}


