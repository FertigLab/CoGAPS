#ifndef GAPSNORM_H_
#define GAPSNORM_H_

#include <iostream>
#include <vector>
#include <boost/tuple/tuple.hpp>

//using namespace gaps;
using namespace std;
//using std::vector;

namespace gaps
{

  class GAPSNorm
  {

  private:
    static void computeMock(double ** M, double const * const * A, double const * const * P,
			    unsigned int nRow, unsigned int nCol, unsigned int nFactor); 

  public:
    GAPSNorm();
    ~GAPSNorm();

    void local_display_matrix(double const * const * Mat, unsigned int n_row,
				     unsigned int n_col);
 
    static double calChi2(double const * const * D, double const * const * S,
			  double const * const * A, double const * const * P,
			  unsigned int nRow, unsigned int nCol,
			  unsigned int nFactor);

    static double calcDeltaLL1E(char matrix_label,
			      double const * const * D, double const * const * S, 
			      double const * const * A, double const * const * P, 
			      const vector<boost::tuple<unsigned int, unsigned int, double> > ElemChange, 
			      unsigned int nRow, unsigned int nCol, unsigned int nFactor); 

    static double calcDeltaLL2E(char matrix_label,
			      double const * const * D, double const * const * S, 
			      double const * const * A, double const * const * P, 
			      const vector<boost::tuple<unsigned int, unsigned int, double> > ElemChange, 
			      unsigned int nRow, unsigned int nCol, unsigned int nFactor);

    static double calcDeltaLLGen(char matrix_label,
			      double const * const * D, double const * const * S, 
			      double const * const * A, double const * const * P, 
			      const vector<boost::tuple<unsigned int, unsigned int, double> > ElemChange, 
			       unsigned int nRow, unsigned int nCol, unsigned int nFactor);
	
	// For use when there are fixed patterns and an entire row is to be changed. 
	static double calcDeltaLLMap(char matrix_label,
				double const * const * D, double const * const * S, 
				double const * const * A, double const * const * P, 
				vector <double> &newPat, unsigned int chPat, unsigned int nRow, 
				unsigned int nCol, unsigned int nFactor);
				
	 // For move exchange / multiple pattern changes
    static double calcDeltaLL2Map(char matrix_label,
				double const * const * D, double const * const * S, 
				double const * const * A, double const * const * P, 
				vector <double> &newPat1, unsigned int chPat1,
				vector <double> &newPat2, unsigned int chPat2, unsigned int nRow, 
				unsigned int nCol, unsigned int nFactor);
	
	/** NEW METHOD 
	*   @short Calculate the parameters involved in exchange action for further use
    *   @return a pair <s, su> of alpha distribution parameters for further calculations. 
    */		
	static pair<double, double> calcAlphaParameters(char the_matrix_label, unsigned int nRow, unsigned int nCol, unsigned int nFactor,
						    double const * const * D, double const * const * S, double ** AOrig,
						    double ** POrig, unsigned int iGene1, unsigned int iPattern1, 
                            unsigned int iGene2, unsigned int iPattern2, unsigned int iSample1, 
							unsigned int iSample2);
      };
}
#endif 
