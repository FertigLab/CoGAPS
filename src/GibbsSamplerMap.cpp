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

// -----------------------------------------------------------------------------
  const double DOUBLE_POSINF = std::numeric_limits<double>::max();
  const double DOUBLE_NEGINF = -std::numeric_limits<double>::max();
  const double epsilon = 1e-10;
// -----------------------------------------------------------------------------



// ******************** CONSTRUCTOR ********************************************
 GibbsSamplerMap::GibbsSamplerMap(unsigned long nEquil, unsigned long nSample, unsigned int nFactor, 
	       double alphaA, double alphaP, double nMaxA, double nMaxP,
	       unsigned long nIterA, unsigned long nIterP, 
	       double max_gibbsmass_paraA, double max_gibbsmass_paraP, 
	       double lambdaA_scale_factor, double lambdaP_scale_factor, 
               unsigned long atomicSize,
	       char label_A,char label_P,char label_D,char label_S,
	       const string & datafile, const string & variancefile,
               const string & simulation_id,  
			   vector <vector <double> >  &parameters, char the_fixed_matrix)
		  // Initialize old parameters
		 : GibbsSampler(nEquil, nSample, nFactor, alphaA, alphaP, nMaxA, nMaxP,
		  nIterA, nIterP, max_gibbsmass_paraA, max_gibbsmass_paraP,
          lambdaA_scale_factor, lambdaP_scale_factor, atomicSize, 
		  label_A, label_P, label_D, label_S, datafile, variancefile,
          simulation_id){
		 //Constructor body
		 _MapValues = parameters;
		 _nFixedMaps = parameters.size();
		 _the_fixed_matrix = the_fixed_matrix;
		}
						
			   
 GibbsSamplerMap::GibbsSamplerMap(unsigned long nEquil, unsigned long nSample, unsigned int nFactor, 
	       double alphaA, double alphaP, double nMaxA, double nMaxP,
	       unsigned long nIterA, unsigned long nIterP, 
	       double max_gibbsmass_paraA, double max_gibbsmass_paraP, 
	       double lambdaA_scale_factor, double lambdaP_scale_factor, 
               unsigned long atomicSize,
	       char label_A,char label_P,char label_D,char label_S,
	       vector<vector<double> > &DVector, vector<vector<double> > &SVector,
               const string & simulation_id, 
			   vector <vector <double> >  &parameters, char the_fixed_matrix)
		 : GibbsSampler(nEquil, nSample, nFactor, 
		  alphaA, alphaP, nMaxA, nMaxP, nIterA, nIterP, max_gibbsmass_paraA, 
		  max_gibbsmass_paraP,lambdaA_scale_factor, lambdaP_scale_factor,
          atomicSize, label_A, label_P, label_D, label_S, DVector, SVector, simulation_id){
		  //Constructor body
		  _MapValues = parameters;
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
 void GibbsSamplerMap::init_Mapped_Matrix(){
  vector <double> thefixedPat;

     switch (_the_fixed_matrix){
    case 'A':{ 
	for (int iCol = 0; iCol < _nFixedMaps; iCol++){
	  thefixedPat = _MapValues.at(iCol);
	  _AMatrix.setCol(thefixedPat, iCol); 
	  }
	 break; }
	case 'P':{
	for (int iRow = 0; iRow < _nFixedMaps; iRow++){
	  thefixedPat = _MapValues.at(iRow);
	  _PMatrix.setRow(thefixedPat, (_nFactor - _nFixedMaps + iRow));
	  }
	 break;}
    } // end switch block
  } //end method
  
 // Keep the atomic domain consistent with the matrix to avoid
 // a mass inconsistency. Get the row or column of the fixed pattern
 // and birth new atoms across it in the middle of the correct bin.   
  void GibbsSamplerMap::initialize_atomic_domain_map(){
    map <unsigned long long, double > getNSync;
	unsigned long long updateloc;
	unsigned int updatebin;
	double updatemass;
    vector <double> fixPat;

    switch (_the_fixed_matrix){
    case 'A':{
	// For A, we start in the 1st column at bin 0. Bins count horizontally
	// so we increment the bin by _nFactor to go down a column. 
	for (int iCol = 0; iCol < _nFixedMaps; iCol++){
	fixPat = _MapValues.at(iCol);
	//reset updatebin
	 updatebin=iCol;
	 for (int iAtom = 0; iAtom < _nRow; iAtom++){
	 // Place the atom in the middle of the bin
	 updateloc = _AAtomicdomain.getMidLocation(updatebin);
	 updatemass = fixPat.at(iAtom);
	 getNSync.insert(pair<unsigned long long,double>(updateloc,updatemass));
	 updatebin+=_nFactor;
	  }
	 }
	 //Proposal is now finished. Accept it to update the atomic domain. 
	 _AAtomicdomain.setProposedAtomMass(getNSync,true);  
	 _AAtomicdomain.acceptProposal(false); 
	 break;
	 } // end case A
	 case 'P':{
	// For P, we start with row nFactor (number of rows in P) - nFixedMaps
    // and go down to the end of the matrix. We increment the bin numbers
	// which here count vertically, by nFactor. 
	for (int iRow = 0; iRow < _nFixedMaps; iRow++){
	fixPat = _MapValues.at(iRow);
	//reset updatebin
	updatebin=_nFactor - _nFixedMaps + iRow;
	 for (int iAtom = 0; iAtom < _nCol; iAtom++){
	 // Place the atom in the middle of the bin
	 updateloc = _PAtomicdomain.getMidLocation(updatebin);
	 updatemass = fixPat.at(iAtom);
	 getNSync.insert(pair<unsigned long long,double>(updateloc,updatemass));
	 updatebin+=_nFactor;
	  }
	 }
	 //Proposal is now finished. Accept it to update the atomic domain. 
	 _PAtomicdomain.setProposedAtomMass(getNSync,true);  
	 _PAtomicdomain.acceptProposal(false); 
	 break;
	 } // end case P
	} // end switch block
	
 } // end initialize_atomic_domain_map
	 
	// Just to see the matrices in the terminal
	 void GibbsSamplerMap:: print_A_and_P(){
	  _AMatrix.display_matrix();
	  _PMatrix.display_matrix();
	  }

  
  //---------------------------------------------------------------------------
  
  
// ************* METHODS FOR COMPUTING LIKELIHOOD FUNCTIONS ******************

// For one pattern change (in birth and death) 
double GibbsSamplerMap::computeDeltaLLBDMap(char the_matrix_label,
				    double const * const * D,
				    double const * const * S,
				    double const * const * A,
				    double const * const * P,
				    vector <double> &newPat, unsigned int chPat){

  double DelLL;

 
  switch(the_matrix_label){
  case 'A':
    { 
	
	DelLL = GAPSNorm::calcDeltaLLMap('A',D,S,A,P,newPat,chPat,_nRow,
					 _nCol,_nFactor); 
					 
      break;} // end of switch-block 'A'
	  
  case 'P':
    { 

	DelLL = GAPSNorm::calcDeltaLLMap('P',D,S,A,P,newPat, chPat,_nRow,
					 _nCol,_nFactor);

      break;} // end of switch-block 'P' 

  } // end of switch block

  return DelLL;

} // end of computeDeltaLLBDMap

// For two pattern changes (in move and exchange)
double GibbsSamplerMap::computeDeltaLLMEMap(char the_matrix_label,
				    double const * const * D,
				    double const * const * S,
				    double const * const * A,
				    double const * const * P,
				    vector <double> &newPat1, unsigned int chPat1,
					vector <double> &newPat2, unsigned int chPat2){

  double DelLL;

 
  switch(the_matrix_label){
  case 'A':
    { 
	
	DelLL = GAPSNorm::calcDeltaLL2Map('A',D,S,A,P,newPat1,chPat1,
					 newPat2, chPat2, _nRow, _nCol,_nFactor); 
					 
      break;} // end of switch-block 'A'
	  
  case 'P':
    { 

	DelLL = GAPSNorm::calcDeltaLL2Map('P',D,S,A,P,newPat1,chPat1,
					 newPat2, chPat2, _nRow, _nCol,_nFactor); 

      break;} // end of switch-block 'P' 

  } // end of switch block

  return DelLL;

} // end of computeDeltaLLMap


// ******************* METHODS FOR THE UPDATE/PROPOSAL *************************
// -----------------------------------------------------------------------------
// For the modified update method, need to check the location of the atomic 
// proposal. If the proposal falls in a fixed area, distribute the mass of the 
// proposal across the row such that the percentage of total mass in each fixed
// bin remains the same. Otherwise continue normally. 
void GibbsSamplerMap::mapUpdate(char the_matrix_label){

  // Send directly to the regular Gibbs and break from this method
  // if the matrix differs from the fixed one
  if (_the_fixed_matrix != the_matrix_label){
   update(the_matrix_label); // (GibbsSampler update method)
   return;
  }

  double rng = 0.1; // no use, filling up the list
  double ** D = _DMatrix.get_matrix();
  double ** S = _SMatrix.get_matrix();
  double ** AOrig = _AMatrix.get_matrix();
  double ** POrig = _PMatrix.get_matrix();
  vector<vector<double> > del_matrix;
  map<unsigned long long, double>::iterator iter;

  bool Q_update;

  switch(the_matrix_label){
  case 'A':
    {
      // ----------- making a proposal from atomic space A:
      _AAtomicdomain.makeProposal(rng);
      get_oper_type('A');
      _atomicProposal = _AAtomicdomain.getProposedAtoms();
      extract_atomicProposal('A');
	  
	  // Check to make sure there are no problems with the proposal
	
       if ( _nChange_atomicProposal> 2){
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
	  if (_nChange_atomicProposal == 1){
	    if (fixed1){
		 if (_oper_type == 'D'){
		  //cout << "Starting mapped death" << endl;
		  mappedDeath('A',D,S,AOrig,POrig);
         }
         else {
		  //cout << "Starting mapped birth" << endl;
	      mappedBirth('A',D,S,AOrig,POrig);
		 } 
		}//end if (fixed pattern) block
		else { // not a fixed pattern, call the methods from Gibbs 
		 if (_oper_type == 'D'){
		  Q_update = death('A',D,S,AOrig,POrig);
         } else {
	      Q_update = birth('A',D,S,AOrig,POrig);
		 }
		 // Modify the proposal in normal Gibbs manner
	     // Update the matrix with improved update, if there are further updates
	     if (Q_update == true){
	      _AMatrix.matrix_Elem_update(_new_matrixElemChange,_oper_type,_new_nChange_matrixElemChange);
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
	   if (fixed1 && fixed2){
	     if (_oper_type == 'M'){
		  //cout << "Starting mapped move" << endl;
	      mappedMove('A',D,S,AOrig,POrig);
         } else {
		  //cout << "Starting mapped exchange" << endl;	     
		  mappedExchange('A',D,S,AOrig,POrig);
		 }	  
		} // end if one or both locations in a fixed pattern 
		else if (!fixed1 && !fixed2) { //both locations not in fixed pattern, normal Gibbs */
	     if (_oper_type == 'M') {
		  //cout << "Starting regular move" << endl;
	      Q_update = move('A',D,S,AOrig,POrig);
         } else {
		  //cout << "Starting regular exchange" << endl;
	      Q_update = exchange('A',D,S,AOrig,POrig);
		 }
		 // Modify the proposal in normal Gibbs manner
	     // Update the matrix with improved update, if there are further updates
	     if (Q_update == true){
	      _AMatrix.matrix_Elem_update(_new_matrixElemChange,_oper_type,_new_nChange_matrixElemChange);
	     }
		}
       } //end Move/Exchange cases
	
      break;
    } // end of case 'A'
  case 'P':
    {
      // ----------- making a proposal from atomic space P:
      _PAtomicdomain.makeProposal(rng);
      get_oper_type('P');
      _atomicProposal = _PAtomicdomain.getProposedAtoms();
      extract_atomicProposal('P');
	  
	  // Check to make sure there are no problems with the proposal
	  if ( _nChange_atomicProposal == 0){}
      if ( _nChange_atomicProposal> 2){
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
	  if (_nChange_atomicProposal == 1){
	    if (fixed1){
		 if (_oper_type == 'D'){
		  //cout << "Starting mapped death" << endl;
		  mappedDeath('P',D,S,AOrig,POrig);
         } else {
		  //cout << "Starting mapped birth" << endl;
	      mappedBirth('P',D,S,AOrig,POrig);
		 }
		}//end if (fixed pattern) block
		else { // not a fixed pattern, call the methods from Gibbs
		 if (_oper_type == 'D'){
		  Q_update = death('P',D,S,AOrig,POrig);
         } else {
	      Q_update = birth('P',D,S,AOrig,POrig);
		 }
		 // Modify the proposal in normal Gibbs manner
	     // Update the matrix with improved update, if there are further updates
	     if (Q_update == true){
	      _PMatrix.matrix_Elem_update(_new_matrixElemChange,_oper_type,_new_nChange_matrixElemChange);
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
	   if (fixed1 && fixed2){
	     if (_oper_type == 'M'){
		  //cout << "Starting mapped move " << endl;
	      mappedMove('P',D,S,AOrig,POrig);
         } else {
		  //cout << "Starting mapped exchange " << endl;
	      mappedExchange('P',D,S,AOrig,POrig);
		 }
		} // end if one or both locations in a fixed pattern 
		else if (!fixed1 && !fixed2) { //both locations not in fixed pattern, normal Gibbs
	     if (_oper_type == 'M') {
	      Q_update = move('P',D,S,AOrig,POrig);
         } else {
	      Q_update = exchange('P',D,S,AOrig,POrig);
		 }
		 // Modify the proposal in normal Gibbs manner
	     // Update the matrix with improved update, if there are further updates
	     if (Q_update == true){
	      _PMatrix.matrix_Elem_update(_new_matrixElemChange,_oper_type,_new_nChange_matrixElemChange);
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
                                       char the_matrix_label){
									   
   unsigned int theBin, theRow, theCol;
   // If the matrix being tested isn't the fixed matrix, stop
   if (the_matrix_label != _the_fixed_matrix){
   return false;
   }
   switch(the_matrix_label){
    case 'A':
	  { 
	    theBin = _AAtomicdomain.getBin(location);
		theCol = getCol('A', theBin);
		// if the change is in rows 0 to nFixedMaps of A,
		// the change is in a fixed row - return true
		if (theCol < _nFixedMaps){
		return true;
		}
		break;
		} //end case A
    case 'P':
	  {
	    theBin = _PAtomicdomain.getBin(location);
		theRow = getRow('P', theBin);
		// if the change is in rows _nFactor - _nFixedMaps to
		// _nFactor of P,  the change is in a fixed row - return true
		if (theRow >= (_nFactor - _nFixedMaps)){
		return true;
		}
		break;
		} //end case P
	} //end switch block
	return false;
}

//----------------------------------------------------------------

// Kill an atom in a fixed pattern. This method
// kills the atom from the original atomic_proposal in the fixed
// pattern fashion. Then it rebirths the atom in the same manner
// if appropriate (determined by MH). 
void GibbsSamplerMap::mappedDeath(char the_matrix_label,
						  double const * const * D,
						  double const * const * S,
						  double ** AOrig,
						  double ** POrig)
{

  double rng = 0.1; // no use, just fill up the list
  double newMass = 0;
  double attemptMass = 0;
  // read in the original _atomicProposal made from the prior
  unsigned long long location = _atomicProposal.begin()->first;
  double origMass = _atomicProposal.begin()->second;
  double delLL = 0;
  double delLLnew = 0;
  bool makeChange;
  unsigned int chBin, chRow, chCol;
  vector <double> proposed_newPat; 

    // Kill the atom in both the atomic domain and matrix based on 
	// the initial proposal
    switch(the_matrix_label){
    case 'A':
      { 
	    chBin = _AAtomicdomain.getBin(location);
		chRow = getRow('A', chBin);
		chCol = getCol('A', chBin);
        makeChange = calc_new_matrix_Pattern('A',proposed_newPat,location, origMass); 
		if (makeChange){
        delLL = computeDeltaLLBDMap('A',D,S,AOrig,POrig,proposed_newPat,chCol); 
        _AAtomicdomain.acceptProposal(false); // "false" only means not to update _iter!
	    update_sysChi2(delLL); 	    // update system Chi2 
	    update_fixed_pattern('A',proposed_newPat,chCol);
		}
	    break;
		}
    case 'P':
      { 
		chBin = _PAtomicdomain.getBin(location);
		chRow = getRow('P', chBin);
		chCol = getCol('P', chBin);
		makeChange = calc_new_matrix_Pattern('P',proposed_newPat,location, origMass);
		if (makeChange){
	    delLL = computeDeltaLLBDMap('P',D,S,AOrig,POrig,proposed_newPat,chRow);
        _PAtomicdomain.acceptProposal(false); // "false" only means not to update _iter!
	    update_sysChi2(delLL);     // update system Chi2
	    update_fixed_pattern('P', proposed_newPat,chRow);
		}
	    break;
		}
    } // end of switch-block
    
    // an attempt to rebirth
	
	// Can't attempt to rebirth if the original atom was never killed
	if (!makeChange){
	 return;
	 }
	 
	vector <double> newestPat;
	bool makeNewChange; // no use. Because there is rebirth, there is nothing to check because
	                    // we will never be lowering any masses 
    attemptMass = -origMass;
    switch(the_matrix_label){
    case 'A':
      {
	    // Check other matrix to see if we can use Gibbs
        if (!checkOtherMatrix('A',chRow, chCol, POrig)) { 
	      newMass = attemptMass;
	    } else {
	      newMass = getMass('A',attemptMass,chRow,chCol,POrig,AOrig,D,S,rng); //Gibbs birth
	      // ------- Q: think about it
	      if (newMass <= epsilon) {
	       newMass = attemptMass;
	      }
	    } // end of if-block 

	    _new_atomicProposal.insert(pair<unsigned long long,double>(location,newMass));
	    extract_new_atomicProposal('A');
	    _AAtomicdomain.setProposedAtomMass(_new_atomicProposal,true);  
		makeNewChange = calc_new_matrix_Pattern('A', newestPat, location, newMass);
        delLLnew = computeDeltaLLBDMap('A',D,S,AOrig,POrig,newestPat,chCol);
	    break;
      } // end of switch-block for A
    case 'P':
      {
	    // Check other matrix to see if we can use Gibbs
	    if (!checkOtherMatrix('P',chRow, chCol, AOrig)) { 
	      newMass = attemptMass;
	     } else {
          newMass = getMass('P',attemptMass,chRow,chCol,AOrig,POrig,D,S,rng); //Gibbs birth
	      // ----- Q: think about it
	      if (newMass <= epsilon) {
	       newMass = attemptMass;
	      }
	    } // end of if-block 

	    _new_atomicProposal.insert(pair<unsigned long long,double>(location,newMass));
	    extract_new_atomicProposal('P');
	    _PAtomicdomain.setProposedAtomMass(_new_atomicProposal,true);
        makeNewChange = calc_new_matrix_Pattern('P', newestPat, location, newMass);
        delLLnew = computeDeltaLLBDMap('P',D,S,AOrig,POrig,newestPat,chRow);
	    break;
      } // end of switch-block for P

    } // end of switch-block

    // M-H sampling to determine whether or not we can accept Gibbs
    if (delLLnew*_annealingTemperature  < log(randgen('U',0,0))) {

     switch(the_matrix_label){
      case 'A':
	  { 
	    _AAtomicdomain.rejectProposal(false);
	    break;
      }
      case 'P':
	  { 
	    _PAtomicdomain.rejectProposal(false);
	    break;
	  }
     } // end of switch-block
	  
    } else {
      switch(the_matrix_label){
      case 'A':
	  { 
          _AAtomicdomain.acceptProposal(false); 
	      update_sysChi2(delLLnew);  // update system Chi2
	      update_fixed_pattern('A',newestPat,chCol);
          break;
	   }
      case 'P':
	  {   
	      _PAtomicdomain.acceptProposal(false);  
	      update_sysChi2(delLLnew);  // update system Chi2
	      update_fixed_pattern('P',newestPat,chRow);
          break;
	   }
      } // end of switch-block
    } // else of if-block for M-H sampling

}  // end of death method

void GibbsSamplerMap::mappedBirth(char the_matrix_label,
						  double const * const * D,
						  double const * const * S,
						  double ** AOrig,
						  double ** POrig)
{
  double rng = 0.1; 
  double newMass = 0;
  // EJF double attemptMass = 0;

  // read in the original _atomicProposal made from the prior
  unsigned long long location = _atomicProposal.begin()->first;
  double origMass = _atomicProposal.begin()->second;
  unsigned int chBin, chRow, chCol;
  vector <double> newestPat;
  double delLL = 0;
  double delLLnew = 0;
  bool makeChange; // again, because we are birthing this just fills up list
   
    switch(the_matrix_label){
    case 'A':
      {
	   // find the row and column of the change
	    chBin = _AAtomicdomain.getBin(location);
		chRow = getRow('A', chBin);
		chCol = getCol('A', chBin);
		
	   // checking conditions for update
	   if (chRow >= _nRow || chCol >= _nFactor) {
	    throw logic_error("Cannot update pattern out of range in A.");
	   }
		
	// Check other matrix to see if we can use Gibbs
	if ( !checkOtherMatrix('A', chRow, chCol, POrig)) {
	
	  _AAtomicdomain.acceptProposal(false); //accept original proposal 
	  makeChange = calc_new_matrix_Pattern('A', newestPat,location, origMass);
      delLL = computeDeltaLLBDMap('A',D,S,AOrig,POrig,newestPat,chCol);
	  update_sysChi2(delLL);  // update system Chi2
	  //Update the matrix:
	  update_fixed_pattern('A',newestPat,chCol);
	  } 
	  else {
	 //Otherwise, do a Gibbs birth
     newMass = getMass('A',origMass,chRow,chCol,POrig,AOrig,D,S,rng);
	 
	 // Update the atomic proposal
	 _new_atomicProposal.insert(pair<unsigned long long,double>(location,newMass));
  	 extract_new_atomicProposal('A');
	 _AAtomicdomain.setProposedAtomMass(_new_atomicProposal,false);
	 
	 // Update the matrix
   	 makeChange = calc_new_matrix_Pattern('A', newestPat, location, newMass);
     delLLnew = computeDeltaLLBDMap('A',D,S,AOrig,POrig,newestPat, chCol);
	 update_fixed_pattern('A',newestPat,chCol);
	 
	 // Update the atomic space and system Chi2
     _AAtomicdomain.acceptProposal(false);
	 update_sysChi2(delLLnew); 
	  }

	break;
     } // end of case-A block
	 
    case 'P':
      {
	  	// find the row and column of the change
	    chBin = _PAtomicdomain.getBin(location);
		chRow = getRow('P', chBin);
		chCol = getCol('P', chBin);
		
	   // checking conditions for update
	   if (chRow >= _nFactor || chCol >= _nCol) {
	    throw logic_error("Cannot update pattern out of range in P.");
	   }
		
	   // Check other matrix to see if we can use Gibbs
	   if ( !checkOtherMatrix('P', chRow, chCol, AOrig)) {
	    
		_PAtomicdomain.acceptProposal(false); // accept original proposal
	    makeChange = calc_new_matrix_Pattern('P', newestPat, location, origMass);
        delLL = computeDeltaLLBDMap('P',D,S,AOrig,POrig,newestPat,chRow);
  	    update_sysChi2(delLL);  // update system Chi2
	    update_fixed_pattern('P',newestPat,chRow); //update matrix

	   } else {
		
	   	 //Otherwise, do a Gibbs birth
        newMass = getMass('P',origMass,chRow,chCol,AOrig,POrig,D,S,rng);
		
		// Update the atomic proposal
	    _new_atomicProposal.insert(pair<unsigned long long,double>(location,newMass));
  	    extract_new_atomicProposal('P');
	    _PAtomicdomain.setProposedAtomMass(_new_atomicProposal,false);
		
		// Update the matrix
		makeChange = calc_new_matrix_Pattern('P', newestPat, location, newMass);
		delLLnew = computeDeltaLLBDMap('P',D,S,AOrig,POrig,newestPat, chRow);
		update_fixed_pattern('P',newestPat,chRow);
	  
	   // Update atomic space and system Chi2
	   _PAtomicdomain.acceptProposal(false);
	   update_sysChi2(delLLnew); 
		}

	break;
      } // end of case-P block
    } // end of switch-block

}  // end of method birth 

// Mapped move when both atoms are in a fixed pattern. 
// Accepts or rejects move using MH and changes both atomic
// space and matrix. 
void GibbsSamplerMap::mappedMove(char the_matrix_label,
						    double const * const * D,
						    double const * const * S,
						    double ** AOrig,
						    double ** POrig)
{
  map<unsigned long long, double>::const_iterator atom;
  double chmass1,chmass2;       
  unsigned long long loc1, loc2;
  unsigned int bin1, bin2;
  unsigned int pat1, pat2;
  double mass1, mass2;
  double newMass1, newMass2;
  atom = _atomicProposal.begin();
  chmass1 = atom->second;
  atom++;
  chmass2 = atom->second;
  if (_atomicProposal.size()==1){
  //cout << "Not doing a move due to update inconsistency."<< endl;
  return;
  }
 
  // extract location, bin #, mass and changed mass corresponding to the 
  // atomic proposal such that "1" refers to a positive mass change and 
  // "2" a negative one.  
  switch(the_matrix_label){
  case 'A':
    {
      if (chmass1 > chmass2) {
	atom--;
	loc1 = atom->first;
	bin1 = _AAtomicdomain.getBin(loc1);
	pat1 = getCol('A',bin1);
	mass1 = _AAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom++;
	loc2 = atom->first;
	bin2 = _AAtomicdomain.getBin(loc2);
	pat2 = getCol('A',bin2);
	mass2 = _AAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      } else {
	loc1 = atom->first;
	bin1 = _AAtomicdomain.getBin(loc1);
    pat1 = getCol('A',bin1);
	mass1 = _AAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom--;
	loc2 = atom->first;
	bin2 = _AAtomicdomain.getBin(loc2);
    pat2 = getCol('A',bin2);
	mass2 = _AAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      }  // end of if-block for comparing chmass1 and chmass2
      break;} // end of case 'A' block
  case 'P':
    {
      if (chmass1 > chmass2) {
	atom--;
	loc1 = atom->first;
	bin1 = _PAtomicdomain.getBin(loc1);
	pat1 = getRow('P',bin1);
	mass1 = _PAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom++;
	loc2 = atom->first;
	bin2 = _PAtomicdomain.getBin(loc2);
	pat2 = getRow('P',bin2);
	mass2 = _PAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      } else {
	loc1 = atom->first;
	bin1 = _PAtomicdomain.getBin(loc1);
	pat1 = getRow('P',bin1);
	mass1 = _PAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom--;
	loc2 = atom->first;
	bin2 = _PAtomicdomain.getBin(loc2);
	pat2 = getRow('P',bin2);
	mass2 = _PAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      }  // end of if-block for comparing chmass1 and chmass2
      break;}  // end of case 'P' block
  } // end of switch-block for extracting the atomic proposal info

  // stop if pattern1 == pattern2
  if (pat1 == pat2){
    return;
  }

     double priorLL = 0.;

  // Metropolis-Hasting move action

  double lambda;
  switch(the_matrix_label){
  case 'A':
    {lambda = _lambdaA;
      break;}
  case 'P':
    {lambda = _lambdaP;
      break;}
  }

  vector <double> newPat1;
  vector <double> newPat2;
  double delLLnew = 0.0; // EJF - MFO check if this changes calc from del LL
  bool makeChange;

   // Calculate the change in the patterns and evaluate corresponding changes
   // in likelihood. 
   _new_atomicProposal.insert(pair<unsigned long long, double>(loc1,newMass1-mass1));  	
   makeChange = calc_new_matrix_Pattern(the_matrix_label, newPat1, loc1, newMass1-mass1);

   _new_atomicProposal.insert(pair<unsigned long long, double>(loc2,newMass2-mass2));
   makeChange = calc_new_matrix_Pattern(the_matrix_label, newPat2, loc2, newMass2-mass2);
   
   delLLnew+= computeDeltaLLMEMap(the_matrix_label, D,S,AOrig,POrig,newPat1,pat1,newPat2,pat2);
  

  // EJF: MFO check double totalLL = priorLL + delLLnew * _annealingTemperature;

  double tmp;
  
    tmp = priorLL + delLLnew*_annealingTemperature;
 

  double rng = log(randgen('U',0,0));
 
    if (tmp < rng) {
	// Proposal rejected by MH. Don't update matrix
      switch(the_matrix_label){
      case 'A':
	  { 
	  _AAtomicdomain.rejectProposal(false);
	  break;
	  }
      case 'P':
	  {
	  _PAtomicdomain.rejectProposal(false);
	  break;
	  }
     } // end of switch-block
    } else if (makeChange){ 
	  // Proposal accepted and change can be made. Need to update the fixed 
	  // pattern matrix with the correct locations that were fixed. 
      switch(the_matrix_label){
      case 'A':
	  { 	   
	   _AAtomicdomain.setProposedAtomMass(_new_atomicProposal,false);
	   _AAtomicdomain.acceptProposal(false); 
	   update_sysChi2(delLLnew);  // update system Chi2
       //update the matrix for loc1
	   update_fixed_pattern('A',newPat1,pat1);
       //update the matrix for loc2
	   update_fixed_pattern('A',newPat2,pat2);
	   break;
	  }
      case 'P':
	  {   
       _PAtomicdomain.setProposedAtomMass(_new_atomicProposal,false);
	   _PAtomicdomain.acceptProposal(false); 
	   update_sysChi2(delLLnew);  // update system Chi2
	   //update the matrix for loc1
	   update_fixed_pattern('P',newPat1,pat1);
	   //update the matrix for loc2
	   update_fixed_pattern('P',newPat2,pat2);
	   break;
	  }
     } // end of switch-block       
    } // end proposal accepted block  
 

  // end of M-H sampling
 
} // end of method move

// Mapped exchange when both atoms are in a fixed pattern
// Determines whether or not proposal can be improved and/or 
// accepted and changes the atomic space and matrix is this is
// true. 
void GibbsSamplerMap::mappedExchange(char the_matrix_label,
						    double const * const * D,
						    double const * const * S,
						    double ** AOrig,
						    double ** POrig)
{
  map<unsigned long long, double>::const_iterator atom;
  double chmass1,chmass2;       
  unsigned long long loc1, loc2;
  unsigned int bin1, bin2;
  unsigned int pat1, pat2;
  double mass1, mass2;
  double newMass1, newMass2;
  atom = _atomicProposal.begin();
  chmass1 = atom->second;
  atom++;
  chmass2 = atom->second;

  // extract location, bin #, mass and changed mass corresponding to the 
  // atomic proposal such that "1" refers to a positive mass change and 
  // "2" a negative one.  
 switch(the_matrix_label){
  case 'A':
    {
      if (chmass1 > chmass2) {
	atom--;
	loc1 = atom->first;
	bin1 = _AAtomicdomain.getBin(loc1);
	pat1 = getCol('A',bin1);
	mass1 = _AAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom++;
	loc2 = atom->first;
	bin2 = _AAtomicdomain.getBin(loc2);
	pat2 = getCol('A',bin2);
	mass2 = _AAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      } else {
	loc1 = atom->first;
	bin1 = _AAtomicdomain.getBin(loc1);
    pat1 = getCol('A',bin1);
	mass1 = _AAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom--;
	loc2 = atom->first;
	bin2 = _AAtomicdomain.getBin(loc2);
    pat2 = getCol('A',bin2);
	mass2 = _AAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      }  // end of if-block for comparing chmass1 and chmass2
      break;} // end of case 'A' block
  case 'P':
    {
      if (chmass1 > chmass2) {
	atom--;
	loc1 = atom->first;
	bin1 = _PAtomicdomain.getBin(loc1);
	pat1 = getRow('P',bin1);
	mass1 = _PAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom++;
	loc2 = atom->first;
	bin2 = _PAtomicdomain.getBin(loc2);
	pat2 = getRow('P',bin2);
	mass2 = _PAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      } else {
	loc1 = atom->first;
	bin1 = _PAtomicdomain.getBin(loc1);
	pat1 = getRow('P',bin1);
	mass1 = _PAtomicdomain.getMass(loc1);
	newMass1 = atom->second + mass1;
	atom--;
	loc2 = atom->first;
	bin2 = _PAtomicdomain.getBin(loc2);
	pat2 = getRow('P',bin2);
	mass2 = _PAtomicdomain.getMass(loc2);
	newMass2 = atom->second + mass2;
      }  // end of if-block for comparing chmass1 and chmass2
      break;}  // end of case 'P' block
  } // end of switch-block for extracting the atomic proposal info

  // Stop if exchanges are in same pattern
  if (pat1 == pat2){
    return;
  }
  

  vector <double> newPat1;
  vector <double> newPat2;

  // preparing quantities for possible Gibbs computation later.
  // EJF  bool exchange = false;
  double priorLL = 0.;

  unsigned int jGene, jSample;
  bool anyNonzero = false;
  bool useGibbs = true;

  unsigned int iGene1, iPattern1, iGene2, iPattern2, iSample1, iSample2;
  switch(the_matrix_label){
  case 'A':
    {
      iGene1 = getRow('A',bin1);
      iPattern1 = getCol('A',bin1);
      iGene2 = getRow('A',bin2);
      iPattern2 = getCol('A',bin2);
      break;}
  case 'P':
    {
      iPattern1 = getRow('P',bin1);
      iSample1 = getCol('P',bin1);
      iPattern2 = getRow('P',bin2);
      iSample2 = getCol('P',bin2);
      break;}
  }

  // ---------------------------------------------
  switch(the_matrix_label){
  case 'A':
    {
      for (jSample = 0; jSample < _nCol; jSample++) {
	if (POrig[iPattern1][jSample] > epsilon) {
	  anyNonzero = true;
	  break;
	}
	if (POrig[iPattern2][jSample] > epsilon) {
	  anyNonzero = true;
	  break;
	}
      }  // end of for-block to determine the existence of corresponding 
      // non-zero elements in P
      if (!anyNonzero)  {  // cannot update in Gibbs way
	useGibbs = false;
      }
      break;}
  case 'P':
    {
      for (jGene = 0; jGene < _nRow; jGene++) {
	if (AOrig[jGene][iPattern1] > epsilon) {
	  anyNonzero = true;
	  break;
	}
	if (AOrig[jGene][iPattern2] > epsilon) {
	  anyNonzero = true;
	  break;
	}
      }  // end of for-block to determine the existence of corresponding 
      // non-zero elements in P
      if (!anyNonzero)  {  // cannot update in Gibbs way
	useGibbs = false;
      }
      break;}
  }  // end of switch-block

  // -------------------------------------------------------------------------
  // EXCHANGE ACTION when initial useGibbs = true and check if Gibbs is usable.

  double s=0., su, mean, sd; // EJF- MFO check
  pair <double, double> alphaparam;

  if (useGibbs == true){

    switch(the_matrix_label){
	
	// ---------- EXCHANGE ACTION WITH A ----------------------------
    case 'A':
	{
	 if (_AAtomicdomain.inDomain(loc1) && _AAtomicdomain.inDomain(loc2))
	 {
      alphaparam = GAPSNorm::calcAlphaParameters('A', _nRow, _nCol, _nFactor, D, S, AOrig, POrig, iGene1, 
	                                              iPattern1, iGene2, iPattern2, iSample1, iSample2);
	   s = alphaparam.first;
	   su = alphaparam.second;

	  s = s * _annealingTemperature;
	  su = su * _annealingTemperature; 
	  mean = su / s;
	  sd = 1./sqrt(s);
	   // end of compute distribution parameters for A	       

	  } // end of if-block for checking whether the changes are in domain (the exchange block)
	
      break;
	  } // end of switch block for EXCHANGE ACTION with A
      
	  // ---------- EXCHANGE ACTION WITH P ----------------------------
    case 'P': 
	{
      if (_PAtomicdomain.inDomain(loc1) && _PAtomicdomain.inDomain(loc2)) 
      {
      alphaparam = GAPSNorm::calcAlphaParameters('P', _nRow, _nCol, _nFactor, D, S, AOrig, POrig, iGene1, 
	                                              iPattern1, iGene2, iPattern2, iSample1, iSample2);
	   s = alphaparam.first;
	   su = alphaparam.second;

	  s = s * _annealingTemperature;
	  su = su * _annealingTemperature; 
	  mean = su / s;
	  sd = 1./sqrt(s);
	  
	   // end of compute distribution parameters for P      

	  } // end of if-block for checking whether the changes are in domain (the exchange block)
	
      break;
	 } // end of switch block for EXCHANGE ACTION with P
 
    }// end of switch block for EXCHANGE ACTION

  } // end of if-block for operations with possibly Gibbs sampling

  if (s == 0. && su == 0.){
    useGibbs = false;
    //cout << "Parameters aren't updated -> useGibbs = false, do M-H" << endl;
  }

  // -------------------------------------------------------------------------
  if(useGibbs == true){

    // set newMass1	 
    // need to retain exponential prior
    double plower = sub_func::pnorm(-mass1, mean, sd, DOUBLE_NEGINF, 0);
    double pupper = sub_func::pnorm(mass2, mean, sd, DOUBLE_NEGINF, 0);
    double u = plower + randgen('U',0,0)*(pupper - plower);
 
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

    if (useGibbs == true){

      gibbsMass1 = sub_func::qnorm(u, mean, sd, DOUBLE_NEGINF, 0);
      if (gibbsMass1 < -mass1) gibbsMass1 = -mass1;
      if (gibbsMass1 > mass2) gibbsMass1 = mass2;
      gibbsMass2 = - gibbsMass1;


      double delLLnew=0.0; // EJF -- MFO check if fixes delLL calc
      bool makeChange1, makeChange2;
	  
      _new_atomicProposal.insert(pair<unsigned long long, double>(loc1,gibbsMass1));
      _new_atomicProposal.insert(pair<unsigned long long, double>(loc2,gibbsMass2));
	  
      switch(the_matrix_label){
      case 'A':
	 {
	  extract_new_atomicProposal('A');
	  
      makeChange1 = calc_new_matrix_Pattern('A',newPat1,loc1,gibbsMass1);
	  makeChange2 = calc_new_matrix_Pattern('A',newPat2,loc2,gibbsMass2);

	  if (makeChange1 && makeChange2){
	  
	  delLLnew+= computeDeltaLLMEMap('A',D,S,AOrig,POrig,newPat1,pat1,newPat2,pat2);
	  
	  update_fixed_pattern('A',newPat1,pat1);
	  update_fixed_pattern('A',newPat2,pat2);
	  
	  _AAtomicdomain.setProposedAtomMass(_new_atomicProposal, false);
	  _AAtomicdomain.acceptProposal(false); //proposal accepted

	  update_sysChi2(delLLnew);  // update system Chi2
	  }
	  break;
	  }
      case 'P':
	 {
	  extract_new_atomicProposal('P');
	  
      makeChange1 = calc_new_matrix_Pattern('P',newPat1,loc1,gibbsMass1);
	  makeChange2 = calc_new_matrix_Pattern('P',newPat2,loc2,gibbsMass2);

	  if (makeChange1 && makeChange2){
	  
	  delLLnew+= computeDeltaLLMEMap('P',D,S,AOrig,POrig,newPat1,pat1,newPat2,pat2);
	  
	  update_fixed_pattern('P',newPat1,pat1);
	  update_fixed_pattern('P',newPat2,pat2);
	  
	  _PAtomicdomain.setProposedAtomMass(_new_atomicProposal, false);
	  _PAtomicdomain.acceptProposal(false); //proposal accepted

	  update_sysChi2(delLLnew);  // update system Chi2
	  }
	  break;
	  }
     }  // end of switch-block       
      
	  return;
	 
    } // end of inner if-block for final updating with Gibbs

  } // end of outer if-block for useGibbs == true && s > epsilon
   

  // ----------------------------
  // We can't use Gibbs, need
  // Metropolis-Hasting exchange action

  double pold = 0.;
  double pnew = 0.;
  double lambda;
  switch(the_matrix_label){
  case 'A':
    {lambda = _lambdaA;
      break;}
  case 'P':
    {lambda = _lambdaP;
      break;}
  }
  
  // Formerly the if(_oper_type == 'E') block in move_exchange
    if (mass1 > mass2) {
      pnew = sub_func::dgamma(newMass1, 2., 1./lambda, false);
      if (newMass1 > newMass2) {
	pold = sub_func::dgamma(mass1, 2., 1./lambda, false);
      } else {
	pold = sub_func::dgamma(mass2, 2., 1./lambda, false);
      }
    } else {
      pnew = sub_func::dgamma(newMass2, 2., 1./lambda, false);
      if (newMass1 > newMass2) {
	pold = sub_func::dgamma(mass1, 2., 1./lambda,  false);
      } else {
	pold = sub_func::dgamma(mass2, 2., 1./lambda, false);
      }
    }  


  
  if (pnew == 0. && pold == 0.) {
    priorLL = 0.0;
  } else if(pnew != 0. && pold == 0.) {
    priorLL = DOUBLE_POSINF;
  } else {
    priorLL = log(pnew / pold);
  }
	 
  double delLLnew=0.0; // EJF - MFO check if fixes delLL
  bool makeNewCh1, makeNewCh2;

   makeNewCh1 = calc_new_matrix_Pattern(the_matrix_label, newPat1, loc1, newMass1-mass1);
   makeNewCh2 = calc_new_matrix_Pattern(the_matrix_label, newPat2, loc2, newMass2-mass2);
   
   delLLnew+= computeDeltaLLMEMap(the_matrix_label, D,S,AOrig,POrig,newPat1,pat1,newPat2,pat2);
  

  // EJF- MFO check double totalLL = priorLL + delLLnew * _annealingTemperature;
 

  double tmp;

    tmp = priorLL + delLLnew*_annealingTemperature;
 

  double rng = log(randgen('U',0,0));

    if (tmp  < rng && priorLL != DOUBLE_POSINF) {
      switch(the_matrix_label){
      case 'A':
	{ _AAtomicdomain.rejectProposal(false);
	  break;}
      case 'P':
	{ _PAtomicdomain.rejectProposal(false);
	  break;}
      } // end of switch-block
    } else if (makeNewCh1 && makeNewCh2){ 
      // Proposal accepted and can be done. Update atomic domain, matrix, and Chi2. 
      switch(the_matrix_label){
      case 'A':
	  { 	   
	   _AAtomicdomain.setProposedAtomMass(_new_atomicProposal,false);
	   _AAtomicdomain.acceptProposal(false); 
	   update_sysChi2(delLLnew);  // update system Chi2
	   // Update the matrix for locs 1 and 2
	   update_fixed_pattern('A',newPat1,pat1);
	   update_fixed_pattern('A',newPat2,pat2);

	   break;
	  }
      case 'P':
	  {   
       _PAtomicdomain.setProposedAtomMass(_new_atomicProposal,false);
	   _PAtomicdomain.acceptProposal(false); 
	   update_sysChi2(delLLnew);  // update system Chi2
	   // Update the matrix for locs 1 and 2
	   update_fixed_pattern('P',newPat1,pat1);
	   update_fixed_pattern('P',newPat2,pat2);

	   break;
	  }
     } // end of switch-block       
    }  


   // end of M-H sampling
  
} // end of method exchange

 // *************** METHODS FOR UPDATING THE MATRIX ***********************************	

  // Calculates the new row in A or P based on the proposed location and mass.
  // If any member of the matrix pattern drops below 0, return false and stop. 
  bool GibbsSamplerMap::calc_new_matrix_Pattern(char the_matrix_label, 
                                     vector <double> &PatternUpdate,
                                     unsigned long long location, double mass){
  unsigned int theBin, thePat;
  // EJF int flag = 0;
  switch (the_matrix_label){
   case 'A':
   {
     // Get the location and corresponding fixed pattern
     theBin = _AAtomicdomain.getBin(location);
     thePat = getCol('A', theBin);
     _AMatrix.get_Col(thePat,PatternUpdate);
	 // Refer to the original vector (sum of pattern[i]=1) to obtain
	 // the distribution of masses
	 for (int iRow=0; iRow < _nRow; iRow++){
	  (PatternUpdate.at(iRow)) += ((_MapValues.at(thePat)).at(iRow))*mass;
	  // Temporary to address negative masses. 
	  // If the change would cause negative masses, return false so
	  // no changes are made. 
	  // THIS IS A TEMPORARY SOLUTION - NEED TO FIND OUT WHY THIS HAPPENS
	  if (PatternUpdate.at(iRow) < 0){
	  	//cout << "Negative mass in A bin # " << theBin << endl;
		return false;
	   }
	  }
	 break;
	 } //end case 'A'
   case 'P':
   {
    // Get the location and corresponding fixed pattern
     theBin = _PAtomicdomain.getBin(location);
     thePat = getRow('P', theBin);
     _PMatrix.get_Row(thePat,PatternUpdate);
	 // Refer to the original vector (sum of pattern[i]=1) to obtain
	 // the distribution of masses
	 int fixedPatNum = (thePat - (_nFactor - _nFixedMaps));
	 for (int iCol=0; iCol < _nCol; iCol++){
	  (PatternUpdate.at(iCol)) += ((_MapValues.at(fixedPatNum)).at(iCol))*mass;
	  // Temporary to address negative masses. 
	  // If the change would cause negative masses, return false so
	  // no changes are made. 
	  // THIS IS A TEMPORARY SOLUTION - NEED TO FIND OUT WHY THIS HAPPENS
	  if (PatternUpdate.at(iCol) < 0){
	  	//cout << "Negative mass in P bin # " << theBin << endl;
		return false;
	   }
	  } 
	 break;
	 } //end case 'P'
  } //end switch block
   return true;
 } //end calc_new_matrix_pattern method
  
  // Updates the matrix A or P with the proposed new pattern. 
  void GibbsSamplerMap::update_fixed_pattern(char the_matrix_label, vector <double> &newPat,
												unsigned int thePat){
  switch (the_matrix_label){
   case 'A':
   {
	 // Update the matrix with the new column:
	 _AMatrix.setCol(newPat, thePat);
	 break;
	 } //end case 'A'
   case 'P':
   {
	 // Update the matrix with the new row:
	 _PMatrix.setRow(newPat, thePat);
	 break;
	 } //end case 'P'
  } //end switch block
 } //end method

