#ifndef __FITHANDLER__H
#define __FITHANDLER__H

#include <iostream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <math.h>
#include "LinAlg.h"
#include "RefHandler.h"
#include <omp.h>
#include "Power.h"
//#include <Eigen/Dense>
#include "OutputHandler.h"
#include <limits>
#include <TMatrixDSym.h>
#include <TDecompBK.h>


class OutputHandler;

/**
* This class handles the fitting procedure.
*/
class FitHandler
{

public:
	//Constructor
	FitHandler();

	//calculates the next iterationstep
	void nextstep(RefHandler& rh, double threshold, double kappa, Power& pow);

	//getter for @_bfperr
	void setfiterrors(RefHandler& rh, double threshold, double kappa);
	void setfiterrors2(ConfigHandler& ch, OutputHandler& oh, size_t num_ipol);

	//getter for the mean of the components of @_ei
	double getDdist();

	//getter for the smoothness
	const double getDsmooth(RefHandler& rh);

	//getter for Chi2
	const double getchi2();

	//getter for @_iterationcounter
	size_t getiterationcounter();

	//setup for a new bin
	void nextbin(RefHandler& rh, size_t num_ipol, double threshold, double kappa, Power& pow);

	//getter for the number of fitparameters
	size_t getnumfitparams();

	//getter fir the fitparameters
	vector<double> getfitparams();

	//getter for @_a
	const vector<double> geta();

	//setter for a specific bin in a specific iterationstep
	void setiteration(RefHandler& rh, size_t num_ipol, double threshold, double kappa, size_t bestiteration, Power& pow);

	//getter for the mean relative error
	const double getmeanerror();

	//calculates the dot product between every gradient vector
	vector<double> getallgraddotproducts(RefHandler& rh, size_t num_bins);

	//calculates the error at a given point
	const double geterrorat(vector<double> point);

	//calculates the reduced Chi^2
	const double getchi2red();

	//adds 0's as fit parameters in order to become usable in Professor 2
	void sortfitparams();

	//getter for @_max
	int getmax();

	//getter for the fit errors
	vector<double> getfiterrors();

	//this function loads an interpolationfile
	void load_fit(ConfigHandler& ch, size_t i, size_t j, size_t num_ipol, Power& pow, RefHandler& rh, OutputHandler& oh);

private:
	/**
	* @_m: storage of the matrix M
	* @_q, @r: storage of the matrices Q & R of the QR-decomposition
	* @_d: identity matrix
	* @_mprime: coloumnwise normalized M
	* @_b: reference values
	* @_a: vector for the fitparameters; is used as indicator of RR constrained parameters
	* @_bprime: normalized @_b
	* @_bfp, @_bfperr: best fit parameters and corresponding errors
	* @_ei: distance measurements between fitfunction and reference values
	* @_sigma: error of the reference values
	* @_max: maximum of powers
	* @_iterationcounter: number of iterations in the fitting process
	* @_power: list that contains the powers of the variables at the terms of the fitting function
	* @_structure: static container for the normalvectors; it stores the derivatives of the anchors points
	*/
	vector<vector<double>> _m, _r, _q, _d, _mprime;
	vector<double> _b, _a, _bprime, _bfp, _bfperr, _ei, _sigma;
	int _max; 
	size_t _iterationcounter;
	vector<vector<int>> _power;
	double _dsmooth_current = 2;
	static vector<vector<vector<double>>> _structure;

	//adding components to @_q & @_r
	void expandqr();

	//check for RR constraints
	void collinearity(double threshold, size_t i, double kappa);

	//calculate @_mprime from @_m
	void makemprime();

	//adding elements to @_m
	void increasem(RefHandler& rh);

	//adding new power terms to @_power
	void windofchange(Power& pow);

	//rescaling @_bfp because of the normalization
	void rescalebestfitparameters();

	//getter for the @_ei
	double getei(RefHandler& rh, size_t i);
	
	//calculates the gradients of the function
	vector<double> normalvecfunction(RefHandler& rh, size_t i);
	
	//this function walks through an iteration without any checks
	void nextstep_walkthrough(RefHandler& rh, Power& pow);	

	double getinvcovmatelement(size_t i, size_t j);
	
};

#endif

