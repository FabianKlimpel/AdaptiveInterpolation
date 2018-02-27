#ifndef __REFHANDLER__H
#define __REFHANDLER__H

#include <iostream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <math.h>
#include <limits>
#include "LinAlg.h"
#include "ConfigHandler.h"
#include <omp.h>

/**
* This class handles the simulated reference data
*/
class RefHandler
{
friend class FitHandler;

public:

	//Constructor
	RefHandler(ConfigHandler& ch);
	
	//calculates the gradient vectors for the smoothness convergence check
	void calculatenormalvectors(size_t num_ipol, double threshold, int power, double kappa);

	//getter for the number of bins
	size_t getnum_analysis();

	//getter for the content of the bins in every run
	vector<vector<double>> getbin_values();

	//getter for the errors of the bins in every run
	vector<vector<double>> getbin_values_err();

	//getter for the sampled parameters
	const vector<vector<double>>& getparameter_values();

	//getter for the number of bins in every observable
	vector<size_t> getnum_bins();

	//deletes the normalvectors
	void clearnormalvectors();

	//calculates the dot product between every gradient vector
	vector<double> getallgraddotproducts();

	//getter for the start of the intervalls of the bins
	vector<vector<double>> getintervall_start();

	//getter for the end of the intervalls of the bins
	vector<vector<double>> getintervall_end();

	//rescaling the parameter sets onto a [0,1] hypercube
	void rescale(vector<double>& minparval, vector<double>& maxparval);
	
	//scales parameters back from a [0,1] hypercube to normal values
	void rerescale(vector<double>& minparval, vector<double>& maxparval);
	
private:
	
	/**
	* @_parameter_values: contains the sampled parameter values
	* @_bin_values: contains the simulated reference data
	* @_bin_values_err: contains the simulated reference data error
	* @_normalvectors: contains the normalvectors of the anchor points
	* @_intervall_start, @_intervall_end: list of the start and the end of the bins respectively
	* @_num_bins: contains the number of bins in each observable
	* @_hypercubes: static variable that stores the points for every anchor point that build their surrounding hypercube
	*/
	vector<vector<double>> _parameter_values, _bin_values, _bin_values_err, _normalvectors, _intervall_start, _intervall_end;
	vector<size_t> _num_bins;
	static vector<vector<size_t>> _hypercubes;

	//reader for the parameter values and the simulated reference data (error)
	void readparametervalues(ConfigHandler& ch);
	void readsimulationdata(ConfigHandler& ch);
	
	//getter for the anchor points that form a hypercube around a certain point
	vector<size_t> gethypercube(size_t i);

	//calculates the normalvector for an anchor point
	vector<double> getnormalvec(vector<size_t>& hypercube, vector<double>& b, size_t i, double threshold, int power, double kappa);

	//calculates the best fit parameters
	vector<double> getfitparams(vector<vector<double> >& points, vector<double>& bpoints, double threshold, int power, double kappa);

	//creates a QR-decomposition
	void expandqrlocal(vector<vector<double> >& mlocal, vector<vector<double> >& qlocal, vector<vector<double> >& rlocal, size_t i);
};

#endif

