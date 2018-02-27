#ifndef __CONFIGHANDLER__H
#define __CONFIGHANDLER__H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <math.h>
#include <dirent.h>
#include <time.h>

using namespace std;

/**
* This class handles the config file provided at the start of the program.
* It handles the reading of the file and stores all necessary information as the source paths of the data.
*/
class ConfigHandler
{
public:
	ConfigHandler(char* configfile);

	/**
	* @_var_names: storage for the names of the variables
	* @_analysis: storage for the names of the anlyses
	* @_observable: storage for the names of the observables
	* @_run: storage for the numbers of the MC runs
	* @_dimension: number of variables
	* @_num_runs: total amount of MC runs
	* @_rcleave_out: # of runs of the total amount of runs that will be left out for a runcombination
	* @_rngseed: seed for the rng
	* @_path: path for the MC data
	* @_paramsfile: name of the file that contains the used parameters
	* @_filename: name of the interpolationfile that will be written to or read from
	* @_summaryflag, @_outdotflag: flags for writing additional summaries
	* @_rescaleflag: flag for rescaling the sampled parameters onto a [0,1] hypercube
	* @_load: flag if an interpolationfile should be read
	*/
	vector<string> _var_names, _analysis, _observable, _run; 
	size_t _dimension = 0, _num_runs = 0, _rcleave_out = 0, _rngseed = 0;
	string _path = "", _paramsfile, _filename = "interpolationresult", _outpath = "";
	int _summaryflag = 0, _outdotflag = 0, _rescaleflag = 0, _runcombs = 0, _load = 0, _covmat = 0;
	/**
	* @_thresholdfit: threshold for the RR constrain in the fitting
	* @_thresholddata: threshold for the RR constrain in the hypercube-fitting of the data
	* @_thresholderr: threshold for the RR constrain in getter of the fitting errors
	* @_ddist: convergence criteria for the distance between fitfunction and data
	* @_chi2mean: number of chi2's to store in order to state a convergence regarding chi2
	* @_dsmoothmean: same as @_chi2mean but for the smoothness dsmooth
	* @_kappa: shifting parameter if the RR constrain is applied
	* @_exponent: exponent for the distance weighting of the hypercube
	* @_smoothness: threshold for a smoothness convergence
	* @_chi2limit: upper limit for the convergence of the chi2 value
	*/
	double _thresholdfit = 0, _thresholddata = 0, _thresholderr = 0, _ddist = 1, _chi2mean = 1, _dsmoothmean = 1, _kappa = 0, _exponent = 1, _smoothness = 0, _chi2limit = 10e10;
	
private:

	//helper for information gathering from strings
	size_t _pos_begin, _pos_end, _length;

	//setter for the member variables
	void read_path(string line);
	void read_varsfile(string line);
	void read_num_runs(string line);
	void set_run_names();
	void read_var_names(string line);
	void set_analyses();

	void read_thresholdfit(string line);
	void read_thresholddata(string line);
	void read_thresholderr(string line);
	void read_ddist(string line);
	void read_chi2mean(string line);
	void read_dsmoothmean(string line);
	void read_kappa(string line);
	void read_exponent(string line);
	void read_smoothness(string line);
	void read_chi2limit(string line);

	void read_summaryflag(string line);
	void read_outdotflag(string line);
	void read_rescaleflag(string line);
	void read_runcombs(string line);
	void read_leave_out(string line);
	void read_rngseed(string line);
	
	void read_load(string line);
	void read_filename(string line);
	void read_covmat(string line);
	void read_outpath(string line);

};

#endif
