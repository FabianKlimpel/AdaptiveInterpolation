#include "RefHandler.h"

using namespace std;




/**
 * Constructor that sets up the anchor points and simulated reference values
 * @ch: object that keeps information about the files that need to be read
 */
RefHandler::RefHandler(ConfigHandler& ch){

	//reading the information
	readparametervalues(ch);
	readsimulationdata(ch);

	
}

/**Reads the used parameters for the simulations
 * @ch: object that keeps information about the files that need to be read
 * @num_runs: this parameter stores the number of runs and therefore the number of different parametersets
 * @value, @exponent: these parameters are for the reconstruction of the actual numbers that are readed from the file
 * @ifile, @line: these are for walking through the file and reading its content linewise
 * @param_values: storage of one parameterset, that will added to the member variable @parameter_values
 **/
void RefHandler::readparametervalues(ConfigHandler& ch){

	int num_runs = ch._num_runs;
	_parameter_values.resize(num_runs);
	
	cout << "reading parameter values ..." << flush;
	
	//walk over every run that was performed
	for(int i = 0; i < num_runs; i++)
	{
		
		double value, exponent;	
		ifstream ifile;
		string line;
		vector<double> param_values;
	
		//opening the specific file
		ifile.open((ch._path + "/" + ch._run[i] + "/" + ch._paramsfile).c_str());
		
		if (ifile.is_open()) 
		{
			//linewise filereading
			while (getline(ifile, line)) 
			{
				//walk over all parameters that are used
				for(size_t j = 0; j < ch._dimension; j++)
				{
					//compare the readed parametername with the name saved in @dr.
					if(line.substr(0, ch._var_names[j].size()) == ch._var_names[j])
					{
						//extract the number via substr() and add to @param_values
						value = atof(line.substr(ch._var_names[j].size() + 1, line.find("e", ch._var_names[j].size() + 1) - ch._var_names[j].size() - 1).c_str());
						exponent = atof(line.substr(line.find("e", ch._var_names[j].size() + 1) + 1).c_str());
						param_values.push_back(value * pow(10, exponent));
							
					}
				}	
					
			}
			
			
			
			//adding the whole parameterset to the overall collection, clear tmp container and take next file
			
			_parameter_values[i] = param_values;
	
			param_values.clear();
			ifile.close();		
		}
		else
			cout << "unable to read from " << (ch._path + "/" + ch._run[i] + "/" + ch._paramsfile).c_str() << endl;
		
	}
	
	cout << "complete" << endl;
}


/**
 * This function reads the data (value and error from each bin) from the simulations, so that a combination of parameterset & binvalue exist
 * @ch: object that keeps information about the files that need to be read
 * @num_runs: this parameter stores the number of runs and therefore the number of different parametersets
 * @histo_flag: 0-1-parameter, that indicates, that the simulated data will be read next in the file (otherwise the referencedata would be read)
 * @start_flag: 0-1-parameter, that indicates, that (given @histo_flag is already set) the binvalues will start in the next line
 * @pos_begin, @pos_end: there parameters are for splitting the readed line into pieces, so that the information can be isolated and extracted
 * @value, @exponent: these parameters are for the reconstruction of the actual numbers that are readed from the file
 * @ifile, @line: these are for walking through the file and reading its content linewise
 * @bin_val, @bin_val_err: these parameters store the binvalues and corresponding binvalue errors of the current analysis and observable and will later added to the overall containers @bin_values & @bin_values_err
**/
void RefHandler::readsimulationdata(ConfigHandler& ch){

	int num_runs = ch._num_runs;
	
	//resizing the vectors for parallel access
	_bin_values.resize(num_runs);
	_bin_values_err.resize(num_runs);
	_intervall_start.resize(ch._analysis.size());
	_intervall_end.resize(ch._analysis.size());
	_num_bins.resize(ch._analysis.size());
	
	cout << "reading simulation data ..." << flush;
	
	//walk over every run that was performed
	#pragma omp parallel for
	for(int i = 0; i < num_runs; i++)
	{
		int histo_flag = 0, start_flag = 0, pos_begin, pos_end;
		double value, exponent;	
		ifstream ifile;
		string line;
		vector<double> bin_val, bin_val_err, start, end;
		size_t tmp = 0;
		
		for(size_t j = 0; j < ch._analysis.size(); j++)
		{
			//opening the specific file
			ifile.open((ch._path + "/" + ch._run[i] + "/plots/" + ch._analysis[j] + "/" + ch._observable[j]).c_str());
				
			if (ifile.is_open()) 
			{
				//linewise filereading
				while (getline(ifile, line)) 
				{		
					//"LineColor" is a parameter in the files, that is only used for the configuration of the plot of the simulated data. Thus, its appearance indicates,
					//that the simulated data can be read in some of the following line in the file. This "soon" is stored in @histo_flag by changing its value
					if(line.substr(0, 9) == "LineColor")
							histo_flag = 1;
					
					
					//@start_flag = 1 means, that the data can be read
					if(start_flag == 1)
					{
						//If the lines of data in the file are over, it is always followed by a line starting with "#".
						//Therefore this is an indicator for the end of the reading and will be stored by resetting the parameters @start_flag and @histo_flag.
						if(line.substr(0, 1) == "#")
						{
							start_flag = 0;
							histo_flag = 0;
							
							continue;
						}
						
						//reading the start of the intervall
						value = atof(line.substr(0, line.find("e", 0)).c_str());	
						exponent = atof(line.substr(line.find("e", 0) + 1, line.find("	") - line.find("e", 0) - 1).c_str());
						start.push_back(value * pow(10, exponent));
						
						pos_begin = line.find("	");
						pos_end = line.find("	", pos_begin + 1);
						
						//reading the end of the intervall
						value = atof(line.substr(pos_begin + 1, line.find("e", pos_begin + 1) - pos_begin - 1).c_str());	
						exponent = atof(line.substr(line.find("e", pos_begin + 1) + 1, pos_end - line.find("e", pos_begin + 1) - 1).c_str());
						end.push_back(value * pow(10, exponent));
						
						pos_begin = pos_end;
						pos_end = line.find("	", pos_begin + 1);
						
						//Isolating the value, formatting it an storing it to @bin_val
						value = atof(line.substr(pos_begin + 1, line.find("e", pos_begin + 1) - pos_begin - 1).c_str());	
						exponent = atof(line.substr(line.find("e", pos_begin + 1) + 1, pos_end - line.find("e", pos_begin + 1) - 1).c_str());
						bin_val.push_back(value * pow(10, exponent));
						
						//same procedure again but for the error of the binentry
						pos_begin = pos_end;
						pos_end = line.find("	", pos_begin + 1);
						
						value = atof(line.substr(pos_begin + 1, line.find("e", pos_begin + 1) - pos_begin - 1).c_str());	
						exponent = atof(line.substr(line.find("e", pos_begin + 1) + 1, pos_end - line.find("e", pos_begin + 1) - 1).c_str());
						
						bin_val_err.push_back(value * pow(10, exponent));
						
						pos_begin = 0;
						pos_end = 0;
											
					}	
					
					//There is always a line starting with "#" that indicates, that in the following lines either the data is written or end.
					//A checkup of @histo_flag gives the safety, that the data is from the simulation. The beginning is stored via setting @start_flag
					if(histo_flag == 1)
						if(line.substr(0, 1) == "#")
							start_flag = 1;	
							

									
				}
				ifile.close();
				
				if(i == 0)
				{
					_num_bins[j] = bin_val.size() - tmp;
					tmp += _num_bins[j];
				}
			}
			else
			{
				cout << "unable to read from " << (ch._path + "/" + ch._run[i] + "/plots/" + ch._analysis[j] + "/" + ch._observable[j]).c_str() << endl;
				//~ break;	
			}
			
			value = 1;
			exponent = 0;
		
			//the intervalls remain the same for every run (assumption!) and therefore they need to be setted only once
			if(_intervall_start[j].empty() || _intervall_end[j].empty())
			{
				_intervall_start[j] = start;
				_intervall_end[j] = end;		
			}
			
			start.clear();
			end.clear();
		}
		
		//adding all values to the containers, clearing the tmp-containers and continue with the next file

		_bin_values[i] = bin_val;
		bin_val.clear();
		_bin_values_err[i] = bin_val_err;
		bin_val_err.clear();
		
		
	}
	
	cout << "complete" << endl;
}

/**
 * This function calculates the gradient vectors of every reference data point
 * @num_ipol: specifies the bin for which the gradient vectors need to determined
 * @threshold: RR constraint parameter for the fitting in the hypercubes
 * @power: specifies the power of the distance weighting
 * @kappa: shift in case of applied RR constraint
 * @b: vector that contains every simulated reference data point of the bin
 * @hypercube: list of indices of the points that contain to a hypercube around one point
 */
void RefHandler::calculatenormalvectors(size_t num_ipol, double threshold, int power, double kappa){

	
	vector<double> b = LinAlg::transpose(_bin_values)[num_ipol];
	vector<size_t> hypercube;
	
	
	_normalvectors.resize(_parameter_values.size());
	hypercube.resize(_parameter_values[0].size() * 2);
		
	//walk over every point
	for(size_t i = 0; i < _parameter_values.size(); i++){
		
		//get the points for the respective hypercube
		hypercube = gethypercube(i);

		//calculate its gradient vector
		_normalvectors[i] = getnormalvec(hypercube, b, i, threshold, power, kappa);
		
	}
	
}

/**
 * This function constructs a hypercube around a certain point
 * @i: number of the point in @_parameter_values around which a hypercube should be constructed
 * @center: vector of the center of the hypercub around which the cube is constructed
 * @result: list of numbers that indicate the taken points for the hypercube
 * @distances: list that contains all the distances between @center and every other point
 * @bin: binary represantation of the location of the possible hypercube point
 */
vector<size_t> RefHandler::gethypercube(size_t i){
	
	if(_hypercubes.size() != _parameter_values.size())
		_hypercubes.resize(_parameter_values.size());
		
	if(!_hypercubes[i].empty())
		return _hypercubes[i];
	
	//setting up @center and the resultvector
	vector<double> center = _parameter_values[i];
	vector<size_t> result;
	
	//setting the resultsize as 2 * #dimensions, so that for every dimension at least 2 points cann be selected
	result.resize(pow(2., _parameter_values[0].size()));
	
	//if a point isn't set, it's value is #number of dimensions + 1, so that these points can be found and won't be used for further calculations
	//this is important for points at the border of the sampleregion
	result.assign(result.size(), _parameter_values.size());
	
	//The same procedure as before is done for the distances, so that a check is possible if the smallest possible hypercube cann be constructed. Therefore it's initial values are set to infinity.
	vector<double> distances;
	distances.resize(pow(2., _parameter_values[0].size()));
	distances.assign(distances.size(), numeric_limits<double>::infinity());	
	size_t bin = 0;
	
	for(size_t j = 0; j < _parameter_values.size(); j++)
		//if a datapoint matches the @center, it will be skipped
		if(i == j)
			continue;
		else
		{
			//component wise check, if the components are bigger than the center point
			for(size_t k = 0; k < _parameter_values[0].size(); k++)
				if(center[k] > _parameter_values[j][k])
					//All 'if' checks were binary questions, therefore a binary representation as summary can be found.
					//@bin is a value that can directly connect to the overall situation of the test point due to its value
					bin += pow(2., k);
		
			//check up, if the new trial point is closer than the stored one in the category specified by @bin
			if(LinAlg::getdistanceofvectors(center, _parameter_values[j]) < distances[bin])
			{
				//if closer than set it as new point
				distances[bin] = LinAlg::getdistanceofvectors(center, _parameter_values[j]);
				result[bin] = j;				
			}
			
			bin = 0;
		}
		
	#pragma omp critical (sethypercube)
	{
		_hypercubes[i] = result;
	}
	
	return _hypercubes[i];
	
}

/**
 * This function calculates the gradient vector of one point
 * @hypercube: list of indices of the points that form the hypercube around a point
 * @b: list of simulated reference data in a bin
 * @i: index of the point of interest
 * @threshold: RR constraint threshold for the fitting
 * @power: specifies the power of the distance weighting
 * @kappa: shift in case of applied RR constraint
 * @fitparams: vector that will store the fitting result and after that the gradient
 * @bpoints: list of simulated reference data that belong to the hypercube
 * @points: list of parameter values that belong to the hypercube
 * @length: storage for the length of a vector
 */
vector<double> RefHandler::getnormalvec(vector<size_t>& hypercube, vector<double>& b, size_t i, double threshold, int power, double kappa){

	vector<double> fitparams, bpoints;
	fitparams.resize(_parameter_values[0].size() + 1);

	vector<vector<double> > points;
	
	//adding the point of interest to the lists
	points.push_back(_parameter_values[i]);
	bpoints.push_back(b[i]);

	//adding the points of the hypercube to the lists
	for(size_t j = 0; j < hypercube.size(); j++)
		if(hypercube[j] < _parameter_values.size())
		{
			points.push_back(_parameter_values[hypercube[j]]);
			bpoints.push_back(b[hypercube[j]]);
		}
		
	//calculate the fit parameters
	fitparams = getfitparams(points, bpoints, threshold, power, kappa);

	//calculating the gradient = the fit parameters of the 1. order only
	fitparams.erase(fitparams.begin());
	double length = LinAlg::getabs(fitparams);
	
	//if the length of the gradient vector is != 0, it will be normalized
	//if it is = 0, then the normalized gradient vector is the vector itself
	if(length != 0)
		for(size_t j = 0; j < fitparams.size(); j++)
			fitparams[j] /= length;
		

	return fitparams;
	
}

/**
 * This function calculates the best fit parameters for the hypercubes
 * @points: matrix that contains all points
 * @bpoints: vector with the functionvalues
 * @threshold: threshold for artificial fitparameter setting
 * @power: list containing the powers for every value at every monomial
 * @result: container for the best fit parameters
 * @mlocal, @qlocal, @rlocal: local equivalent to @m, @q, @r respectively
 * @distances: vector for weighting the hypercubepoints with its distance to the center
 * @riihat, @bihat: regularization as in @collinearity()
 * @dlocal: local equivalent to @d
 */
vector<double> RefHandler::getfitparams(vector<vector<double> >& points, vector<double>& bpoints, double threshold, int power, double kappa){

	//setting up the parametervector, setting nans for finding paramaters to fit
	vector<double> result;
	result.resize(points[0].size() + 1);
	result.assign(result.size(), nan("1")); //muss gesetzt werden fuer fitfunktion, sonst annahme, dass wert bereits gesetzt ist
	vector<vector<double> > mlocal, qlocal, rlocal;
	mlocal.resize(points.size());	
	
	//calculating the distances of the points of the hypercube to the center
	vector<double> distances;
	
	distances.resize(points.size());
	
	//the zeroth component is set to 1 and therefore it won't be weighted otherwise
	distances[0] = 1.;
	double sum = 0;
	//calculate every distance and increase its "importance decrease by distance" by the power of @power
	for(size_t i = 1; i < distances.size(); i++)
	{
		distances[i] = pow(LinAlg::getdistanceofvectors(points[0], points[i]), power);
		sum += distances[i];
	}
		
	vector<double> weights;
	weights.push_back(1);
	for(size_t i = 1; i < distances.size(); i++)
		weights.push_back(1 - distances[i] / sum);
	
	
	//rescale the functionvalues by the new weights
	//~ bpoints[0] /= distances[0];
	//~ for(size_t i = 1; i < bpoints.size(); i++)
		//~ bpoints[i] /= distances[i];
	
	for(size_t i = 0; i < bpoints.size(); i++)
		bpoints[i] *= weights[i];
		
		
	//build an rescale @mlocal as in @increasem() etc.
	//~ mlocal[0].push_back(distances[0]);
//~ 
	//~ for(size_t i = 1; i < points.size(); i++)
	//~ {
		//~ mlocal[i].push_back(1. / distances[i]);		
	//~ }
	//~ expandqrlocal(mlocal, qlocal, rlocal, 0);
//~ 
	//~ for(size_t i = 0; i < points[0].size(); i++)
	//~ {
		//~ for(size_t j = 0; j < points.size(); j++)
			//~ mlocal[j].push_back(points[j][i] / distances[j]);
			//~ 
		//~ expandqrlocal(mlocal, qlocal, rlocal, i);
		//~ 
	//~ }
	
	for(size_t i = 0; i < points.size(); i++)
		mlocal[i].push_back(weights[i]);		
	
	expandqrlocal(mlocal, qlocal, rlocal, 0);

	for(size_t i = 0; i < points[0].size(); i++)
	{
		for(size_t j = 0; j < points.size(); j++)
			mlocal[j].push_back(points[j][i] * weights[j]);
			
		expandqrlocal(mlocal, qlocal, rlocal, i);
		
	}
	
	//regularize the QR-decomposition as in @collinearity()
	double riihat, bihat;
	vector<vector<double> > dlocal;
	dlocal.resize(rlocal.size());
	for(size_t i = 0; i < dlocal.size(); i++)
		dlocal[i].assign(rlocal.size(), 0.);
		
	for(size_t i = 0; i < dlocal.size(); i++)
		dlocal[i][i] = 1.;
		
		
	for(size_t i = 0; i < rlocal.size(); i++)
		if(rlocal[i][i] < threshold)
		{

			riihat = sqrt(rlocal[i][i] * rlocal[i][i] + kappa * 1.);
			bihat = (rlocal[i][i] / riihat) * bpoints[i];
			result[i] = bihat / riihat;
		}
	
	//return the fit
	return LinAlg::getbestfitparameters(rlocal, result, LinAlg::multmatvec(LinAlg::transpose(qlocal), bpoints));
	
	
}

/**
 * This function does the same as @FitHandler::expandqr() but takes the matrices as arguments and is therefore more functionally
 * @mlocal, @qlocal, @rlocal: analog to the FitHandler member variables @_m, @_q, @_r respectively
 * @i: counter of the respective iteration
 */
void RefHandler::expandqrlocal(vector<vector<double> >& mlocal, vector<vector<double> >& qlocal, vector<vector<double> >& rlocal, size_t i){
	if(qlocal.empty() || rlocal.empty())
	{
	
		qlocal.resize(mlocal.size());
		rlocal.resize(mlocal[0].size());
			
		for(size_t k = 0; k < qlocal.size(); k++)
		{
			qlocal[k].resize(mlocal[0].size());
		}
		
		for(size_t k = 0; k < rlocal.size(); k++)
		{
			rlocal[k].resize(mlocal[0].size());
		}
		
		for(size_t k = 0; k < qlocal.size(); k++)
			qlocal[k][0] = LinAlg::getcol(mlocal,0)[k] / LinAlg::getabs(LinAlg::getcol(mlocal,0));

		rlocal[0][0] = LinAlg::getabs(LinAlg::getcol(mlocal,0));
		return;
	}


	qlocal.resize(mlocal.size());
	rlocal.resize(mlocal[0].size());


	for(size_t k = 0; k < qlocal.size(); k++)
	{
	qlocal[k].resize(mlocal[0].size());
	}
	
	for(size_t k = 0; k < rlocal.size(); k++)
	{	
	rlocal[k].resize(mlocal[0].size());
	}
	

	double tmp;
	for(size_t l = 0; l <= i; l++)
	{
		tmp = 0;
		for(size_t k = 0 ; k < qlocal.size(); k++)
		{
			tmp += qlocal[k][l] * mlocal[k][i + 1];

		}
		rlocal[l][i + 1] = tmp;
	}


	vector<double> sum;
	sum.assign(qlocal.size(), 0);
	for(size_t l = 0; l < sum.size(); l++)
		for(size_t k = 0; k <= i; k++)
			sum[l] += rlocal[k][i + 1] * qlocal[l][k];


	for(size_t k = 0; k < qlocal.size(); k++)
		qlocal[k][i + 1] = mlocal[k][i + 1] - sum[k];
		
	

	rlocal[i + 1][i + 1] = LinAlg::getabs(LinAlg::getcol(qlocal,i + 1));
	
	tmp = LinAlg::getabs(LinAlg::getcol(qlocal, i + 1));
	if(tmp != 0)
		for(size_t k = 0; k < qlocal.size(); k++)
			qlocal[k][i + 1] /= tmp;
			
}

/**
 * This function is a getter for the total amount of bins
 */
size_t RefHandler::getnum_analysis(){

	return LinAlg::transpose(_bin_values).size();
	
}

/**
 * This function is a getter for the simulated reference data
 */
vector<vector<double>> RefHandler::getbin_values(){

	return _bin_values;
	
}

/**
 * This function is a getter for the simulated reference data errors
 */
vector<vector<double>> RefHandler::getbin_values_err(){

	return _bin_values_err;
	
}

/**
 * This function clears the normal vectors
 */
void RefHandler::clearnormalvectors(){

		_normalvectors.clear();
	
}

/**
 * This function is a getter for the parameter values
 */
const vector<vector<double>>& RefHandler::getparameter_values(){

		return _parameter_values;
	
}

/**
 * This function is a getter for the number of bins in each histogram
 */
vector<size_t> RefHandler::getnum_bins(){
	
	return _num_bins;

}

/**
 * This function calculates every dot product of every anchor point with every anchor point
 * @result: summary of all dot products
 */
vector<double> RefHandler::getallgraddotproducts(){
	
	vector<double> result;
	
	//calculate the dot product of every vector with every other
	for(size_t i = 0; i < _normalvectors.size(); i++)
		for(size_t j = 0; j < _normalvectors.size(); j++)
			//if they are the same, the calculation will be skipped
			if(i != j)
				result.push_back(LinAlg::dotproduct(_normalvectors[i], _normalvectors[j]));
	
	return result;
	
}

/**
 * This function returns the list of the start of all bins
 */
vector<vector<double>> RefHandler::getintervall_start(){
	
	return _intervall_start;
	
}

/**
 * This function returns the list of the end of all bins
 */
vector<vector<double>> RefHandler::getintervall_end(){
	
	return _intervall_end;
	
}

/**
 * This function scales the parameter sets onto a [0,1] hypercube
 * @minparval: the minimum of all parameter values
 * @maxparval: the maximum of all parameter values
 * @diff: the difference between the minimum and maximum value of each dimension
 */
void RefHandler::rescale(vector<double>& minparval, vector<double>& maxparval){

	//walk over every dimension and calculate the difference between the maximum and the minimum value
	vector<double> diff;
	for(size_t i = 0; i < minparval.size(); i++)
		diff.push_back(maxparval[i] - minparval[i]);
	
	//walk over every parameter set dimensionwise and rescale it onto a [0,1] hypercube
	for(size_t i = 0; i < _parameter_values.size(); i++)
		for(size_t j = 0; j < _parameter_values[i].size(); j++)
			_parameter_values[i][j] = (_parameter_values[i][j] - minparval[j]) / diff[j];
}

/**
 * This function scales parameter set from a [0,1] hypercube to regular values
 * @minparval: the minimum of all parameter values
 * @maxparval: the maximum of all parameter values
 * @diff: the difference between the minimum and maximum value of each dimension
 */
void RefHandler::rerescale(vector<double>& minparval, vector<double>& maxparval){
	
	//walk over every dimension and calculate the difference between the maximum and the minimum value
	vector<double> diff;
	for(size_t i = 0; i < minparval.size(); i++)
		diff.push_back(maxparval[i] - minparval[i]);
		
	//calculate the reverse mapping as in @RefHandler::rescale()
	for(size_t i = 0; i < _parameter_values.size(); i++)
		for(size_t j = 0; j < _parameter_values[i].size(); j++)
			_parameter_values[i][j] = _parameter_values[i][j] * diff[j] + minparval[j];	
		
}
