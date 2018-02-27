#include "FitHandler.h"

using namespace std;
//~ using namespace Eigen;

/**
 * Constructor
 */
FitHandler::FitHandler(){
}

/**
 * This function builds the QR-decomposition or increase the matrices, if they already exist
 * @i: represents the iterationstep
 * @_q, @_r: the respective matrices of the QR-decomposition
 * @_mprime: coloumnwise rescaled @_m, c.f. @FitHandler::makemprime() 
 * @tmp: temporary storage for the increase of @_r
 * @sum: temporary storage for the increase of @_q
 */
void FitHandler::expandqr(){
	//in the first iteration, @q & @r need to be resized
	if(_q.empty() || _r.empty())
	{
		
		_q.resize(_mprime.size());
		_r.resize(_mprime[0].size());
			
		for(size_t k = 0; k < _q.size(); k++)
			_q[k].resize(_mprime[0].size());
		
		
		for(size_t k = 0; k < _r.size(); k++)
			_r[k].resize(_mprime[0].size());
		
		//calculating the first elements by using http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5342430&tag=1 on top of eq. 5
		for(size_t k = 0; k < _q.size(); k++)
			_q[k][0] = LinAlg::getcol(_mprime,0)[k] / LinAlg::getabs(LinAlg::getcol(_mprime,0));

		_r[0][0] = LinAlg::getabs(LinAlg::getcol(_mprime,0));
		return;
	}

	//in a later iteration, the size should be adapted according to the current size of @_mprime 
	_q.resize(_mprime.size());
	_r.resize(_mprime[0].size());


	for(size_t k = 0; k < _q.size(); k++)
		_q[k].resize(_mprime[0].size());
	
	
	for(size_t k = 0; k < _r.size(); k++)
		_r[k].resize(_mprime[0].size());
	
	
	//the following calculation add the new elements to @_q & @_r
	double tmp;
	for(size_t l = 0; l <= _iterationcounter - 1; l++)
	{
		tmp = 0;
		for(size_t k = 0 ; k < _q.size(); k++)
			tmp += _q[k][l] * _mprime[k][_iterationcounter];
			
		_r[l][_iterationcounter] = tmp;
	}


	vector<double> sum;
	sum.assign(_q.size(), 0);
	for(size_t l = 0; l < sum.size(); l++)
		for(size_t k = 0; k <= _iterationcounter - 1; k++)
			sum[l] += _r[k][_iterationcounter] * _q[l][k];
			

	for(size_t k = 0; k < _q.size(); k++)
		_q[k][_iterationcounter] = _mprime[k][_iterationcounter] - sum[k];


	_r[_iterationcounter][_iterationcounter] = LinAlg::getabs(LinAlg::getcol(_q, _iterationcounter));

	tmp = LinAlg::getabs(LinAlg::getcol(_q, _iterationcounter));
	for(size_t k = 0; k < _q.size(); k++)
	{
		_q[k][_iterationcounter] /= tmp;
		
		//If the absolut value of @q[k][i + 1] is 0, the result would become nan. Setting it instead to 0 stabilizes the calculations.
		if(std::isnan(_q[k][_iterationcounter])) 
			_q[k][_iterationcounter] = 0;
	}
}

/**
 * This function checks for the collinearity of a coloumn of the matrix @_m in comparison with the previously coloumns.
 * This leads to a small values at the corresponding diagonalterm in @_r. The "smallness" is checked and, if it's true, the according fitparameter will be artificial set in order to
 * keep the solution from diverging.
 * @threshold: this value represents the threshold, below which a diagonalterm of @r is considered as small
 * @i: this values represents the bin that is concerned
 * @riihat, @bihat: helpervariables for setting of the fitparameter in @_a
 */
void FitHandler::collinearity(double threshold, size_t i, double kappa){

	//check if the value in @r is small
	if(fabs(_r[i][i]) < threshold)
	{
		//setting the helpervariables for the creation of a smaller @a than one would get if the fit would run normally.
		//the shift has it's source in @kappa and the according term of @_d.
		double riihat = sqrt(_r[i][i] * _r[i][i] + kappa * _d[i][i]);
		double bihat = (_r[i][i] / riihat) * _b[i];
		_a[i] = bihat / riihat;
	}
}

/**
 * This function normalizes every coloumn of @_m by using the absolut value of the respective coloumn.
 * @coloumn: represents of the coloumn of interest
 * @abs: absolut value of the coloumn
 */
void FitHandler::makemprime(){

	//setting @_mprime to the same size as @_m
	_mprime.resize(_m.size());
	for(size_t i = 0; i < _mprime.size(); i++)
		_mprime[i].resize(_m[i].size());

	//setting the absolut value of the coloumn, because while assigning it row wise, it changes after every assigned value
	double abs = LinAlg::getabs(LinAlg::getcol(_m, _iterationcounter));
	
	//setting the normalization
	for(size_t i = 0; i < _m.size(); i++)
		_mprime[i][_iterationcounter] = _m[i][_iterationcounter] / abs;

}

/**
 * This function adds a new coloumn to @_m
 * @rh: container for the anchor points
 * @_power: matrix, that contains all the potencies of the values
 * @size: old number of coloumns in @_m
 * @tmp: temporary storage of the product of the different powers of the values
 */
void FitHandler::increasem(RefHandler& rh){

	//set the size of @_m, if not done yet
	if(_m.empty())
		_m.resize(rh._parameter_values.size());

	size_t size = _m[0].size();
	for(size_t i = 0; i < _m.size(); i++)
		_m[i].resize(size + 1);

	double tmp = 1;
	//walking over every row of @_m
	for(size_t j = 0; j < _m.size(); j++)
	{	
		//in every row the product of the values and the respective powers is calculated
		for(size_t i = 0; i < rh._parameter_values[0].size(); i++)
			tmp *= pow(rh._parameter_values[j][i], _power[_m[0].size() - 1][i]);
			
		//the new terms will be added to the last coloumn in every row and @tmp is resetted
		_m[j][_m[0].size() - 1] = tmp;
		tmp = 1;
	}

}

/**
 * This function sets the list of powers needed wrt @_max
 * @pow: This object calculates the list of powers of a certain order
 */
void FitHandler::windofchange(Power& pow){
	
		#pragma omp critical (windofchange)
		_power = pow.getpoweroforder(_max);
	
}

/**
 * This function rescales the best fit parameters due to normalization of @_m and @_b
 * @vec: vector, containing the best fit parameters
 */
void FitHandler::rescalebestfitparameters(){

	//backwards calculating of the norms that influenzed @_m and @_b in order to get the regular best fit parameters
	for(size_t i = 0; i < _bfp.size(); i++)
		_bfp[i] *= LinAlg::getabs(_b) / LinAlg::getabs(LinAlg::getcol(_m, i));
		
	//~ for(size_t i = 0; i < _b.size(); i++)
		//~ if(_b[i] > 10)
			//~ cout << i << "\t" << _b[i] << endl;
}

/**
 * This function calculates the normalized distancevector between the fitfunction and the datapoints at a certain point.
 * @rh: container for all anchor points
 * @i: # of datapoint at which its distance will be calculated
 * @functionvalue: functionvalue at the given point, represented by @i
 * @functiongradient: gradient of the fitfunction at the given point; used for normalization
 * @tmp: temporary storage while calculating the parts of a monomial
 */
double FitHandler::getei(RefHandler& rh, size_t i){

	double functionvalue = 0;
	
	//taking the monomials, saved in @m and multiply it by the respective fitparameter
	for(size_t j = 0; j < _bfp.size(); j++)
		functionvalue += _m[i][j] * _bfp[j];
	
	//use the absolute value of the difference of the functionvalue and the datapoint
	functionvalue = fabs(functionvalue - _b[i]);
	
	//setting up the gradient of the function
	vector<double> functiongradient;
	functiongradient.assign(rh._parameter_values[0].size(), 0);
	double tmp = 1;
	
	//walking over every component of the gradient
	for(size_t j = 0; j < functiongradient.size(); j++)
	{
		//walk over every fitparameter
		for(size_t k = 1; k < _bfp.size(); k++) // leave out constant offset
		{
			//walk over every value
			for(size_t l = 0; l < rh._parameter_values[0].size(); l++)
				//if the component of @functiongradient and @l match, the respective term will be derived
				if(j == l)
				{
					//setting a zero term to zero and therefore the whole product too leads to stability against the function @pow() itself
					if(rh._parameter_values[i][l] == 0)
						tmp = 0;
					else
						//derive the factor x^n in the monomial by simply calculate n*x^(n - 1)
						tmp *= pow(rh._parameter_values[i][l], _power[k][l] - 1) * _power[k][l];
				}
				else
					//every other factor will be left as it is
					tmp *= pow(rh._parameter_values[i][l], _power[k][l]);
		
			//finally multiply the prefactor to @tmp and add it as one monomial to the respective component in the gradient
			tmp *= _bfp[k];
			functiongradient[j] += tmp;
			tmp = 1;
		}
	}

	//if the length of the gradient is zero, @functionvalue itself will be returned, otherwise it will be normalized
	if(LinAlg::getabs(functiongradient) == 0)
		return functionvalue;
	else
		if(functionvalue > 0)
			return functionvalue / LinAlg::getabs(functiongradient);
		else
			return -functionvalue / LinAlg::getabs(functiongradient);
}

/**
 * This function calculates the mean of the overall distance between the data and the fit function as one breakoff criteria
 * @sum: sum of all differences
 */
double FitHandler::getDdist(){

	double sum = 0;
	
	//sum up all differences
	for(size_t i = 0; i < _ei.size(); i++)
		sum += _ei[i];

	//return the mean
	return sum / _ei.size();

}

/**
 * This function calculates the Chi^2 value of the fit.
 * @result: resulting Chi^2
 * @functionvalue: functionvalue of the fit at a certain point
 */
const double FitHandler::getchi2(){

	double result = 0, functionvalue = 0;
	//~ ofstream off, off2, off3, off4;
	//~ off.open("read");
	//~ off2.open("readm");
	//~ off3.open("readbfp");
	//~ off4.open("readmono");
	//~ int doneflag = 1;
	
	//walk over every polynomial
	for(size_t i = 0; i < _m.size(); i++)
	{
		//if the bin entry is 0, it can be a simple misscalculation in Rivet etc.
		if(_b[i] == 0)
			continue;
			
		//calculate the monomials and add them to the result
		//~ for(size_t j = 0; j < _m[i].size(); j++){off4 << _m[i][j] * _bfp[j] << " ";
			//~ functionvalue += _m[i][j] * _bfp[j];}
		//~ off4 << endl;
		
		for(size_t j = 0; j < _m[i].size(); j++)
			functionvalue += _m[i][j] * _bfp[j];

		//get the difference between the functionvalue and the datapoint		
		functionvalue -= _b[i];
		
		//if the uncertainty is 0, the calculation would rise a nan
		//this is prevented by setting it to the arbitrary value of 1e-10
		if(_sigma[i] == 0)
			_sigma[i] = 1e-10;
		
		//~ if(_iterationcounter == 334)	
		//~ {
			//~ off << i << " " << functionvalue << " " << _b[i] << " " << _sigma[i] << endl;
		//~ 
			//~ for(size_t j = 0; j < _m[i].size(); j++)
				//~ off2 << _m[i][j] << " ";
			//~ off2 << endl;
			//~ 
			//~ if(doneflag)
			//~ {
				//~ for(size_t j = 0; j < _bfp.size(); j++)
					//~ off3 << _bfp[j] << endl;
				//~ doneflag = 0;
			//~ }
			//~ 
		//~ }
		//calculating the squares of @functionvalue and @sigma, divide them and add this new summand to new overall Chi^2 value @result
		result += (functionvalue * functionvalue) / (_sigma[i] * _sigma[i]);
		//~ if ((functionvalue * functionvalue) / (_sigma[i] * _sigma[i]) > 10e3)
			//~ cout << i << "\t" << (functionvalue * functionvalue) / (_sigma[i] * _sigma[i]) << endl;
		functionvalue = 0;

	}
	//~ off.close();
	//~ off2.close();
	//~ off3.close();
	//~ off4.close();
		
	return result;
}

/**
 * This function calculates the Chi^2 value of the fit.
 * @result: resulting Chi^2
 * @functionvalue: functionvalue of the fit at a certain point
 * @skips: number of bins that were skipped due to 0 as entry, this needs to be tracked in order to get the ndof
 */
const double FitHandler::getchi2red(){

	double result = 0, functionvalue = 0;
	size_t skips = 0;
	
	//walk over every polynomial
	for(size_t i = 0; i < _m.size(); i++)
	{
		//if the bin entry is 0, it can be a simple misscalculation in Rivet etc.
		if(_b[i] == 0)
		{
			skips++;
			continue;
		}
		
		//calculate the monomials and add them to the result
		for(size_t j = 0; j < _m[0].size(); j++)
			functionvalue += _m[i][j] * _bfp[j];
			
		
		//get the difference between the functionvalue and the datapoint		
		functionvalue -= _b[i];
		
		//if the uncertainty is 0, the calculation would rise a nan
		//this is prevented by setting it to the arbitrary value of 1e-10
		if(_sigma[i] == 0)
			_sigma[i] = 1e-10;
			
		
		//calculating the squares of @functionvalue and @sigma, divide them and add this new summand to new overall Chi^2 value @result
		result += (functionvalue * functionvalue) / (_sigma[i] * _sigma[i]);
		
		functionvalue = 0;
		
	}

	//return Chi^2 / ndof
	return result / (_m.size() - skips - _m[0].size());
}

/**
 * This function calculates the normalvectors
 * @rh: handles the reference data
 * @i: # of the anchor point
 * @functiongradient: noramlized gradient vector
 * @tmp: temporary storage for a monomial the gradient structure
 * @tmpstrc: storage for a new @_structure
 */
vector<double> FitHandler::normalvecfunction(RefHandler& rh, size_t i){
	
	//if the static variable @_structure is not initialized yet, its size will be set
	if(_structure.size() == 0)
	{
		//the critical pragma prevents multiple accesses to the size at the same time
		#pragma omp critical (structurebuild)
		{
			//set the size to the # of anchor points
			_structure.resize(rh._parameter_values.size());
			
			//set for every anchor point the size to the dimension, so there is one for every derivative
			for(size_t j = 0; j < _structure.size(); j++)
				_structure[j].resize(rh._parameter_values[0].size());
		}
	}
	
	//set up the gradient and resize to the number of dimensions, initialized as 0
	vector<double> functiongradient;
	functiongradient.assign(rh._parameter_values[0].size(), 0);
	double tmp = 1;
	
	//calculating the components of the gradient
	for(size_t j = 0; j < functiongradient.size(); j++)
	{
		//ifenough monomials were already calculated, they can be taken directly
		if(_structure[i][j].size() > _bfp.size()) 
			//walking over every monomial and calculate the contribution to the gradient vector
			for(size_t k = 1; k < _bfp.size(); k++)
				functiongradient[j] += _bfp[k] * _structure[i][j][k];
		else
		{
			//if not enough monomials were already calculated in structure, they need to be calculated
			//a critical pragma prevents multiple accesses at the same time and therefore potentially problems
			#pragma omp critical (structureadd)
			{
			//walking over every monomial, calculate what is already available
			for(size_t k = 1; k < _structure[i][j].size(); k++)
				functiongradient[j] += _bfp[k] * _structure[i][j][k];
		
			//use a new variable and calculate the new elements for @_structure there, replace them later
			vector<double> tmpstruc = _structure[i][j];
			
			//set the size to the number of fit parameters as the new maximum needed
			tmpstruc.resize(_bfp.size());
			
			//start the calculation after the last calculated component and walk up to the new maximum needed
			for(size_t k = _structure[i][j].size(); k < tmpstruc.size(); k++)
			{
				//extracting the respective parameter values
				for(size_t l = 0; l < rh._parameter_values[0].size(); l++)
					//if the parameter is the one that is derived in the component of the gradient, then its derivative is used
					if(j == l)
						//if the value of the parameter is 0 & the power of the parameter to derive is != 1, the whole monomial will be 0 but the special case in the derivative of 0^0 = 1
						//if the power is 0, the whole monomial will be 0 after derived 
						if((rh._parameter_values[i][l] == 0 && _power[k][l] != 1) || _power[k][l] == 0)
						{	
							tmp = 0;
							continue;
						}
														
						else			
							//multiply the derived factor
							tmp *= pow(rh._parameter_values[i][l], _power[k][l] - 1) * _power[k][l];	

					//multiply the rest to @tmp
					else
						//if at least one of the not derived parameters is 0 and its power is !=0, the whole monomial become 0
						if(rh._parameter_values[i][l] == 0 && _power[k][l] != 0)
						{
							tmp = 0;
							continue;
						}
						//else the rest will be multiplied
						else
							tmp *= pow(rh._parameter_values[i][l], _power[k][l]);
					
				//put the new monomial part to the @tmpstruc at the specific point in the list
				tmpstruc[k] = tmp;

				//adding the monomial to the gradient component
				functiongradient[j] += tmp * _bfp[k];
				tmp = 1;
			}
			
				//replace the new structure list as the static list for all threads
				//another if-condition concerning the length of the list is meant to prevent possible waiting threads at the beginning of this critical pragma to enter afterwards and replace
				//a longer list by a shorter one
				//therefore this calculation will appear very often in the beginning of the program but fewer and fewer afterwards if more was already calculated. Therefore it is kept static.
				if(tmpstruc.size() > _structure[i][j].size())
					_structure[i][j] = tmpstruc;
			
			}
		}
	}
	
	//normalizing the gradient
	tmp = LinAlg::getabs(functiongradient);
	
	//in order to avoid NaN's, the gradient is only normalized if it isn't a zerovector, else it's normalized vector is the vector itself
	if(tmp != 0)
		for(size_t i = 0; i < functiongradient.size(); i++)
			functiongradient[i] /= tmp;
	
	return functiongradient;
}



/**
 * This function returns Dsmooth
 * @rh: container for reference gradient vectors
 * @result: result of the calculation that will be returned
 */
const double FitHandler::getDsmooth(RefHandler& rh){
	
	if(_dsmooth_current == 2)
	{
		double result = 0;
			
		//walking over all anchor points and adding the dotproduct of the normalvectors
		for(size_t i = 0; i < rh._parameter_values.size(); i++)
			result += LinAlg::dotproduct(rh._normalvectors[i], normalvecfunction(rh, i));
			
		
		_dsmooth_current = result / (double) rh._parameter_values.size();
		
	}
	
	//returning the mean
	return _dsmooth_current;

}

/**
 * This function calculates the errors of the fit
 * @rh: contains the QR-decomposition
 * @threshold: threshold for the RR constraint
 * @kappa: shift in case of applied RR constraint
 * @result: vector containing the fitparameters and errors
 * @err: vector containing the errors of the reference data
 * @merr, @rerr, @qerr: matrices M, Q & R for the errors
 */
void FitHandler::setfiterrors(RefHandler& rh, double threshold, double kappa){

	
	vector<double> result, err;
	result.resize(_m[0].size());
	
	//setting the vector to NaN in order to identify RR constraints later on
	result.assign(result.size(), nan("1"));
	err.resize(_m[0].size());

	vector<vector<double> > merr, rerr, qerr;
	merr.resize(_m.size());

	//constructing @merr, @qerr, @rerr as in the fit itself by using @rh's method @rh::expandqrlocal()
	for(size_t i = 0; i < _m.size(); i++)
	{	
		//~ merr[i].push_back(_m[i][0] * _m[i][0]);	
		merr[i].push_back(_m[i][0]);					
	}
	rh.expandqrlocal(merr, qerr, rerr, 0);

	for(size_t i = 1; i < _m[0].size(); i++)
	{
		for(size_t j = 0; j < _m.size(); j++)
			//~ merr[j].push_back(_m[j][i] * _m[j][i]);
			merr[j].push_back(_m[j][i]);
			
		rh.expandqrlocal(merr, qerr, rerr, i - 1);
		
	}
	
	//calculating the variances
	for(size_t i = 0; i < err.size(); i++)
		//~ err[i] = _sigma[i] * _sigma[i];
		err[i] = _sigma[i];
	
	//calculating the RR constraint
	double riihat, bihat;
	for(size_t i = 0; i < rerr.size(); i++)
		if(rerr[i][i] < threshold)
		{
			riihat = sqrt(rerr[i][i] * rerr[i][i] + kappa * _d[i][i]);
			bihat = (rerr[i][i] / riihat) * err[i];
			result[i] = bihat / riihat;
		}
		
	//getting the fit parameters (= variances of the fitparameters)
	result = LinAlg::getbestfitparameters(rerr, result, LinAlg::multmatvec(LinAlg::transpose(qerr), err));
	
	//extracting the errors
	//~ for(size_t i = 0; i < result.size(); i++)
		//~ if(result[i] >= 0)
			//~ result[i] = sqrt(result[i]);
		//~ else
			//~ result[i] = sqrt(-result[i]);
	
	//create the absolute value of the fit errors as the real errors
	for(size_t i = 0; i < result.size(); i++)
		result[i] = fabs(result[i]);
		
	_bfperr = result;
}


void FitHandler::setfiterrors2(ConfigHandler& ch, OutputHandler& oh, size_t num_ipol){

	TMatrixDSym tmds(_bfp.size());

	double tmp;
	
	for(size_t row = 0; row < _bfp.size(); row++)
		for(size_t col = row; col < _bfp.size(); col++)
		{
			
			tmp = getinvcovmatelement(row, col);
		
			if(row == col)
				tmds[row][col] = tmp;

			else
			{
				tmds[row][col] = tmp;
				tmds[col][row] = tmp;
			}
		}
	
	//~ for(size_t i = 0; i < 10; i++)
	//~ {
		//~ for(size_t j = 0; j < 10; j++)
			//~ cout << tmds[i][j] << "\t";
		//~ cout << endl;
	//~ }
	//~ cout << endl;
	
	TDecompBK tdbk(tmds);
	tdbk.Invert(tmds);

	//~ for(size_t i = 0; i < 10; i++)
	//~ {
		//~ for(size_t j = 0; j < 10; j++)
			//~ cout << tmds[i][j] << "\t";
		//~ cout << endl;
	//~ }
	
	for(size_t i = 0; i < _bfp.size(); i++)
		if(tmds[i][i] == std::numeric_limits<double>::infinity())
		{
			_bfperr.push_back(std::numeric_limits<double>::max());
			cout << _bfperr.size() << endl;
			cout << _bfperr.back() << endl;
		}
		else
			if(-tmds(i, i) == std::numeric_limits<double>::infinity())
				_bfperr.push_back(-std::numeric_limits<double>::max());
			else
				_bfperr.push_back(sqrt(tmds[i][i]));

	if(ch._covmat)
		oh.write_covmat(tmds, num_ipol, ch);
}

double FitHandler::getinvcovmatelement(size_t i, size_t j)
{
	double result = 0;
	
	
	for(size_t k = 0; k < _m.size(); k++)
	{ 
		if(_sigma[k] == 0)
			continue;
			//~ return std::numeric_limits<double>::infinity();
		
		result += _m[k][i] * _m[k][j] / (_sigma[k] * _sigma[k]);
		
	}

	return result;
}


/**
 * This function calculates the next iterationstep in a bin
 * @rh: container of the reference data
 * @threshold: threshold for the RR constraint of the fit
 * @kappa: shift in case of applied RR constraint
 * @pow: getter for the list of powers of a certain order
 */
void FitHandler::nextstep(RefHandler& rh, double threshold, double kappa, Power& pow){
	
	_dsmooth_current = 2;
	
	//if every power mentioned in @_power is already in use in @_m, new components of a higher power needs to be calculated
	if(_m[0].size() == _power.size())
	{
		//@_max is the order of the polynom. It will be increased and the new powers will be calculated
		_max++;
		windofchange(pow);
	}

	//setting up the new member variables
	increasem(rh);
	_iterationcounter++;
	makemprime();

	_a.resize(_mprime[0].size());

	expandqr();


	
	_bprime = LinAlg::normalizevec(_b);
	
	_bprime = LinAlg::multmatvec(LinAlg::transpose(_q), _bprime);
		
	
	_d.resize(_m[0].size());
	for(size_t i = 0; i < _d.size(); i++)
		_d[i].assign(_d.size(), 0);
	for(size_t i = 0; i < _d.size(); i++)
		_d[i][i] = 1;
	for(size_t i = 0; i < _a.size(); i++)
	{
		_a[i] = nan("1");
		collinearity(threshold, i, kappa);
	}
	
	//getting the fit and rescale it
	_bfp = LinAlg::getbestfitparameters(_r, _a, _bprime);
	rescalebestfitparameters();
	
	//clearing the old @_ei and calculating the new ones
	//~ _ei.clear();
	//~ for(size_t i = 0; i < rh._parameter_values.size(); i++)
		//~ _ei.push_back(getei(rh, i));
		

	
}

/**
 * This function walks through an iteration without any checks, so that the matrices are created only
 * @rh: handles the reference data
 * @pow: getter for a list of powers
 */
void FitHandler::nextstep_walkthrough(RefHandler& rh, Power& pow){

	
	//if every power mentioned in @_power is already in use in @_m, new components of a higher power needs to be calculated
	if(_m[0].size() == _power.size())
	{
		//@_max is the order of the polynom. It will be increased and the new powers will be calculated
		_max++;
		windofchange(pow);
	}

	//setting up the new member variables
	increasem(rh);
	_iterationcounter++;
	makemprime();
	expandqr();
	
}


/**
 * This function is a getter for @_iterationcounter
 */
size_t FitHandler::getiterationcounter(){

	return _iterationcounter;
	
}

/**
 * This function sets up the object for a new bin
 * @rh: container of the reference data
 * @num_ipol: number of the bin
 * @threshold: threshold for the RR constraint of the fit
 * @kappa: shift in case of applied RR constraint
 * @pow: getter for the list of powers of a certain order
 */
void FitHandler::nextbin(RefHandler& rh, size_t num_ipol, double threshold, double kappa, Power& pow){

	//deleting old information
	_iterationcounter = 0;
	_max = 0;
	
	_power.clear();
	_bfp.clear();
	_ei.clear();

	_bfperr.clear();
	_m.clear();
	_r.clear();
	_q.clear();
	_d.clear();
	_mprime.clear();
	_a.clear();


	//setting up the new reference vector and the corresponding error
	_b = LinAlg::transpose(rh._bin_values)[num_ipol];
	_sigma = LinAlg::transpose(rh._bin_values_err)[num_ipol];
		
	//calculating the first some powers, and setting up the matrices and @_a

	windofchange(pow);
	increasem(rh);
	makemprime();
	_a.resize(_mprime[0].size());
	expandqr();
	
	//normalize @_b and preparing it as new righthand side of the LinAlg problem
	_bprime = LinAlg::normalizevec(_b);
	_bprime = LinAlg::multmatvec(LinAlg::transpose(_q), _bprime);

	//setting up the identity matrix
	_d.resize(_m.size());
	for(size_t i = 0; i < _d.size(); i++)
		_d[i].assign(_d.size(), 0);
	for(size_t i = 0; i < _d.size(); i++)
		_d[i][i] = 1;
		
	//checking for RR constraints
	for(size_t i = 0; i < _a.size(); i++)
	{
		_a[i] = nan("1");
		collinearity(threshold, i, kappa);
	}

	//getting the fit and rescaling it
	_bfp = LinAlg::getbestfitparameters(_r, _a, _bprime);
	
	rescalebestfitparameters();

	//getting @_ei
	//~ for(size_t i = 0; i < rh._parameter_values.size(); i++)
		//~ _ei.push_back(getei(rh, i));
}

/**
 * This function is a getter for the size of @_bfp
 */
size_t FitHandler::getnumfitparams(){

	return _bfp.size();
	
}

/**
 * This function is getter for @_bfp
 */
vector<double> FitHandler::getfitparams(){

	return _bfp;
	
}

/**
 * This function is a getter for @_a
 */
const vector<double> FitHandler::geta(){
	
	return _a;
	
}

/**
 * This function sets the object up for a certain bin and calculates a given number of iterations
 * @rh: container of the reference data
 * @num_ipol: number of the bin
 * @threshold: threshold for the RR constraint
 * @kappa: shift in case of applied RR constraint
 * @bestiteration: number of iterations that will be performed
 */
void FitHandler::setiteration(RefHandler& rh, size_t num_ipol, double threshold, double kappa, size_t bestiteration, Power& pow){

	//setting up the bin
	nextbin(rh, num_ipol, threshold, kappa, pow);
	
	//Walking the number of iterations along. The final state can be addressed through the interfaces of the object itself.
	if(bestiteration > 0)
	{
		for(size_t i = 0; i < bestiteration - 1; i++)
			nextstep_walkthrough(rh, pow);
		
		nextstep(rh, threshold, kappa, pow);
	}
}

/**
 * This function calculates the mean of the relative errors of the fit function
 * @relerr: storage for the relative error
 */
const double FitHandler::getmeanerror(){
	
	//check, if the fit was done and the errors were calculated
	if(_bfp.size() != 0 && _bfperr.size() != 0)
	{
		
		double relerr = 0;
		int skips = 0;
		
		//walk over every fit parameter
		for(size_t i = 0; i < _bfp.size(); i++)
		{
			
			//if a fit parameter is 0, a relative error would be inf
			if(_bfp[i] == 0)
			{
				skips++;
				//~ cout << "error: fitparameter is zero" << endl;
				continue;
			}
			
			//calculating the relative error of a fit parameter and adding it up
			relerr += _bfperr[i] / _bfp[i];
		}
			
		//return the mean
		return relerr / (_bfp.size() - skips);
	}
	
	//if there was no fit or fit error, this function returns 0
	return 0;
	
}

/**
 * This function calculates the dot product of all normalized gradients at every anchor point with every other
 * @rh: container for the anchor points
 * @num_runs: number of MC runs
 * @result: vector containing all dot products
 */
vector<double> FitHandler::getallgraddotproducts(RefHandler& rh, size_t num_runs){
	
	vector<double> result;
	vector<vector<double>> tmp;
	tmp.resize(num_runs);
	
	//walk over every possible combination of anchor points
	for(size_t i = 0; i < num_runs; i++)
	{
		if(tmp[i].empty())
			tmp[i] = normalvecfunction(rh, i);
		for(size_t j = 0; j < num_runs; j++)
		{
			if(tmp[j].empty())
				tmp[j] = normalvecfunction(rh, j);
			//only calculate the dot product, if the points aren't the same
			if(i != j)
				result.push_back(LinAlg::dotproduct(tmp[i], tmp[j]));
		}
	}
	
	return result;
	
	
}

/**
 * This function calculates the relative error of the fit function at a given point
 * @point: the point of interest
 * @result: stores the error of the fit function
 * @tmp: stores the function value of the fit function
 */
const double FitHandler::geterrorat(vector<double> point){

	double result = 0, tmp = 0;
	
	//calculate the error as Gaussian error propagation
	for(size_t i = 0; i < _bfperr.size(); i++)
	{
		for(size_t j = 0; j < _power[0].size(); j++)
			result += pow(point[j], _power[i][j]) * pow(point[j], _power[i][j]);
		result *= _bfperr[i] * _bfperr[i];
	}
	
	//calculate the function value
	for(size_t i = 0; i < _bfp.size(); i++)
	{
		for(size_t j = 0; j < _power[0].size(); j++)
			tmp += pow(point[j], _power[i][j]);
		tmp *= _bfp[i];
	}
	
	//return the relative error
	return sqrt(result) / tmp;
}

/**
 * This function adds 0's to the list of fit parameters and their errors in order to fill up a certain order.
 * That is needed in order to become usable for Professor 2. 
 */
void FitHandler::sortfitparams(){

    for(size_t i = _bfp.size(); i < _power.size(); i++)
    {
		_bfp.push_back(0);
		_bfperr.push_back(0);
	}
	
}

/**
 * This function is getter for @_max
 */
int FitHandler::getmax(){
	
	return _max;
	
}

/**
 * This function is getter for @_bfperr
 */
vector<double> FitHandler::getfiterrors(){
	
	return _bfperr;
	
}

/**
 * This function loads an interpolationfile
 * @ch: handles the names of analyses, observables, and the interpolationfilename
 * @i, @j, @num_ipol: indices for finding the data from the right bin
 * @pow: getter for list of powers
 * @rh: handles the reference data
 * @ifile: handles the data input
 * @pos_begin, @pos_end: helper for slicing the string
 * @line: keeps the line read by @ifile
 */
void FitHandler::load_fit(ConfigHandler& ch, size_t i, size_t j, size_t num_ipol, Power& pow, RefHandler& rh, OutputHandler& oh){

	//if the smoothness wasn't calculated yet, @_dsmooth_current = 2 is set due to _dsmooth_current is in [-1,1]
	_dsmooth_current = 2;
		
	//set up the file reading
	ifstream ifile;
	ifile.open(ch._filename);
	size_t pos_begin, pos_end;
		
	if(ifile.is_open())
	{
		string line;
		
		//The header information are not important for the loaded interpolationfile,
		//therefore this loop is a walkthrough up to the end of it. This is indicated by "---"
		while(getline(ifile, line))
			if(line == "---")
				break;
			
		//walk through the interpolation information
		while(getline(ifile, line))
			//if a line contains the name of the name of the analysis and the observable ...
			if(line.find(("/" + ch._analysis[i] + "/" + ch._observable[i].substr(0, ch._observable[i].find(".")).c_str())) != string::npos)
			{
				pos_begin = line.find("#");
				pos_end = line.find(" ", pos_begin);
				// ... check if it is also the right bin number aka @j
				//if so, then quit the loop
				if(line.substr(pos_begin + 1, pos_end - pos_begin - 1) == to_string(j))
					break;
			}
		
		//At that point, the right bin was found if it is in the file. The fit parameters and their errors are in the lines below.
			
		//The next line contains the fit parameters.
		getline(ifile, line);
		
		//The polynomial order will be extracted
		pos_begin = 9;
		pos_end = line.find(" ", pos_begin);
		_max = atoi(line.substr(pos_begin, pos_end - pos_begin).c_str());
		
		//The list of powers is created
		windofchange(pow);
			
		pos_begin = pos_end;
		
		//the information about the powers provides information about the number of best fit parameters => resize
		_bfp.resize(_power.size());
		
		//walk through the line and extract all parameter values and add them to the list of best fit parameters
		for(size_t k = 0; k < _power.size(); k++)
		{
			pos_end = line.find(" ", pos_begin + 1);
			_bfp[k] = atof(line.substr(pos_begin + 1, pos_end - pos_begin - 1).c_str());
			pos_begin = pos_end;
		}

		if(!ch._covmat)
		{
			//the next line contains the errors
			getline(ifile, line);

			//#errors = #fit parameters = _power.size()
			_bfperr.resize(_power.size());
			pos_begin = 11;

			//same walkthrough as for the best fit parameters
			for(size_t k = 0; k < _power.size(); k++)
			{
				pos_end = line.find(" ", pos_begin + 1);
				_bfperr[k] = atof(line.substr(pos_begin + 1, pos_end - pos_begin - 1).c_str());
					pos_begin = pos_end;
			}
		}

		
		//finish reading
		ifile.close();

		//For Professor, the orders need to be filled up. This would be waste of calculation time, therefore the useless added 0's will be deleted by pop_back()'s until something useful appears
		while(_bfp.back() == 0 && _bfp.size() > 1)
		{
			_bfp.pop_back();
			_bfperr.pop_back();
		}
			
		//@_iterationcounter serves as a reference point for the # of parameters for the ndof, so it will be adjusted
		//@_iterationcounter will increase after every iteration while fitting but the first one. So it will be shifted by -1 but for bins that were only once iterated.
		if(_bfp.size() == 0)
			_iterationcounter = 0;
		else
			_iterationcounter = _bfp.size() - 1;
	
		//set up the function values and their errors
		_b = LinAlg::transpose(rh._bin_values)[num_ipol];
		_sigma = LinAlg::transpose(rh._bin_values_err)[num_ipol];
		
		//build @_m
		for(size_t k = 0; k <= _iterationcounter; k++)
			increasem(rh);
		
		if(ch._covmat)
		{
			_bfperr.clear();
			setfiterrors2(ch, oh, num_ipol);
		}
		
	}
	else
		cout << "unable to read " << ch._filename << endl;
	
}

