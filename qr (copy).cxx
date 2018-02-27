#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>
#include <limits>
#include "datareading.h"
#include "refreading.h"
#include <fstream>
#include "linalg.h"

using namespace std;

// example: http://www.mia.uni-saarland.de/Teaching/NAVC-SS11/sol_c8.pdf
// page 5
//http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5342430&tag=1

//~ vector<vector<double> > m, r, q, d, mprime;
//~ vector<double> b, a, bprime;
double kappa;
linalg la();

//~ /**
 //~ * This function returns a coloumn out of a matrix
 //~ * @mat: matrix, out of which the coloumn will be extracted
 //~ * @j: number of coloumn
 //~ * @tmp: coloumn that will be returned
 //~ */
//~ vector<double> getcol(vector<vector<double> > mat, size_t j){
//~ 
	//~ vector<double> tmp;
//~ 
	//~ for(size_t i = 0; i < mat.size(); i++)
		//~ tmp.push_back(mat[i][j]);
//~ 
	//~ return tmp;
//~ 
//~ }
//~ 
//~ /**
 //~ * This function delivers the absolut value of a vector
 //~ * @vec: vector of interest
 //~ * @abs: absolut value
 //~ */
//~ double getabs(vector<double> vec){
//~ 
	//~ double abs = 0.;
//~ 
	//~ //summing up the sqaure of all components of the vector
	//~ for(size_t i = 0; i < vec.size(); i++)
		//~ abs += vec[i] * vec[i];
//~ 
	//~ //return the sqrt of the sum
	//~ return sqrt(abs);
//~ 
//~ }

//~ /**
 //~ * This function builds the QR-decomposition or increase the matrices, if they already exist
 //~ * @i: represents the iterationstep
 //~ * @q, @r: the respective matrices of the QR-decomposition
 //~ * @mprime: coloumnwise rescaled @m, c.f. @makemprime() 
 //~ * @tmp: temporary storage for the increase of @r
 //~ * @sum: temporary storage for the increase of @q
 //~ */
//~ void expandqr(size_t i){
	//~ //in the first iteration, @q & @r need to be resized
	//~ if(q.empty() || r.empty())
	//~ {
		//~ 
		//~ q.resize(mprime.size());
		//~ r.resize(mprime[0].size());
			//~ 
		//~ for(size_t k = 0; k < q.size(); k++)
		//~ {
			//~ q[k].resize(mprime[0].size());
		//~ }
		//~ 
		//~ for(size_t k = 0; k < r.size(); k++)
		//~ {
			//~ r[k].resize(mprime[0].size());
		//~ }
		//~ //calculating the first elements by using http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5342430&tag=1 on top of eq. 5
		//~ for(size_t k = 0; k < q.size(); k++)
			//~ q[k][0] = la.getcol(mprime,0)[k] / la.getabs(la.getcol(mprime,0));
//~ 
		//~ r[0][0] = la.getabs(la.getcol(mprime,0));
		//~ return;
	//~ }
//~ 
	//~ //in a later iteration, the size should be adapted according to the current size of mprime 
	//~ q.resize(mprime.size());
	//~ r.resize(mprime[0].size());
//~ 
//~ 
	//~ for(size_t k = 0; k < q.size(); k++)
	//~ {
	//~ q[k].resize(mprime[0].size());
	//~ }
	//~ 
	//~ for(size_t k = 0; k < r.size(); k++)
	//~ {
	//~ r[k].resize(mprime[0].size());
	//~ }
	//~ 
	//~ //the following calculation add the new elements to @q & @r
	//~ double tmp;
	//~ for(size_t l = 0; l <= i; l++)
	//~ {
		//~ tmp = 0;
		//~ for(size_t k = 0 ; k < q.size(); k++)
		//~ {
			//~ tmp += q[k][l] * mprime[k][i + 1];
			//~ 
		//~ }
		//~ 
		//~ r[l][i + 1] = tmp;
	//~ }
//~ 
//~ 
	//~ vector<double> sum;
	//~ sum.assign(q.size(), 0);
	//~ for(size_t l = 0; l < sum.size(); l++)
		//~ for(size_t k = 0; k <= i; k++)
			//~ sum[l] += r[k][i + 1] * q[l][k];
			//~ 
//~ 
	//~ for(size_t k = 0; k < q.size(); k++)
	//~ {
		//~ q[k][i + 1] = mprime[k][i + 1] - sum[k];
		//~ 
	//~ }
//~ 
	//~ r[i + 1][i + 1] = la.getabs(la.getcol(q,i + 1));
//~ 
	//~ tmp = la.getabs(la.getcol(q, i + 1));
	//~ for(size_t k = 0; k < q.size(); k++){
		//~ q[k][i + 1] /= tmp;
		//~ 
		//~ //If the absolut value of @q[k][i + 1] is 0, the result would become nan. Setting it instead to 0 stabilizes the calculations.
		//~ if(std::isnan(q[k][i + 1])) 
			//~ q[k][i + 1] = 0;
	//~ }
//~ }

//~ /**
 //~ * This function transposes a given matrix
 //~ * @mat: matrix, that will be transposed
 //~ * @tmp: temporary storage of the transposed matrix
 //~ */
//~ vector<vector<double> > transpose(vector<vector<double> > mat){
//~ 
	//~ vector<vector<double> > tmp;
	//~ 
	//~ //setting @tmp's size to the size of @mat
	//~ tmp.resize(mat[0].size());
	//~ for(size_t i = 0; i < tmp.size(); i++)
		//~ tmp[i].resize(mat.size());
//~ 
	//~ //transposing by switching the according indices
	//~ for(size_t i = 0; i < mat.size(); i++)
		//~ for(size_t j = 0; j < mat[0].size(); j++)
			//~ tmp[j][i] = mat[i][j];
//~ 
	//~ return tmp;
//~ }
//~ 
//~ /**
 //~ * This function multiplicates a matrix and a vector
 //~ * @mat: matrix of the product
 //~ * @vec: vector of the product
 //~ * @result: returned vector, represents the product of @mat and @vec
 //~ * @tmp: helper of summing up a row of @mat and the components of @vec
 //~ */
//~ vector<double> multmatvec(vector<vector<double> > mat, vector<double> vec){
//~ 
		//~ vector<double> result;
		//~ 
		//~ //setting @result's size to the number of rows in @mat
		//~ result.resize(mat.size());
		//~ double tmp = 0;
//~ 
		//~ //walking over every row of @mat
		//~ for(size_t i = 0; i < mat.size(); i++)
		//~ {
			//~ //Walking over every coloumn of @mat at a given row and walking over every component of @vec. The product is added to @tmp
			//~ for(size_t j = 0; j < vec.size(); j++)
				//~ tmp += mat[i][j] * vec[j];
//~ 
			//~ //the sum is assigned to one component of the result and tmp is resetted
			//~ result[i] = tmp;
			//~ tmp = 0;
		//~ }
//~ 
		//~ return result;
//~ 
//~ }

//~ /**
 //~ * This function calculates the product of two matrices
 //~ * @mat1, @mat2: the matrices that need to be multiplicated
 //~ * @result: result of the product of @mat1 & @mat2
 //~ * @tmp: temporary storage while calulating a component of @result
 //~ */
//~ vector<vector<double> > multmatmat(vector<vector<double> > mat1, vector<vector<double> > mat2){
//~ 
		//~ vector<vector<double> > result;
		//~ 
		//~ //setting the size of result by reading the #rows of @mat1 and the #coloumn of @mat2
		//~ result.resize(mat1.size());
//~ 
		//~ for(size_t i = 0; i < result.size(); i++)
			//~ result[i].resize(mat2.size());
//~ 
		//~ double tmp = 0;
//~ 
		//~ //walking over every row of @mat1
		//~ for(size_t i = 0; i < mat1.size(); i++)
		//~ {
			//~ //walking over every coloumn of @mat1
			//~ for(size_t j = 0; j < mat1[0].size(); j++)
			//~ {
				//~ //walking over every row of @mat2
				//~ for(size_t k = 0; k < mat2.size(); k++)
					//~ //summing up all the terms that result in a component of @result
					//~ tmp += mat1[i][k] * mat2[k][j];
//~ 
			//~ //the sum is assigned to one component of the result and tmp is resetted
			//~ result[i][j] = tmp;
			//~ tmp = 0;
			//~ }
		//~ }
//~ 
		//~ return result;
//~ }

//~ /**
 //~ * This function checks for the collinearity of a coloumn of the matrix @m in comparison with the previously coloumns.
 //~ * This leads to a small values at the corresponding diagonalterm in @r. The "smallness" is checked and, if it's true, the according fitparameter will be artificial set in order to
 //~ * keep the solution from diverging.
 //~ * @threshold: this value represents the threshold, below which a diagonalterm of @r is considered as small
 //~ * @i: this values represents the bin that is concerned
 //~ * @riihat, @bihat: helpervariables for setting of the fitparameter in @a
 //~ */
//~ void collinearity(double threshold, size_t i){
//~ 
	//~ //check if the value in @r is small
	//~ if(r[i][i] < threshold)
	//~ {
		//~ //setting the helpervariables for the creation of a smaller @a than one would get if the fit would run normally.
		//~ //the shift has it's source in @kappa and the according term of @d.
		//~ double riihat = sqrt(r[i][i] * r[i][i] + kappa * d[i][i]);
		//~ double bihat = (r[i][i] / riihat) * b[i];
		//~ a[i] = bihat / riihat;
	//~ }
//~ }
//~ 
//~ /**
 //~ * This function normalizes every coloumn of @m by using the absolut value of the respective coloumn.
 //~ * @coloumn: represents of the coloumn of interest
 //~ * @abs: absolut value of the coloumn
 //~ */
//~ void makemprime(size_t coloumn){
//~ 
	//~ //setting @mprime to the same size as @m
	//~ mprime.resize(m.size());
	//~ for(size_t i = 0; i < mprime.size(); i++)
		//~ mprime[i].resize(m[i].size());
//~ 
	//~ //setting the absolut value of the coloumn, because while assigning it rowwise, it changes after every assigned value
	//~ double abs = la.getabs(la.getcol(m, coloumn));
	//~ 
	//~ //setting the normalization
	//~ for(size_t i = 0; i < m.size(); i++)
		//~ mprime[i][coloumn] = m[i][coloumn] / abs;
//~ 
//~ }
//~ 
//~ /**
 //~ * This function adds a new coloumn to @m
 //~ * @power: matrix, that contains all the potencies of the values
 //~ * @x: matrix of the values
 //~ * @size: old number of coloumns in @m
 //~ * @tmp: temporary storage of the product of the different powers of the values
 //~ */
//~ void increasem(vector<vector<int> > power, vector<vector<double> > x){
//~ 
	//~ //set the size of @m, if not done yet
	//~ if(m.empty())
		//~ m.resize(x.size());
//~ 
	//~ size_t size = m[0].size();
	//~ for(size_t i = 0; i < m.size(); i++)
		//~ m[i].resize(size + 1);
//~ 
	//~ double tmp = 1;
	//~ //walking over every row of @m
	//~ for(size_t j = 0; j < m.size(); j++)
	//~ {	
		//~ //in every row the product of the values and the respective powers is calculated
		//~ for(size_t i = 0; i < x[0].size(); i++)
			//~ tmp *= pow(x[j][i], power[m[0].size() - 1][i]);
			//~ 
		//~ //the new terms will be added to the last coloumn in every row and @tmp is resetted
		//~ m[j][m[0].size() - 1] = tmp;
		//~ tmp = 1;
	//~ }
//~ 
//~ }

//~ /**
 //~ * This function calculates the new powers for every value
 //~ * @power: list of all powers for the respective values
 //~ * @size: number of values
 //~ * @max: the maximum sum of elements in a row of @power
 //~ * @sum: sum of all elements in a row of @power
 //~ * @start: represents the first row of @power that will be expanded
 //~ */
//~ void windofchange(vector<vector<int> >& power, size_t size, int max){
//~ 
	//~ //if the list is empty, the first entrys will be set 0 for all values
	//~ if(power.empty())
	//~ {
		//~ vector<int> tmp;
		//~ tmp.assign(size, 0);
		//~ power.push_back(tmp);
		//~ return;
	//~ }
//~ 
	//~ vector<vector<int> > tmp = power;
	//~ int sum = 0;
	//~ size_t start = 0;
	//~ 
	//~ //if there is already at least the first order polynomials in @power, the calculation can start later
	//~ //this construction exist only due to the fact, that the zeroth order is explicitly implemented
	//~ if(tmp.size() > 1)
		//~ start++;
//~ 
	//~ //walk over every value
	//~ for(size_t j = 0; j < tmp[0].size(); j++)
	//~ {
		//~ //walk down the rows, beginning at @start
		//~ for(size_t i = start; i < tmp.size(); i++)
		//~ {
			//~ //the respective power of a value that is already inside the list will be increased
			//~ tmp[i][j]++;
			//~ 
			//~ //the sum of the new powers will be calculated in order to check, if it sums up to the new order
			//~ for(size_t k = 0; k < tmp[0].size(); k++)
				//~ sum += tmp[i][k];
			//~ if(sum == max)
				//~ power.push_back(tmp[i]);
			//~ sum = 0;
			//~ 
			//~ //decreasing the element, so that there doesn't remain any changes in earlier part of @power
			//~ tmp[i][j]--;
//~ 
		//~ }
//~ 
		//~ //if the next value will get its power set, there would be redundandent entries in @power, therefore the calculation will start later
		//~ if(tmp.size() > 1)
			//~ start++;
//~ 
	//~ }
//~ 
//~ }

//~ /**
 //~ * This function normalizes a vector. If the length of the vector is 0, the norm will be the vector itself in order to prevent nan's.
 //~ * @vec: vector, that will be normalized
 //~ * @abs: absolute value of the vector
 //~ * @result: if the length is != 0, the normalized vector will be stored in this variable
 //~ */
//~ vector<double> normalizevec(vector<double> vec){
//~ 
	//~ //calculate the absolute value of @vec
	//~ double abs = getabs(vec);
	//~ 
	//~ //if the length of @vec is != 0, it will be normalized, else the vector itself will be returned 
	//~ if(abs != 0)
	//~ {
		//~ vector<double> result;
		//~ result.resize(vec.size());
//~ 
		//~ //calculating the normalized components of the vector
		//~ for(size_t i = 0; i < vec.size(); i++)
			//~ result[i] = vec[i] / abs; 
			//~ 
		//~ return result;
	//~ }
	//~ else
		//~ return vec;
	//~ 
//~ 
//~ }

//~ /**
 //~ * This function rescales the best fit parameters due to normalization of @m and @b
 //~ * @vec: vector, containing the best fit parameters
 //~ */
//~ void rescalebestfitparameters(vector<double>& vec){
//~ 
	//~ //backwards calculating of the norms that influenzed @m and @b in order to get the regular best fit parameters
	//~ for(size_t i = 0; i < vec.size(); i++)
		//~ vec[i] *= la.getabs(b) / la.getabs(la.getcol(m,i));
	//~ 
//~ }

//~ /**
 //~ * This function solves a problem of the type A*x=b where x is a vector containing the parameters to fit
 //~ * @a: matrix A
 //~ * @x: vector of parameters to fit
 //~ * @b: vector b
 //~ * @result: resulting parameters
 //~ */
//~ vector<double> getbestfitparameters(vector<vector<double> > a, vector<double> x, vector<double> b){
//~ 
	//~ //setting the size of @result
	//~ vector<double> result;
	//~ result.resize(x.size());
//~ 
	//~ //starting at the bottom line, walking upwards and calculating the parameters
	//~ for(size_t i = x.size(); i > 0; i--)
	//~ {
		//~ 
		//~ //nan is an indicator in a component of @x for a component, that needs to be calculated
		//~ if(std::isnan(x[i - 1])){
			//~ //calculating the respective parameter
			//~ result[i - 1] = b[i - 1] / a[i - 1][i - 1];
			//~ 
			//~ //forward the solution to every row above the regarded one
			//~ for(size_t j = 0; j < i - 1; j++)
				//~ b[j] -= a[j][i - 1] * result[i - 1];
//~ 
		//~ }
		//~ else{
			//~ //If a values isn't nan, it was set by @collinearity(). It's value will be forwarded as above.
			//~ result[i - 1] = x[i - 1];
			//~ for(size_t j = 0; j < i - 1; j++)
				//~ b[j] -= a[j][i - 1] * result[i - 1];
		//~ }
//~ 
	//~ }
//~ 
	//~ 
	//~ return result;
//~ 
//~ }

//~ /**
 //~ * This function calculates the normalized distancevector between the fitfunction and the datapoints at a certain point.
 //~ * @x: list of all parametervalues
 //~ * @power: list of all powers for the respective monomial
 //~ * @bfm: best fit parameters
 //~ * @i: # of datapoint at which its distance will be calculated
 //~ * @functionvalue: functionvalue at the given point, represented by @i
 //~ * @functiongradient: gradient of the fitfunction at the given point; used for normalization
 //~ * @tmp: temporary storage while calculating the parts of a monomial
 //~ */
//~ double getei(vector<vector<double> > x, vector<vector<int> > power, vector<double> bfm, size_t i){
//~ 
	//~ double functionvalue = 0;
	//~ 
	//~ //taking the monomials, saved in @m and multiply it by the respective fitparameter
	//~ for(size_t j = 0; j < bfm.size(); j++)
		//~ functionvalue += m[i][j] * bfm[j];
		//~ 
	//~ //use the absolute value of the difference of the functionvalue and the datapoint
	//~ functionvalue = fabs(functionvalue - b[i]);
//~ 
	//~ //setting up the gradient of the function
	//~ vector<double> functiongradient;
	//~ functiongradient.assign(x[0].size(), 0);
	//~ double tmp = 1;
//~ 
	//~ //walking over every component of the gradient
	//~ for(size_t j = 0; j < functiongradient.size(); j++)
	//~ {
		//~ //walk over every fitparameter
		//~ for(size_t k = 1; k < bfm.size(); k++) // leave out constant offset
		//~ {
			//~ //walk over every value
			//~ for(size_t l = 0; l < x[0].size(); l++)
				//~ //if the component of @functiongradient and @l match, the respective term will be derived
				//~ if(j == l)
				//~ {
					//~ //setting a zero term to zero and therefore the whole product too leads to stability against the function @pow() itself
					//~ if(x[i][l] == 0)
						//~ tmp = 0;
					//~ else
						//~ //derive the factor x^n in the monomial by simply calculate n*x^(n - 1)
						//~ tmp *= pow(x[i][l], power[k][l] - 1) * power[k][l];
				//~ }
				//~ else
					//~ //every other factor will be left as it is
					//~ tmp *= pow(x[i][l], power[k][l]);
		//~ 
			//~ //finally multiply the prefactor to @tmp and add it as one monomial to the respective component in the gradient
			//~ tmp *= bfm[k];
			//~ functiongradient[j] += tmp;
			//~ tmp = 1;
		//~ }
	//~ }
	//~ 
	//~ //if the length of the gradient is zero, @functionvalue itself will be returned, otherwise it will be normalized
	//~ if(la.getabs(functiongradient) == 0)
		//~ return functionvalue;
	//~ else
		//~ if(functionvalue > 0)
			//~ return functionvalue / la.getabs(functiongradient);
		//~ else
			//~ return -functionvalue / la.getabs(functiongradient);
//~ }
//~ 
//~ /**
 //~ * This function calculates the mean of the overall distance between the data and the fitfunction as one breakoff criteria
 //~ * @ei: vector, containing every distance vector
 //~ * @sum: sum of all differences
 //~ */
//~ double getDdist(vector<double> ei){
//~ 
	//~ double sum = 0;
	//~ 
	//~ //sum up all differences
	//~ for(size_t i = 0; i < ei.size(); i++)
		//~ sum += ei[i];
//~ 
	//~ //return the mean
	//~ return sum / ei.size();
//~ 
//~ }

//~ double getchi2(vector<double> bfm, vector<double> sigma, vector<double> bfmerr){
//~ 
	//~ double result = 0, functionvalue = 0, fiterr = 0;
	//~ 
	//~ cout << "chi2 berechnung:" << endl;
	//~ for(size_t i = 0; i < m.size(); i++)
	//~ {
		//~ for(size_t j = 0; j < m[0].size(); j++)
			//~ functionvalue += m[i][j] * bfm[j];
		//~ 
		//~ for(size_t k = 0; k < bfmerr.size(); k++)
			//~ if(std::isnan(bfmerr[k]))
				//~ continue;
			//~ else
				//~ fiterr += m[i][k] * m[i][k] * bfmerr[k] * bfmerr[k];
				//~ 
		//~ functionvalue -= b[i];
		//~ cout << functionvalue * functionvalue << "\t" << sigma[i] * sigma[i] << endl;
		//~ result += (functionvalue * functionvalue) / (sigma[i] * sigma[i] + fiterr);
		//~ 
		//~ functionvalue = 0;
		//~ fiterr = 0;
	//~ }
	//~ cout << endl;
	//~ 
	//~ return result;
//~ }


//~ /**
 //~ * This function calculates the Chi^2 value of the fit.
 //~ * @bfm: best fit parameters
 //~ * @sigma: errors of the datapoints
 //~ * @result: resulting Chi^2
 //~ * @functionvalue: functionvalue of the fit at a certain point
 //~ */
//~ double getchi2(vector<double> bfm, vector<double> sigma){
//~ 
	//~ double result = 0, functionvalue = 0;
	//~ 
	//~ //walk over every polynomial
	//~ for(size_t i = 0; i < m.size(); i++)
	//~ {
		//~ //calculate the monomials and add them to the result
		//~ for(size_t j = 0; j < m[0].size(); j++)
			//~ functionvalue += m[i][j] * bfm[j];
		//~ 
		//~ //get the difference between the functionvalue and the datapoint		
		//~ functionvalue -= b[i];
		//~ 
		//~ //if the uncertainty is 0, the calculation would rise a nan
		//~ //this is prevented by setting it to the arbitrary value of 1e-10
		//~ if(sigma[i] == 0)
			//~ sigma[i] = 1e-10;
		//~ 
		//~ //calculating the squares of @functionvalue and @sigma, divide them and add this new summand to new overall Chi^2 value @result
		//~ result += (functionvalue * functionvalue) / (sigma[i] * sigma[i]);
		//~ 
		//~ functionvalue = 0;
		//~ 
	//~ }
//~ 
	//~ 
	//~ return result;
//~ }

//~ /**
 //~ * This function calculates the distance between two vectors
 //~ * @a, @b: vectors of interest
 //~ */
//~ double getdistanceofvectors(vector<double> a, vector<double> b){
//~ 
	//~ //calculate the difference componentwise
	//~ for(size_t i = 0; i < a.size(); i++)
		//~ a[i] = a[i] - b[i];
		//~ 
	//~ //return its absolut value
	//~ return getabs(a);
	//~ 
//~ }

//~ /**
 //~ * This function constructs a hypercube around a certain point
 //~ * @x: matrix that contains every point
 //~ * @i: number of the point in @x around which a hypercube should be constructed
 //~ * @center: vector of the center of the hypercub around which the cube is constructed
 //~ * @result: list of numbers that indicate the taken points for the hypercube
 //~ * @distances: list that contains all the distances between @center and every other point
 //~ * @div: quotient of the differences between points
 //~ * @flag: flag that will be set if 2 datapoints carry the same information for a hypercube
 //~ */
//~ vector<size_t> gethypercube(vector<vector<double> > x, size_t i){
//~ 
	//~ //setting up @center and the resultvector
	//~ vector<double> center = x[i];
	//~ vector<size_t> result;
	//~ //setting the resultsize as 2 * #dimensions, so that for every dimension at least 2 points cann be selected
	//~ result.resize(x[0].size() * 2);
	//~ //if a point isn't set, it's value is #number of dimensions + 1, so that these points can be found and won't be used for further calculations
	//~ //this is important for points at the border of the sampleregion
	//~ result.assign(result.size(), x.size());
	//~ 
	//~ //The same procedure as before is done for the distances, so that a check is possible if the smallest possible hypercube cann be constructed. Therefore it's initialvalues is set to infinity.
	//~ vector<double> distances;
	//~ distances.resize(x[0].size() * 2);
	//~ distances.assign(distances.size(), numeric_limits<double>::infinity());
	//~ 
	//~ //walk over every datapoint
	//~ for(size_t j = 0; j < x.size(); j++)
		//~ //if a datapoint matches the @center, it will be skipped
		//~ if(i == j)
			//~ continue;
		//~ else
			//~ //walk over every component of a point
			//~ for(size_t k = 0; k < x[0].size(); k++)
			//~ {
				//~ //if the component is bigger than the @center and closer than the closest at that point, it will be set instead
				//~ //the bigger points will be set here as the even indices, the smaller as the odd
				//~ //the break statement leads to not using the same datapoint in several components (worst case: hypercube with 2 datapoints)
				//~ if(center[k] > x[j][k] && getdistanceofvectors(center, x[j]) < distances[2 * k])
				//~ {
					//~ 
					//~ distances[2 * k] = getdistanceofvectors(center, x[j]);
					//~ result[2 * k] = j;
					//~ break;
				//~ }
				//~ 
				//~ //same as above, but for smaller values
				//~ if(center[k] < x[j][k] && getdistanceofvectors(center, x[j]) < distances[2 * k + 1])
				//~ {
					//~ 
					//~ distances[2 * k + 1] = getdistanceofvectors(center, x[j]);
					//~ result[2 * k + 1] = j;
					//~ break;
				//~ }
			//~ 
			//~ }
//~ 
	//~ //filtering selected points that are in a linear relationship with another selected point and therefore don't carry any new informations
	//~ //this will only be considered if there is more than one dimension
	//~ if(x[0].size() > 1)
	//~ {
		//~ double div;
		//~ int flag = 0;
		//~ //check every combination of points
		//~ for(size_t j = 0; j < result.size(); j++)
			//~ for(size_t k = 0; k < result.size(); k++)
				//~ //if they are the same points of no point was selected, these values will be skipped
				//~ if(j == k || result[j] == x.size() || result[k] == x.size())
					//~ continue;
				//~ else
				//~ {
					//~ //calculate the ratio of the first component of the differences of the points to the center
					//~ div = (x[result[j]][0] - center[0]) / (x[result[k]][0] - center[0]);
					//~ 
					//~ //if both points lie on the same side of the center, @div is >0, else it isn't and the point remains of interest
					//~ if(div > 0)
					//~ {
						//~ //set up a "ready to delete this point"-flag
						//~ flag = 1;
						//~ //walk over every other component of the points
						//~ for(size_t l = 1; l < x[0].size(); l++)
						//~ {
							//~ //if these points have at least another quotient in one component, they remain in the list
							//~ if(div != (x[result[j]][l] - center[l]) / (x[result[k]][l] - center[l])) //strahlensatz
							//~ {
								//~ flag = 0;
								//~ break;
							//~ }
						//~ 
						//~ }
						//~ //if both points lie on a straight line, one is deleted
						//~ if(flag)
							//~ result.erase(result.begin() + k);
						//~ flag = 0;
					//~ }
				//~ }
				//~ 
	//~ }
//~ 
	//~ 
	//~ return result;
//~ }

//~ /**
 //~ * This function does the same as @expandqr() but for a fit of a hypercube.
 //~ * @mlocal, @qlocal, @rlocal: analog to @m, @q, @r respectively
 //~ * @i: counter of the respective iteration
 //~ */
//~ void expandqrlocal(vector<vector<double> >& mlocal, vector<vector<double> >& qlocal, vector<vector<double> >& rlocal, size_t i){
	//~ if(qlocal.empty() || rlocal.empty())
	//~ {
	//~ 
		//~ qlocal.resize(mlocal.size());
		//~ rlocal.resize(mlocal[0].size());
			//~ 
		//~ for(size_t k = 0; k < qlocal.size(); k++)
		//~ {
			//~ qlocal[k].resize(mlocal[0].size());
		//~ }
		//~ 
		//~ for(size_t k = 0; k < rlocal.size(); k++)
		//~ {
			//~ rlocal[k].resize(mlocal[0].size());
		//~ }
		//~ 
		//~ for(size_t k = 0; k < qlocal.size(); k++)
			//~ qlocal[k][0] = la.getcol(mlocal,0)[k] / la.getabs(la.getcol(mlocal,0));
//~ 
		//~ rlocal[0][0] = la.getabs(la.getcol(mlocal,0));
		//~ return;
	//~ }
//~ 
//~ 
	//~ qlocal.resize(mlocal.size());
	//~ rlocal.resize(mlocal[0].size());
//~ 
//~ 
	//~ for(size_t k = 0; k < qlocal.size(); k++)
	//~ {
	//~ qlocal[k].resize(mlocal[0].size());
	//~ }
	//~ 
	//~ for(size_t k = 0; k < rlocal.size(); k++)
	//~ {	
	//~ rlocal[k].resize(mlocal[0].size());
	//~ }
	//~ 
//~ 
	//~ double tmp;
	//~ for(size_t l = 0; l <= i; l++)
	//~ {
		//~ tmp = 0;
		//~ for(size_t k = 0 ; k < qlocal.size(); k++)
		//~ {
			//~ tmp += qlocal[k][l] * mlocal[k][i + 1];
//~ 
		//~ }
		//~ rlocal[l][i + 1] = tmp;
	//~ }
//~ 
//~ 
	//~ vector<double> sum;
	//~ sum.assign(qlocal.size(), 0);
	//~ for(size_t l = 0; l < sum.size(); l++)
		//~ for(size_t k = 0; k <= i; k++)
			//~ sum[l] += rlocal[k][i + 1] * qlocal[l][k];
//~ 
//~ 
	//~ for(size_t k = 0; k < qlocal.size(); k++)
		//~ qlocal[k][i + 1] = mlocal[k][i + 1] - sum[k];
		//~ 
	//~ 
//~ 
	//~ rlocal[i + 1][i + 1] = la.getabs(la.getcol(qlocal,i + 1));
	//~ 
	//~ tmp = la.getabs(la.getcol(qlocal, i + 1));
	//~ if(tmp != 0)
		//~ for(size_t k = 0; k < qlocal.size(); k++)
			//~ qlocal[k][i + 1] /= tmp;
			//~ 
//~ }

//~ /**
 //~ * This function calculates the best fit parameters for the hypercubes
 //~ * @points: matrix that contains all points
 //~ * @bpoints: vector with the functionvalues
 //~ * @threshold: threshold for artificial fitparameter setting
 //~ * @power: list containing the powers for every value at every monomial
 //~ * @result: container for the best fit parameters
 //~ * @mlocal, @qlocal, @rlocal: local equivalent to @m, @q, @r respectively
 //~ * @distances: vector for weighting the hypercubepoints with its distance to the center
 //~ * @riihat, @bihat: regularization as in @collinearity()
 //~ * @dlocal: local equivalent to @d
 //~ */
//~ vector<double> getfitparams(vector<vector<double> > points, vector<double> bpoints, double threshold, int power){
//~ 
	//~ //setting up the parametervector, setting nans for finding paramaters to fit
	//~ vector<double> result;
	//~ result.resize(points[0].size() + 1);
	//~ result.assign(result.size(), nan("1")); //muss gesetzt werden fuer fitfunktion, sonst annahme, dass wert bereits gesetzt ist
	//~ vector<vector<double> > mlocal, qlocal, rlocal;
	//~ mlocal.resize(points.size());	
	//~ 
	//~ //calculating the distances of the points of the hypercube to the center
	//~ vector<double> distances;
	//~ 
	//~ distances.resize(points.size());
	//~ 
	//~ //the zeroth component is set to 1 and therefore it won't be weighted otherwise
	//~ distances[0] = 1.;
	//~ 
	//~ //calculate every distance and increase its "importance decrease by distance" by the power of @power
	//~ for(size_t i = 1; i < distances.size(); i++)
		//~ distances[i] = pow(getdistanceofvectors(points[0], points[i]), power);
	//~ 
	//~ //rescale the functionvalues by the new weights
	//~ bpoints[0] /= distances[0];
	//~ for(size_t i = 1; i < bpoints.size(); i++)
		//~ bpoints[i] /= distances[i];
	//~ 
		//~ 
	//~ //build an rescale @mlocal as in @increasem() etc.
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
	//~ 
	//~ //regularize the QR-decomposition as in @collinearity()
	//~ double riihat, bihat;
	//~ vector<vector<double> > dlocal;
	//~ dlocal.resize(rlocal.size());
	//~ for(size_t i = 0; i < dlocal.size(); i++)
		//~ dlocal[i].assign(rlocal.size(), 0.);
		//~ 
	//~ for(size_t i = 0; i < dlocal.size(); i++)
		//~ dlocal[i][i] = 1.;
		//~ 
		//~ 
	//~ for(size_t i = 0; i < rlocal.size(); i++)
		//~ if(rlocal[i][i] < threshold)
		//~ {
//~ 
			//~ riihat = sqrt(rlocal[i][i] * rlocal[i][i] + kappa * 1.);
			//~ bihat = (rlocal[i][i] / riihat) * bpoints[i];
			//~ result[i] = bihat / riihat;
		//~ }
	//~ 
	//~ //return the fit
	//~ return getbestfitparameters(rlocal, result, multmatvec(transpose(qlocal), bpoints));
	//~ 
	//~ 
//~ }

//~ vector<double> getnormalvec(vector<size_t>& hypercube, vector<vector<double> >& x, vector<double> b, size_t i, double threshold, int power){
//~ 
	//~ vector<double> fitparams, bpoints;
	//~ fitparams.resize(x[0].size() + 1);
//~ 
	//~ vector<vector<double> > points;
	//~ points.push_back(x[i]);
	//~ bpoints.push_back(b[i]);
//~ 
	//~ for(size_t j = 0; j < hypercube.size(); j++)
		//~ if(hypercube[j] < x.size())
		//~ {
			//~ points.push_back(x[hypercube[j]]);
			//~ bpoints.push_back(b[hypercube[j]]);
		//~ }
		//~ 
//~ 
	//~ fitparams = getfitparams(points, bpoints, threshold, power);
//~ 
	//~ 
	//~ fitparams.erase(fitparams.begin());
	//~ double length = la.getabs(fitparams);
	//~ if(length != 0)
		//~ for(size_t j = 0; j < fitparams.size(); j++)
			//~ fitparams[j] /= length;
		//~ 
	//~ cout << "normal vector:" << endl;
	//~ for(size_t j = 0; j < fitparams.size(); j++)
		//~ cout << fitparams[j] << "\t";
	//~ cout << endl;
	//~ 
	//~ return fitparams;
	//~ 
//~ }


//~ vector<vector<double> > calculatenormalvectors(vector<vector<double> > x, vector<double> b, vector<double> sigma, double threshold, int power){
//~ 
	//~ vector<vector<double> > normalvec;
	//~ normalvec.resize(x.size());
	//~ for(size_t i = 0; i < normalvec.size(); i++)
		//~ normalvec[i].resize(x[i].size() + 1);
	//~ 
		//~ 
	//~ vector<size_t> hypercube;
	//~ hypercube.resize(x[0].size() * 2);
		//~ 
	//~ for(size_t i = 0; i < x.size(); i++){
		//~ 
		//~ hypercube = gethypercube(x, i);
//~ 
		//~ normalvec[i] = getnormalvec(hypercube, x, b, i, threshold, power);
		//~ 
	//~ }
	//~ 
	//~ 
	//~ return normalvec;
//~ }

//~ vector<double> normalvecfunction(vector<vector<double> > x, vector<double> bfm, vector<vector<int> > power, size_t i){
	//~ 
	//~ vector<double> functiongradient;
	//~ functiongradient.assign(x[0].size(), 0);
	//~ double tmp = 1;
	//~ 
	//~ for(size_t j = 0; j < functiongradient.size(); j++)
	//~ {
		//~ for(size_t k = 1; k < bfm.size(); k++)
		//~ {
			//~ for(size_t l = 0; l < x[0].size(); l++)
				//~ if(j == l)
				//~ {
					//~ if(x[i][l] == 0)
						//~ tmp = 0;
					//~ else
						//~ tmp *= pow(x[i][l], power[k][l] - 1) * power[k][l];	
					//~ 
				//~ }
				//~ else
					//~ tmp *= pow(x[i][l], power[k][l]);
					//~ 
			//~ tmp *= bfm[k];
			//~ functiongradient[j] += tmp;
			//~ 
			//~ tmp = 1;
		//~ }
	//~ }
	//~ 
	//~ tmp = la.getabs(functiongradient);
	//~ 
	//~ if(tmp != 0)
		//~ for(size_t i = 0; i < functiongradient.size(); i++)
			//~ functiongradient[i] /= tmp;
	//~ 
	//~ return functiongradient;
//~ }

//~ double dotproduct(vector<double> a, vector<double> b){
	//~ 
	//~ if(getabs(a) == 0 && getabs(b) == 0)
		//~ return 1.;
		//~ 
	//~ double result = 0;
	//~ 
	//~ for(size_t i = 0; i < a.size(); i++)
		//~ result += a[i] * b[i];
	//~ 
	//~ return result;
//~ 
	//~ 
//~ }

//~ double getDsmooth(vector<vector<double> > normalvectors, vector<vector<double> > x, vector<double> bfm, vector<vector<int> > power){
	//~ 
	//~ double result = 0;
		//~ 
	//~ for(size_t i = 0; i < x.size(); i++)
		//~ result += la.dotproduct(normalvectors[i], normalvecfunction(x, bfm, power, i));
		//~ 
//~ 
	//~ return result / (double) x.size();
//~ 
//~ }
//~ 
//~ vector<double> getfiterrors(vector<double> sigma, double threshold){
//~ 
	//~ vector<double> result, err;
	//~ result.resize(m[0].size());
	//~ result.assign(result.size(), nan("1")); //muss gesetzt werden fuer fitfunktion, sonst annahme, dass wert bereits gesetzt ist
	//~ err.resize(m[0].size());
//~ 
	//~ vector<vector<double> > merr, rerr, qerr;
	//~ merr.resize(m.size());
//~ 
	//~ 
	//~ for(size_t i = 0; i < m.size(); i++)
	//~ {	
		//~ merr[i].push_back(m[i][0] * m[i][0]);			
	//~ }
	//~ expandqrlocal(merr, qerr, rerr, 0);
//~ 
	//~ for(size_t i = 1; i < m[0].size(); i++)
	//~ {
		//~ for(size_t j = 0; j < m.size(); j++)
			//~ merr[j].push_back(m[j][i] * m[j][i]);
			//~ 
		//~ expandqrlocal(merr, qerr, rerr, i - 1);
		//~ 
	//~ }
	//~ 
	//~ for(size_t i = 0; i < err.size(); i++)
		//~ err[i] = sigma[i] * sigma[i];
	//~ 
	//~ double riihat, bihat;
	//~ for(size_t i = 0; i < rerr.size(); i++)
		//~ if(rerr[i][i] < threshold)
		//~ {
			//~ riihat = sqrt(rerr[i][i] * rerr[i][i] + kappa * d[i][i]);
			//~ bihat = (rerr[i][i] / riihat) * err[i];
			//~ result[i] = bihat / riihat;
		//~ }
		//~ 
			//~ 
	//~ 
	//~ result = la.getbestfitparameters(rerr, result, multmatvec(transpose(qerr), err));
	//~ 

	//~ 
	//~ for(size_t i = 0; i < result.size(); i++)
		//~ if(result[i] >= 0)
			//~ result[i] = sqrt(result[i]);
		//~ else
			//~ result[i] = sqrt(-result[i]);
		//~ 
	//~ return result; 
//~ }

int main() {

	//parameters:
	double threshold, thresholddata, thresholderr, ddist, chi2mean, dsmoothmean;
	int exponent;
	//double dsmooth;
	
	threshold = 1e-10;
	thresholddata = 1e-10;
	thresholderr = 1e-10;
	ddist = 0.05;
	//dsmooth = 0.5;
	chi2mean = 15;
	dsmoothmean = 15;
	kappa = 0.05;
	exponent = 3;
	

	vector<vector<double> > x;
	datareading dr("interpolation0");
	refreading rr(dr);
	size_t num_analysis = transpose(rr.bin_values).size();
	vector<double> sigma;
	x = rr.parameter_values;
	
	
	
	vector<vector<int> > power;
	int max, iterationcounter, quitflagchi2, quitflagdsmooth;
	vector<double> bfm, ei;
	vector<vector<double> > normalvectors;
	vector<double> bfmerr;	
	vector<double> chi2values, dsmoothvalues;
	double sum;
	ofstream output, outputchi2;
	outputchi2.open("summary");
	
	for(size_t num_ipol = 0; num_ipol < num_analysis; num_ipol++)
	{
		
	power.clear();
	bfm.clear();
	ei.clear();
	normalvectors.clear();
	bfmerr.clear();
	chi2values.clear();
	dsmoothvalues.clear();
	m.clear();
	r.clear();
	q.clear();
	d.clear();
	mprime.clear();
		
	b = transpose(rr.bin_values)[num_ipol];

	sigma = transpose(rr.bin_values_err)[num_ipol];
	


	//~ vector<double> tmp;
	//~ tmp.push_back(-1.);
	//~ tmp.push_back(-1.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(-2.);
	//~ tmp.push_back(-2.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(-3.);
	//~ tmp.push_back(-3.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(-4.);
	//~ tmp.push_back(-4.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(-5.);
	//~ tmp.push_back(-5.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(-6.);
	//~ tmp.push_back(-6.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(1.);
	//~ tmp.push_back(1.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(2.);
	//~ tmp.push_back(2.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(3.);
	//~ tmp.push_back(3.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(4.);
	//~ tmp.push_back(4.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(5.);
	//~ tmp.push_back(5.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(6.);
	//~ tmp.push_back(6.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(-1.);
	//~ tmp.push_back(1.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(-2.);
	//~ tmp.push_back(2.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(-3.);
	//~ tmp.push_back(3.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(-4.);
	//~ tmp.push_back(4.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(-5.);
	//~ tmp.push_back(5.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(-6.);
	//~ tmp.push_back(6.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(1.);
	//~ tmp.push_back(-1.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(2.);
	//~ tmp.push_back(-2.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(3.);
	//~ tmp.push_back(-3.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(4.);
	//~ tmp.push_back(-4.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(5.);
	//~ tmp.push_back(-5.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(6.);
	//~ tmp.push_back(-6.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(0.);
	//~ tmp.push_back(1.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(0.);
	//~ tmp.push_back(2.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(0.);
	//~ tmp.push_back(3.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(0.);
	//~ tmp.push_back(4.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(0.);
	//~ tmp.push_back(5.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(0.);
	//~ tmp.push_back(6.);
	//~ x.push_back(tmp);
	//~ tmp.clear();
	//~ tmp.push_back(0.);
	//~ tmp.push_back(0.);
	//~ x.push_back(tmp);
	//~ tmp.clear();


	//~ b.push_back(2.);
	//~ b.push_back(8.);
	//~ b.push_back(18.);
	//~ b.push_back(32.);
	//~ b.push_back(50.);
	//~ b.push_back(72.);
	//~ b.push_back(2.);
	//~ b.push_back(8.);
	//~ b.push_back(18.);
	//~ b.push_back(32.);
	//~ b.push_back(50.);
	//~ b.push_back(72.);
	//~ b.push_back(2.);
	//~ b.push_back(8.);
	//~ b.push_back(18.);
	//~ b.push_back(32.);
	//~ b.push_back(50.);
	//~ b.push_back(72.);
	//~ b.push_back(2.);
	//~ b.push_back(8.);
	//~ b.push_back(18.);
	//~ b.push_back(32.);
	//~ b.push_back(50.);
	//~ b.push_back(72.);
	//~ b.push_back(1.);
	//~ b.push_back(4.);
	//~ b.push_back(9.);
	//~ b.push_back(16.);
	//~ b.push_back(25.);
	//~ b.push_back(36.);
	//~ b.push_back(0.);

	//~ sigma.push_back(1.);
	//~ sigma.push_back(1.4);
	//~ sigma.push_back(1.7);
	//~ sigma.push_back(2.);
	//~ sigma.push_back(2.2);
	//~ sigma.push_back(2.6);
	//~ sigma.push_back(1.);
	//~ sigma.push_back(1.4);
	//~ sigma.push_back(1.7);
	//~ sigma.push_back(2.);
	//~ sigma.push_back(2.2);
	//~ sigma.push_back(2.6);
	//~ sigma.push_back(1.);
	//~ sigma.push_back(1.4);
	//~ sigma.push_back(1.7);
	//~ sigma.push_back(2.);
	//~ sigma.push_back(2.2);
	//~ sigma.push_back(2.6);
	//~ sigma.push_back(1.);
	//~ sigma.push_back(1.4);
	//~ sigma.push_back(1.7);
	//~ sigma.push_back(2.);
	//~ sigma.push_back(2.2);
	//~ sigma.push_back(2.6);
	//~ sigma.push_back(1.);
	//~ sigma.push_back(1.4);
	//~ sigma.push_back(1.7);
	//~ sigma.push_back(2.);
	//~ sigma.push_back(2.2);
	//~ sigma.push_back(2.6);
	//~ sigma.push_back(1.);
	
	
	max = 0;
	iterationcounter = 0;
	
	
	windofchange(power, x[0].size(), max);
	increasem(power, x);
	makemprime(iterationcounter);
	a.resize(mprime[0].size());
	expandqr(0);
	bprime = la.normalizevec(b);
	bprime = la.multmatvec(la.transpose(q),bprime);


	d.resize(m.size());
	for(size_t i = 0; i < d.size(); i++)
		d[i].assign(d.size(), 0);
	for(size_t i = 0; i < d.size(); i++)
		d[i][i] = 1;
		
		
	for(size_t i = 0; i < a.size(); i++)
	{
		a[i] = nan("1");
		collinearity(threshold, i);
	}


	bfm = la.getbestfitparameters(r, a, bprime);
	rescalebestfitparameters(bfm);

	
	for(size_t i = 0; i < x.size(); i++)
		ei.push_back(getei(x, power, bfm, i));

	bfmerr = getfiterrors(sigma, thresholderr);

	normalvectors = calculatenormalvectors(x, b, sigma, thresholddata, exponent);


	
	//convergence check
	chi2values.push_back(getchi2(bfm, sigma));
	sum = 0.;
	quitflagchi2 = 1;
	quitflagdsmooth = 1;

	cout << "setted a:\t";
	for(size_t i = 0; i < a.size(); i++)
		if(!std::isnan(a[i])) cout << i << "\t";
	cout << endl;
	

	
	cout << "Ddist:\t" << getDdist(ei) << endl;
	cout << "Dsmooth:\t" << getDsmooth(normalvectors, x, bfm, power) << endl;
	cout << "chi2:\t" << getchi2(bfm, sigma) << endl;
	cout << "iterationcounter:\t" << iterationcounter << endl;
	cout << endl;

	
/////////////////////////////////////////////////////////////////////////////////////////////////////
	
	while(getDdist(ei) > ddist || quitflagdsmooth || quitflagchi2)
	//~ while(getDsmooth(normalvectors,x,bfm,power) < dsmooth || getchi2(bfm,sigma, bfmerr) > chi2)
	{
	if(m[0].size() == power.size())
	{
		max++;
		windofchange(power, x[0].size() ,max);
	}

	increasem(power, x);
	iterationcounter++;
	makemprime(iterationcounter);
	

	a.resize(mprime[0].size());

	expandqr(iterationcounter - 1);


	
	bprime = normalizevec(b);

	//~ cout << "brime:" << endl;
	//~ for(size_t i = 0; i < bprime.size(); i++)
		//~ cout << bprime[i] << "\t";
	//~ cout << endl;
	
		//~ cout << "Q:" << endl;
	//~ for(size_t i = 0; i < q.size(); i++)
	//~ {
		//~ for(size_t j = 0; j < q[0].size(); j++)
			//~ cout << q[i][j] << "\t";
		//~ cout << endl;
	//~ }
	
	bprime = multmatvec(transpose(q),bprime);
		
	
	d.resize(m[0].size());
	for(size_t i = 0; i < d.size(); i++)
		d[i].assign(d.size(), 0);
	for(size_t i = 0; i < d.size(); i++)
		d[i][i] = 1;
	for(size_t i = 0; i < a.size(); i++)
	{
		a[i] = nan("1");
		collinearity(threshold, i);
	}
	
	cout << "setted a\t";
	for(size_t i = 0; i < a.size(); i++)
		if(!std::isnan(a[i])) cout << i << "\t";
	cout << endl;
	
	bfm = la.getbestfitparameters(r, a, bprime);
	rescalebestfitparameters(bfm);
	
	//~ cout << "bfm" << endl;
	//~ for(size_t i = 0; i < bfm.size(); i++)
		//~ cout << bfm[i] << "\t";
	//~ cout << endl;

	ei.clear();
	for(size_t i = 0; i < x.size(); i++)
		ei.push_back(getei(x, power, bfm, i));


	chi2values.push_back(getchi2(bfm, sigma));
	dsmoothvalues.push_back(getDsmooth(normalvectors, x, bfm, power));
	
	quitflagchi2 = 1;
	quitflagdsmooth = 1;
	
	if(chi2values.size() >= chi2mean)
	{
		for(size_t i = 0; i < chi2values.size(); i++)
			sum += chi2values[i];
			
		sum /= chi2values.size();
	
		if(chi2values[chi2values.size() - 1] >= sum)
			quitflagchi2 = 0;
	
	
		chi2values.erase(chi2values.begin());
	}
	
	sum = 0;
	
	if(dsmoothvalues.size() >= dsmoothmean)
	{
		for(size_t i = 0; i < dsmoothvalues.size(); i++)
			sum += dsmoothvalues[i];
			
		sum /= dsmoothvalues.size();
	
		if(dsmoothvalues[dsmoothvalues.size() - 1] <= sum)
			quitflagdsmooth = 0;
	
	
		dsmoothvalues.erase(dsmoothvalues.begin());
	}
	
	sum = 0;
	
	bfmerr = getfiterrors(sigma, thresholderr);
	
	cout << "Ddist:\t" << getDdist(ei) << endl;
	cout << "Dsmooth:\t" << getDsmooth(normalvectors, x, bfm, power) << endl;
	cout << "chi2:\t" << getchi2(bfm, sigma) << endl;
	cout << "iterationcounter\t" << iterationcounter << endl;
	cout << endl;

	if(iterationcounter > 10)
	{
		//break;
		ddist += ddist * 0.01;
		//dsmooth -= dsmooth * 0.01;
	}
		
	
	}
	
	//cout << "ddist:\t" << ddist << "\tdsmooth:\t" << dsmoothmean << endl << endl;
	
	
	cout << "bfm:\t";
	for(size_t i = 0; i < bfm.size(); i++)
		cout << bfm[i] << "\t";
	cout << endl << endl;
	
	cout << "errors:\t";
	for(size_t i = 0; i < bfmerr.size(); i++)
		cout << bfmerr[i] << "\t";
	cout << endl;
	
	
	output.open(("interpolation/" + dr.analysis[num_ipol] + "_" + dr.observable[num_ipol] + "_" + to_string(num_ipol)).c_str());
	
	output << "setted a:\t";
	for(size_t i = 0; i < a.size(); i++)
		if(!std::isnan(a[i]))
			output << i << "\t";
	output << endl;
	
	output << "Ddist:\t" << getDdist(ei) << endl;
	output << "Dsmooth:\t" << getDsmooth(normalvectors, x, bfm, power) << endl;
	output << "chi2:\t" << getchi2(bfm, sigma) << endl;
	output << "# iteration:\t" << iterationcounter << endl;
	
	output << "params:\t";
	for(size_t i = 0; i < bfm.size(); i++)
		output <<  bfm[i] << "\t";
	output << endl;
	output << "errors:\t";
	for(size_t i = 0; i < bfmerr.size(); i++)
		output << bfmerr[i] << "\t";
	output.close();
	
	
	
	outputchi2 << getchi2(bfm, sigma) << "\t" << iterationcounter << endl;
}

	outputchi2.close();
}
