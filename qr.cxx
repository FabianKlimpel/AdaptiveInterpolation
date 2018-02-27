#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>
#include <limits>
#include "ConfigHandler.h"
#include <fstream>
#include "LinAlg.h"
#include "RefHandler.h"
#include "FitHandler.h"
#include <omp.h>
#include "WeightsHandler.h"
#include "Power.h"
#include "OutputHandler.h"


using namespace std;

// example: http://www.mia.uni-saarland.de/Teaching/NAVC-SS11/sol_c8.pdf
// page 5
//http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5342430&tag=1

//calling the program: ./qr configfile weightsfile


//the hypercubes of the data remains always the same once it is calculated as well as the power structure of the monomials, therefore those are static variables
vector<vector<size_t>> RefHandler::_hypercubes;
vector<vector<vector<double>>> FitHandler::_structure;

int main(int argc, char** argv) {
	
	
	//setting up the parallel environment
	omp_set_dynamic(0);
	omp_set_nested(1);
	omp_set_num_threads(8);
	
	//check if a config file is given
	if(argc > 1)
	{
	//read the config file
	ConfigHandler ch(argv[1]);
	
	if(argc > 2)
	{
		//read the weights file and sort it as readed in the config file
		WeightsHandler wh(argv[2]);
		wh.sort(ch);
	}
		
	//load the reference data
	RefHandler rh(ch);
	size_t num_analysis = rh.getnum_analysis();
	
	cout << "setting up summary variables ..." << flush;

	
	//calculating the maximum and minimum of every parameter value
	vector<double> maxparval, minparval;
	
	maxparval.assign(rh.getparameter_values()[0].size(), -numeric_limits<double>::infinity());
	minparval.assign(rh.getparameter_values()[0].size(), numeric_limits<double>::infinity());

	//walking over every parametervalue and compare it to the current minimum/maximum
	//if the value is smaller/bigger than the current, it will be set
	//the parallel environment just gives a small boost in the runtime due to the writes to the vector entry
	#pragma omp parallel for
	for(size_t i = 0; i < rh.getparameter_values().size(); i++)
		for(size_t j = 0; j < rh.getparameter_values()[0].size(); j++)
		{
			if(rh.getparameter_values()[i][j] < minparval[j])
				#pragma omp critical (min)
				minparval[j] = rh.getparameter_values()[i][j];
			if(rh.getparameter_values()[i][j] > maxparval[j])
				#pragma omp critical (max)
				maxparval[j] = rh.getparameter_values()[i][j];
		}	
	
	//@pow is created and serves as power calculator
	Power pow(rh.getparameter_values()[0].size());
		
	//rescaling the parameters if the respective flag is set	
	if(ch._rescaleflag)
		#pragma omp single
		rh.rescale(minparval, maxparval);
	
	//create @oh that handles the output
	OutputHandler oh(ch, rh, minparval, maxparval);

	//if the program should not load an interpolationfile (= calculate it here), the interpolation file will be set up
	if(!ch._load)
		oh.setup_outputfile(ch, minparval, maxparval);
		
	//if a summary should be written, the respective file will be created
	if(ch._summaryflag)
		oh.setup_summary(ch);
		
	cout << "complete" << endl << "start fitting ..." << endl;
		
		//if-condition for loading an interpolation or calculating it
		if(ch._load)
		{
	
			//loop over every histogram in a parallel way
			//the reference data contains every data point of every histogram and is the same for all threads,
			//but the normalvectors depend on the respective bin and they are stored in the object as well, so it must be firstprivate for all threads
			#pragma omp parallel for firstprivate(rh) shared(pow, ch, oh) schedule(dynamic)
			//~ #pragma omp single
				for(size_t i = 0; i < ch._analysis.size(); i++)
				{
					/**
					 * In @rh, the histograms and their bins are stored in a matrix, but the fits in @fh are stored in a single list bin by bin. The second dimension of the fits are the fitparameters.
					 * Therefore the number of loaded fits need to be tracked in order to get the mapping between the storage constructions. This is done by @num_ipol
					 */
					size_t num_ipol = 0;
					
					//Due to the parallel environment, skips for @num_ipol need be possible so that every thread handles the right data. This for-loop simply skips through to the right histogram.
					for(size_t j = 0; j < i; j++)
						num_ipol += rh.getnum_bins()[j];
					
					//walking over every bin in a histogram
					for(size_t j = 0; j < rh.getnum_bins()[i]; j++)
					{
						
							
						//create @fh and load the needed fit
						FitHandler fh;			
						fh.load_fit(ch, i, j, num_ipol, pow, rh, oh);
						
						//get the normalvectors of the bin
						rh.clearnormalvectors();
						rh.calculatenormalvectors(num_ipol, ch._thresholddata, ch._exponent, ch._kappa);
					
						//write the result
						#pragma omp critical (cout)
						oh.write_binresult(num_ipol, ch, rh, fh);
						
						//write out a dotproduct-summary if the flag is set; every bin got its own filename, so it is threadsafe
						if(ch._outdotflag)
							oh.write_dotproduct(num_ipol, ch, rh, fh);
						
						//write a sumary if the flag is set; every thread writes into the same file, so it needs a critical pragma
						if(ch._summaryflag)
							#pragma omp critical (summaryout)
							oh.write_summary(fh, rh, ch);
						
						//if a bin is done, increase @num_ipol and skip to the next bin
						num_ipol++;
					}
				}
			}
			
		else
		{
			
			//the else-statement is the calculation of a new interpolationfunction
			#pragma omp parallel for firstprivate(rh) shared(pow, ch, oh) schedule(dynamic)
			//~ #pragma omp single
			for(size_t num_ipol = 0; num_ipol < num_analysis; num_ipol++)
			//~ for(size_t num_ipol = 17; num_ipol < 18; num_ipol++)
			{
				cout << "start with bin " << num_ipol << " in thread #" << omp_get_thread_num() << endl;
				
				//setting up @fh and several summary variables that serve as convergence checks
				FitHandler fh;
				
				int quitflagchi2, quitflagdsmooth;
				vector<double> chi2values, dsmoothvalues;
				double sum;
	
			
				//set @fh up on the the right bin
				fh.nextbin(rh, num_ipol, ch._thresholdfit, ch._kappa, pow);
				
				//clear potentially stored normal vectors in @rh
				rh.clearnormalvectors();

				//clear the convergence criteria
				chi2values.clear();
				dsmoothvalues.clear();
			
				//calculate the new normal vectors of the data
				rh.calculatenormalvectors(num_ipol, ch._thresholddata, ch._exponent, ch._kappa);
			
				//after the first iteration, this is the currently best iteration
				//the number of the iteration will always be stored, so that it is possible to get to that later
				//the respective quality parameter is stored, too
				size_t bestiteration = fh.getiterationcounter();
				double bestchi2 = fh.getchi2();
			
				//convergence check setup
				chi2values.push_back(fh.getchi2());
				sum = 0.;
				quitflagchi2 = 1;
				quitflagdsmooth = 1;
				
		/////////////////////////////////////////////////////////////////////////////////////////////////////
				
				cout << "fit done once for bin " << num_ipol << " by thread #" << omp_get_thread_num() << endl;
			
			//if-condition for rescaling or not	
			if(ch._rescaleflag)
			{
				//set @quitflag to 1 is for entering the while-loop
				int quitflag = 1;
				//@iteration_results serves as comparison of an iteration to earlier ones, so that a better/worse statement can be made
				vector<double> iteration_results;
				
				//while a fit becomes better with more iterations, the following code will loop
				while(quitflag)
				{					
					//set up new monomials for the bin
					fh.nextstep(rh, ch._thresholdfit, ch._kappa, pow);
					
					//store a product of the current Chi^2 value and a mapped smoothness
					iteration_results.push_back(fh.getchi2() * (1 - fh.getDsmooth(rh)) / (1 + fh.getDsmooth(rh)));
					
					//if enough iterations were made, the convergence check will be performed
					if(iteration_results.size() >= ch._chi2mean)
					{
						//calculate the mean of @iteration_results
						for(size_t i = 0; i < iteration_results.size(); i++)
							sum += iteration_results[i];
					
						sum /= iteration_results.size();
			
						//if the current value is equal or worse the sliding mean, the loop ends
						if(iteration_results[iteration_results.size() - 1] >= sum)
							quitflag = 0;
			
						//delete the first entry in @iteration_results = keeping the list always at the same length
						iteration_results.erase(iteration_results.begin());
					}
			
					//reset @sum
					sum = 0.;
			
					//if an iteration is better than the best iteration performed, this iteration will be stored as the best iteration
					if(iteration_results.back() < bestchi2)
					{
						bestchi2 = iteration_results.back();
						bestiteration = fh.getiterationcounter();
					}
					
					if(fh.getiterationcounter() % 100 == 0)
						cout << "Bin: " << num_ipol << ", Iteration: " << fh.getiterationcounter() << ", Chi^2_{mod}: " << iteration_results.back() << endl;
					
				
				}	

				cout << "fit complete for bin " << num_ipol << " by thread #" << omp_get_thread_num() << endl;
				
				//recalculate the best iteration
				fh.setiteration(rh, num_ipol, ch._thresholdfit, ch._kappa, bestiteration, pow);			
								
				//calculate the errors of the fitparameters
				//~ fh.setfiterrors(rh, ch._thresholderr, ch._kappa);
				fh.setfiterrors2(ch, oh, num_ipol);
			
			}
			else
			{
				while(quitflagdsmooth || quitflagchi2)
				{

					fh.nextstep(rh, ch._thresholdfit, ch._kappa, pow);
								
					chi2values.push_back(fh.getchi2());
					dsmoothvalues.push_back(fh.getDsmooth(rh));
			
					quitflagchi2 = 1;
					quitflagdsmooth = 1;
			
					if(chi2values.size() >= ch._chi2mean)
					{
						for(size_t i = 0; i < chi2values.size(); i++)
							sum += chi2values[i];
					
						sum /= chi2values.size();
			
						if(chi2values[chi2values.size() - 1] >= sum)
							quitflagchi2 = 0;
			
						chi2values.erase(chi2values.begin());
					}
			
					sum = 0.;
			
					if(dsmoothvalues.size() >= ch._dsmoothmean)
					{
						for(size_t i = 0; i < dsmoothvalues.size(); i++)
							sum += dsmoothvalues[i];
					
						sum /= dsmoothvalues.size();
			
						if(dsmoothvalues[dsmoothvalues.size() - 1] <= sum)
							quitflagdsmooth = 0;
			
			
						dsmoothvalues.erase(dsmoothvalues.begin());
					}
			
					sum = 0.;
					
					if(fh.getchi2() < bestchi2)
					{
						bestchi2 = fh.getchi2();
						bestiteration = fh.getiterationcounter();
					}
					
				
				}
				
				rh.rescale(minparval, maxparval);
				
				fh.setiteration(rh, num_ipol, ch._thresholdfit, ch._kappa, bestiteration, pow);			
								
				//~ fh.setfiterrors(rh, ch._thresholderr, ch._kappa);
				fh.setfiterrors2(ch, oh, num_ipol);
				
				rh.rerescale(minparval, maxparval);
			}	
				
				cout << "prepare output for bin " << num_ipol << " by thread #" << omp_get_thread_num() << endl;
				
				//write result to the terminal; the critical pragma prevents formatting problems
				#pragma omp critical (cout)
				oh.write_binresult(num_ipol, ch, rh, fh);
					
				//write the dotproduct-summary, if the flag is set
				if(ch._outdotflag)
					oh.write_dotproduct(num_ipol, ch, rh, fh);
				
				//write the summary, if the flag is set; every thread writes to the file, so it needs a critical pragma
				if(ch._summaryflag)
					#pragma omp critical (summaryout)
					oh.write_summary(fh, rh, ch);
					
				//store the result for the interpolationfile
				oh.update(num_ipol, fh);

				//write the result to the interpolationfile
				#pragma omp critical (ipolresultoutput)
				oh.write_ipolfile(ch, rh);
			}
		}
	cout << "fitting complete" << endl;
}
else
	cout << "error: no config file given" << endl;
}
