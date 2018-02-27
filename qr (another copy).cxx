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
//~ #include <ctime>


using namespace std;

// example: http://www.mia.uni-saarland.de/Teaching/NAVC-SS11/sol_c8.pdf
// page 5
//http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=5342430&tag=1

vector<vector<size_t>> RefHandler::_hypercubes;
vector<vector<vector<double>>> FitHandler::_structure;

int main(int argc, char** argv) {
	
	omp_set_dynamic(0);
	omp_set_nested(0);
	omp_set_num_threads(8);
	
	if(argc > 1)
	{
		
	ConfigHandler ch(argv[1]);
	
	if(argc > 2)
	{
		WeightsHandler wh(argv[2]);
		wh.sort(ch);
	}
		
	RefHandler rh(ch);
	size_t num_analysis = rh.getnum_analysis();
	
	cout << "setting up summary variables ..." << flush;
	
	vector<vector<double> > resultval, resulterr;
	vector<int> resultorder;
	resultval.resize(num_analysis);
	resulterr.resize(num_analysis);
	resultorder.resize(num_analysis);
	
	//~ vector<double> resultchi2, resultchi2red, resultiteration, resultdsmooth;
	
	//~ if(ch._summaryflag)
	//~ {
		//~ resultchi2.resize(num_analysis);
		//~ resultchi2red.resize(num_analysis);
		//~ resultiteration.resize(num_analysis);
		//~ resultdsmooth.resize(num_analysis);
	//~ }
	
	
	//~ vector<double> errcenter;
	//~ errcenter.resize(num_analysis);
	


	vector<double> maxparval, minparval;
	
	maxparval.assign(rh.getparameter_values()[0].size(), -numeric_limits<double>::infinity());
	minparval.assign(rh.getparameter_values()[0].size(), numeric_limits<double>::infinity());
	
	vector<double> center;
	center.resize(rh.getparameter_values()[0].size());
	

	
	//~ vector<size_t> unsmooth;
	//~ vector<double> meanerrs;
	//~ meanerrs.resize(num_analysis);
	
	vector<double> distances;
	
	
	if(ch._outdotflag)
		distances.resize(ch._num_runs * ch._num_runs - ch._num_runs);

	
	Power pow(rh.getparameter_values()[0].size());
	
	#pragma omp parallel
	{
		#pragma omp for
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
			
		if(ch._rescaleflag)
		#pragma omp single
		{
			rh.rescale(minparval, maxparval);
		}
	
		if(ch._outdotflag)
		{
			#pragma omp for schedule(dynamic)
			for(size_t i = 0; i < ch._num_runs; i++)
			{		
				size_t skips = 0;
				for(size_t j = 0; j < ch._num_runs; j++)
					if(i == j)
						skips--;				
					else
						distances[i * ch._num_runs - i - 1 + j + skips] = LinAlg::getdistanceofvectors(rh.getparameter_values()[i], rh.getparameter_values()[j]);
			}
		}
	
		#pragma omp for
		for(size_t i = 0; i < maxparval.size(); i++)
				center[i] = (maxparval[i] + minparval[i]) / 2.;
				
		#pragma omp single nowait
		{
			cout << "complete" << endl << "start fitting ..." << endl;
		}
			
		if(ch._load)
		{
			#pragma omp single
			{
				
				size_t num_ipol = 0;
	
			

				for(size_t i = 0; i < ch._analysis.size(); i++)
				{
					for(size_t j = 0; j < rh.getnum_bins()[i]; j++)
					{
						
				//~ for(size_t i = 0; i < 1; i++)
				//~ {
					//~ for(size_t j = 0; j < 1; j++)
					//~ {
						//~ if(num_ipol < 559)
						//~ {
							//~ num_ipol++;
							//~ continue;
						//~ }
						FitHandler fh;
						fh.load_fit(ch, i, j, num_ipol, pow, rh);
						
						rh.clearnormalvectors();
						rh.calculatenormalvectors(num_ipol, ch._thresholddata, ch._exponent, ch._kappa);
					
						
						cout << endl << "Result for bin " << num_ipol << " in ";
								
						size_t store = num_ipol;
						for(size_t k = 0; k < rh.getnum_bins().size(); k++)
						{
							store -= rh.getnum_bins()[k];
							if(store > num_ipol)
							{
								cout << ch._analysis[k] << "/" << ch._observable[k].substr(0, ch._observable[k].find(".")) << ":" << endl;
								break;
							}
						}
							

							cout << "TID:\t\t\t" << omp_get_thread_num() << endl;
							cout << "Dsmooth:\t\t" << fh.getDsmooth(rh) << endl;
							cout << "chi2:\t\t\t" << fh.getchi2() << endl;
							cout << "iterationcounter:\t" << fh.getiterationcounter() << endl;
							cout << "rel. error at center:\t" << fh.geterrorat(center) << endl;
							cout << "mean rel. error:\t" << fh.getmeanerror() << endl;
							cout << "-------------------------" << endl;
						
							
						if(ch._outdotflag)
						{
							ofstream outdot;
							outdot.open(("dotproduct" + to_string(num_ipol)).c_str());
							vector<double> rhdot = rh.getallgraddotproducts(), fhdot = fh.getallgraddotproducts(rh, ch._num_runs);
							outdot << fh.geterrorat(center) << "\t" << fh.getmeanerror() << "\t" << 0.0 << "\t" << fh.getDsmooth(rh) << "\t" << fh.getchi2() << "\t" << fh.getchi2red() << "\t" << fh.getiterationcounter() << endl;

							for(size_t k = 0; k < rhdot.size(); k++)
								outdot << rhdot[k] << "\t" << fhdot[k] << "\t" << distances[k] << "\t";
							outdot << endl;
						
							outdot.close();
								
							rhdot.clear();
							fhdot.clear();
						}

						/**
						 * Wenn das unten klappt, dann muesste das hier editiert werden
						 */
						if(ch._summaryflag)
						{
							ofstream outsummary;
							outsummary.open("summary", ofstream::out | ofstream::app);
							outsummary << fh.getchi2() << "\t" << fh.getchi2red() << "\t" << fh.getiterationcounter() << "\t" << fh.getDsmooth(rh) << endl;
							outsummary.close();
							//~ resultchi2[num_ipol] = fh.getchi2();
							//~ resultchi2red[num_ipol] = fh.getchi2red();
							//~ resultiteration[num_ipol] = fh.getiterationcounter();
							//~ resultdsmooth[num_ipol] = fh.getDsmooth(rh);
						}
							
						num_ipol++;
					}
				}
			}
			
		}
		else
		{
			#pragma omp for firstprivate(rh) schedule(dynamic)
			//~ #pragma omp single
			for(size_t num_ipol = 0; num_ipol < num_analysis; num_ipol++)
			//~ for(size_t num_ipol = 0; num_ipol < 8; num_ipol++)
			{
				FitHandler fh;
				
				int quitflagchi2, quitflagdsmooth;
				vector<double> chi2values, dsmoothvalues;
				double sum;
	
			
					
				fh.nextbin(rh, num_ipol, ch._thresholdfit, ch._kappa, pow);
				
				rh.clearnormalvectors();

				chi2values.clear();
				dsmoothvalues.clear();
			
				rh.calculatenormalvectors(num_ipol, ch._thresholddata, ch._exponent, ch._kappa);
			

				size_t bestiteration = fh.getiterationcounter();
				double bestchi2 = fh.getchi2();

			
				//convergence check
				chi2values.push_back(fh.getchi2());
				sum = 0.;
				quitflagchi2 = 1;
				quitflagdsmooth = 1;
				
			
		/////////////////////////////////////////////////////////////////////////////////////////////////////
			
				//~ while(fh.getDdist() > ch._ddist || quitflagdsmooth || quitflagchi2)
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

				//~ meanerrs[num_ipol] = fh.getmeanerror();
			
				
				fh.setfiterrors(rh, ch._thresholderr, ch._kappa);
				
				
				
				#pragma omp critical (cout)
				{
					cout << endl << "Result for bin " << num_ipol << " in ";
					
					size_t store = num_ipol;
					for(size_t i = 0; i < rh.getnum_bins().size(); i++)
					{
						store -= rh.getnum_bins()[i];
						if(store > num_ipol)
						{
							cout << ch._analysis[i] << "/" << ch._observable[i].substr(0, ch._observable[i].find(".")) << ":" << endl;
							break;
						}
					}
					
					cout << "RR constraint:\t";
					for(size_t i = 0; i < fh.getnumfitparams(); i++)
						if(!isnan(fh.geta()[i]))
							cout << i << "\t";
					cout << endl;
					cout << "TID:\t\t\t" << omp_get_thread_num() << endl;
					//~ cout << "Ddist:\t\t\t" << fh.getDdist() << endl;
					cout << "Dsmooth:\t\t" << fh.getDsmooth(rh) << endl;
					cout << "chi2:\t\t\t" << fh.getchi2() << endl;
					cout << "iterationcounter:\t" << fh.getiterationcounter() << endl;
					cout << "rel. error at center:\t" << fh.geterrorat(center) << endl;
					cout << "mean rel. error:\t" << fh.getmeanerror() << endl;
					cout << "-------------------------" << endl;
				}
				
				if(ch._outdotflag)
				{
					ofstream outdot;
					outdot.open(("dotproduct" + to_string(num_ipol)).c_str());
					vector<double> rhdot = rh.getallgraddotproducts(), fhdot = fh.getallgraddotproducts(rh, ch._num_runs);
					//~ outdot << fh.geterrorat(center) << "\t" << fh.getmeanerror() << "\t" << fh.getDdist() << "\t" << fh.getDsmooth(rh) << "\t" << fh.getchi2() << "\t" << fh.getchi2red() << "\t" << fh.getiterationcounter() << endl;
					outdot << fh.geterrorat(center) << "\t" << fh.getmeanerror() << "\t" << 0.0 << "\t" << fh.getDsmooth(rh) << "\t" << fh.getchi2() << "\t" << fh.getchi2red() << "\t" << fh.getiterationcounter() << endl;

					for(size_t i = 0; i < rhdot.size(); i++)
						outdot << rhdot[i] << "\t" << fhdot[i] << "\t" << distances[i] << "\t";
					outdot << endl;
			
					outdot.close();
					
					rhdot.clear();
					fhdot.clear();
				}
				
				//~ errcenter[num_ipol] = fh.geterrorat(center);

				if(ch._summaryflag)
				{
					#pragma omp critical (summaryout)
					{
						ofstream outsummary;
						outsummary.open("summary", ofstream::out | ofstream::app);
						outsummary << fh.getchi2() << "\t" << fh.getchi2red() << "\t" << fh.getiterationcounter() << "\t" << fh.getDsmooth(rh) << endl;
						outsummary.close();
						//~ resultchi2[num_ipol] = fh.getchi2();
						//~ resultchi2red[num_ipol] = fh.getchi2red();
						//~ resultiteration[num_ipol] = fh.getiterationcounter();
						//~ resultdsmooth[num_ipol] = fh.getDsmooth(rh);
					}
				}
				
				
				resultorder[num_ipol] = fh.getmax();
				
				fh.sortfitparams();
				resultval[num_ipol] = fh.getfitparams();
				resulterr[num_ipol] = fh.getfiterrors();

				//~ rh.rerescale(minparval, maxparval);
			}
		}
	}
	cout << "fitting complete" << endl << "writing output ..." << endl;
	
	ofstream output;//, outputchi2;
	
	//~ if(ch._summaryflag)
		//~ outputchi2.open("summary");
	
	if(!ch._load)
	{
		output.open(ch._filename);
		output << "DataFormat: binned 1" << endl;
		output << "ParamNames: ";
		for(size_t i = 0; i < ch._var_names.size(); i++)
			output << ch._var_names[i] << " ";
		output << endl;
		
		output << "Dimension: " << ch._var_names.size() << endl;
		
		
		output << "MinParamVals: ";
		for(size_t i = 0; i < minparval.size(); i++)
			output << minparval[i] << " ";
		output << endl << "MaxParamVals: ";
		for(size_t i = 0; i < maxparval.size(); i++)
			output << maxparval[i] << " ";
			
		if(ch._rescaleflag)
			output << endl << "DoParamScaling: 1";
		else
			output << endl << "DoParamScaling: 0";
		output << endl << "NumInputs: " << ch._num_runs << endl;
			
		maxparval.clear();
		minparval.clear();
		
		output << "Runs: ";
		for(size_t i = 0; i < ch._run.size(); i++)
			output << ch._run[i] << " ";
		output << endl << "---" << endl;
	}

	size_t offset = 0;
	
	for(size_t i = 0; i < ch._analysis.size(); i++)
	{
		for(size_t j = 0; j < rh.getnum_bins()[i]; j++)
		{
			//~ if(j + offset >= num_analysis)
				//~ goto endloop;
				
			if(!ch._load)
			{	
				output << "/" << ch._analysis[i] << "/" << ch._observable[i].substr(0, ch._observable[i].find(".")) << "#" << j << " " << rh.getintervall_start()[i][j] 
					<< " " << rh.getintervall_end()[i][j] << endl;
		
				output << "  val: " << ch._var_names.size() << " " << resultorder[j + offset] << " ";
				for(size_t k = 0; k < resultval[j + offset].size(); k++)
					output << resultval[j + offset][k] << " ";
				output << endl;


				output << "  err: " << ch._var_names.size() << " " << resultorder[j + offset] << " ";
				for(size_t k = 0; k < resulterr[j + offset].size(); k++)
					output << resulterr[j + offset][k] << " ";
				output << endl;
			}
			//~ if(ch._summaryflag)
				//~ outputchi2 << resultchi2[j + offset] << "\t" << resultchi2red[j + offset] << "\t" << resultiteration[j + offset] << "\t" << resultdsmooth[j + offset] << endl;
		}
		offset += rh.getnum_bins()[i];
	}

	//~ endloop:
	
		//~ if(ch._summaryflag)
			//~ outputchi2.close();
			
		if(!ch._load)
			output.close();
		
		
		//~ offset = 0;
		//~ ofstream outweights;
		//~ outweights.open("weights_errcenter");
		//~ for(size_t j = 0; j < ch._analysis.size(); j++)
		//~ {
			//~ for(size_t i = 0; i < rh.getnum_bins()[j]; i++)
				//~ outweights << "/" << ch._analysis[j] << "/" << ch._observable[j] << "\t" << fabs(errcenter[offset + i]) << endl;
			//~ offset += rh.getnum_bins()[j];
		//~ }
		//~ outweights.close();

	
}
else
	cout << "error: no config file given" << endl;
}
