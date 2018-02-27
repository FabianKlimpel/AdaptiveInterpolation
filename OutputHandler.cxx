#include "OutputHandler.h"

using namespace std;

/**
 * Constructor
 * Here, the sizes of the vectors that will later contain fitparameters etc. will be set. Additionally, the distances and centers will be calculated
 * @ch: carries the flag for the dotproduct-summary
 * @rh: handler of the reference data
 * @minparval: minimum of the parametervalues
 * @maxparval: maximum of the parametervalues
 */
OutputHandler::OutputHandler(ConfigHandler& ch, RefHandler& rh, vector<double>& minparval, vector<double>& maxparval)
{
	//resize the container of the fit results
	//the order is set to -1, so that an identification of a fit, that wasn't done yet can be performed
	resultval.resize(rh.getnum_analysis());
	resulterr.resize(rh.getnum_analysis());
	resultorder.assign(rh.getnum_analysis(), -1);
	
	//if the dotproduct-summary should be written, the distances will be calculated
	if(ch._outdotflag)
	{
		//Resizing the @distances vectors. No distances between an anchor point and itself will be calculated, therefore the size is adapted accordingly.
		distances.resize(ch._num_runs * ch._num_runs - ch._num_runs);
		//Calculate the distances
		setDistances(ch, rh);
	}
	
	//resize and calculate the centers
	center.resize(rh.getparameter_values()[0].size());	
	setCenter(minparval, maxparval);
}	

/**
 * This function calculates the distances between the anchors points
 * @ch: carries the number of MC runs
 * @rh: handles the reference data
 * @skips: no distances between an anchor point and itself is calculated, so the index shift needs to be registered by this variable
 */
void OutputHandler::setDistances(ConfigHandler& ch, RefHandler& rh)
{
	//walk over every point in parallel
	#pragma omp parallel for schedule(dynamic)
	for(size_t i = 0; i < ch._num_runs; i++)
	{		
		size_t skips = 0;
		
		//compare every point with every other
		for(size_t j = 0; j < ch._num_runs; j++)
			//if both points are the same, the iteration is skipped
			if(i == j)
				skips--;				
			else
				//calculate the distance between the vectors and store them
				distances[i * ch._num_runs - i - 1 + j + skips] = LinAlg::getdistanceofvectors(rh.getparameter_values()[i], rh.getparameter_values()[j]);
	}
}

/**
 * This function calculates the center of the hypercube of the parametervalues
 * @minparval: minimum of the parametervalues
 * @maxparval: maximum of the parametervalues
 */
void OutputHandler::setCenter(vector<double>& minparval, vector<double>& maxparval)
{
	//walk over every dimension in parallel
	#pragma omp parallel for
	for(size_t i = 0; i < maxparval.size(); i++)
		//calculate the mean in every dimension
		center[i] = (maxparval[i] + minparval[i]) / 2.;
}	

/**
 * This function sets up the interpolationfile and writes its header
 * @ch: container for information that will be written to the header
 * @minparval: minimum of the parametervalues
 * @maxparval: maximum of the parametervalues
 * @output: delivers the output to file functionality
 */
void OutputHandler::setup_outputfile(ConfigHandler& ch, vector<double>& minparval, vector<double>& maxparval)
{
	//set up the outputstream
	ofstream output;
	output.open((ch._outpath + ch._filename).c_str());
	
	//if runcombinations were performed, the seed for the rng will be written in the first line
	// starting the line with "#" prevents reading that line in "Professor"
	if(ch._runcombs)
		output << "# " << ch._rngseed << endl;
		
	//write out that the histograms are binned
	output << "DataFormat: binned 1" << endl;
	
	//write the parameter names
	output << "ParamNames: ";
	for(size_t i = 0; i < ch._var_names.size(); i++)
		output << ch._var_names[i] << " ";
	output << endl;
			
	//write the dimension of the hypercube	
	output << "Dimension: " << ch._var_names.size() << endl;
				
	//write out the minimum and maximum parameter values		
	output << "MinParamVals: ";
	for(size_t i = 0; i < minparval.size(); i++)
		output << minparval[i] << " ";
	output << endl << "MaxParamVals: ";
	for(size_t i = 0; i < maxparval.size(); i++)
		output << maxparval[i] << " ";
		
	//write out if the parameters were mapped onto a [0,1] intervall or not			
	if(ch._rescaleflag)
		output << endl << "DoParamScaling: 1";
	else
		output << endl << "DoParamScaling: 0";
		
	//write the number of samples used for the interpolation
	output << endl << "NumInputs: " << ch._num_runs << endl;
				
	//write the number of the used runs; that is important for future loading of the interpolation
	output << "Runs: ";
	for(size_t i = 0; i < ch._run.size(); i++)
		output << ch._run[i] << " ";
		
	//"---" indicates the end of the header
	output << endl << "---" << endl;
	output.close();
}

/**
 * This function creates the summary file and writes its header
 * @outsummary: delivers the output to file functionality
 */
void OutputHandler::setup_summary(ConfigHandler& ch)
{	
		//set up the file and write the header
		ofstream outsummary;
		outsummary.open((ch._outpath + "summary").c_str());
		outsummary << "Chi^2" << "\t" << "Chi2^2,red" << "\t" << "Iterations" << "\t" << "Dsmooth" << endl;
		outsummary.close();
}

/**
 * This function write the result of the fit for a bin to the console
 * @num_ipol: # of the bin
 * @ch: container for the observable name and the analysis name
 * @rh: handler of the reference data
 * @fh: handler of the fit
 * @store: counter to identify the right analysis/observable
 */
void OutputHandler::write_binresult(size_t num_ipol, ConfigHandler& ch, RefHandler& rh, FitHandler& fh)
{
	cout << endl << "Result for bin " << num_ipol<< " in ";
					
	//set store to the bin number				
	size_t store = num_ipol;
	
	//walk over the analyses
	for(size_t k = 0; k < rh.getnum_bins().size(); k++)
	{
		//substract the number of bins in the histograms
		store -= rh.getnum_bins()[k];
		
		//store is an unsigned int, so if @rh.getnum_bins()[k] > @store, the value of @store will be very large and definetly bigger than @num_ipol
		//this shows that the interpolation is performed in the analysis number @k
		if(store > num_ipol)
		{
			//write out the analysis name and the observable name
			cout << ch._analysis[k] << "/" << ch._observable[k].substr(0, ch._observable[k].find(".")) << ":" << endl;
			break;
		}
	}

	//if the interpolation is loaded, no RR constraints are performed, else the constrained monomial numbers will be printed to the terminal
	if(!ch._load)
	{
		cout << "RR constraint:\t";
		for(size_t i = 0; i < fh.getnumfitparams(); i++)
			if(!std::isnan(fh.geta()[i]))
				cout << i << "\t";
		cout << endl;
	}
	
	//write further summary variables will be printed to the terminal
	cout << "TID:\t\t\t" << omp_get_thread_num() << endl;
	cout << "Dsmooth:\t\t" << fh.getDsmooth(rh) << endl;
	cout << "chi2:\t\t\t" << fh.getchi2() << endl;
	cout << "iterationcounter:\t" << fh.getiterationcounter() << endl;
	cout << "rel. error at center:\t" << fh.geterrorat(center) << endl;
	cout << "mean rel. error:\t" << fh.getmeanerror() << endl;
	cout << "-------------------------" << endl;	
}

/**
 * This function writes the dotproduct-summary
 * @num_ipol: # of the bin
 * @ch: contains the number of runs
 * @rh: handles the reference data
 * @fh: handles the fit
 * @outdot: delivers the output to file functionality
 * @rhdot, @fhdot: stores the dot products of the reference data and the fit at every anchor point
 */
void OutputHandler::write_dotproduct(size_t num_ipol, ConfigHandler& ch, RefHandler& rh, FitHandler& fh)
{
	//set up the output
	ofstream outdot;
	outdot.open((ch._outpath + "dotproduct" + to_string(num_ipol)).c_str());
	
	//calculate all dot products
	vector<double> rhdot = rh.getallgraddotproducts(), fhdot = fh.getallgraddotproducts(rh, ch._num_runs);
	
	//write some summary parameters
	outdot << fh.geterrorat(center) << "\t" << fh.getmeanerror() << "\t" << 0.0 << "\t" << fh.getDsmooth(rh) << "\t" << fh.getchi2() << "\t" 
		<< fh.getchi2red() << "\t" << fh.getiterationcounter() << endl;
		
	//write the dot products of the reference data, the fit and the distances between the anchor points used for the dot product
	for(size_t k = 0; k < rhdot.size(); k++)
		outdot << rhdot[k] << "\t" << fhdot[k] << "\t" << distances[k] << "\t";
	outdot << endl;
						
	outdot.close();
}
	
/**
 * This function writes a summary of a fit to the summary file
 * @fh: handles the fit
 * @rh: handles the reference data
 * @outsummary: delivers the output to file functionality
 */
void OutputHandler::write_summary(FitHandler& fh, RefHandler& rh, ConfigHandler& ch)
{
	//open the file and continue writing at its end
	ofstream outsummary;
	outsummary.open((ch._outpath + "summary").c_str(), ofstream::out | ofstream::app);
	//write out the resulting Chi^2, the reduced Chi^2, the number of iterations needed and the smoothness
	outsummary << fh.getchi2() << "\t" << fh.getchi2red() << "\t" << fh.getiterationcounter() << "\t" << fh.getDsmooth(rh) << endl;
	outsummary.close();
}

/**
 * This function updates the list of fits
 * @num_ipol: # of the bin
 * @fh: handles the fit
 */
void OutputHandler::update(size_t num_ipol, FitHandler& fh)
{
	//fills up the fitparameters with 0's until the polynomial order is full
	fh.sortfitparams();
	
	//save the order, the parameters and their errors to the storage variables
	resultorder[num_ipol] = fh.getmax();
	resultval[num_ipol] = fh.getfitparams();
	resulterr[num_ipol] = fh.getfiterrors();
}

/**
 * This function writes fits to the interpolationfile
 * @ch: handles the filename, the analysis name, the observable name and the variable names
 * @rh: handles the reference data
 * @output: delivers the output to file functionality
 */
void OutputHandler::write_ipolfile(ConfigHandler& ch, RefHandler& rh)
{
	//set up the output and continue writing at the end of the file
	ofstream output;
	output.open((ch._outpath + ch._filename).c_str(), ofstream::out | ofstream::app);
	
	/**
	 * These for-loops write fit results to the interpolationfile. All the results are stored in a list. While a block of sequential results are available, those will be written to the file.
	 */
	//walk over every analysis
	//@start_analysis serves as an offset to skip already written histograms
	for(size_t i = start_analysis; i < ch._analysis.size(); i++)
	{
		//walk over every bin
		//@start_bins serves as an offset to skip already written bins in an histogram
		for(size_t j = start_bins; j < rh.getnum_bins()[i]; j++)
		{
			//check if the fit at the current position was already performed
			if(!resultval[j + offset].empty() && !resulterr[j + offset].empty() && resultorder[j + offset] != -1)
			{
				//write out the fit result
				cout << "writing: " << ch._analysis[i] << "/" << ch._observable[i].substr(0, ch._observable[i].find(".")) << "#" << j << " from thread #" 
						<< omp_get_thread_num() << endl;
								
								
				output << "/" << ch._analysis[i] << "/" << ch._observable[i].substr(0, ch._observable[i].find(".")) << "#" << j << " " << rh.getintervall_start()[i][j] 
						<< " " << rh.getintervall_end()[i][j] << endl;
							
				output << "  val: " << ch._var_names.size() << " " << resultorder[j + offset] << " ";
				for(size_t k = 0; k < resultval[j + offset].size(); k++)
					output << resultval[j + offset][k] << " ";
				output << endl;

				output << "  err: " << ch._var_names.size() << " " << resultorder[j + offset] << " ";
				for(size_t k = 0; k < resulterr[j + offset].size(); k++)
					if(resulterr[j + offset][k] == std::numeric_limits<double>::infinity())
						output << std::numeric_limits<double>::max() << " ";
					else
						output << resulterr[j + offset][k] << " ";
				output << endl;
									
				//if the fit is written, the entries in the summary list is cleared
				//the order is resetted to -1 as additional indicator that there is no fit stored anymore
				resultval[j + offset].clear();
				resulterr[j + offset].clear();
				resultorder[j + offset] = -1;
									
				//if there are are more bins in a histogram, the offset will be increased for future for loops
				if(start_bins < rh.getnum_bins()[i] - 1)
					start_bins++;
				else
				{
					//if the histogram is fully written to file, the bin # is set to 0 and the next histogram will be written in future for loops
					start_bins = 0;
					start_analysis++;
				}
			}
			else
				//leave the loop if a fit isn't performed
				goto endloop;
						
		}
		//offset represents the number of bins of completely written histograms
		offset += rh.getnum_bins()[i];
	}
	endloop:
	output.close();
}

//~ void OutputHandler::write_covmat(MatrixXd& mat, size_t num_ipol, ConfigHandler& ch)
//~ {
	//~ ofstream outcovmat;
	//~ outcovmat.open((ch._outpath + "covmat_" + to_string(num_ipol)).c_str());
	//~ 
	//~ for(size_t row = 0; row < (size_t) mat.rows(); row++)
	//~ {
		//~ for(size_t col = 0; col < (size_t) mat.cols(); col++)
			//~ if(mat(row, col) == std::numeric_limits<double>::infinity())
				//~ outcovmat << std::numeric_limits<double>::max() << " ";
			//~ else
				//~ if(-mat(row, col) == std::numeric_limits<double>::infinity())
					//~ outcovmat << -std::numeric_limits<double>::max() << " ";
				//~ else
					//~ outcovmat << mat(row, col) << " ";
		//~ outcovmat << "\n";
	//~ }
	//~ 
	//~ outcovmat.close();
	//~ 
//~ }

void OutputHandler::write_covmat(TMatrixDSym& tmds, size_t num_ipol, ConfigHandler& ch)
{
	ofstream outcovmat;
	outcovmat.open((ch._outpath + "covmat_" + to_string(num_ipol)).c_str());
	
	for(size_t row = 0; row < (size_t) tmds.GetNrows(); row++)
	{
		for(size_t col = 0; col < (size_t) tmds.GetNcols(); col++)
			if(tmds[row][col] == std::numeric_limits<double>::infinity())
				outcovmat << std::numeric_limits<double>::max() << " ";
			else
				if(-tmds[row][col] == std::numeric_limits<double>::infinity())
					outcovmat << -std::numeric_limits<double>::max() << " ";
				else
					outcovmat << tmds[row][col] << " ";
		outcovmat << "\n";
	}
	
	outcovmat.close();
	
}
