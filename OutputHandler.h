#ifndef __OUTPUTHANDLER__H
#define __OUTPUTHANDLER__H

#include <iostream>
#include <vector>
#include "FitHandler.h"
#include "RefHandler.h"
#include "ConfigHandler.h"
#include <fstream>
#include <omp.h>
#include "LinAlg.h"
//#include <Eigen/Dense>
#include <TMatrixDSym.h>

using namespace std;
//~ using namespace Eigen;

class FitHandler;

/**
* This Class is a container for output functions
*/
class OutputHandler
{
public:

	//Constructor that sets up some necessary parameters
	OutputHandler(ConfigHandler& ch, RefHandler& rh, vector<double>& minparval, vector<double>& maxparval);

	//setup for the interpolationfile
	void setup_outputfile(ConfigHandler& ch, vector<double>& minparval, vector<double>& maxparval);

	//writes a summary of the fit to the terminal
	void write_binresult(size_t num_ipol, ConfigHandler& ch, RefHandler& rh, FitHandler& fh);

	//writes a dotproduct-summary
	void write_dotproduct(size_t num_ipol, ConfigHandler& ch, RefHandler& rh, FitHandler& fh);

	//writes to the summaryfile
	void write_summary(FitHandler& fh, RefHandler& rh, ConfigHandler& ch);

	//stores performed fits into lists
	void update(size_t num_ipol, FitHandler& fh);
	
	//writes performed fits to the interpolationfile
	void write_ipolfile(ConfigHandler& ch, RefHandler& rh);

	//setup for the summary file
	void setup_summary(ConfigHandler& ch);

	//~ void write_covmat(MatrixXd& mat, size_t num_ipol, ConfigHandler& ch);
	void write_covmat(TMatrixDSym& tmds, size_t num_ipol, ConfigHandler& ch);

private:

	/**
	 * @center: center of the hypercube of the parametervalues
	 * @distances: distances of every anchor point to another
	 * @start_bins, @start_analysis, @offset: variables that are used to identify what fit will be/was written to the interpolationfile
	 * @resultval, @resulterr, @resultorder: storage for performed fits
	 */
	vector<double> center, distances;

	size_t start_bins = 0, start_analysis = 0, offset = 0;
	vector<vector<double> > resultval, resulterr;
	vector<int> resultorder;

	//calculates the center of the hypercube of the parametervalues
	void setCenter(vector<double>& minparval, vector<double>& maxparval);
	
	//calculates the distances between the anchor points
	void setDistances(ConfigHandler& ch, RefHandler& rh);
};

#endif
