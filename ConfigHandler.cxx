#include "ConfigHandler.h"

/**
 * Constructor that reads a config file and sets all member variables according to the setted ones
 * @configfile: name of the config file
 * @ifile: input stream for the config file
 * @line: string that gets the content of a line in the configfile
 */
ConfigHandler::ConfigHandler(char* configfile){
	
	cout << "reading config file ...";
	ifstream ifile;
	ifile.open(configfile);
	string line;
	
	if(ifile.is_open())
	{
		//linewise file reading
		while(getline(ifile, line))
		{
			//checking for signal words and calling the respective function
			if(line.substr(0, 4) == "path")
				read_path(line);
		
			if(line.substr(0, 8) == "varsfile")
				read_varsfile(line);
				
			if(line.substr(0, 7) == "numruns")
				read_num_runs(line);
				
			if(line.substr(0, 8) == "varnames")
				read_var_names(line);
			
			if(line.substr(0, 12) == "thresholdfit")
				read_thresholdfit(line);
				
			if(line.substr(0, 13) == "thresholddata")
				read_thresholddata(line);
				
			if(line.substr(0, 12) == "thresholderr")
				read_thresholderr(line);
				
			if(line.substr(0, 5) == "ddist")
				read_ddist(line);
				
			if(line.substr(0, 8) == "chi2mean")
				read_chi2mean(line);
				
			if(line.substr(0, 9) == "chi2limit")
				read_chi2limit(line);
				
			if(line.substr(0, 11) == "dsmoothmean")
				read_dsmoothmean(line);
				
			if(line.substr(0, 10) == "smoothness")
				read_smoothness(line);
				
			if(line.substr(0, 5) == "kappa")
				read_kappa(line);
				
			if(line.substr(0, 8) == "exponent")
				read_exponent(line);
				
			if(line.substr(0, 7) == "summary")
				read_summaryflag(line);
				
			if(line.substr(0, 6) == "outdot")
				read_outdotflag(line);
				
			if(line.substr(0, 7) == "rescale")
				read_rescaleflag(line);
				
			if(line.substr(0, 8) == "runcombs")
				read_runcombs(line);
				
			if(line.substr(0, 8) == "leaveout")
				read_leave_out(line);
				
			if(line.substr(0, 4) == "load")
				read_load(line);
				
			if(line.substr(0, 8) == "filename")
				read_filename(line);
				
			if(line.substr(0, 7) == "rngseed")
				read_rngseed(line);
				
			if(line.substr(0, 6) == "covmat")
				read_covmat(line);
				
			if(line.substr(0, 7) == "outpath")
				read_outpath(line);
		}
	
		//creating the folder names of the runs based of @_num_runs
		set_run_names();
		
		//setting all analyses and observables that are in the first run
		set_analyses();
		
	}

	cout << "complete" << endl;
	
}

/**
 * This function reads the path of the data
 * @line: contains the path
 */
void ConfigHandler::read_path(string line){

	_path = line.substr(5);
	
}

/**
 * This function reads the name of the file that contains the anchor points
 * @line: contains the file name
 */
void ConfigHandler::read_varsfile(string line){
	
	_paramsfile = line.substr(9);
	
}

/**
 * This function reads the total amount of runs
 * @line: contains the total amount of runs
 */
void ConfigHandler::read_num_runs(string line){

	_num_runs = atoi(line.substr(8).c_str());
	
}

/**
 * This function reads the names of the variables
 * @line: contains the names of the variables
 */
void ConfigHandler::read_var_names(string line){

	_pos_begin = 8;
	
	//stepping through the whole line
	while(true)
	{
		//searching for spaces as indicator of the end of a variable name
		_pos_end = line.find(" ", _pos_begin + 1);
		
		//if a space is located at a position bigger than the size of the line, the line is finished
		//the last name is set and the loop ends
		if(_pos_end >= line.size())
		{
			_var_names.push_back(line.substr(_pos_begin + 1));
			break;
		}
		
		//extracting the name and setting the found space as the begin for the next search
		_var_names.push_back(line.substr(_pos_begin + 1, _pos_end - _pos_begin - 1));
		_pos_begin = _pos_end;
	}
	
	//setting the dimension
	_dimension = _var_names.size();
}

/**
 * This function sets the names of the folders
 */
void ConfigHandler::set_run_names(){

	//if no runs were performed/set, the function will call an error and return
	if(_num_runs == 0)
	{
		cout << "error: no runs set" << endl;
		return;
	}
		
	//writing the names
	for(size_t i = 0; i < _num_runs; i++)
	{	
		//setting a pure enumeration of the run names
		_run.push_back(to_string(i));
		
		//adding 0's, so that it fits the folder names
		if(i < 1000)
			_run[i].insert(_run[i].begin(), '0');
		
		if(i < 100)
			_run[i].insert(_run[i].begin(), '0');
			
		if(i < 10)
			_run[i].insert(_run[i].begin(), '0');
	}
	
	//if runcombs should be constructed ...
	if(_runcombs)
	{
		//init the rng
		srand(_rngseed);
		
		//dice the runs that will be left out and delete them from the @_run list
		for(size_t i = 0; i < _rcleave_out; i++)
		{
			_run.erase(_run.begin() + (rand() % _run.size()));
		}
		
		//correct _num_runs to the new number of used runs
		_num_runs -= _rcleave_out;
	}
			
	
}
//source: http://www.linuxquestions.org/questions/programming-9/c-list-files-in-directory-379323/
/**
 * This function reads the analyses that were performed
 * @dirbasis: base directory
 * @dir: directory in @dirbasis that contains the observables of an anlysis
 * @dirreader: reader for the observables
 * @diff: helper for making @_analysis and @_observable of equal length
 */
void ConfigHandler::set_analyses(){
		
	//if no path is given, an error is called and the function returns
	if(_path.empty())
	{
		cout << "error: no path set" << endl;
		return;
	}
	
	DIR* dirbasis;
	DIR* dir;
	struct dirent* dirreader;
	size_t diff;
	
	//moving into the directory that contains all the analysis folders
	if((dirbasis = opendir((_path + "/0000/plots/").c_str())) != NULL)
	
		//reading all names that are in the folder
		while((dirreader = readdir(dirbasis)) != NULL)
		{
			//skipping things that aren't analysis folders
			if(string(dirreader->d_name) == "index.html" || string(dirreader->d_name) == "." || string(dirreader->d_name) == "..")
				continue;
			
			//adding the analysis name
			_analysis.push_back(string(dirreader->d_name));
			
			//moving into the analysis folder
			if((dir = opendir((_path + "/0000/plots/" + string(dirreader->d_name)).c_str())) != NULL)
				
				//reading every observable of an analysis
				while((dirreader = readdir(dir)) != NULL)
				{
					//skipping things that aren't observables
					if(string(dirreader->d_name) == "index.html" || string(dirreader->d_name) == "." || string(dirreader->d_name) == "..")
						continue;
				
					//adding observable names
					_observable.push_back(string(dirreader->d_name));
		
				}
			else
				//cout an error if the observable couldn't be read
				cout << "error: could not open " << (_path + "/0000/plots/" + string(dirreader->d_name)).c_str() << endl;
		
			closedir(dir);	
		
			//getting the difference between the number of observables in the analysis & filling @_analysis with the same names
			//so that in both lists are the same amount of entries
			diff = _observable.size() - _analysis.size();
			for(size_t i = 0; i < diff; i++)
				_analysis.push_back(_analysis.back());
			
		
		}
	else
		//cout an error if the analysis couldn't be read
		cout << "error: could not open " << (_path + "/0000/plots/").c_str() << endl;
	
	closedir(dirbasis);
	
}

/**
 * This function reads the threshold for the RR constrain of the fit
 * @line: contains the threshold value
 */
void ConfigHandler::read_thresholdfit(string line){

		_thresholdfit = atof(line.substr(13).c_str());
	
}

/**
 * This function reads the threshold for the RR constrain of the hypercube fitting
 * @line: contains the threshold value
 */
void ConfigHandler::read_thresholddata(string line){

		_thresholddata = atof(line.substr(14).c_str());
	
}

/**
 * This function reads the threshold for the RR constrain of the error of the fit
 * @line: contains the threshold value
 */
void ConfigHandler::read_thresholderr(string line){

		_thresholderr = atof(line.substr(13).c_str());
	
}

/**
 * This function reads the convergence criteria for Ddist
 * @line: contains the value for the convergence
 */
void ConfigHandler::read_ddist(string line){

		_ddist = atof(line.substr(6).c_str());
	
}

/**
 * This function reads the number of Chi2 values for the mean
 * @line: contains the number of points
 */
void ConfigHandler::read_chi2mean(string line){

		_chi2mean = atof(line.substr(9).c_str());
	
}

/**
 * This function reads the number of Dsmooth values for the mean
 * @line: contains the number of points
 */
void ConfigHandler::read_dsmoothmean(string line){

		_dsmoothmean = atof(line.substr(12).c_str());
	
}

/**
 * This function reads the kappa for the RR constrain
 * @line: contains the value
 */
void ConfigHandler::read_kappa(string line){

		_kappa = atof(line.substr(6).c_str());
	
}

/**
 * This function reads the exponent for the distance weighting of the hypercube fitting
 * @line: contains the exponent
 */
void ConfigHandler::read_exponent(string line){

		_exponent= atof(line.substr(9).c_str());
	
}

/**
 * This function reads the smoothness convergence
 * @line: contains the value
 */
void ConfigHandler::read_smoothness(string line){

		_smoothness = atof(line.substr(11).c_str());
	
}

/**
 * This function reads the upper Chi2 limit
 * @line: contains the value
 */
void ConfigHandler::read_chi2limit(string line){
	
	_chi2limit = atof(line.substr(10).c_str());
	
}

/**
 * This function reads the flag for writing the summary output
 * @line: contains the value
 */
void ConfigHandler::read_summaryflag(string line){
	
	_summaryflag = atoi(line.substr(8).c_str());
	
}

/**
 * This function reads the flag for writing the dot product summary
 * @line: contains the value
 */
void ConfigHandler::read_outdotflag(string line){
	
	_outdotflag = atoi(line.substr(7).c_str());
	
}

/**
 * This function reads the rescaleflag
 * @line: contains the value
 */
void ConfigHandler::read_rescaleflag(string line){
	
	_rescaleflag = atoi(line.substr(8).c_str());
	
}

/**
 * This function reads the runcombinations flag
 * @line: contains the value
 */
void ConfigHandler::read_runcombs(string line){
	
	_runcombs = atoi(line.substr(9).c_str());
	
}

/**
 * This function reads how many runs should be left out for a runcombination
 * @line: contains the value
 */
void ConfigHandler::read_leave_out(string line){
	
	_rcleave_out = atoi(line.substr(9).c_str());
	
}

void ConfigHandler::read_load(string line){
	
	_load = atoi(line.substr(5).c_str());
	
}

/**
 * This function reads the name of the interpolationfile to which the result will be written or from which will be read
 * @line: contains the value
 */
void ConfigHandler::read_filename(string line){
	
	_filename = line.substr(9);
	
}

/**
 * This function reads the rng seed, that will be used
 * @line: contains the value
 */
void ConfigHandler::read_rngseed(string line){
	
	_rngseed = atoi(line.substr(8).c_str());
	
}

void ConfigHandler::read_covmat(string line){
	
	_covmat = atoi(line.substr(7).c_str());
	
}

void ConfigHandler::read_outpath(string line){
	
	_outpath = line.substr(8);
	
}
