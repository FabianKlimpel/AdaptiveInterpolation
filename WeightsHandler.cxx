#include "WeightsHandler.h"

using namespace std;

/**
 * Constructor, that reads a given weightfile
 * @weightsfile: name of the file, that contains the informations about the weights
 * @pos_begin, @pos_end: begin und end positions in a string for data extraction
 * @ifile: reader for the file
**/ 
WeightsHandler::WeightsHandler(char* weightsfile){

	size_t pos_begin, pos_end;
	string line;
	ifstream ifile(weightsfile);

	if (ifile.is_open()) 
	{
		cout << "reading weights file ...";

		//linewise stepping through the file
		while (getline(ifile, line)) 
		{
			if(line.empty()) continue;
			if(line[0] == '#') continue;
			
			//"/" is an indicator for content in the line
			if (line.substr(0, 1) == "/")
			{
				//extracting the analysisname
				pos_begin = 1;
				pos_end = line.find("/", pos_begin + 1);
		
				_analysis.push_back(line.substr(pos_begin, pos_end - pos_begin));

				//extracting the observablename
				pos_begin = pos_end + 1;
				pos_end = line.find("\t", pos_begin + 1);

				_observable.push_back(line.substr(pos_begin, pos_end - pos_begin));

				//extracting the actual weight of it
				pos_begin = pos_end;

				_weights.push_back(atoi(line.substr(pos_begin, line.size() - pos_begin).c_str()));
					
			}
		}
			

		ifile.close();
		cout << "complete" << endl;
	}
	else
		cout << "unable to read weights file" << endl;
}

/**
 * This function sorts the weights so, that it can be easier combined with the order in the ConfigHandler object.
 * Furthermore, the ConfigHandler object will be truncated if a weight is not mentioned in @_weights
 * @ch: this object contains the analyses and observables in the data folders
 * @tmp_weights: temporary storage of the weights
 * @to_delete: list that stores all observables that will be deleted
 * @i, @j: countingvariables
 * @notfound: flag for adding an index to @to_delete
**/ 
void WeightsHandler::sort(ConfigHandler& ch){

	cout << "sorting weights ...";
	vector<int> tmp_weights;
	vector<size_t> to_delete;
	size_t i, j;
	
	//flag is initialized as active to delete
	int notfound = 1;

	//walk over all bins (aka @ch._analysis.size())
	for(i = 0; i < ch._analysis.size(); i++)
	{
		//walk over all weights in this object
		for(j = 0; j < _analysis.size(); j++)	
			//If the observable & analysis of ch and this object matches, add it to @tmp_weights. 
			//Therefore the weights will be in a list in the same order as the bins mentioned in the interpolationfile
			if(ch._analysis[i] == _analysis[j] && ch._observable[i] == (_observable[j] + ".dat").c_str())
			{
				tmp_weights.push_back(_weights[j]);
				
				//deactivate the flag for this loop iteration
				notfound = 0;
				break;
			}

				
		//if the flag remains active (= an analysis in @ch is not in @_weights) after every weight is checked, the respective observable in @ch will be deleted later
		if(notfound)
			to_delete.push_back(i);
			
		//activate the flag for the next iteration
		notfound = 1;
	}
	
	//Walking from back to front over all observables that will be deleted. That way, the indices remain the same.
	for(i = to_delete.size(); i > 0; i--)
	{
		ch._analysis.erase(ch._analysis.begin() + to_delete[i - 1]);
		ch._observable.erase(ch._observable.begin() + to_delete[i - 1]);
	}
		
	_weights = tmp_weights;
	
	cout << "complete" << endl;			
}

/**
 * This function return @_weights
 */
vector<int> WeightsHandler::getweights(){
	
	return _weights;
	
}
