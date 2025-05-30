#include <string>
#include <vector>
#include <fstream>
#include "TTree.h"
#include "TFile.h"
#include "RecozorTypes.h"

// run using $> root 'ampVsTime.cc("[datafile name]")'


void GetRecozorList(std::string filename)
{

#ifdef CINT
#pragma link C++ class pe_data+;
#pragma link C++ class sipm_data+;
#pragma link C++ class std::vector<sipm_data>+;
#pragma link C++ class event_data_v2+;
#endif


//	gInterpreter->GenerateDictionary("std::vector<pe_data>", "vector");
//	gInterpreter->GenerateDictionary("std::vector<sipm_data>", "vector");
	
	TFile* file = TFile::Open(filename.c_str());
	TTree* tree = (TTree*)file->Get("reco_tree");
	
	unsigned int* eventID;
	double* time;
	double* deltaTime;
	std::vector<sipm_data>* sipmData;
//	pmt_data* pmtData;
//	pmt_data* topVeto;
//	pmt_data* bottomVeto;
	
	tree->SetBranchAddress("event_id", &eventID);
	tree->SetBranchAddress("time", &time);
	tree->SetBranchAddress("delta_time", &deltaTime);
	tree->SetBranchAddress("sipm", &sipmData);
//	tree->SetBranchAddress("pmt", &pmtData);
//	tree->SetBranchAddress("top_veto", &topVeto);
//	tree->SetBranchAddress("bottom_veto", &bottomVeto);
	
	std::vector<double> amps;
	std::vector<double> times;
	
	for (int iEntry = 0; tree->LoadTree(iEntry) >= 0; ++iEntry) {
		// Load the data for the given tree entry
		tree->GetEntry(iEntry);
		
		for(int i = 0; i<(*sipmData).size(); i++){
			for(int j = 0; j<((*sipmData)[i]).pes.size(); j++){
				amps.emplace_back(((*sipmData)[i]).pes[j].amplitude);
				times.emplace_back(((*sipmData)[i]).pes[j].time);
			}
		}
	}
	
	std::fstream ampFile("amp.dat", std::ios::out | std::ios::binary);
	ampFile.write((char*)&amps[0], amps.size() * sizeof(double));
	ampFile.close();
	
	std::fstream timeFile("time.dat", std::ios::out | std::ios::binary);
	timeFile.write((char*)&amps[0], amps.size() * sizeof(double));
	timeFile.close();
	
	file->Close();
}