#include <string>
#include <vector>

// run using $> root 'ampVsTime.cc("[datafile name]")'

typedef struct {
	float amplitude;
	float amplitude_error;
	// float npe;
	// float npe_error;
	float time;
	float time_error;
	
	// for debugging purpose
	float found_amplitude;
	float found_time;
} pe_data;

typedef struct {
	unsigned short ch; // channel_id
	unsigned short id; // sipm_id
	float chi2ndf;
	float baseline;
	std::vector<pe_data> pes;
} sipm_data ;

typedef struct {
	float amplitude;
	float time;
} pmt_data ;

typedef struct {
	unsigned int event_id;
	double time;
	double delta_time;
	std::vector<sipm_data> sipm;
	pmt_data pmt;
	pmt_data top_veto;
	pmt_data bottom_veto;
} event_data_v2 ;

void ampVsTime(std::string filename)
{
	TFile* file = TFile::Open(filename.c_str());
	TTree* tree = (TTree*)file->Get("reco_tree");
	
	unsigned int eventID;
	double time;
	double deltaTime;
	std::vector<sipm_data> sipmData;
	pmt_data pmtData;
	pmt_data topVeto;
	pmt_data bottomVeto;
	
	tree->SetBranchAddress("event_id", &eventID);
	tree->SetBranchAddress("time", &time);
	tree->SetBranchAddress("delta_time", &deltaTime);
	tree->SetBranchAddress("sipm", &sipmData);
	tree->SetBranchAddress("pmt", &pmtData);
	tree->SetBranchAddress("top_veto", &topVeto);
	tree->SetBranchAddress("bottom_veto", &bottomVeto);
	
	for (int iEntry = 0; tree->LoadTree(iEntry) >= 0; ++iEntry) {
		// Load the data for the given tree entry
		tree->GetEntry(iEntry);
		
		// Now, `variable` is set to the value of the branch
		// "branchName" in tree entry `iEntry`
		printf("%d\n", eventID);
	}
	
}