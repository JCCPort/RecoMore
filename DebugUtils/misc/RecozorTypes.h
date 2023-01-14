#ifndef RECOMORE_RECOZORTYPES_H
#define RECOMORE_RECOZORTYPES_H
#include <vector>


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

#ifdef __CLING__
#pragma link C++ class pe_data;
//#pragma link C++ class std::vector<pe_data>;
#pragma link C++ class sipm_data;
//#pragma link C++ class pmt_data;
#pragma link C++ class std::vector<sipm_data>;
//#pragma link C++ class event_data;
#pragma link C++ class event_data_v2;
#endif

#endif //RECOMORE_RECOZORTYPES_H
