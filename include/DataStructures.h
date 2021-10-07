#ifndef RECOMORE_DATASTRUCTURES_H
#define RECOMORE_DATASTRUCTURES_H

#include <utility>
#include <vector>

struct Waveform{
	unsigned short channel_;
	unsigned int event_;
	std::vector<float> waveform_{};
};


class WCData {
public:
	void addRow(const Waveform&);
private:
	std::vector<Waveform> waveforms_{};
};


typedef struct {
	float amplitude;
	float amplitudeError;
	float time;
	float timeError;

	// For debugging purpose, initial estimates of parameters.
	float foundAmplitude;
	float foundTime;
} PEData;

typedef struct {
	unsigned short ch; // channel_id
	unsigned short id; // sipm_id
	float chi2ndf;
	float baseline;
	std::vector<PEData> pes;
} waveformData;

typedef struct {
	float amplitude;
	float time;
} PMTData;

typedef struct {
	unsigned int eventID;
	double time;
	double deltaTime;
	std::vector<waveformData> sipm;
	PMTData pmt;
	PMTData topVeto;
	PMTData bottomVeto;
} eventDataV2;

#endif //RECOMORE_DATASTRUCTURES_H
