#ifndef RECOMORE_DATASTRUCTURES_H
#define RECOMORE_DATASTRUCTURES_H

#include <utility>
#include <vector>

struct Waveform {
	unsigned short channel_;
	unsigned int event_;
	std::vector<float> waveform_{};
};


class WCData {
public:
	void addRow(const Waveform &);

	std::vector<Waveform> getWaveforms() { return waveforms_; };
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
	unsigned short ch;
	float chi2ndf;
	float baseline;
	std::vector<PEData> pes;
} WaveformData;


typedef struct {
	unsigned int eventID;
	std::vector<WaveformData> sipm;
} EventData;

#endif //RECOMORE_DATASTRUCTURES_H
