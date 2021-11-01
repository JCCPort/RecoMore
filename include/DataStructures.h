#ifndef RECOMORE_DATASTRUCTURES_H
#define RECOMORE_DATASTRUCTURES_H

#include <utility>
#include <vector>
#include <string>

struct WaveformData {
	unsigned short channel;
	std::vector<float> waveform{};
};

typedef struct {
	unsigned int eventID;
	std::string TDCCorrTime;
	std::string date;
	std::vector<WaveformData> chData;
} EventData;


class WCData {
public:
	void addRow(const EventData &);

	std::vector<EventData> getEvents() { return events_; };
private:
	std::vector<EventData> events_{};
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
} ChannelFitData;


typedef struct {
	unsigned int eventID;
	std::string TDCCorrTime;
	std::string date;
	std::vector<ChannelFitData> sipm;
} EventFitData;

#endif //RECOMORE_DATASTRUCTURES_H
