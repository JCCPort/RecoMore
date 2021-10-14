#ifndef RECOMORE_DATASTRUCTURES_H
#define RECOMORE_DATASTRUCTURES_H

#include <utility>
#include <vector>

struct WaveformData {
	unsigned short channel;
	std::vector<float> waveform{};
};

typedef struct {
	unsigned int eventID;
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
	double amplitude;
	double amplitudeError;
	double time;
	double timeError;

	// For debugging purpose, initial estimates of parameters.
	double foundAmplitude;
	double foundTime;
} PEData;

typedef struct {
	unsigned short ch;
	float chi2ndf;
	float baseline;
	std::vector<PEData> pes;
} ChannelFitData;


typedef struct {
	unsigned int eventID;
	std::vector<ChannelFitData> sipm;
} EventFitData;

#endif //RECOMORE_DATASTRUCTURES_H
