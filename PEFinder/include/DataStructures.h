#ifndef RECOMORE_DATASTRUCTURES_H
#define RECOMORE_DATASTRUCTURES_H

#include <utility>
#include <vector>
#include <string>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/vector.hpp>


struct ChannelData {
	unsigned short channel;
	std::vector<float> waveform{};
};

typedef struct {
	unsigned int             eventID;
	std::string              TDCCorrTime;
	std::string              date;
	std::vector<ChannelData> chData;
} EventData;


class WCData {
public:
	void addRow(const EventData &);

	std::vector<EventData> getEvents() { return events_; };
private:
	std::vector<EventData> events_{};
};


struct PEData {
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & amplitude;
		ar & amplitudeError;
		ar & time;
		ar & timeError;
	}
	float amplitude;
	float amplitudeError;
	float time;
	float timeError;

	// For debugging purpose, initial estimates of parameters.
	float foundAmplitude;
	float foundTime;
};

struct ChannelFitData {
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & ch;
		ar & redChiSq;
		ar & baseline;
		ar & pes;
	}
	unsigned short ch;
	float          redChiSq;
	float          baseline;
	std::vector<PEData> pes;
};


struct EventFitData{
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & eventID;
		ar & TDCCorrTime;
		ar & date;
		ar & SiPM;
	}
	unsigned int eventID{};
	std::string TDCCorrTime;
	std::string date;
	std::vector<ChannelFitData> SiPM;
};


class [[maybe_unused]] FitParams {
public:
	[[maybe_unused]] FitParams(unsigned int, double, std::vector<double>, std::vector<double>);

	[[maybe_unused]] FitParams(double, const std::vector<PEData> &);

	[[maybe_unused]] std::vector<double *> makeFitterParams();

	[[maybe_unused]] std::vector<float> makeGuesserParams();

	unsigned int numPEs_;
	std::vector<double> params;
};

#endif //RECOMORE_DATASTRUCTURES_H
