#ifndef RECOMORE_DATASTRUCTURES_H
#define RECOMORE_DATASTRUCTURES_H

#include <utility>
#include <vector>
#include <string>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/vector.hpp>


// Input types

struct DigitiserChannel {
	unsigned int channel;
	std::vector<float> waveform{};
};

struct DigitiserEvent{
	DigitiserChannel getDigitiserChannel(unsigned int channelNumber){
		for(auto & channelWF : chData){
			if(channelWF.channel == channelNumber){
				return channelWF;
			}
		}
		throw std::runtime_error("Channel " + std::to_string(channelNumber) + " not found in event " + std::to_string(eventID) + ".");
	}
	
	std::vector<unsigned int> getChannels(){
		std::vector<unsigned int> channels;
		for(const auto& channel: chData){
			channels.push_back(channel.channel);
		}
		return channels;
	}
	
	unsigned int                  eventID;
	std::string                   TDCCorrTime;
	std::string                   date;
	std::vector<DigitiserChannel> chData;
};


class DigitiserRun {
public:
	void addEvent(const DigitiserEvent &);
	
	 std::vector<DigitiserEvent> getEvents() { return events_; };
	 DigitiserEvent getEvent(int eventNumber);
	 DigitiserChannel getChannelWaveform(int eventNumber, int channelNumber);
private:
	std::vector<DigitiserEvent> events_{};
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


// Output types

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
	
	ChannelFitData getChannel(int channelNumber){
		for(auto & channelWF : SiPM){
			if(channelWF.ch == channelNumber){
				return channelWF;
			}
		}
		throw std::runtime_error("Channel " + std::to_string(channelNumber) + " not found in event " + std::to_string(eventID) + ".");
	}
	
	std::vector<unsigned short> getChannels(){
		std::vector<unsigned short> channels;
		for(const auto& channel: SiPM){
			channels.push_back(channel.ch);
		}
		return channels;
	}
	
	unsigned int eventID{};
	std::string TDCCorrTime;
	std::string date;
	std::vector<ChannelFitData> SiPM;
};


class FitData {
public:
	void addRow(const EventFitData &);
	void setRows(const std::vector<EventFitData> &);
	
	std::vector<EventFitData> getFitEvents() { return fitEvents_; };
	EventFitData getEventFit(int eventNumber);
	ChannelFitData getChannelFit(int eventNumber, int channelNumber);
private:
	std::vector<EventFitData> fitEvents_{};
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
