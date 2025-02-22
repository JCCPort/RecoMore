#ifndef RECOMORE_DATASTRUCTURES_H
#define RECOMORE_DATASTRUCTURES_H

#include <utility>
#include <vector>
#include <string>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/vector.hpp>


// Input types

struct DigitiserChannel {
    unsigned int ID;
    std::vector<float> waveform{};
};

struct DigitiserEvent {
    DigitiserChannel getChannel(const unsigned int channelNumber) {
        for (auto &channel: channels) {
            if (channel.ID == channelNumber) {
                return channel;
            }
        }
        throw std::runtime_error(
                "Channel " + std::to_string(channelNumber) + " not found in event " + std::to_string(ID) + ".");
    }

    std::vector<unsigned int> getChannelIDs() {
        std::vector<unsigned int> channels_;
        channels_.reserve(channels.size());
        for (const auto &channel: channels) {
            channels_.push_back(channel.ID);
        }
        return channels_;
    }

    unsigned int ID{};
    std::string correctedTime;
    std::string date;
    std::vector<DigitiserChannel> channels;
};


class DigitiserRun {
public:
    void addEvent(const DigitiserEvent &);

    std::vector<DigitiserEvent> getEvents() { return events; };

    DigitiserEvent getEvent(unsigned int eventNumber);

    DigitiserChannel getEventChannel(unsigned int eventNumber, unsigned int channelNumber);

private:
    std::vector<DigitiserEvent> events{};
};


struct Photoelectron {
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version) {
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
    float initialAmplitude;
    float initialTime;
};


// Output types

struct FitChannel {
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & ID;
        ar & reducedChiSq;
        ar & baseline;
        ar & PEs;
    }

    unsigned int ID;
    float reducedChiSq;
    float baseline;
    std::vector<Photoelectron> PEs;
};


struct FitEvent {
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & ID;
        ar & correctedTime;
        ar & date;
        ar & channels;
    }

    FitChannel getChannel(const unsigned int channelNumber) {
        for (auto &channelWF: channels) {
            if (channelWF.ID == channelNumber) {
                return channelWF;
            }
        }
        throw std::runtime_error(
                "Channel " + std::to_string(channelNumber) + " not found in event " + std::to_string(ID) + ".");
    }

    std::vector<unsigned int> getChannelIDs() {
        std::vector<unsigned int> channels_;
        for (const auto &channel: channels) {
            channels_.push_back(channel.ID);
        }
        return channels_;
    }

    unsigned int ID{};
    std::string correctedTime;
    std::string date;
    std::vector<FitChannel> channels;
};


class FitRun {
public:
    void addEvent(const FitEvent &);

    void setEvents(const std::vector<FitEvent> &);

    std::vector<FitEvent> getEvents() { return events; };

    FitEvent getEvent(unsigned int eventNumber);

    FitChannel getEventChannel(unsigned int eventNumber, unsigned int channelNumber);

    std::vector<unsigned int> getEventIDs();

private:
    std::vector<FitEvent> events{};
};


class FitParams {
    public:

    FitParams() = default;

    // Constructor that creates Photoelectron objects from amplitude/time vectors.
    FitParams(unsigned int, double, const std::vector<double>&, const std::vector<double>&);

    // Constructor that accepts a vector of Photoelectron directly.
    FitParams(double, const std::vector<Photoelectron> &);

    void makeSolverParams(std::vector<double*>* solverParams, std::vector<double>* times, std::vector<double>* amplitudes, double* baseline);

    std::vector<float> makeGuesserParams();

    int getNumParams() const;

    void sortPEsByTime();

    void addPE(const Photoelectron &pe)
    {
        PEParams_.push_back(pe);
        numPEs_++;
    }

    unsigned int numPEs_ = 0;
    double baseline_ = 0.;
    std::vector<Photoelectron> PEParams_;
};

#endif //RECOMORE_DATASTRUCTURES_H
