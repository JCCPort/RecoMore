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


#endif //RECOMORE_DATASTRUCTURES_H
