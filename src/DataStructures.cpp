#include "../include/DataStructures.h"

void WCData::addRow(const Waveform& wf) {
	waveforms_.emplace_back(wf);
}
