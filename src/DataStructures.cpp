#include "../include/DataStructures.h"

void WCData::addRow(const EventData& wf) {
	events_.emplace_back(wf);
}
