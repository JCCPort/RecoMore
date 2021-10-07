#ifndef RECOMORE_DATAREADING_H
#define RECOMORE_DATAREADING_H

#include "DataStructures.h"

WCData ReadWCDataFile(const std::string&);

std::vector<float> readIdealWFs(unsigned int, int, const std::string&, unsigned int);

#endif //RECOMORE_DATAREADING_H
