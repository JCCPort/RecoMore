#ifndef RECOMORE_DATAREADING_H
#define RECOMORE_DATAREADING_H

#include "DataStructures.h"

bool startsWith(const char *a, const char *b);

WCData ReadWCDataFile(const std::string& fileName);

#endif //RECOMORE_DATAREADING_H
