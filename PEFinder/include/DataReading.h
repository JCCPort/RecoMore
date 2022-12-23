#ifndef RECOMORE_DATAREADING_H
#define RECOMORE_DATAREADING_H

#include "DataStructures.h"

WCData ReadWCDataFile(const std::string &);

WCData ReadWCDataFileDat(const std::string &);

WCData ReadWCDataFileBinary(const std::string &fileName);

std::vector<double> readIdealWFs(unsigned int, int, const std::string &, unsigned int);

#endif //RECOMORE_DATAREADING_H
