#ifndef RECOMORE_DATAREADING_H
#define RECOMORE_DATAREADING_H

#include "DataStructures.h"
#include "../../CReader_/include/DataStructures.h"
#include "DataReading.h"

WCData ReadWCDataFile(const std::string &);

WCData ReadWCDataFileDat(const std::string &);

WCData ReadWCDataFileBinary(const std::string &);

std::vector<double> readIdealWFs(unsigned int, int, const std::string &, unsigned int);

FitData ReadRecoMoreFile(const std::string &fileName);

FitData ReadRecoMoreOutput(const std::string &fileName);

FitData ReadRecoMoreBinaryOutput(const std::string &);

#endif //RECOMORE_DATAREADING_H
