#ifndef RECOMORE_DATAREADING_H
#define RECOMORE_DATAREADING_H

#include "DataStructures.h"

DigitiserRun ReadWCDataFile(const std::string &);

DigitiserRun ReadWCDataFileDat(const std::string &);

DigitiserRun ReadWCDataFileBinary(const std::string &);

std::vector<double> readIdealWFs(unsigned int, int, const std::string &, unsigned int);

FitRun ReadRecoMoreOutput(const std::string &fileName);

FitRun ReadRecoMoreTextOutput(const std::string &fileName);

FitRun ReadRecoMoreBinaryOutput(const std::string &);

#endif //RECOMORE_DATAREADING_H
