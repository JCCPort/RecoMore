#ifndef RECOMORE_DATAREADING_H
#define RECOMORE_DATAREADING_H

#include "DataStructures.h"

DigitiserRun ReadWCDataFile(const std::string &, bool positivePulse);

DigitiserRun ReadWCDataFileDat(const std::string &, bool positivePulse);

DigitiserRun ReadWCDataFileBinary(const std::string &, bool positivePulse);

std::vector<double> readIdealWFs(unsigned int ch, unsigned int interpFactor, const std::string &idealWFDir, unsigned int expectedSize, bool positivePulse);

FitRun ReadRecoMoreOutput(const std::string &fileName);

FitRun ReadRecoMoreTextOutput(const std::string &fileName);

FitRun ReadRecoMoreBinaryOutput(const std::string &);

#endif //RECOMORE_DATAREADING_H
