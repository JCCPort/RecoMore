import typing

import re as re
import datetime
from copy import copy

# RECOMORE REGEX SEARCHES ================================================================================

eventRegexRM = re.compile(r"EVENT=(?P<evt>\d*), DATE=(?P<date>\d*\.\d*\.\d*), TDCCorrTime=(?P<corrTime>\d*h\d*m\d*\.\d*s)\n")

channelRegexRM = re.compile(r"Ch=(?P<ch>\d*), RedChisq=(?P<redChisq>\d*\.\d*), Baseline=(?P<baseline>[+-]?\d*\.\d*)\n")

peRegexRM = re.compile(
    r"(?P<amp>\d*\.\d*),(?P<time>\d*\.\d*)\n")

dateTimeRegex = re.compile(r"(?P<year>\d*)\.(?P<month>\d*)\.(?P<day>\d*) (?P<hour>\d*)h(?P<min>\d*)m(?P<sec>\d*\.\d*)s")


class PE:
    def __init__(self, amplitude, amplitudeError, time, timeError):
        self.amplitude: float = amplitude
        self.amplitudeError: float = amplitudeError
        self.time: float = time
        self.timeError: float = timeError


class ChannelWF:
    def __init__(self):
        self.eventTime: datetime = datetime.datetime(1970, 1, 1)
        self.extraTimePrecision: float = 0
        self.eventNumber: int = -1
        self.channelNumber: int = -1

        self.PEData: typing.List[PE] = []


def readRecoMore(RMDataFile) -> typing.List[ChannelWF]:
    line = RMDataFile.readline()
    WFs = []
    while line:
        tempWF = ChannelWF()
        while "EVENT=" in line:
            lineSearch = eventRegexRM.search(line)
            tempWF.eventNumber = int(lineSearch.group("evt"))
            if lineSearch.group("date") == "0.0.0":
                print("Date was not correctly recorded for event {}".format(tempWF.eventNumber))
                tempWF.date = datetime.datetime(1970, 1, 1)
                dateTimeString = "1970.01.01 " + "00h00m00.000000000s"
            else:
                tempWF.date = datetime.datetime(*map(int, lineSearch.group("date").split('.')))
                dateTimeString = lineSearch.group("date") + " " + lineSearch.group("corrTime")

            testString = dateTimeString.replace("h", ".")
            testString = testString.replace("m", ".")
            testString = testString.replace(" ", ".")
            testString = testString[:-11]

            extraTimePrecision = dateTimeString.split(".")[-1].split("s")[0]

            time = datetime.datetime(*map(int, testString.split('.')))
            tempWF.eventTime = time + datetime.timedelta(seconds=float("." + extraTimePrecision))
            if len(extraTimePrecision) > 6:
                tempWF.extraTimePrecision = int(extraTimePrecision[6:])
            else:
                tempWF.extraTimePrecision = 0
            line = RMDataFile.readline()
            while "Ch=" in line:
                PEs = []
                chRegSearch = channelRegexRM.search(line)
                tempWF.channelNumber = int(chRegSearch.group("ch"))
                line = RMDataFile.readline()
                while line != "\n":
                    peRegSearch = peRegexRM.search(line)
                    if peRegSearch is None:
                        line = RMDataFile.readline()
                        continue
                    PEs.append(PE(float(peRegSearch.group("amp")), 0, float(peRegSearch.group("time")), 0))
                    line = RMDataFile.readline()
                tempWF.PEData = PEs
                line = RMDataFile.readline()
                if tempWF.date != datetime.datetime(1970, 1, 1):
                    WFs.append(copy(tempWF))
        line = RMDataFile.readline()
        line = RMDataFile.readline()
    return WFs
