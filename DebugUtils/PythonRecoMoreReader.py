import dataclasses
import typing

import re as re
import datetime
from copy import copy

import numpy as np

# WAVECATCHER REGEX SEARCHES ================================================================================


versionRegex = re.compile(r"=== DATA FILE SAVED WITH SOFTWARE VERSION: V(?P<ver>\d*\.\d*\.\d*) ===")

sysInfoRegex = re.compile(r"=== WAVECATCHER SYSTEM OF TYPE (?P<type>\d*) WITH (?P<numCh>\d*) CHANNELS AND GAIN: ("
                          r"?P<gain>\d*\.\d*) ===\n")

coincidenceRegex = re.compile(r"=== Rate coincidence masks (4 bits) for ch 0 to 15: 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 ==="
                              r" Posttrig in ns for SamBlock 0 to 8: 210 210 210 210 210 210 210 210 210 ===\n")

dataSampleRegex = re.compile(r"=== DATA SAMPLES \[(?P<numSamples>\d*)] in Volts =="
                             r" NB OF CHANNELS ACQUIRED: (?P<numCh>\d*) =="
                             r" Sampling Period: (?P<samplePeriod>\d*\.\d*) ps =="
                             r" INL Correction: (?P<INLCorr>\d*) =="
                             r" MEASUREMENTS: (?P<measurs>\d*) ===\n")

eventRegex = re.compile(r"=== EVENT (?P<evt>\d*) ===\n")

timeRegex = re.compile(r"=== UnixTime = (?P<unixTime>\d*\.\d*)"
                       r" date = (?P<date>\d*\.\d*\.\d*)"
                       r" time = (?P<time>\d*h\.\d*m\.\d*s\.\d*ms) =="
                       r" TDC = (?P<TDC>\d*) =="
                       r" TDC corrected time = (?P<corrTime>\d*h\d*m\d*s,\d*\.\d*\.\d*ns) =="
                       r" Nb of channels = (?P<numCh>\d*) ===\n")

channelRegex = re.compile(r"=== CH: (?P<ch>\d*) EVENTID: (?P<eid>\d*) FCR: (?P<fcr>\d*) ===\n")


# RECOMORE REGEX SEARCHES ================================================================================

eventRegexRM = re.compile(r"EVENT=(?P<evt>\d*), DATE=(?P<date>\d*\.\d*\.\d*), TDCCorrTime=(?P<corrTime>\d*h\d*m\d*\.\d*s)\n")

channelRegexRM = re.compile(r"Ch=(?P<ch>\d*), RedChiSq=(?P<redChiSq>\d*\.\d*), Baseline=(?P<baseline>[+-]?\d*\.\d*)\n")

peRegexRM = re.compile(
    r"(?P<amp>-?\d*\.\d*),(?P<time>-?\d*\.\d*)\n")

dateTimeRegex = re.compile(r"(?P<year>\d*)\.(?P<month>\d*)\.(?P<day>\d*) (?P<hour>\d*)h(?P<min>\d*)m(?P<sec>\d*\.\d*)s")


@dataclasses.dataclass
class PE:
    def __init__(self, amplitude, amplitudeError, time, timeError):
        self.amplitude: float = amplitude
        self.amplitudeError: float = amplitudeError
        self.time: float = time
        self.timeError: float = timeError


@dataclasses.dataclass
class RMChannelEvent:
    def __init__(self):
        self.channel: int = 0
        self.baseline: float = 0
        self.PEs: typing.List[PE] = []


@dataclasses.dataclass
class RMEvent:
    def __init__(self):
        self.eventNumber: int = 0
        self.TDCCorrTime: str
        self.date: str
        self.channels: typing.List[RMChannelEvent] = []


@dataclasses.dataclass
class Run:
    def __init__(self):
        self.Events: typing.List[RMEvent] = []


@dataclasses.dataclass
class RawChannelEvent:
    def __init__(self):
        self.channelNumber: int = -1
        self.rawData: np.array = None


@dataclasses.dataclass
class RawEvent:
    def __init__(self):
        self.eventNumber: int = 0
        self.TDCCorrTime: str
        self.date: str
        self.channels: typing.List[RawChannelEvent] = []


class ChannelWF:
    def __init__(self):
        self.eventTime: datetime = datetime.datetime(1970, 1, 1)
        self.extraTimePrecision: float = 0
        self.eventNumber: int = -1
        self.channelNumber: int = -1
        self.xVals = None

        self.rawData: np.array = None


class RecoMoreEvent:
    def __init__(self):
        self.eventTime: datetime = datetime.datetime(1970, 1, 1)
        self.extraTimePrecision: float = 0
        self.eventNumber: int = -1
        self.channelNumber: int = -1
        self.redChiSq: float = 0
        self.baseline: float = 0

        self.PEData: typing.List[PE] = []


def plotWF(channelWF, *args, **kwargs):
    import matplotlib.pyplot as plt
    xVals = np.array(range(0, len(channelWF.rawData))) * channelWF.parentRun.WCConfig.samplePeriod
    plt.plot(xVals, channelWF.rawData,
             label="Event: {}\n"
                   "Time: {}\n"
                   "Channel: {}".format(channelWF.eventNumber, channelWF.eventTime, channelWF.channelNumber),
             *args, **kwargs)
    plt.xlabel("Time (ns)")
    plt.ylabel("Voltage (V)")
    plt.legend()


def readWCWaveforms(WCDataFile) -> typing.List[ChannelWF]:
    line: str = WCDataFile.readline()
    WFs: typing.List[ChannelWF] = []

    tempWF = ChannelWF()
    while line:
        if "=== EVENT" in line:
            eventSearch = eventRegex.search(line)

            tempWF.eventNumber = int(eventSearch.group("evt"))
            line = WCDataFile.readline()

            timeSearch = timeRegex.search(line)
            dateString = timeSearch.group("date")

            if dateString == "0.0.0":
                print("Time failed to be written to file for event {}".format(eventSearch.group("evt")))
                tempWF.eventTime = datetime.datetime(1970, 1, 1)
            else:
                tempWF.unixTime = float(timeSearch.group("unixTime"))
                tempWF.numChannels = int(timeSearch.group("numCh"))

                testString = timeSearch.group("corrTime").replace("h", ".")
                testString = testString.replace("m", ".")
                testString = testString.replace(" ", ".")
                testString = testString.replace("s,", ".")
                testString2 = testString.replace(".", "")

                extraTimePrecision = testString2[6:-2]
                testString = timeSearch.group("date") + '.' + testString[:8]
                time = datetime.datetime(*map(int, testString.split('.')))
                tempWF.eventTime = time + datetime.timedelta(seconds=float("." + extraTimePrecision))
                tempWF.extraTimePrecision = int(testString2[12:-2])

            line = WCDataFile.readline()

        if "=== CH:" in line:
            channelSearch = channelRegex.search(line)
            tempWF.channelNumber = int(channelSearch.group("ch"))
            line = WCDataFile.readline()

            tempWF.rawData = np.asarray(line.split(' ')[:-1], dtype=np.float32)
            tempWF.xVals = [i * 0.3125 for i in range(len(tempWF.rawData))]
            WFs.append(copy(tempWF))
            # tempWF = ChannelWF()
        line = WCDataFile.readline()
    return WFs


def readRecoMore(RMDataFile) -> typing.List[RecoMoreEvent]:
    line = RMDataFile.readline()
    WFs = []
    while line:
        tempWF = RecoMoreEvent()
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
                tempWF.redChiSq = float(chRegSearch.group("redChiSq"))
                tempWF.baseline = float(chRegSearch.group("baseline"))
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
