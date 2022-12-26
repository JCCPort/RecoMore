"""C++ readers for wavecatcher data files."""
from __future__ import annotations
import CReader
import typing

__all__ = [
    "ChannelData",
    "ChannelFitData",
    "EventData",
    "EventFitData",
    "PE",
    "ReadWCDataFile",
    "ReadWCDataFileBinary",
    "ReadWCDataFileDat",
    "WCData"
]


class ChannelData:
    def __init__(self, arg0: int, arg1: typing.List[float]) -> None: ...

    @property
    def channel(self) -> int:
        """
        :type: int
        """

    @channel.setter
    def channel(self, arg0: int) -> None:
        pass

    @property
    def waveform(self) -> typing.List[float]:
        """
        :type: typing.List[float]
        """

    @waveform.setter
    def waveform(self, arg0: typing.List[float]) -> None:
        pass

    pass


class ChannelFitData:
    def __init__(self, arg0: int, arg1: float, arg2: float, arg3: typing.List[PE]) -> None: ...

    @property
    def baseline(self) -> float:
        """
        :type: float
        """

    @baseline.setter
    def baseline(self, arg0: float) -> None:
        pass

    @property
    def ch(self) -> int:
        """
        :type: int
        """

    @ch.setter
    def ch(self, arg0: int) -> None:
        pass

    @property
    def pes(self) -> typing.List[PE]:
        """
        :type: typing.List[PE]
        """

    @pes.setter
    def pes(self, arg0: typing.List[PE]) -> None:
        pass

    @property
    def redChiSq(self) -> float:
        """
        :type: float
        """

    @redChiSq.setter
    def redChiSq(self, arg0: float) -> None:
        pass

    pass


class EventData:
    def __init__(self, arg0: int, arg1: str, arg2: str, arg3: typing.List[ChannelData]) -> None: ...

    @property
    def TDCCorrTime(self) -> str:
        """
        :type: str
        """

    @TDCCorrTime.setter
    def TDCCorrTime(self, arg0: str) -> None:
        pass

    @property
    def chData(self) -> typing.List[ChannelData]:
        """
        :type: typing.List[ChannelData]
        """

    @chData.setter
    def chData(self, arg0: typing.List[ChannelData]) -> None:
        pass

    @property
    def date(self) -> str:
        """
        :type: str
        """

    @date.setter
    def date(self, arg0: str) -> None:
        pass

    @property
    def eventID(self) -> int:
        """
        :type: int
        """

    @eventID.setter
    def eventID(self, arg0: int) -> None:
        pass

    pass


class EventFitData:
    def __init__(self, arg0: int, arg1: str, arg2: str, arg3: typing.List[ChannelFitData]) -> None: ...

    @property
    def SiPM(self) -> typing.List[ChannelFitData]:
        """
        :type: typing.List[ChannelFitData]
        """

    @SiPM.setter
    def SiPM(self, arg0: typing.List[ChannelFitData]) -> None:
        pass

    @property
    def TDCCorrTime(self) -> str:
        """
        :type: str
        """

    @TDCCorrTime.setter
    def TDCCorrTime(self, arg0: str) -> None:
        pass

    @property
    def date(self) -> str:
        """
        :type: str
        """

    @date.setter
    def date(self, arg0: str) -> None:
        pass

    @property
    def eventID(self) -> int:
        """
        :type: int
        """

    @eventID.setter
    def eventID(self, arg0: int) -> None:
        pass

    pass


class PE:
    def __init__(self, arg0: float, arg1: float, arg2: float, arg3: float, arg4: float, arg5: float) -> None: ...

    @property
    def amplitude(self) -> float:
        """
        :type: float
        """

    @amplitude.setter
    def amplitude(self, arg0: float) -> None:
        pass

    @property
    def amplitudeError(self) -> float:
        """
        :type: float
        """

    @amplitudeError.setter
    def amplitudeError(self, arg0: float) -> None:
        pass

    @property
    def foundAmplitude(self) -> float:
        """
        :type: float
        """

    @foundAmplitude.setter
    def foundAmplitude(self, arg0: float) -> None:
        pass

    @property
    def foundTime(self) -> float:
        """
        :type: float
        """

    @foundTime.setter
    def foundTime(self, arg0: float) -> None:
        pass

    @property
    def time(self) -> float:
        """
        :type: float
        """

    @time.setter
    def time(self, arg0: float) -> None:
        pass

    @property
    def timeError(self) -> float:
        """
        :type: float
        """

    @timeError.setter
    def timeError(self, arg0: float) -> None:
        pass

    pass


class WCData:
    def addRow(self, arg0: EventData) -> None:
        """
        Add entry row to WCData
        """

    def getEvents(self) -> typing.List[EventData]:
        """
        Get rows
        """

    pass


def ReadWCDataFile(arg0: str) -> WCData:
    """
    Read WaveCatcher data files.
    """


def ReadWCDataFileBinary(arg0: str) -> WCData:
    """
    Read binary WaveCatcher data files.
    """


def ReadWCDataFileDat(arg0: str) -> WCData:
    """
    Read plain text WaveCatcher data files.
    """


__version__ = 'dev'
