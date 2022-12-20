import dataclasses
import typing

import numpy as np


"""Raw data types"""


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


"""RecoMore data types"""


@dataclasses.dataclass
class PE:
    def __init__(self):
        self.amplitude: float
        self.amplitudeError: float
        self.time: float
        self.timeError: float

        self.foundAmplitude: float
        self.foundTime: float


@dataclasses.dataclass
class RMChannelEvent:
    def __init__(self):
        self.channel: int = 0
        self.baseline: float = 0
        self.redChiSq: float = 0
        self.PEs: typing.List[PE] = []


@dataclasses.dataclass
class RMEvent:
    def __init__(self):
        self.eventNumber: int = 0
        self.TDCCorrTime: str
        self.date: str
        self.channels: typing.List[RMChannelEvent] = []