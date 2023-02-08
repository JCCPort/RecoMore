#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "DataReading.cpp"
#include "DataStructures.cpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)


namespace py = pybind11;

PYBIND11_MODULE(CReader, m) {
	m.doc() = "C++ readers for wavecatcher data files.";
	py::class_<PEData>(m, "PE")
			.def(py::init<const float&, const float&, const float&, const float&, const float&, const float&>())
			.def_readwrite("amplitude", &PEData::amplitude)
			.def_readwrite("amplitudeError", &PEData::amplitudeError)
			.def_readwrite("time", &PEData::time)
			.def_readwrite("timeError", &PEData::timeError)
			.def_readwrite("foundAmplitude", &PEData::foundAmplitude)
			.def_readwrite("foundTime", &PEData::foundTime);
	
	py::class_<ChannelFitData>(m, "ChannelFitData", py::dynamic_attr())
			.def(py::init<const unsigned short&, const float&, const float&, const std::vector<PEData>&>())
			.def_readwrite("ch", &ChannelFitData::ch)
			.def_readwrite("redChiSq", &ChannelFitData::redChiSq)
			.def_readwrite("baseline", &ChannelFitData::baseline)
			.def_readwrite("pes", &ChannelFitData::pes);
	
	py::class_<EventFitData>(m, "EventFitData", py::dynamic_attr())
			.def(py::init<const unsigned int&, const std::string&, const std::string&, const std::vector<ChannelFitData>>())
			.def_readwrite("eventID", &EventFitData::eventID)
			.def_readwrite("TDCCorrTime", &EventFitData::TDCCorrTime)
			.def_readwrite("date", &EventFitData::date)
			.def_readwrite("SiPM", &EventFitData::SiPM);
	
	py::class_<DigitiserChannel>(m, "ChannelData", py::dynamic_attr())
			.def(py::init<const unsigned short&, const std::vector<float>&>())
			.def_readwrite("channel", &DigitiserChannel::channel)
			.def_readwrite("waveform", &DigitiserChannel::waveform);
	
	py::class_<DigitiserEvent>(m, "EventData", py::dynamic_attr())
			.def(py::init<const unsigned int&, const std::string&, const std::string&, const std::vector<DigitiserChannel>>())
			.def_readwrite("eventID", &DigitiserEvent::eventID)
			.def_readwrite("TDCCorrTime", &DigitiserEvent::TDCCorrTime)
			.def_readwrite("date", &DigitiserEvent::date)
			.def_readwrite("chData", &DigitiserEvent::chData);
	
	py::class_<DigitiserRun>(m, "WCData")
			.def("addEvent", &DigitiserRun::addEvent, "Add entry row to WCData")
			.def("getEvents", &DigitiserRun::getEvents, "Get all events")
			.def("getEvent", &DigitiserRun::getEvent, "Get event by event number")
			.def("getChannelWaveform", &DigitiserRun::getChannelWaveform, "Get channel waveform by event and channel number");
	
	py::class_<FitData>(m, "FitData")
			.def("addEvent", &FitData::addRow, "Add entry row to FitData")
			.def("setRows", &FitData::setRows, "Set whole event fit vector at once.")
			.def("getFitEvents", &FitData::getFitEvents, "Get all fit events")
			.def("getEventFit", &FitData::getEventFit, "Get event fit by event number")
			.def("getChannelFit", &FitData::getChannelFit, "Get channel fit by event and channel number");
	
	m.def("ReadWCDataFileDat", &ReadWCDataFileDat, "Read plain text WaveCatcher data files.");
	m.def("ReadWCDataFileBinary", &ReadWCDataFileBinary, "Read binary WaveCatcher data files.");
	m.def("ReadWCDataFile", &ReadWCDataFile, "Read WaveCatcher data files.");
	m.def("ReadRecoMoreTextOutput", &ReadRecoMoreTextOutput, "Read RecoMore output text files.");
	m.def("ReadRecoMoreBinaryOutput", &ReadRecoMoreBinaryOutput, "Read RecoMore output binary files.");
	m.def("ReadRecoMoreOutput", &ReadRecoMoreOutput, "Read RecoMore output files.");

#ifdef VERSION_INFO
	m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
	m.attr("__version__") = "dev";
#endif
}
