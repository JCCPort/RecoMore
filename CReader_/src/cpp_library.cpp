#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "DataReading.cpp"
#include "DataStructures.cpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)


namespace py = pybind11;

PYBIND11_MODULE(CReader, m) {
	m.doc() = "C++ readers for wavecatcher data files.";
// TODO(josh): Change python class names to match C++ class names... a bit unnecessarily confusing
	py::class_<Photoelectron>(m, "PE")
			.def(py::init<const float&, const float&, const float&, const float&, const float&, const float&>())
			.def_readwrite("amplitude", &Photoelectron::amplitude)
			.def_readwrite("amplitudeError", &Photoelectron::amplitudeError)
			.def_readwrite("time", &Photoelectron::time)
			.def_readwrite("timeError", &Photoelectron::timeError)
			.def_readwrite("foundAmplitude", &Photoelectron::initialAmplitude)
			.def_readwrite("foundTime", &Photoelectron::initialTime);
	
	py::class_<FitChannel>(m, "ChannelFitData", py::dynamic_attr())
			.def(py::init<const unsigned short&, const float&, const float&, const std::vector<Photoelectron>&>())
			.def_readwrite("ch", &FitChannel::ID)
			.def_readwrite("redChiSq", &FitChannel::reducedChiSq)
			.def_readwrite("baseline", &FitChannel::baseline)
			.def_readwrite("pes", &FitChannel::PEs);
	
	py::class_<FitEvent>(m, "EventFitData", py::dynamic_attr())
			.def(py::init<const unsigned int&, const std::string&, const std::string&, const std::vector<FitChannel>>())
			.def_readwrite("eventID", &FitEvent::ID)
			.def_readwrite("TDCCorrTime", &FitEvent::correctedTime)
			.def_readwrite("date", &FitEvent::date)
			.def_readwrite("SiPM", &FitEvent::channels);
	
	py::class_<DigitiserChannel>(m, "ChannelData", py::dynamic_attr())
			.def(py::init<const unsigned short&, const std::vector<float>&>())
			.def_readwrite("channel", &DigitiserChannel::ID)
			.def_readwrite("waveform", &DigitiserChannel::waveform);
	
	py::class_<DigitiserEvent>(m, "EventData", py::dynamic_attr())
			.def(py::init<const unsigned int&, const std::string&, const std::string&, const std::vector<DigitiserChannel>>())
			.def_readwrite("eventID", &DigitiserEvent::ID)
			.def_readwrite("TDCCorrTime", &DigitiserEvent::correctedTime)
			.def_readwrite("date", &DigitiserEvent::date)
			.def_readwrite("chData", &DigitiserEvent::channels);
	
	py::class_<DigitiserRun>(m, "WCData")
			.def("addEvent", &DigitiserRun::addEvent, "Add entry row to WCData")
			.def("getEvents", &DigitiserRun::getEvents, "Get all events")
			.def("getEvent", &DigitiserRun::getEvent, "Get event by event number")
			.def("getEventChannel", &DigitiserRun::getEventChannel, "Get channel waveform by event and channel number");
	
	py::class_<FitRun>(m, "FitData")
            .def_property("events", &FitRun::getEvents, &FitRun::setEvents)
			.def("addEvent", &FitRun::addEvent, "Add entry row to FitData")
			.def("setEvents", &FitRun::setEvents, "Set whole event fit vector at once.")
			.def("getEvents", &FitRun::getEvents, "Get all fit events")
			.def("getEvent", &FitRun::getEvent, "Get event fit by event number")
			.def("getEventChannel", &FitRun::getEventChannel, "Get channel fit by event and channel number");
	
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
