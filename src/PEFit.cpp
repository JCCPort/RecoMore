#include "../include/PEFit.h"
#include "../Globals.h"
#include <cmath>
#include "ceres/ceres.h"
#include "ceres/cubic_interpolation.h"
#include "../Utils.h"
#include <memory>
#include <utility>

const double samplingRate2Inv = 1 / (0.01 * pdfSamplingRate);
const double pdfT0SampleConv = (double) pdfT0Sample;

struct npe_pdf_functor {
	npe_pdf_functor(double x, double y, ceres::CubicInterpolator<ceres::Grid1D<double, true>> *PDFInterpolator,
	                unsigned int numPES) : x_(
			x), y_(y), PDFInterpolator_(PDFInterpolator), numPES_(numPES) {};

	template<typename T>
	bool operator()(T const *const *params, T *residual) const {
		T f;
		T X_(x_);
		residual[0] = params[0][0];
		for (unsigned int i = 0; i < numPES_; i++) {
			unsigned int i2 = 2 * i;
			PDFInterpolator_->Evaluate((X_ - params[i2 + 2][0]), &f);
			residual[0] += (params[i2 + 1][0] * f);
		}
		residual[0] -= y_;
		return true;
	}

private:
	const double x_;
	const double y_;
	ceres::CubicInterpolator<ceres::Grid1D<double> > *PDFInterpolator_;
	const unsigned int numPES_;
};


float npe_pdf_func(float X, const std::vector<float> &p, std::vector<double> *idealWaveform) {
	// This way of passing x as a list then choosing the 0th index is from ROOT's syntax for fitting where you can fit
	//  with an arbitrary number of dimensions (N input values -> single output, an N dimensional function)

	// parameter of the function:
	// p[0] : number N of PE (fixed)
	// p[1] : baseline
	// p[2] : ampl of PE #0
	// p[3] : time of PE #0
	// [...]
	// p[2+N*2] : ampl of PE #N
	// p[3+N*2] : time of PE #N

	const int NPE = (int) (p[0]);
	const float BASELINE = p[1];

	float value = BASELINE;

	// This is adding up the contribution of each fit PE to a specific bin
	for (int PE = 0; PE < NPE; ++PE) {
		const float PE_CHARGE = p[2 + PE * 2];
		const float PE_TIME = p[3 + PE * 2];
		// TODO(josh): Start using interpolation for this too instead of using floor
		const int PE_PDF_BIN =
				pdfT0Sample + (int) std::floor(0.5 + (X - PE_TIME) * samplingRate2Inv);
		if ((PE_PDF_BIN >= 0) && (PE_PDF_BIN < pdfNSamples)) {
			float thing = idealWaveform->at(PE_PDF_BIN);
			value += PE_CHARGE * thing;
		}
	}
	return value;
}


float ampDiff = 0;
float timeDiff = 0;
float baselineDiff = 0;

int N = 0;
int N2 = 0;

void
fitPE(const EventData *event, const std::vector<std::vector<double>> *idealWaveforms, std::shared_ptr<SyncFile> file) {
	EventFitData evFitDat;
	evFitDat.eventID = event->eventID;
	evFitDat.date = event->date;
	evFitDat.TDCCorrTime = event->TDCCorrTime;
	for (const auto &waveformData: event->chData) {
		auto residualWaveform = waveformData; // This will be the variable that is modified to be the residual distribution after each iteration

		ChannelFitData chFit;
		unsigned int ch = waveformData.channel;
		if (ch == 15) {
			return;
		}
		chFit.ch = ch;

		// Making a pointer to the ideal waveform for this channel to improve speed of passing.
		std::vector<double> tmp = (*idealWaveforms)[ch];
		std::vector<double> *chIdealWaveform = &tmp;

		// Baseline calculation
		float initBaseline = averageVector(waveformData.waveform, 0, 150, 0.0015);
		chFit.baseline = initBaseline; // Will want to replace this with the fit baseline

		// Start loop that will break when no more PEs are present
		std::vector<PEData> pesFound;
		unsigned int numPEsFound = 0;

		// This is getting estimate PEs that will then be passed as initial guesses to the minimiser.
		PEData guessPE;
		while (true) {

			std::sort(pesFound.begin(), pesFound.end(), comparePETime()); // Need to sort for amplitude adjustment.

			// TODO(josh): Create a struct that corresponds to the params (number PEs, baseline, PE info...)
			std::vector<float> params;
			params.push_back((float) pesFound.size());
			params.push_back(initBaseline);
			for (auto pe: pesFound) {
				params.push_back(pe.amplitude);
				params.push_back(pe.time);
			}

			// Amplitude adjustment: if the latest PE found is before other one(s),
			//  its tail is going to add some amplitude to the following one. Compares
			//  real and fit amplitude at the time bin corresponding to the PE time.
			for (auto &pe: pesFound) {
				// TODO(josh): Improve the adjustment by averaging the shift based off of a few bins around the PE time
				unsigned int peTimeBinPos = std::floor(pe.time / pdfSamplingRate); // This should use a variable
				// corresponding to the input sampling rate, however for now they're the same
				float fitVal = npe_pdf_func(pe.time, params, chIdealWaveform);
				float extraAmplitude = fitVal - residualWaveform.waveform[peTimeBinPos];
				float newAmplitude = pe.amplitude + extraAmplitude;
				if (newAmplitude > 0) { // TODO(josh): We need to consider the situations that this would ever be true
					pe.amplitude = newAmplitude;
				}
			}

			params = std::vector<float>{};
			params.push_back((float) pesFound.size());
			params.push_back(initBaseline);
			for (auto pe: pesFound) {
				params.push_back(pe.amplitude);
				params.push_back(pe.time);
			}

//			writeVector("rawWaveform.csv", waveformData.waveform);

			// Compute residual
			// TODO(josh): I suspect the residual is running into issues for short waveforms.
//			std::vector<float> fitVecForPlot;
			for (unsigned int k = 0; k < residualWaveform.waveform.size(); ++k) {
				// TODO(josh): Should it be k or k + 0.5?
				float val = npe_pdf_func(float(k) * pdfSamplingRate, params, chIdealWaveform);
				residualWaveform.waveform[k] = residualWaveform.waveform[k] - val;
//				fitVecForPlot.emplace_back(val);
			}

//			writeVector("fit.csv", fitVecForPlot);

//			writeVector("residual.csv", residualWaveform.waveform);


			// Get initial guesses for the next PE
			auto minPosIt = std::min_element(residualWaveform.waveform.begin(), residualWaveform.waveform.end());
			unsigned int minTimePos = std::distance(residualWaveform.waveform.begin(), minPosIt);

			if (-residualWaveform.waveform[minTimePos] < 0.009) {
				break;
			}

			guessPE.amplitude = -residualWaveform.waveform[minTimePos];

			guessPE.time = float(minTimePos) * pdfSamplingRate;

			if ((minTimePos > 1) && (minTimePos < residualWaveform.waveform.size())) {
				// improve initial time for a new PE based on average time
				// over 3 consecutive sample ponderated by the amplitude
				// of each sample... help a lot to resolve PEs very close!
				double timeSum = 0;
				double ponderationSum = 0;
				for (unsigned int b = minTimePos - 1; b <= minTimePos + 1; ++b) {
					double binCenter = (minTimePos - 0.5) * (pdfSamplingRate);
					double binVal = residualWaveform.waveform[minTimePos];
					timeSum += binCenter * binVal;
					ponderationSum += binVal;
				}
				guessPE.time = float(PEFinderTimeOffset * 0.1) + timeSum / ponderationSum;
			}

			numPEsFound += 1;
			pesFound.push_back(guessPE);
			residualWaveform = waveformData;

			if (numPEsFound > 100) {
				break;
			}
		}

		using namespace ceres;
		Problem problem;

		std::vector<double> initialTimes;
		std::vector<double> initialAmplitudes;
		for (auto &k: pesFound) {
			initialAmplitudes.push_back(k.amplitude);
			initialTimes.push_back(k.time);
		}

		std::vector<double> times = initialTimes;
		std::vector<double> amplitudes = initialAmplitudes;

		// Applying running correction to parameter estimates.
		for (int i = 0; i < times.size(); i++) {
			times[i] += timeDiff;
			amplitudes[i] += ampDiff;
		}

		// Converting the time into an ideal waveform PDF index to simplify npe_pdf_functor method call
		for (double &time: times) {
			time = time * samplingRate2Inv;
		}

		// Creating interpolator to allow for the ideal waveform to be used as a continuous (and differentiable) function.
		ceres::Grid1D grid = ceres::Grid1D<double>(chIdealWaveform->data(), 0, (int) chIdealWaveform->size());
		auto idealPDFInterpolator = new ceres::CubicInterpolator<ceres::Grid1D<double> >(grid);


		// Creating vector of references to parameter values that the fitter will use and modify. Note this means that
		// the references are the initial values before the fit and the final values after the fit.
		double baseline = initBaseline;
		std::vector<double *> params;
		params.push_back(&baseline);
		for (int i = 0; i < pesFound.size(); i++) {
			params.push_back(&amplitudes[i]);
			params.push_back(&times[i]);
		}

		// Creating the x values that the solver will use
		std::vector<double> xValues;
		for (unsigned int j = 0; j < waveformData.waveform.size(); j++) {
			xValues.push_back(
					((double) j * 100) + pdfT0SampleConv);  // Multiplying index to match position on ideal PDF
		}

		// Set up the only cost function (also known as residual). This uses
		// auto-differentiation to obtain the derivative (Jacobian).
		for (unsigned int j = 0; j < waveformData.waveform.size(); ++j) {
			auto costFunction = new DynamicAutoDiffCostFunction<npe_pdf_functor, 4>(new npe_pdf_functor(xValues[j],
			                                                                                            waveformData.waveform[j],
			                                                                                            idealPDFInterpolator,
			                                                                                            pesFound.size()));
			costFunction->AddParameterBlock(1);
			for (auto pe: pesFound) {
				costFunction->AddParameterBlock(1);
				costFunction->AddParameterBlock(1);
			}
			costFunction->SetNumResiduals(1);
			problem.AddResidualBlock(costFunction, nullptr, params);
		}

		// Run the solver!
		Solver::Options options;
		options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
		options.minimizer_progress_to_stdout = false;
		Solver::Summary summary;
		Solve(options, &problem, &summary);

		// Going back from ideal waveform PDF index to time
		for (double &time: times) {
			time = time / samplingRate2Inv;
		}

		std::vector<float> finalParams;
		finalParams.push_back(pesFound.size());
		finalParams.push_back(initBaseline);
		for (int k = 0; k < pesFound.size(); k++) {
			finalParams.push_back((float) amplitudes[k]);
			finalParams.push_back((float) times[k]);
		}

		double chiSq = 0;
		for (int j = 0; j < waveformData.waveform.size(); j++) {
			double observed = waveformData.waveform[j];
			double expected = npe_pdf_func(float(j) * pdfSamplingRate, finalParams, chIdealWaveform);
			chiSq += std::pow(observed - expected, 2) /
			         (pdfResidualRMS / 1000); //TODO(josh): Need to calculate the residual RMS on a per waveform basis?
		}
		double redChiSq = chiSq / (finalParams.size() - 1);
		N2++;
		meanReducedChisq = meanReducedChisq + (redChiSq - meanReducedChisq) / N2;

		delete idealPDFInterpolator;

//		std::cout << summary.FullReport() << "\n";


//
//		std::vector<float> fullFitVecForPlot;
//		for (unsigned int k = 0; k < waveformData.waveform.size(); k++) {
//			auto thing = npe_pdf_func(k * pdfSamplingRate, finalParams, chIdealWaveform);
//			fullFitVecForPlot.emplace_back(thing);
//		}
//
//		writeVector("fullFit.csv", fullFitVecForPlot);

//		for (int k = 0; k < pesFound.size(); k++) {
//			std::cout << "Amplitude:\t" << initialAmplitudes[k] << " -> " << amplitudes[k] << std::endl;
//			std::cout << "Time:\t\t" << initialTimes[k] << " -> " << times[k] << std::endl;
//			std::cout << std::endl;
//		}

		std::vector<PEData> FitPEs;

		for (int k = 0; k < pesFound.size(); k++) {
			PEData pe;
			pe.amplitude = float(amplitudes[k]);
			pe.time = float(times[k]);
			FitPEs.emplace_back(pe);
			N++;
			ampDiff = ampDiff + ((pe.amplitude - initialAmplitudes[k]) - ampDiff) / N;
			timeDiff = timeDiff + ((pe.time - initialTimes[k]) - timeDiff) / N;
			baselineDiff = baselineDiff + ((baseline - initBaseline) - baselineDiff) / N;
		}

		chFit.pes = FitPEs;
		chFit.baseline = float(baseline);
		evFitDat.sipm.push_back(chFit);


		// TODO(josh): No magic numbers for different lengths of data
	}
	Writer writer(std::move(file));
	writer.writeEventInfo(evFitDat);
}