#include "../include/PEFit.h"
#include "../Utils.h"
#include <mutex>
#include <atomic>
#include <cmath>
#include "ceres/ceres.h"
#include "ceres/cubic_interpolation.h"
#include <memory>
#include <utility>


struct npe_pdf_functor {
	npe_pdf_functor(double x, double y, ceres::CubicInterpolator<ceres::Grid1D<double, true>> *PDFInterpolator,
	                unsigned int numPES) : x_(
			x), y_(y), PDFInterpolator_(PDFInterpolator), numPES_(numPES) {};

	template<typename T>
	bool operator()(T const *const *params, T *__restrict__ residual) const {
		T f;
		T X_(x_);
		residual[0] = params[0][0];
		for (unsigned int i = 0; i < numPES_; ++i) {
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
		const int PE_PDF_BIN = pdfT0Sample + (int) std::floor(0.5 + (X - PE_TIME) * samplingRate2Inv);
		if ((PE_PDF_BIN >= 0) && (PE_PDF_BIN < pdfNSamples)) {
			value += PE_CHARGE * idealWaveform->at(PE_PDF_BIN);
		}
	}
	return value;
}


void
fitPE(const EventData *event, const std::vector<std::vector<double>> *idealWaveforms, std::shared_ptr<SyncFile> file) {
	EventFitData evFitDat;
	evFitDat.eventID = event->eventID;
	evFitDat.date = event->date;
	evFitDat.TDCCorrTime = event->TDCCorrTime;

	for (const auto &WFData: event->chData) { // Looping through all channels for a given event
		auto residualWF = WFData; // This will be the variable that is modified to be the residual distribution after each iteration

		ChannelFitData chFit;
		unsigned int ch = WFData.channel;
		chFit.ch = ch;

		// Making a pointer to the ideal waveform for this channel to improve speed of passing.
		std::vector<double> tmp = (*idealWaveforms)[ch];
		std::vector<double> *chIdealWF = &tmp;

		// Baseline calculation
		float initBaseline = averageVector(WFData.waveform, 0, 50, 0.01);
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
//			auto tempRes = WFData.waveform;
			for (int i = 0; i < pesFound.size(); i++) {
				// TODO(josh): Improve the adjustment by averaging the shift based off of a few bins around the PE time
				unsigned int peTimeBinPos = std::floor(
						pesFound[i].time / pdfSamplingRate); // This should use a variable
				// corresponding to the input sampling rate, however for now they're the same
				float fitVal = npe_pdf_func(pesFound[i].time, params, chIdealWF);
				float extraAmplitude = fitVal - WFData.waveform[peTimeBinPos];
				float newAmplitude = pesFound[i].amplitude + extraAmplitude;
				if (newAmplitude >
				    WFSigThresh) { // TODO(josh): We need to consider the situations that this would ever be true
					params[2 + (2 * i)] = newAmplitude;
					pesFound[i].amplitude = newAmplitude;
				}
			}

			// Compute residual
//			std::vector<float> fitVecForPlot; // Debug line
			for (unsigned int k = 0; k < residualWF.waveform.size(); ++k) {
				// TODO(josh): Should it be k or k + 0.5?
				float fitVal = npe_pdf_func(float(k) * pdfSamplingRate, params, chIdealWF);
				residualWF.waveform[k] = residualWF.waveform[k] - fitVal;
//				fitVecForPlot.emplace_back(fitVal); // Debug line
			}

			/** Debugging code - writing fits to CSV for plotting **/
//			writeVector("rawWaveform.csv", WFData.waveform);
//			writeVector("fit.csv", fitVecForPlot);
//			writeVector("residual.csv", residualWF.waveform);
			/** ================================================= **/

			// Get initial guesses for the next PE
			auto minPosIt = std::min_element(residualWF.waveform.begin(), residualWF.waveform.end());
			unsigned int minTimePos = std::distance(residualWF.waveform.begin(), minPosIt);

			// If lowest point in waveform isn't below threshold there are no more PEs
			if (-residualWF.waveform[minTimePos] < WFSigThresh) {
				break;
			}

			guessPE.amplitude = -residualWF.waveform[minTimePos];
			guessPE.time = float(minTimePos) * pdfSamplingRate;


			// This is effectively checking in what direction the residual is skewed.
			// (t1*a1)/(a1*a2*a3) + (t2*a2)/(a1*a2*a3) + (t3*a3)/(a1*a2*a3)
			// If the estimated PE time is larger than truth the residual will be negative on the left,
			// and positive on the right (of the PE time), this means a1/(a1*a2*a3) will be less than one,
			// and a3/(a1+a2+a3) will be greater than one, shifting the time to the right...
			if ((minTimePos > 1) && (minTimePos < residualWF.waveform.size() - 1)) {
				// improve initial time for a new PE based on average time
				// over 3 consecutive sample ponderated by the amplitude
				// of each sample... help a lot to resolve PEs very close!
				double timeSum = 0;
				double ponderationSum = 0;
//				for (unsigned int b = minTimePos - 1; b <= minTimePos + 1; ++b) {
				double binCenter = (minTimePos - 0.5) * (pdfSamplingRate);
				double binVal = residualWF.waveform[minTimePos];
				timeSum += binCenter * binVal;
				ponderationSum += binVal;
//				}
				guessPE.time = float(PEFinderTimeOffset * 0.1) + timeSum / ponderationSum;
			}
			// Changing the multiple of PEFinerTimeOffset affects speed and chisq significantly

			numPEsFound += 1;
			pesFound.push_back(guessPE);
			residualWF = WFData;

			if (numPEsFound > maxPEs) {  // To handle the possibility of the algorithm being overly keen.
				break;
			}
		} // End of PE find loop

		using namespace ceres;
		Problem problem;

		std::vector<double> initialTimes;
		std::vector<double> initialAmplitudes;
		for (auto &k: pesFound) {
			initialAmplitudes.push_back(k.amplitude);
			initialTimes.push_back(k.time);
		}

		// Want to create a copy of the initial estimates to modify in below running correction.
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
		ceres::Grid1D grid = ceres::Grid1D<double>(chIdealWF->data(), 0, (int) chIdealWF->size());
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
		for (unsigned int j = 0; j < WFData.waveform.size(); j++) {
			xValues.push_back(
					((double) j * 100) + pdfT0SampleConv);  // Multiplying index to match position on ideal PDF
		}

		// Set up the only cost function (also known as residual). This uses
		// auto-differentiation to obtain the derivative (Jacobian).
		for (unsigned int j = 0; j < WFData.waveform.size(); ++j) {
			auto costFunction = new DynamicAutoDiffCostFunction<npe_pdf_functor, 4>(new npe_pdf_functor(xValues[j],
			                                                                                            WFData.waveform[j],
			                                                                                            idealPDFInterpolator,
			                                                                                            pesFound.size()));
			costFunction->AddParameterBlock(1); // Baseline param
			for (auto pe: pesFound) {
				costFunction->AddParameterBlock(1); // Amplitude param for PE
				costFunction->AddParameterBlock(1); // Time param for PE
			}
			costFunction->SetNumResiduals(1);
			problem.AddResidualBlock(costFunction, nullptr, params);
		}

		// Run the solver!
		Solver::Options options;
		options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
		options.parameter_tolerance = 1e-5; // default is 1e-8, check if this is tolerance for any or all params
		options.minimizer_progress_to_stdout = false;
		Solver::Summary summary;
		Solve(options, &problem, &summary);

//        std::cout << summary.FullReport() << "\n";


		// Going back from ideal waveform PDF index to time
		for (double &time: times) {
			time = time / samplingRate2Inv;
		}

		std::vector<float> finalParams;
		finalParams.push_back(pesFound.size());
		finalParams.push_back(baseline);
		for (int k = 0; k < pesFound.size(); k++) {
			finalParams.push_back((float) amplitudes[k]);
			finalParams.push_back((float) times[k]);
		}

		// TODO(josh): ChiSq calculation needs improvement/validating.
		double chiSq = 0;
		for (int j = 0; j < WFData.waveform.size(); j++) {
			double observed = WFData.waveform[j];
			double expected = npe_pdf_func(float(j) * pdfSamplingRate, finalParams, chIdealWF);
			chiSq += std::pow(observed - expected, 2) /
			         (pdfResidualRMS / 1000); //TODO(josh): Need to calculate the residual RMS on a per waveform basis?
		}
		double redChiSq = chiSq / (finalParams.size() - 1);
		sysProcWFCount++;
		meanReducedChisq = meanReducedChisq + (redChiSq - meanReducedChisq) / sysProcWFCount;

		delete idealPDFInterpolator;

		/** Debugging code - writing fits to CSV for plotting **/
//		std::vector<float> fullFitVecForPlot;
//		for (unsigned int k = 0; k < WFData.waveform.size(); k++) {
//			auto thing = npe_pdf_func(k * pdfSamplingRate, finalParams, chIdealWF);
//			fullFitVecForPlot.emplace_back(thing);
//		}
//		writeVector("fullFit.csv", fullFitVecForPlot);
		/** ================================================= **/

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
			FitPEs.push_back(pe);
			sysProcPECount++;
			ampDiff = ampDiff + ((pe.amplitude - initialAmplitudes[k]) - ampDiff) / sysProcPECount;
			timeDiff = timeDiff + ((pe.time - initialTimes[k]) - timeDiff) / sysProcPECount;
			baselineDiff = baselineDiff + ((baseline - initBaseline) - baselineDiff) / sysProcPECount;
		}

		chFit.pes = FitPEs;
		chFit.baseline = float(baseline);
		evFitDat.SiPM.push_back(chFit);
	}
	Writer writer(std::move(file));
	writer.writeEventInfo(evFitDat);
}

/**
 *
 * @param events
 * @param count
 * @param m
 * @param idealWaveforms
 * @param file
 * @return
 */
bool fitBatchPEs(const std::vector<EventData> &events, std::atomic<unsigned long> &count, std::mutex &m,
                 const std::vector<std::vector<double>> *idealWaveforms, const std::shared_ptr<SyncFile> &file) {
	for (const auto &event: events) {
		m.lock();
		++count;
		m.unlock();
		fitPE(&event, idealWaveforms, file);
	}
	return true;
}