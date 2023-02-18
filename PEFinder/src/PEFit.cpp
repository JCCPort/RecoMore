#include "../include/PEFit.h"
#include "../include/Utils.h"
#include <mutex>
#include <atomic>
#include <cmath>
#include "ceres/ceres.h"
#include "ceres/cubic_interpolation.h"
#include <memory>


struct NPEPDFFunctor {
	NPEPDFFunctor(std::vector<float> waveformTimes, std::vector<float> waveformAmplitudes,
	              ceres::CubicInterpolator<ceres::Grid1D<float, true>> *PDFInterpolator, unsigned int numPES)
			: waveformTimes_(std::move(waveformTimes)), waveformAmplitudes_(std::move(waveformAmplitudes)), PDFInterpolator_(PDFInterpolator), numPES_(numPES) {};
	
	template<typename T>
	inline bool operator()(T const *const *params, T *__restrict__ residual) const {
		for (int j = 0; j < waveformTimes_.size(); j++) {
			T f;
			T X_(waveformTimes_[j]);
			residual[j] = params[0][0];
			for (unsigned int i = 0; i < numPES_; ++i) {
				auto peParams = params[i];
				PDFInterpolator_->Evaluate((X_ - peParams[0]), &f);
				residual[j] += (peParams[1] * f);
			}
			residual[j] -= waveformAmplitudes_[j];
		}
		return true;
	}

private:
	const std::vector<float>                        waveformTimes_;
	const std::vector<float>                        waveformAmplitudes_;
	ceres::CubicInterpolator<ceres::Grid1D<float> > *PDFInterpolator_;
	const unsigned int                              numPES_;
};


inline float NPEPDFFunc(float X, const std::vector<float> &p, const std::vector<double> *idealWaveform) {
	// TODO(josh): Replace p with an instance of fit parameters type
	
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
	
	int   NPE      = (int) (p[0]);
	float BASELINE = p[1];
	
	float value = BASELINE;
	
	// This is adding up the contribution of each fit PE to a specific bin
	for (int PE = 0; PE < NPE; ++PE) {
		float PE_CHARGE  = p[2 + PE * 2];
		float PE_TIME    = p[3 + PE * 2];
		// TODO(josh): Start using interpolation for this too instead of using floor
		int   PE_PDF_BIN = pdfT0Sample + (int) std::floor(0.5 + (X - PE_TIME) * samplingRate2Inv);
		if ((PE_PDF_BIN >= 0) && (PE_PDF_BIN < pdfNSamples)) {
			value += PE_CHARGE * (float)idealWaveform->at(PE_PDF_BIN);
		}
	}
	return value;
}


inline void updateGuessCorrector(const std::vector<double>& amps, const std::vector<double>& times,
                                 const std::vector<double>& initialAmps, const std::vector<double>& initialTimes,
                                 double baseline, double initBaseline, const std::vector<Photoelectron>& pesFound){
	for (int k   = 0; k < pesFound.size(); k++) {
		sysProcPECount++;
		ampDiff  = ampDiff + ((amps[k] - initialAmps[k]) - ampDiff) / sysProcPECount;
		timeDiff = timeDiff + ((times[k] - initialTimes[k]) - timeDiff) / sysProcPECount;
	}
	baselineDiff = baselineDiff + ((baseline - initBaseline) - baselineDiff) / sysProcPECount;
}


inline bool getNextPEGuess(DigitiserChannel* residualWF, Photoelectron *guessPE){
	// Get initial guesses for the next PE
	const auto         minPosIt   = std::min_element(residualWF->waveform.begin(), residualWF->waveform.end());
	const unsigned int minTimePos = std::distance(residualWF->waveform.begin(), minPosIt);
	
	// If lowest point in waveform isn't below threshold there are no more PEs
	if (-residualWF->waveform[minTimePos] < WFSigThresh) {
		return false;
	}
	
	guessPE->amplitude = -residualWF->waveform[minTimePos];
	guessPE->time      = float(minTimePos) * pdfSamplingRate;
	return true;
}


inline void amplitudeCorrection(std::vector<Photoelectron> *pesFound, std::vector<float> *params, std::vector<float> waveform, const std::vector<double> *chIdealWF){
	for (int i = 0; i < pesFound->size(); i++) {
		// TODO(josh): Improve the adjustment by averaging the shift based off of a few bins around the PE time
		const unsigned int peTimeBinPos   = std::floor(pesFound->at(i).time / pdfSamplingRate);
		const float        fitVal         = NPEPDFFunc(pesFound->at(i).time, *params, chIdealWF);
		const float        extraAmplitude = fitVal - waveform[peTimeBinPos];
		const float        newAmplitude   = pesFound->at(i).amplitude + extraAmplitude;
		if (newAmplitude > WFSigThresh) { // TODO(josh): We need to consider the situations that this would ever be true
			params->at(2 + (2 * i)) = newAmplitude;
			pesFound->at(i).amplitude = newAmplitude;
		}
	}
}


void
fitEvent(const DigitiserEvent *event, const std::vector<std::vector<double>> *idealWaveforms, std::shared_ptr<SyncFile> outputFile, std::mutex &lock) {
	std::vector<FitChannel> chFits;
	std::cout << "Fitting to event" << std::endl;
	for (const auto &channel: event->channels) { // Looping through all channels for a given event
		auto residualWF = channel; // This will be the variable that is modified to be the residual distribution after each iteration
		
		FitChannel         chFit{};
		const unsigned int ch = channel.ID;
		if (std::count(skipChannels.begin(), skipChannels.end(), ch)) {
			continue;
		}
		chFit.ID = ch;
		
		// Making a pointer to the ideal waveform for this channel to improve speed of passing.
		const std::vector<double> tmp        = (*idealWaveforms)[ch];
		const std::vector<double> *chIdealWF = &tmp;
		
		// Baseline calculation
		float initBaseline = averageVector(channel.waveform, 2);
		chFit.baseline = initBaseline;

		// Start loop that will break when no more PEs are present
		std::vector<Photoelectron> pesFound;
		unsigned int               numPEsFound = 0;
		
		// This is getting estimate PEs that will then be passed as initial guesses to the minimiser.
		Photoelectron guessPE{};

		// Initial baseline removal.
		for (float &k: residualWF.waveform) {
			k = k - initBaseline;
		}
		
		while (true) {
			
			std::sort(pesFound.begin(), pesFound.end(), comparePETime()); // Need to sort for amplitude adjustment.
			
			std::vector<float> params;
			params.push_back((float) pesFound.size());
			params.push_back(initBaseline);
			for (const auto &pe: pesFound) {
				params.push_back(pe.amplitude);
				params.push_back(pe.time);
			}
			
			// Amplitude adjustment: if the latest PE found is before other one(s),
			//  its tail is going to add some amplitude to the following one. Compares
			//  real and fit amplitude at the time bin corresponding to the PE time.
			amplitudeCorrection(&pesFound, &params, channel.waveform, chIdealWF);
			
			// Compute residual
			for (unsigned int  k = 0; k < residualWF.waveform.size(); ++k) {
				// TODO(josh): Should it be k or k + 0.5?
				const float fitVal = NPEPDFFunc((float)k * pdfSamplingRate, params, chIdealWF);
				residualWF.waveform[k] = residualWF.waveform[k] - fitVal + initBaseline;
			}

			// Keep correcting baseline as new PEs are found.
			if(!pesFound.empty()){
				initBaseline = averageVector(residualWF.waveform, 2);
				for (float &k: residualWF.waveform) {
					k = k - initBaseline;
				}
			}
			
			if(!getNextPEGuess(&residualWF, &guessPE)){
				break;
			}
			
			numPEsFound += 1;
			pesFound.push_back(guessPE);
			residualWF = channel;
            std::cout << "Initial guess loop" << std::endl;
			
			if (numPEsFound > maxPEs) {  // To handle the possibility of the algorithm being overly keen.
				break;
			}
		} // End of PE find loop
		
		if(numPEsFound == 0){
			continue;
		}
		
		ceres::Problem problem{};
		
		std::vector<double> initialTimes;
		std::vector<double> initialAmplitudes;
		for (auto &k: pesFound) {
			initialAmplitudes.push_back(k.amplitude);
			initialTimes.push_back(k.time);
		}
		
		// Want to create a copy of the initial estimates to modify in below running correction.
		std::vector<double> times      = initialTimes;
		std::vector<double> amplitudes = initialAmplitudes;
		for (int i = 0; i < times.size(); i++) {
			times[i] += timeDiff;
			amplitudes[i] += ampDiff;
		}
		
		// Converting the time into an ideal waveform PDF index to simplify NPEPDFFunctor method call
		for (double &time: times) {
			time = time * samplingRate2Inv;
		}
		
		// Creating interpolator to allow for the ideal waveform to be used as a continuous (and differentiable) function.
		std::vector<float> temp(chIdealWF->begin(), chIdealWF->end());
		ceres::Grid1D      grid                 = ceres::Grid1D<float>(temp.data(), 0, (int) chIdealWF->size());
		auto               idealPDFInterpolator = new ceres::CubicInterpolator<ceres::Grid1D<float> >(grid);
		
		
		// Creating vector of references to parameter values that the fitter will use and modify. Note this means that
		// the references are the initial values before the fit and the final values after the fit.
		double              baseline = initBaseline;
		std::vector<double> params   = {};
		params.push_back(baseline);
		for (int i = 0; i < pesFound.size(); i++) {
			params.push_back(amplitudes[i]);
			params.push_back(times[i]);
		}
		
		// Creating the x values that the solver will use
		std::vector<float> xValues;
		for (unsigned int  j = 0; j < channel.waveform.size(); j++) {
			xValues.push_back(((float)j * 100.0f) + pdfT0SampleConv);  // Multiplying index to match position on ideal PDF
		}
		
		// Set up the only cost function (also known as residual). This uses
		// auto-differentiation to obtain the derivative (Jacobian).
		auto functor = new NPEPDFFunctor(xValues,
		                                 channel.waveform,
		                                 idealPDFInterpolator,
		                                 pesFound.size());
		auto costFunction = new ceres::DynamicAutoDiffCostFunction<NPEPDFFunctor>(functor);
		
		costFunction->SetNumResiduals((int) channel.waveform.size());
		
		costFunction->AddParameterBlock(1); // Baseline param
		for ([[maybe_unused]]const auto &pe: pesFound) {
			costFunction->AddParameterBlock(2); // Params for one PE
		}

        std::cout << "About to make parameter blocks" << std::endl;

		// Formatting parameters to allow grouping of params for one PE
		std::vector<double *> parameterBlocks;
		double                x1[] = {params[0]};
		parameterBlocks.push_back(x1);
		double x2[pesFound.size()][2];
		
		for (int i = 0; i < pesFound.size(); i++) {
			x2[i][0] = params[(2 * i) + 1];
			x2[i][1] = params[(2 * i) + 2];
			parameterBlocks.push_back(x2[i]);
		}

        std::cout << "Made parameter blocks" << std::endl;

		auto lossFunction(new ceres::ArctanLoss(WFSigThresh));

        std::cout << "Made loss function" << std::endl;
		
		problem.AddResidualBlock(costFunction, lossFunction, parameterBlocks);
		
		for (int i = 1; i < pesFound.size() + 1; i++) {
			problem.SetParameterLowerBound(parameterBlocks[i], 0, 0);
		}
		
		// Run the solver!
		ceres::Solver::Options options;
//		options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
		options.linear_solver_type           = ceres::DENSE_QR;
		options.parameter_tolerance          = 1e-5; // default is 1e-8, check if this is tolerance for any or all params
		options.minimizer_progress_to_stdout = false;
		ceres::Solver::Summary summary;
		Solve(options, &problem, &summary);

        std::cout << summary.FullReport() << "\n";
		
		// Going back from ideal waveform PDF index to time
		for (double &time: times) {
			time = time / samplingRate2Inv;
		}
		
		std::vector<Photoelectron> FitPEs;
		for (int            k = 0; k < pesFound.size(); k++) {
			Photoelectron pe{};
			pe.amplitude = float(amplitudes[k]);
			pe.time      = float(times[k]);
			FitPEs.push_back(pe);
		}
		chFit.PEs      = FitPEs;
		chFit.baseline = float(baseline);
		
		
		std::vector<float> finalParams;
		finalParams.push_back((float) FitPEs.size());
		finalParams.push_back((float) baseline);
		for (const auto &PE: FitPEs) {
			finalParams.push_back((float) PE.amplitude);
			finalParams.push_back((float) PE.time);
		}

		float            chiSq    = 0;
		for (unsigned int j        = 0; j < channel.waveform.size(); j++) {
			const float observed = channel.waveform[j];
			const float expected = NPEPDFFunc((float) j * pdfSamplingRate, finalParams, chIdealWF);
			chiSq += (float)std::pow(observed - expected, 2) / (pdfResidualRMS / 1000);
			//TODO(josh): Need to calculate the residual RMS on a per waveform basis?
		}
		float            redChiSq = chiSq / ((float) channel.waveform.size() - ((float) finalParams.size() - 1));
		
		sysProcWFCount++;
		meanReducedChiSq = meanReducedChiSq + (redChiSq - meanReducedChiSq) / sysProcWFCount;
		
		chFit.reducedChiSq = redChiSq;
		
		lock.lock();
		reducedChiSqs.emplace_back(redChiSq);
		lock.unlock();
		
		chFits.push_back(chFit);
		
		// Update correction values for initial guesses
		updateGuessCorrector(amplitudes, times, initialAmplitudes, initialTimes, baseline, initBaseline, pesFound);
		
		delete idealPDFInterpolator;
	}
	
	FitEvent evFitDat{event->ID, event->correctedTime, event->date, chFits};
	
	Writer writer(std::move(outputFile));
	writer.writeEventInfo(evFitDat);
}

/**
 *
 * @param events
 * @param count
 * @param lock
 * @param idealWaveforms
 * @param file
 * @return
 */
bool batchFitEvents(const std::vector<DigitiserEvent> &events, std::atomic<unsigned long> &count, std::mutex &lock,
                    const std::vector<std::vector<double>> *idealWaveforms, const std::shared_ptr<SyncFile> &file) {
	for (const auto &event: events) {
		lock.lock();
		++count;
		lock.unlock();
		fitEvent(&event, idealWaveforms, file, lock);
	}
	return true;
}