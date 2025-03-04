#include "../include/PEFit.h"
#include "../include/Utils.h"
#include "../include/PETemplate.h"
#include "../Globals.h"
#include <mutex>
#include <atomic>
#include <cmath>
#include "ceres/ceres.h"
#include "ceres/cubic_interpolation.h"
#include <memory>
#include <utility>


/**
 * @class TemplateFunctor
 * @brief A functor (callable) for computing residuals of a model consisting of multiple photoelectrons plus a baseline.
 *
 * TemplateFunctor uses a single-photoelectron (PE) template (represented by a ceres::CubicInterpolator) and sums the
 * contributions from multiple PEs at each sample point in the waveform. The difference between the summed model
 * and the actual waveform value is returned as the residual for fitting in Ceres.
 */
struct TemplateFunctor {
	/**
	 * @brief Constructor.
	 *
	 * @param waveformTimes         Indices (on the PETemplate) corresponding to the samples of the waveform.
	 * @param waveformAmplitudes    The observed waveform amplitudes.
	 * @param templateInterpolator       A pointer to a ceres::CubicInterpolator for the single PE template.
	 * @param numTemplateSamples         Number of samples in the single PE template.
	 * @param numPES                Number of photoelectrons (PEs) to consider in the model.
	 * @param dataSampleSpacing     The time spacing between samples in the observed waveform.
	 * @param templateSampleSpacing  The time spacing between samples in the single PE template.
	 * @param templateMaxIntensityIdx The index of the maximum intensity in the single PE template.
	 */
	TemplateFunctor(std::vector<float> waveformTimes, std::vector<float> waveformAmplitudes, ceres::CubicInterpolator<ceres::Grid1D<float, true>>* templateInterpolator,
		const unsigned int numTemplateSamples, const unsigned int numPES, const float dataSampleSpacing, const float templateSampleSpacing, const unsigned int templateMaxIntensityIdx)
		: waveformTimes_(std::move(waveformTimes)), waveformAmplitudes_(std::move(waveformAmplitudes)), templateInterpolator_(templateInterpolator), numTemplateSamples_(numTemplateSamples), numPES_(numPES),
		dataSampleSpacing_(dataSampleSpacing), templateSampleSpacing_(templateSampleSpacing), templateMaxIntensityIdx_(templateMaxIntensityIdx) {
		waveformSize_ = waveformTimes_.size();
		sampleSpacingRatio_ = dataSampleSpacing_ / templateSampleSpacing_;

		lowerOffset_ = (templateMaxIntensityIdx_ + 1) / sampleSpacingRatio_;
		upperOffset_ = (numTemplateSamples_ - (templateMaxIntensityIdx_ + 1))/ sampleSpacingRatio_;
	}

	/**
	 * @brief The operator that computes the residuals for Ceres.
	 *
	 * @tparam T       The scalar type for autodifferentiation.
	 * @param params   Array of parameter blocks. The first element contains the baseline;
	 *                 each subsequent pair contains [amplitude, time] for one PE.
	 * @param residual The array of residuals (difference between model and observed data).
	 * @return True if successful, otherwise false.
	 */
	template <typename T>
	inline bool operator()(T const* const* params, T* __restrict__ residual) const {
		const T baseline = params[0][0];

		// 1) Initialize residual with (baseline - observed) for all j
		for (unsigned int j = 0; j < waveformSize_; j++) {
			residual[j] = baseline - T(waveformAmplitudes_[j]);
		}

		// 2) Accumulate each PE’s contribution
		for (unsigned int i = 0; i < numPES_; i++) {
			const unsigned int i2 = 2 * i;
			const T amplitude     = params[i2 + 1][0];
			const T time          = params[i2 + 2][0];

			// Precompute the wave time window where this PE is valid
			// we want 0 <= (X - time) < numTemplateSamples_
			// => time <= X < time + numTemplateSamples_

			// Extract scalar value from time: if T is a Jet, use time.a, else time itself, divide by sampleSpacingRatio_ here to reduce the number of divisions.
			double ratioAdjustedTime;
			if constexpr (std::is_same_v<T, double>) {
				ratioAdjustedTime = time/sampleSpacingRatio_;
			} else {
				ratioAdjustedTime = time.a/sampleSpacingRatio_;
			}

			// Compute indices as ints
			const unsigned int j1 = std::max(0, static_cast<int>(std::ceil(ratioAdjustedTime - lowerOffset_)));
			const unsigned int j2 = std::min(waveformSize_, static_cast<int>(std::ceil(ratioAdjustedTime + upperOffset_)));

			// 3) Loop only over valid j’s for this PE
			for (unsigned int j = j1; j < j2; j++) {
				const T pos = T(waveformTimes_[j]) - time;
				// pos should now be guaranteed to be in [0, numTemplateSamples_)
				T f;
				templateInterpolator_->Evaluate(pos, &f);
				residual[j] += amplitude * f;
			}
		}

		return true;
	}


	/**
	 * @brief Update the number of photoelectrons in the model.
	 * @param numPEs  The new number of PEs.
	 */
	void updateNumPEs(const unsigned int numPEs) {
		if (numPEs != numPES_) {
			numPES_ = numPEs;
		}
	}

private:
	const std::vector<float>                              waveformTimes_;      ///< Index positions on the template.
	const std::vector<float>                              waveformAmplitudes_; ///< Measured waveform values.
	ceres::CubicInterpolator<ceres::Grid1D<float, true>>* templateInterpolator_; ///< Interpolator for the single PE template.
	unsigned int numTemplateSamples_; ///< Number of samples in the single PE template.
	unsigned int numPES_; ///< Number of photoelectrons in the model.
	const float dataSampleSpacing_;
	const float templateSampleSpacing_;
	const unsigned int templateMaxIntensityIdx_;

	int waveformSize_;
	double sampleSpacingRatio_;
	double lowerOffset_;
	double upperOffset_;
};


/**
 * @brief Compute the photoelectron current fit value at a given index @p X using cubic interpolation.
 *
 * @param X                The index (in template sample space).
 * @param fitParams        A collection of fit parameters (baseline, amplitudes, times).
 * @param templateInterpolator  Pointer to a ceres::CubicInterpolator for the single PE template.
 * @param ChPETemplate     Pointer to the PETemplate object describing the single PE waveform.
 * @return                 The computed model value at index @p X.
 */
inline float currentFitAmp(const float X, const FitParams& fitParams, const ceres::CubicInterpolator<ceres::Grid1D<float, true>>* templateInterpolator, const PETemplate* ChPETemplate) {
	float value = fitParams.baseline;

	for (int PE = 0; PE < fitParams.numPEs; ++PE) {
		const float amplitude = fitParams.PEs[PE].amplitude;
		const float time      = fitParams.PEs[PE].time / ChPETemplate->getTimeSpacing();

		// Now we want to interpolate if we're in range. pos is the position on the PE template waveform relative to the PE time in terms of index position
		if (const float pos = X - time; pos >= 0.0f && pos < static_cast<float>(ChPETemplate->getNumberOfSamples())) {
			// Call the double overload
			auto   posDouble         = static_cast<double>(pos);
			double templateValDouble = 0.0;
			templateInterpolator->Evaluate(posDouble, &templateValDouble);

			// Accumulate
			const auto templateVal = static_cast<float>(templateValDouble);
			value += amplitude * templateVal;
		}
	}

	return value;
}


/**
 * @brief Update global correction factors for guesses (amplitude, time, baseline) based on a post-fit result.
 *
 * This function updates global running averages (ampDiff, timeDiff, baselineDiff) for future initial guesses.
 *
 * @param amps          The fitted PE amplitudes.
 * @param times         The fitted PE times.
 * @param initialAmps   The initial guessed PE amplitudes.
 * @param initialTimes  The initial guessed PE times.
 * @param baseline      The fitted baseline.
 * @param initBaseline  The initial guessed baseline.
 * @param numPEs        Number of PEs in the current fit.
 */
inline void updateGuessCorrector(const std::vector<double>& amps,
                                 const std::vector<double>& times,
                                 const std::vector<double>& initialAmps,
                                 const std::vector<double>& initialTimes,
                                 const double               baseline,
                                 const double               initBaseline,
                                 const unsigned int         numPEs) {
	for (int k = 0; k < numPEs; k++) {
		sysProcPECount++;
		ampDiff  = ampDiff + ((amps[k] - initialAmps[k]) - ampDiff) / sysProcPECount;
		timeDiff = timeDiff + ((times[k] - initialTimes[k]) - timeDiff) / sysProcPECount;
		// TODO(josh): Keep the below but only for debug builds
#ifdef IS_DEBUG
		if (ampDiff != ampDiff) {
			throw std::runtime_error("AmpDiff is NaN");
		}
		if (timeDiff != timeDiff) {
			throw std::runtime_error("TimeDiff is NaN");
		}
#endif
	}
	baselineDiff = baselineDiff + ((baseline - initBaseline) - baselineDiff) / sysProcPECount;
}


/**
 * @brief Locate the next candidate photoelectron from the current residual waveform.
 *
 * This function finds the most negative point in the residual waveform, and if it is below a threshold
 * (i.e., above a certain absolute amplitude), it is deemed a new PE candidate.
 * Note that positive pulses are handled previously by multiplying data and templates by -1.
 *
 * @param residualWF   Pointer to the digitizer channel representing the current residual waveform.
 * @param guessPE      Output parameter to store the guessed PE amplitude/time.
 * @param sampleSpacing
 * @return             True if a PE guess was found, false otherwise.
 */
inline bool getNextPEGuess(const std::vector<float>& residualWF, Photoelectron* guessPE, const float sampleSpacing) {
	// Get initial guesses for the next PE
	const auto         minPosIt   = std::ranges::min_element(residualWF);
	const unsigned int minTimePos = std::distance(residualWF.begin(), minPosIt);

	// If lowest point in waveform isn't below threshold there are no more PEs
	if (-residualWF[minTimePos] < WFSigThresh) {
		return false;
	}

	guessPE->amplitude = -residualWF[minTimePos];
	guessPE->time      = static_cast<float>(minTimePos) * sampleSpacing;

	return true;
}


/**
 * @brief Adjust the amplitudes of the currently guessed PEs by comparing the current fit to the observed waveform.
 *
 * For each PE, it calculates the current fit value at its time and compares with the observed waveform amplitude,
 * and applies a shift if necessary.
 *
 * @param fitParams         The current set of fit parameters (including guessed PEs).
 * @param waveform          The observed waveform.
 * @param templateInterpolator   Pointer to the single PE template interpolator.
 * @param ChPETemplate      The single PE template.
 * @param sampleSpacing
 */
inline void amplitudeCorrection(FitParams fitParams, const std::vector<float>& waveform, const ceres::CubicInterpolator<ceres::Grid1D<float>>* templateInterpolator, const PETemplate* ChPETemplate, const float sampleSpacing) {
	for (int i = 0; i < fitParams.numPEs; ++i) {
		// TODO(josh): Improve the adjustment by averaging the shift based off of a few bins around the PE time
		const unsigned int peTimeBinPos = std::floor(fitParams.PEs[i].time / sampleSpacing);

		const float fitVal = currentFitAmp(ChPETemplate->getFractionalIndex(fitParams.PEs[i].time),
		                                   // time in index position
		                                   fitParams,            // vector of parameters
		                                   templateInterpolator, // your cubic interpolator
		                                   ChPETemplate);
		const float extraAmplitude = fitVal - waveform[peTimeBinPos];
		if (const float newAmplitude = fitParams.PEs[i].amplitude + extraAmplitude; newAmplitude > WFSigThresh) {
			// TODO(josh): We need to consider the situations that this would ever be true
			fitParams.PEs[i].amplitude = newAmplitude;
		}
	}
}


/**
 * @brief Perform the fitting of a single event, updating channel-by-channel, writing results to output.
 *
 * The procedure is roughly:
 *   1. Estimate baseline.
 *   2. Iteratively find possible PEs and refine their amplitudes/times by subtracting the model from the residual.
 *   3. Use Ceres solver to refine the final guess.
 *   4. Compute goodness-of-fit and store results.
 *
 * @param event         Pointer to the DigitiserEvent to fit.
 * @param PETemplates   A map of channel ID to PETemplate, used to get the single PE template for each channel.
 * @param outputFile    Shared pointer to the output file structure (thread-safe).
 * @param lock          Mutex for thread-safety in file writing.
 * @param sampleSpacing
 */
void fitEvent(const DigitiserEvent* event, const std::unordered_map<unsigned int, PETemplate>& PETemplates, std::shared_ptr<SyncFile> outputFile, std::mutex& lock, const float sampleSpacing) {
	std::vector<FitChannel> chFits;
	for (const auto& channel : event->channels) { // Looping through all channels for a given event

		// This will be the variable that is modified to be the residual distribution after each iteration
		auto residualWF = channel.waveform;

		FitChannel         chFit{};
		const unsigned int ch = channel.ID;
		chFit.ID = ch;

		constexpr float    preSigWindowTime = 18.f; // ns. This is the time window before the signal that is used to calculate the noise level.
		std::vector<float> preSignalWF(residualWF.begin(), residualWF.begin() + static_cast<int>(preSigWindowTime / sampleSpacing));
		//        auto stdDevNoise = calculateVariance(preSignalWF, calculateMean(preSignalWF));
		auto stdDevNoise    = calculateStandardDeviation(preSignalWF);
		templateResidualRMS = static_cast<float>(stdDevNoise);

		// Making a pointer to the ideal waveform for this channel to improve speed of passing.
		if (ch > PETemplates.size() - 1) {
			throw std::runtime_error("Channel number exceeds the number of ideal waveforms");
		}

		const PETemplate*          ChPETemplate = &PETemplates.at(ch);
		const std::vector<double>* chIdealWF    = ChPETemplate->getVoltagesRef();

		// Baseline calculation
		float initBaseline = averageVector(residualWF, 0.01);
		chFit.baseline     = initBaseline;

		// This is getting estimate PEs that will then be passed as initial guesses to the minimiser.
		Photoelectron guessPE{};

		// Creating the x values that the solver will use. These are the index positions on the ideal WF for the positions on the real WF.
		std::vector<float> xValues;
		for (unsigned int j = 0; j < channel.waveform.size(); j++) {
			xValues.push_back(ChPETemplate->getFractionalIndex(static_cast<float>(j) * sampleSpacing) + 1.f);
		}

		// Creating interpolator to allow for the ideal waveform to be used as a continuous (and differentiable) function.
		std::vector<float> temp(chIdealWF->begin(), chIdealWF->end());
		auto               grid                 = ceres::Grid1D<float>(temp.data(), 0, static_cast<int>(chIdealWF->size()));
		auto               templateInterpolator = new ceres::CubicInterpolator<ceres::Grid1D<float>>(grid);

		// Set up the only cost function (also known as residual). This uses auto-differentiation to obtain the derivative (Jacobian).
		// Creating them here to be able to re-use them.
		auto functor = new TemplateFunctor(xValues, channel.waveform, templateInterpolator, ChPETemplate->getNumberOfSamples(), 0, /* will set numPEs later */
		                                  sampleSpacing, ChPETemplate->getTimeSpacing(), ChPETemplate->getTMaxIntensity());
		auto costFunction = new ceres::DynamicAutoDiffCostFunction<TemplateFunctor, 16>(functor, ceres::Ownership::DO_NOT_TAKE_OWNERSHIP);

		// Remove initial baseline
		for (float& k : residualWF) {
			k -= initBaseline;
		}

		// Container for fit results
		FitParams fitParams;

		unsigned int numPEsFound = 0;

		// Iteratively find PEs until none can be found, or we exceed a max
		while (true) {
			fitParams.sortPEsByTime();
			fitParams.baseline = initBaseline;

			// Amplitude adjustment: if the latest PE found is before other one(s),
			//  its tail is going to add some amplitude to the following one. Compares
			//  real and fit amplitude at the time bin corresponding to the PE time.
			amplitudeCorrection(fitParams, channel.waveform, templateInterpolator, ChPETemplate, sampleSpacing);

			// Update residual waveform
			for (unsigned int k = 0; k < xValues.size(); ++k) {
				float fitVal = currentFitAmp(xValues[k],           // time in index position
				                             fitParams,            // vector of parameters
				                             templateInterpolator, // your cubic interpolator
				                             ChPETemplate);
				residualWF[k] = residualWF[k] - fitVal + initBaseline;
			}

			// Recompute baseline if PEs exist
			if (fitParams.numPEs != 0) {
				initBaseline = averageVector(residualWF, 0.01);
				for (float& k : residualWF) {
					k = k - initBaseline;
				}
			}

			if (!getNextPEGuess(residualWF, &guessPE, sampleSpacing)) {
				break;
			}
			numPEsFound += 1;
			fitParams.addPE(guessPE);

			// Reset residual to original waveform.
			residualWF = channel.waveform;

			if (numPEsFound >= maxPEs) {
				// To handle the possibility of the algorithm being overly keen.
				break;
			}
		} // End of PE find loop

		// Carry out one final adjustment to the amplitude of the PEs found and recalculate baseline, residual.
		amplitudeCorrection(fitParams, channel.waveform, templateInterpolator, ChPETemplate, sampleSpacing);

		initBaseline = averageVector(residualWF, 0.01);
		// =========================================

		// TODO(JOSH): ADD A CHECKER TO SEE IF FIT IS GOOD ENOUGH AT THIS POINT AND THEN SKIP USING CERES IF IT IS.

		if (numPEsFound == 0) {
			continue;
		}

		// Want to create a copy of the initial estimates for later improvement of correction factors.
		std::vector<double> initialTimes;
		std::vector<double> initialAmplitudes;
		for (auto& k : fitParams.PEs) {
			initialAmplitudes.push_back(k.amplitude);
			initialTimes.push_back(k.time);
		}

		// Creating vector of references to parameter values that the fitter will use and modify. Note this means that
		// the references are the initial values before the fit and the final values after the fit.
		std::vector<double>  times;
		std::vector<double>  amplitudes;
		double               baseline;
		std::vector<double*> params = {};
		fitParams.makeSolverParams(&params, &times, &amplitudes, &baseline);

		// Make adjustment to initial guesses with running correction factors.
		for (int i = 0; i < times.size(); i++) {
			times[i] += timeDiff;
			amplitudes[i] += ampDiff;
		}

		// Convert times to template index space
		for (double& time : times) {
			time = time / ChPETemplate->getTimeSpacing();
		}

		costFunction->SetNumResiduals(static_cast<int>(channel.waveform.size()));

		costFunction->AddParameterBlock(1); // Baseline param
		for ([[maybe_unused]] const auto& pe : fitParams.PEs) {
			costFunction->AddParameterBlock(1); // Params for one PE amplitude
			costFunction->AddParameterBlock(1); // Params for one PE time
		}

		functor->updateNumPEs(fitParams.numPEs);

		// Configure Ceres problem
		ceres::Problem::Options option;
		option.loss_function_ownership = ceres::Ownership::DO_NOT_TAKE_OWNERSHIP;
		// Necessary to allow reuse whilst preventing double deletion.
		option.cost_function_ownership = ceres::Ownership::DO_NOT_TAKE_OWNERSHIP;
		option.manifold_ownership      = ceres::Ownership::DO_NOT_TAKE_OWNERSHIP;
		ceres::Problem problem(option);
		problem.AddResidualBlock(costFunction, nullptr, params);

		// Lower bound parameters except for baseline
		for (unsigned int i = 1; i < params.size(); i++) {
			problem.SetParameterLowerBound(params[i], 0, 0);
		}

		// Run the solver!
		if (std::abs(parameterTolerance) >= 1e-10f) {
			ceres::Solver::Options options;
			options.minimizer_progress_to_stdout = false;
			options.parameter_tolerance          = parameterTolerance;

			// (Optional) Parameter block ordering, improves performance.
			auto ordering = std::make_unique<ceres::ParameterBlockOrdering>();
			ordering->AddElementToGroup(params[0], 0); // baseline
			for (int i = 0; i < fitParams.numPEs; i++) {
				ordering->AddElementToGroup(params[2 * i + 1], i); // amplitude
				ordering->AddElementToGroup(params[2 * i + 2], i); // time
			}
			options.linear_solver_ordering = std::move(ordering);

			ceres::Solver::Summary summary;
			Solve(options, &problem, &summary);

			//        std::cout << summary.FullReport() << "\n";
		}

		// Updating the amplitudes and times using params which is a vector of pointers to doubles
		//		std::cout << "\n\n=====================" << std::endl;
		//		std::cout << "Fit values" << std::endl;
		//		std::cout << "Baseline : " << (postfitParams[0]) << std::endl;
		//		baseline = (postfitParams[0]);
		//		for (int i = 0; i < pesFound.size(); i++) {
		//			std::cout << "Amplitude " << i << " : " << (postfitParams[2*i + 1]) << std::endl;
		//			std::cout << "Time " << i << " : " << (postfitParams[2*i + 2]) << std::endl;
		//			amplitudes[i] = (postfitParams[2*i + 1]);
		//			times[i] = (postfitParams[2*i + 2]);
		//		}
		//		std::cout << "=====================" << std::endl;

		// Convert times back from index space
		for (double & time : times) {
			time = time * ChPETemplate->getTimeSpacing();
		}

		// Build final FitParams object
		FitParams finalFitParams{fitParams.numPEs, *params[0], amplitudes, times};
		chFit.PEs      = finalFitParams.PEs;
		chFit.baseline = finalFitParams.baseline;

		// Compute chi-squared
		float               chiSq = 0.f;
		std::vector<double> observedValues;
		std::vector<double> predictedValues;
		observedValues.reserve(channel.waveform.size());
		predictedValues.reserve(channel.waveform.size());

		for (unsigned int j = 0; j < channel.waveform.size(); j++) {
			const float observed = channel.waveform[j];
			const float expected = currentFitAmp(xValues[j], finalFitParams, templateInterpolator, ChPETemplate);
			observedValues.push_back(observed);
			predictedValues.push_back(expected);

			const float diff = observed - expected;
			chiSq += (diff * diff) / (templateResidualRMS * templateResidualRMS);
		}
		float redChiSq = chiSq / (static_cast<float>(channel.waveform.size()) - static_cast<float>(finalFitParams.getNumParams()));

		sysProcWFCount++;
		meanReducedChiSq   = meanReducedChiSq + (redChiSq - meanReducedChiSq) / sysProcWFCount;
		chFit.reducedChiSq = redChiSq;

		{
			// Lock only for critical region (push_back on global list)
			std::lock_guard<std::mutex> guard(lock);
			reducedChiSqs.emplace_back(redChiSq);
		}

		chFits.push_back(chFit);

		// Update correction factors for future initial guesses
		updateGuessCorrector(amplitudes, times, initialAmplitudes, initialTimes, finalFitParams.baseline, initBaseline, finalFitParams.numPEs);

		// Cleanup
		delete templateInterpolator;
		delete functor;
		delete costFunction;
	}

	FitEvent evFitDat{event->ID, event->correctedTime, event->date, chFits};

	Writer writer(std::move(outputFile));
	writer.writeEventInfo(evFitDat);
}

/**
 * @brief Batch process a range of events, fitting them and optionally writing output.
 *
 * @param events       A vector of DigitiserEvent objects to process.
 * @param count        An atomic counter for the number of events processed.
 * @param lock         A mutex to serialize access to shared resources.
 * @param PETemplates  A map of channel ID to PETemplate objects for single PE shapes.
 * @param file         A shared pointer to a thread-safe output file structure.
 * @param sampleSpacing
 * @return             True if processing is successful for all events, otherwise false.
 */
bool batchFitEvents(const std::vector<DigitiserEvent>& events, std::atomic<unsigned long>& count, std::mutex& lock, const std::unordered_map<unsigned int, PETemplate>& PETemplates, const std::shared_ptr<SyncFile>& file, const float sampleSpacing) {
	for (const auto& event : events) {
		lock.lock();
		++count;
		lock.unlock();
		fitEvent(&event, PETemplates, file, lock, sampleSpacing);
	}
	return true;
}
