#ifndef PETEMPLATE_H
#define PETEMPLATE_H

#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cmath>
#include <iostream>

class PETemplate {
public:
    /**
     * Constructor that reads and interpolates the ideal waveform from a file.
     *
     * - We infer the raw time spacing from consecutive lines in the file
     *   and confirm it is uniform (within a floating-point tolerance).
     * - The final interpolated time spacing is `rawDt / interpFactor`.
     * - We apply linear interpolation between each pair of raw points,
     *   multiplying *only* the intermediate steps by a sign factor for
     *   positive or negative pulse (following original snippet).
     * - We find the index of the minimum amplitude (if `positivePulse == false`)
     *   or maximum amplitude (if `positivePulse == true`) in the **interpolated**
     *   waveform. We store that in `tMaxIntensity`.
     * - Finally, we check that the amplitude at `tMaxIntensity` is near −1 (for
     *   negative pulses) or near +1 (for positive pulses).
     *
     * @param channelNumber   The channel index (used for finding the file).
     * @param interpFactor    Number of points to insert between consecutive
     *                        raw points during interpolation.
     * @param idealWFDir      Path to directory containing the ideal waveform
     *                        files (e.g., "/path/to/ideal/"). This constructor
     *                        looks for "ch<channel>.txt" in that directory.
     * @param positivePulse   If `true`, this is a positive pulse; otherwise negative.
     */
    PETemplate(unsigned int channelNumber,
               unsigned int interpFactor,
               const std::string &idealWFDir,
               bool positivePulse)
        : channelNumber(channelNumber)
    {
        // Build the file path for the channel.
        std::string idealWFPath = idealWFDir + "ch" + std::to_string(channelNumber) + ".txt";
        std::ifstream idealWFFile(idealWFPath);

        // Basic file checks:
        if (!idealWFFile.is_open()) {
            throw std::runtime_error("Ideal PE PDF file: " + idealWFPath + " not found.");
        }
        if (idealWFFile.peek() == std::ifstream::traits_type::eof()) {
            throw std::runtime_error("Ideal PE PDF file: " + idealWFPath + " is empty.");
        }
        if (idealWFFile.fail()) {
            throw std::runtime_error("Ideal PE PDF file: " + idealWFPath + " could not be opened.");
        }

        // Read all times and amplitudes from the file.
        std::vector<double> rawTimes;
        std::vector<double> rawAmps;
        {
            double timeVal, ampVal;
            while (idealWFFile >> timeVal >> ampVal) {
                rawTimes.push_back(timeVal);
                rawAmps.push_back(ampVal);
            }
        }

        // Need at least two points to infer a time step
        if (rawTimes.size() < 2) {
            throw std::runtime_error(
                "Not enough samples in " + idealWFPath + " to determine time spacing."
            );
        }

        // Check that raw time steps are constant within tolerance
        const double referenceDt = rawTimes[1] - rawTimes[0];
        if (referenceDt <= 0.0) {
            throw std::runtime_error(
                "Non-positive time spacing detected in " + idealWFPath +
                " (the first time difference is <= 0)."
            );
        }

        constexpr double timeTolerance = 1e-3; // Adjust as needed
        for (size_t i = 2; i < rawTimes.size(); ++i) {
            double dt = rawTimes[i] - rawTimes[i - 1];
            if (dt <= 0.0) {
                throw std::runtime_error(
                    "Non-positive or reversed time spacing detected in " + idealWFPath +
                    " at index " + std::to_string(i) +
                    ". rawTimes[" + std::to_string(i - 1) + "] = " + std::to_string(rawTimes[i - 1]) +
                    ", rawTimes[" + std::to_string(i)     + "] = " + std::to_string(rawTimes[i]) +
                    ", dt = " + std::to_string(dt)
                );
            }

            const double diff = std::fabs(dt - referenceDt);
            if (diff > timeTolerance) {
                throw std::runtime_error(
                    "Non-uniform time spacing detected in " + idealWFPath +
                    " at index " + std::to_string(i) + ".\n" +
                    "  referenceDt = " + std::to_string(referenceDt) + "\n" +
                    "  currentDt   = " + std::to_string(dt) + "\n" +
                    "  difference  = " + std::to_string(diff) + "\n" +
                    "  tolerance   = " + std::to_string(timeTolerance) + "\n" +
                    "  rawTimes[" + std::to_string(i - 1) + "] = " + std::to_string(rawTimes[i - 1]) + "\n" +
                    "  rawTimes[" + std::to_string(i)     + "] = " + std::to_string(rawTimes[i])
                );
            }
        }


        // Final interpolation spacing
        double interpolatedDt = referenceDt / static_cast<double>(interpFactor);
        timeSpacing = static_cast<float>(interpolatedDt);

        // signFactor as in the original snippet
        // (Recall: The original code multiplied only the intermediate steps by signFactor,
        //  while the final raw sample each segment was unmodified.)
        const float signFactor = positivePulse ? -1.f : 1.f;

        // Interpolate
        voltages.reserve((rawTimes.size() - 1)*interpFactor + 1);
        // Start with the first amplitude as-is
        voltages.push_back(rawAmps[0]);

        for (size_t i = 1; i < rawTimes.size(); ++i) {
            double deltaAmp = (rawAmps[i] - rawAmps[i - 1]) / static_cast<double>(interpFactor);

            // Insert linearly interpolated points (intermediate ones use signFactor)
            for (unsigned int step = 1; step < interpFactor; ++step) {
                double interpAmp = rawAmps[i - 1] + step * deltaAmp * signFactor;
                voltages.push_back(interpAmp);
            }

            // Push back the final raw amplitude unmodified
            voltages.push_back(rawAmps[i]);
        }

        // Find the index of the amplitude extremum
        //   - min amplitude if positivePulse == false
        //   - max amplitude if positivePulse == true
        unsigned int extremumIndex = 0;
        double extremumValue = voltages[0];

        if (positivePulse) {
            // Looking for maximum
            for (size_t i = 1; i < voltages.size(); ++i) {
                if (voltages[i] > extremumValue) {
                    extremumValue = voltages[i];
                    extremumIndex = i;
                }
            }
        } else {
            // Looking for minimum
            for (size_t i = 1; i < voltages.size(); ++i) {
                if (voltages[i] < extremumValue) {
                    extremumValue = voltages[i];
                    extremumIndex = i;
                }
            }
        }

        // Store that index in tMaxIntensity
        tMaxIntensity = extremumIndex;

        // Now check the amplitude at that index
        // For negative pulses: amplitude should be near -1
        // For positive pulses: amplitude should be near +1
        const double amplitudeAtExtremum = voltages[extremumIndex];
        const double expectedValue = positivePulse ? 1.0 : -1.0;
        const double checkTolerance = 0.00001; // for example: ±0.00001 from target
        const double diff = std::fabs(amplitudeAtExtremum - expectedValue);

        if (diff > checkTolerance) {
            throw std::runtime_error(
                "Amplitude check failed at index " + std::to_string(extremumIndex) +
                ". Expected ~" + std::to_string(expectedValue) + " but got " +
                std::to_string(amplitudeAtExtremum) + " (diff=" + std::to_string(diff) + ")"
            );
        }
    }

    // Accessors
    inline std::vector<double> getVoltages() const { return voltages; }
    inline std::vector<double> const *getVoltagesRef() const { return &voltages; }
    inline unsigned int getChannelNumber() const { return channelNumber; }

    /**
     * Time spacing between consecutive interpolated points
     * (in the same units as in the input file).
     */
    inline float getTimeSpacing() const { return timeSpacing; }

    /**
     * The index in 'voltages' where the extremum amplitude was found.
     * - If `positivePulse == true`, it is the index of the maximum amplitude.
     * - If `positivePulse == false`, it is the index of the minimum amplitude.
     */
    inline float getTMaxIntensity() const { return tMaxIntensity; }

private:
    std::vector<double> voltages;      ///< The interpolated waveform amplitudes
    unsigned int channelNumber;        ///< E.g., 0, 1, 2, ...
    float timeSpacing = 0.0f;          ///< Interpolated spacing between consecutive samples
    unsigned int tMaxIntensity = 0;    ///< Index of the extremum amplitude
};



#endif //PETEMPLATE_H
