import unittest
import os

from CReader import ReadWCDataFile, ReadRecoMoreOutput

class TestCReader(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        setUpClass is called once for this test class. Here, we attempt to
        read in known reference files (one WaveCatcher data file and one
        RecoMore file) so we can inspect them in the actual tests.
        """
        cls.ref_wc_file = "TestData/R185.bin"
        cls.ref_rm_file = "TestData/R185PES_Reference.dat"

        # Make sure reference files exist (optional)
        # If they don't exist, you could skip tests or raise an error
        if not os.path.isfile(cls.ref_wc_file):
            raise FileNotFoundError(f"Missing reference WC data file: {cls.ref_wc_file}")
        if not os.path.isfile(cls.ref_rm_file):
            raise FileNotFoundError(f"Missing reference RecoMore file: {cls.ref_rm_file}")

        # Read the WaveCatcher data
        cls.wc_data = ReadWCDataFile(cls.ref_wc_file)

        # Read the RecoMore data
        cls.rm_data = ReadRecoMoreOutput(cls.ref_rm_file)

    def testWaveCatcherDataLoad(self):
        """
        Ensure the WaveCatcher reference file loads and
        that the contents match the expectations.
        """
        self.assertIsNotNone(self.wc_data, "ReadWCDataFile returned None.")

        # Retrieve the list of events
        events = self.wc_data.getEvents()
        self.assertEqual(len(events), 1044, "No events found in the WaveCatcher data file.")

        # Check the first event for known/expected fields
        first_event = events[0]

        self.assertEqual(first_event.eventID, 1, "First event ID is not the expected value.")
        self.assertGreater(len(first_event.chData), 0, "No channels in the first event.")
        self.assertEqual(first_event.date, "2022.4.4", "First channel has unexpected date.")
        self.assertEqual(first_event.TDCCorrTime, "16h11m57s,234.853.185ns", "First channel has unexpected time.")
        self.assertEqual(first_event.chData[0].channel, 0, "First event first channel has unexpected ID.")
        self.assertEqual(len(first_event.chData), 14, "First event has unexpected number of channels.")
        self.assertEqual(len(first_event.chData[0].waveform), 1024, "First event waveform data has unexpected length.")
        self.assertAlmostEqual(first_event.chData[3].waveform[260], -0.01776123046875, places=4, msg="Waveform data has unexpected value at pos 260 for event 1, channel 3.")


        last_event = events[-1]
        self.assertEqual(last_event.eventID, 1044, "Last event ID is not the expected value.")
        self.assertGreater(len(last_event.chData), 0, "No channels in the first event.")
        self.assertEqual(last_event.date, "2022.4.4", "Last event has unexpected date.")
        self.assertEqual(last_event.TDCCorrTime, "16h28m29s,781.000.000ns", "Last event has unexpected time.")
        self.assertEqual(last_event.chData[0].channel, 0, "Last event first channel has unexpected ID.")
        self.assertEqual(len(last_event.chData), 14, "Last event has unexpected number of channels.")
        self.assertEqual(len(last_event.chData[0].waveform), 1024, "Last event waveform data has unexpected length.")
        self.assertAlmostEqual(last_event.chData[1].waveform[338], -0.0348510742, places=4, msg="Waveform data has unexpected value at pos 338 for event 1044, channel 1.")


    def testRecoMoreDataLoad(self):
        """
        Ensure the RecoMore reference file loads and
        that the contents match the expectations.
        """
        self.assertIsNotNone(self.rm_data, "ReadRecoMoreOutput returned None.")

        # Retrieve the list of events
        events = self.rm_data.getEvents()
        self.assertGreater(len(events), 0, "No events found in the RecoMore data file.")

        # Check the first event for known/expected fields
        first_event = events[0]

        self.assertEqual(first_event.eventID, 1, "First event ID is not the expected value.")
        self.assertEqual(len(first_event.SiPM), 5, "First event has unexpected number of SiPM channels.")
        self.assertEqual(first_event.date, "2022.4.4", "First event has unexpected date.")
        self.assertEqual(first_event.TDCCorrTime, "16h11m57.234853185s", "First event has unexpected time.")


        first_channel = first_event.SiPM[0]
        self.assertEqual(first_channel.ch, 3, "First SiPM channel has unexpected ID.")
        self.assertEqual(len(first_channel.pes), 1, "No photoelectrons found where 1 was expected.")
        self.assertAlmostEqual(first_channel.redChiSq, 1.100520, places=4, msg="Unexpected channel reduced chi-squared.")
        self.assertAlmostEqual(first_channel.baseline, 6.70000008540228e-05, places=4, msg="Unexpected channel baseline.")

        first_pe = first_event.SiPM[0].pes[0]
        self.assertAlmostEqual(first_pe.amplitude, 0.040065, places=6, msg="Unexpected PE amplitude.")
        self.assertAlmostEqual(first_pe.time, 81.906364, places=4, msg="Unexpected PE time.")


        # Check the last event for known/expected fields
        last_event = events[-1]

        self.assertEqual(last_event.eventID, 1044, "Last event ID is not the expected value.")
        self.assertEqual(last_event.date, "2022.4.4", "Last event has unexpected date.")
        self.assertEqual(last_event.TDCCorrTime, "16h28m29.781000000s", "Last event has unexpected time.")
        self.assertEqual(len(last_event.SiPM), 4, "Last event has an unexpected number of SiPM channels.")

        # --- Channel 1 ---
        channel1 = last_event.SiPM[0]
        self.assertEqual(channel1.ch, 1, "First SiPM channel in last event has an unexpected ID.")
        self.assertAlmostEqual(channel1.redChiSq, 0.814246, places=6, msg="Unexpected reduced chi-squared for channel 1.")
        self.assertAlmostEqual(channel1.baseline, -0.000029, places=6, msg="Unexpected baseline for channel 1.")
        self.assertEqual(len(channel1.pes), 2, "Channel 1 should have 2 photoelectrons.")

        pe = channel1.pes[0]
        self.assertAlmostEqual(pe.amplitude, 0.020311, places=6, msg="Unexpected amplitude for first PE in channel 1.")
        self.assertAlmostEqual(pe.time, 105.202232, places=6, msg="Unexpected time for first PE in channel 1.")

        pe = channel1.pes[1]
        self.assertAlmostEqual(pe.amplitude, 0.021380, places=6, msg="Unexpected amplitude for second PE in channel 1.")
        self.assertAlmostEqual(pe.time, 105.939987, places=6, msg="Unexpected time for second PE in channel 1.")

        # --- Channel 3 ---
        channel3 = last_event.SiPM[1]
        self.assertEqual(channel3.ch, 3, "Second SiPM channel in last event has an unexpected ID.")
        self.assertAlmostEqual(channel3.redChiSq, 1.443000, places=6, msg="Unexpected reduced chi-squared for channel 3.")
        self.assertAlmostEqual(channel3.baseline, 0.000054, places=6, msg="Unexpected baseline for channel 3.")
        self.assertEqual(len(channel3.pes), 1, "Channel 3 should have 1 photoelectron.")

        pe = channel3.pes[0]
        self.assertAlmostEqual(pe.amplitude, 0.040654, places=6, msg="Unexpected amplitude for PE in channel 3.")
        self.assertAlmostEqual(pe.time, 93.132629, places=6, msg="Unexpected time for PE in channel 3.")

        # --- Channel 11 ---
        channel11 = last_event.SiPM[2]
        self.assertEqual(channel11.ch, 11, "Third SiPM channel in last event has an unexpected ID.")
        self.assertAlmostEqual(channel11.redChiSq, 0.903420, places=6, msg="Unexpected reduced chi-squared for channel 11.")
        self.assertAlmostEqual(channel11.baseline, 0.000051, places=6, msg="Unexpected baseline for channel 11.")
        self.assertEqual(len(channel11.pes), 4, "Channel 11 should have 4 photoelectrons.")

        pe = channel11.pes[0]
        self.assertAlmostEqual(pe.amplitude, 0.022264, places=6, msg="Unexpected amplitude for first PE in channel 11.")
        self.assertAlmostEqual(pe.time, 105.721817, places=6, msg="Unexpected time for first PE in channel 11.")

        pe = channel11.pes[1]
        self.assertAlmostEqual(pe.amplitude, 0.016481, places=6, msg="Unexpected amplitude for second PE in channel 11.")
        self.assertAlmostEqual(pe.time, 111.178963, places=6, msg="Unexpected time for second PE in channel 11.")

        pe = channel11.pes[2]
        self.assertAlmostEqual(pe.amplitude, 0.026122, places=6, msg="Unexpected amplitude for third PE in channel 11.")
        self.assertAlmostEqual(pe.time, 111.596588, places=6, msg="Unexpected time for third PE in channel 11.")

        pe = channel11.pes[3]
        self.assertAlmostEqual(pe.amplitude, 0.008681, places=6, msg="Unexpected amplitude for fourth PE in channel 11.")
        self.assertAlmostEqual(pe.time, 128.190704, places=6, msg="Unexpected time for fourth PE in channel 11.")

        # --- Channel 13 ---
        channel13 = last_event.SiPM[3]
        self.assertEqual(channel13.ch, 13, "Fourth SiPM channel in last event has an unexpected ID.")
        self.assertAlmostEqual(channel13.redChiSq, 1.541667, places=6, msg="Unexpected reduced chi-squared for channel 13.")
        self.assertAlmostEqual(channel13.baseline, -0.000064, places=6, msg="Unexpected baseline for channel 13.")
        self.assertEqual(len(channel13.pes), 1, "Channel 13 should have 1 photoelectron.")

        pe = channel13.pes[0]
        self.assertAlmostEqual(pe.amplitude, 0.042507, places=6, msg="Unexpected amplitude for PE in channel 13.")
        self.assertAlmostEqual(pe.time, 97.757668, places=6, msg="Unexpected time for PE in channel 13.")


    def testInvalidFileHandling(self):
        """
        (Optional) Example: test that an invalid or non-existent file
        raises an appropriate exception or handles errors gracefully.
        """
        non_existent = "no_such_file.dat"
        with self.assertRaises(RuntimeError):
            # Or FileNotFoundError, or your custom error, depending on how your C++ function handles it
            ReadWCDataFile(non_existent)


if __name__ == "__main__":
    unittest.main()
