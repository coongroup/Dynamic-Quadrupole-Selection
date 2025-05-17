from unittest import TestCase
from trianglelib import Spectrum, ExperimentSpectra, TriangleProfile


class TestSpectrum(TestCase):
    def setUp(self):
        self.spectrum = Spectrum()
        self.spectrum.add_ms_data([0, 201.1, 202.5, 204, 207, 210], [60, 0, 100, 200, 150, 300])

    def test_add_msdata(self):
        self.setUp()
        self.assertEqual(len(self.spectrum.get_ms_data()), 2)
        self.assertEqual(len(self.spectrum.get_ms_data()[1]), 6)

    def test_filter_data(self):
        self.setUp()
        self.spectrum.filter_data(50)
        self.assertEqual(len(self.spectrum.get_ms_data()[1]), len(self.spectrum.get_ms_data()[0]),
                         "filterData: Different length mass and intensity lists")


class TestExperimentSpectra(TestCase):
    def setUp(self):
        """s1 = Spectrum()
        s1.isolationCenter = 202.1
        s1.isolation_width = 5.5
        s1.addMSData([201, 202, 203, 204, 205, 206],      [100, 150, 0, 100, 100, 0])
        s1.filterData()
        s2 = Spectrum()
        s2.isolationCenter = 203.1
        s2.isolation_width = 5.5
        s2.addMSData([201, 202, 203, 204, 205, 206, 207], [50, 100, 100, 101, 0, 0, 100])
        s2.filterData()
        s3 = Spectrum()
        s3.isolationCenter = 204.1
        s3.isolation_width = 5.5
        s3.addMSData([201, 202, 203, 204, 205, 206, 207], [40, 50, 200, 102, 100, 50, 100])
        s3.filterData()
        s4 = Spectrum()
        s4.isolationCenter = 205.1
        s4.isolation_width = 5.5
        s4.addMSData([201, 202, 203, 204, 205, 206, 207], [30, 40, 100, 103, 100, 100, 100])
        s4.filterData()
        s5 = Spectrum()
        s5.isolationCenter = 206.1
        s5.isolation_width = 5.5
        s5.addMSData([202, 203, 204, 205, 206, 207, 208],     [10, 0, 100, 104, 0, 200, 0])
        s5.filterData()
        s6 = Spectrum()
        s6.isolationCenter = 207.1
        s6.isolation_width = 5.5
        s6.addMSData([1], [100])
        s6.filterData()"""
        s1 = Spectrum()
        s1.isolationCenter = 521
        s1.isolation_width = 5.5
        s1.add_ms_data([202], [100])
        s2 = Spectrum()
        s2.isolationCenter = 522
        s2.isolation_width = 5.5
        s2.add_ms_data([202], [200])
        s3 = Spectrum()
        s3.isolationCenter = 523
        s3.isolation_width = 5.5
        s3.add_ms_data([202], [200])
        s4 = Spectrum()
        s4.isolationCenter = 524
        s4.isolation_width = 5.5
        s4.add_ms_data([202], [100])
        s5 = Spectrum()
        s5.isolationCenter = 524
        s5.isolation_width = 5.5
        s5.add_ms_data([204], [100])
        self.experiment = ExperimentSpectra(1)
        self.experiment.addSpectrum(s1)
        self.experiment.addSpectrum(s2)
        self.experiment.addSpectrum(s3)
        self.experiment.addSpectrum(s4)
        self.experiment.addSpectrum(s5)

    def test_extractProfilesFromExperiment(self):
        self.setUp()
        self.triangleProfiles = self.experiment.extract_profiles_from_experiment()
        assert False
