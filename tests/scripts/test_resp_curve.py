import unittest
from Observation import Observation

import numpy as np
import os

class TestObservation(unittest.TestCase):

    def setUp(self):
        # set up a fake observation object
        self.obs = Observation('../objects/OB0011/')

        if not os.path.exists(self.obs.abba_dir+'/ABBA_subtracted_HK.fits'):
            self.obs.initialize('object', 'HK')
            self.obs.rectify_and_analyze('object','HK')
            self.obs.ABBA_subtract('HK')

    def test_make_response_curve(self):
        # Test response curve generation with a test A0V spectrum
        test_spectrum_path = "../objects/uka0v.fits"
        filter_type = "HK"
        self.obs.make_response_curve(test_spectrum_path, filter_type)

        # Check that the response curve file was created
        response_curve_path = self.obs.response_curve_dir + f"/{filter_type}_response_curve.fits"
        self.assertTrue(os.path.exists(response_curve_path))

        # Check that the response curve has the expected dimensions
        wave, flux = self.obs.get_response_curve(filter_type)
        self.assertEqual(wave.shape, flux.shape)

    def test_get_response_curve(self):
        '''Test the get_response_curve method'''

        # Test case with valid input file
        filter_type = "HK"
        test_spectrum_path = "../objects/uka0v.fits"
        self.obs.make_response_curve(test_spectrum_path, filter_type)
        wave, flux = self.obs.get_response_curve(filter_type)
        self.assertIsInstance(wave, np.ndarray)
        self.assertIsInstance(flux, np.ndarray)
        self.assertGreater(len(wave), 0)
        self.assertGreater(len(flux), 0)

        # Test case with missing input file
        with self.assertRaises(FileNotFoundError):
            self.obs.get_response_curve("nonexistent_filter")
    
if __name__ == '__main__':
    unittest.main()



