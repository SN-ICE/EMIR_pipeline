import os
from test_library import *
from Observation import Observation

err_dict = {}

try:
	o = Observation('../objects/OB0011/')
	o.initialize('object', 'HK')
	o.rectify_and_analyze('object', 'HK')
	o.ABBA_subtract('HK')
except Exception as e:
	tb = traceback.format_exc()
	err_dict['INSTANTIATION_ERR'] = tb
	print_test_summary("FLUX CALIBRATION TESTS", err_dict)
	exit()



#o._clean_files()
print_test_summary("FLUX CALIBRATION TESTS", err_dict)



