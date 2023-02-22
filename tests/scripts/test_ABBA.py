import os
from test_library import *
from Observation import Observation

err_dict = {}

try:
	o = Observation('../objects/OB0011/')
except Exception as e:
	tb = traceback.format_exc()
	err_dict['INSTANTIATION_ERR'] = tb
	print_test_summary("ABBA TESTS", err_dict)
	exit()


# directory not initialized test
try:
	o.ABBA_subtract('HK')
	FILE_INIT_NOT_CAUGHT = "No file initialization was not caught"
except FileNotFoundError:
	FILE_INIT_NOT_CAUGHT = None
except Exception as e:
	tb = traceback.format_exc()
	FILE_INIT_NOT_CAUGHT = tb

err_dict['FILE_INIT_NOT_CAUGHT'] = FILE_INIT_NOT_CAUGHT

# initialze EMIR object files
try:
	o.initialize('object', 'HK')
except Exception as e:
	tb = traceback.format_exc()
	err_dict['FILE_INIT_ERR'] = tb
	print_test_summary("ABBA TESTS", err_dict)
	exit()

# analysis not completed test
try:

	o.ABBA_subtract('HK')
	PREANALYSIS_NOT_CAUGHT = "no pre-analyzed files not caught"

except FileNotFoundError:
	PREANALYSIS_NOT_CAUGHT = None
except Exception as e:
	tb = traceback.format_exc()
	PREANALYSIS_NOT_CAUGHT = tb

err_dict['ABBA ERROR FULL'] = PREANALYSIS_NOT_CAUGHT


# initialze EMIR object files
try:
	o.rectify_and_analyze('object', 'HK')
except Exception as e:
	tb = traceback.format_exc()
	err_dict['RECT_AND_CAL'] = tb
	print_test_summary("ABBA TESTS", err_dict)
	exit()

# now THIS should work
try:
	o.ABBA_subtract('HK')
	if os.path.exists(o.abba_dir+"/ABBA_subtracted_HK.fits"):
		FULL_ERR = None
	else:
		FULL_ERR = "Did not generate ABBA subtracted fits file"
except Exception as e:
	tb = traceback.format_exc()
	FULL_ERR = tb

err_dict['FULL_ERR'] = FULL_ERR
	
#o._clean_files()
print_test_summary("ABBA TEST", err_dict)



