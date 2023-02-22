
import os
from test_library import *
from Observation import Observation

err_dict = {}

try:
	o = Observation('../objects/OB0011/')
except Exception as e:
	tb = traceback.format_exc()
	err_dict['INIT_ERR'] = tb
	print_test_summary("RECTIFICATION AND CALIBRATION TESTS", err_dict)
	exit()

# initialize components test

try:
	o._copy_files('arc', 'HK')
	if len(o.arc_files['HK']) == len(os.listdir(o.emir_path+"/data"))-10:
		for f in o.arc_files['HK']:
			if f not in os.listdir(o.emir_path+"/data"):
				COPY_ERR = "files not copied over correctly"
				break
			COPY_ERR = None
	else:
		COPY_ERR = "File system not cleaned"
				
	o._generate_obs_res('arc', 'HK')
	o.rectify_and_analyze('arc','HK')

	if os.path.exists(o.emir_path+"/obsidOB0011_arc_results/rectwv_coeff.json"):
		VANILLA_ERR = None
	else:
		VANILLA_ERR = "Error in numina run, check log files"

	o._clean_emir_directory()
except Exception as e:
	tb = traceback.format_exc()
	VANILLA_ERR = tb

err_dict['COPY_ERR'] = COPY_ERR
err_dict['ARC_REC_AND_CAL_COMP'] = VANILLA_ERR


# initialize black box test
try:
	o.initialize('arc', 'HK')
	o.rectify_and_analyze('arc','HK')

	if len(o.arc_files['HK']) == len(os.listdir(o.emir_path+"/data"))-10:
		for f in o.arc_files['HK']:
			if f not in os.listdir(o.emir_path+"/data"):
				COPY_ERR = "files not copied over correctly"
				break
			COPY_ERR = None
	else:
		COPY_ERR = "File system not cleaned"

	if os.path.exists(o.emir_path+"/obsidOB0011_arc_results/rectwv_coeff.json"):
		VANILLA_ERR = None
	else:
		VANILLA_ERR = "Error in numina run, check log files"

except Exception as e:
	tb = traceback.format_exc()
	VANILLA_ERR = tb
	COPY_ERR = tb

err_dict['BB_COPY_ERR'] = COPY_ERR
err_dict['ARC_REC_AND_CAL_FULL'] = VANILLA_ERR



print_test_summary("RECTIFY AND CALIBRATE STEPS TESTS", err_dict)



