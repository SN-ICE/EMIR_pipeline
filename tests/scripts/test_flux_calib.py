import os
from test_library import *
from Observation import Observation

err_dict = {}

try:
	o = Observation('../objects/OB0011/')
	if not os.path.exists(o.abba_dir+'/ABBA_subtracted_YJ.fits'):
		o.initialize('object', 'YJ')
		o.rectify_and_analyze('object','YJ')
		o.ABBA_subtract('YJ')
except Exception as e:
	tb = traceback.format_exc()
	err_dict['INIT_ERR'] = tb
	print_test_summary("FLUX CALIBRATION", err_dict)
	exit()

try:
	o.make_response_curve('../objects/uka0v.fits', 'YJ')
	err_dict['RESP_CURVE_NO_KWARGS'] = None
except Exception as e:
	tb = traceback.format_exc()
	err_dict['RESP_CURVE_NO_KWARGS'] = tb

try:
	o.make_response_curve('../objects/uka0v.fits', 'YJ', SN_position=(20,30))
	err_dict['RESP_CURVE_SN_POS'] = None
except Exception as e:
	tb = traceback.format_exc()
	err_dict['RESP_CURVE_SN_POS'] = tb

try:
	o.make_response_curve('../objects/uka0v.fits', 'YJ', prelim_smooth_radius=300)
	err_dict['RESP_CURVE_PRELIM_SMOOTH'] = None
except Exception as e:
	tb = traceback.format_exc()
	err_dict['RESP_CURVE_PRELIM_SMOOTH'] = tb

try:
	o.make_response_curve('../objects/uka0v.fits', 'YJ', correct_telluric_lines=False)
	err_dict['RESP_CURVE_NO_TELL'] = None
except Exception as e:
	tb = traceback.format_exc()
	err_dict['RESP_CURVE_NO_TELL'] = tb





print_test_summary("FLUX CALIBRATION", err_dict)




