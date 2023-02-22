
import os
import shutil

from test_library import *
from Observation import Observation


error_dict = {}

# vanilla test
try:
	o = Observation("../objects/OB0011", "../output/OB0011")
	o._clean_files()
	VANILLA_ERR = None
except Exception as e:
	tb = traceback.format_exc()
	VANILLA_ERR = tb

error_dict['VANILLA_ERR'] = VANILLA_ERR

# vanilla test with results in same place as observation
try:
	o = Observation("../objects/OB0011")
	if os.path.exists("../objects/OB0011/results"):
		VANILLA_PLUS_ERR = None
	else:
		VANILLA_PLUS_ERR = "Did not create results directory"
	o._clean_files()
except Exception as e:
    tb = traceback.format_exc()
    VANILLA_PLUS_ERR = tb
        
error_dict['VANILLA_PLUS_ERR'] = VANILLA_PLUS_ERR

# observation directory does not exist
try:
	o = Observation("nonexistant_dir")
	INFILE_ERR = "Did not throw an error"
except FileNotFoundError as e:
    INFILE_ERR = None
except Exception as e:
	tb = traceback.format_exc()
	INFILE_ERR = tb

error_dict['INFILE_ERR'] = INFILE_ERR


# observation directory does not have an 'object' directory

dirname = '../objects/fake_obj'
os.mkdir(dirname)
try:
	o = Observation(dirname)
	OBJ_DIR_ERR = "Did not throw an error"
except FileNotFoundError as e:
	OBJ_DIR_ERR = None
except Exception as e:
	tb = traceback.format_exc()
	OBJ_DIR_ERR = tb

shutil.rmtree(dirname)

error_dict['OBJ_DIR_ERR'] = OBJ_DIR_ERR

# test name harvest 1

try:
	o = Observation("../objects/OB0011", "../output/OB0011")
	o._clean_files()

	if o.name == "OB0011":
		NAME_ERR1 = None
	else:
		NAME_ERR1 = "%s is not the correct name"%o.name
   
except Exception as e:
	tb = traceback.format_exc()
	NAME_ERR1 = tb

error_dict['NAME_ERR1'] = NAME_ERR1

# test name harvest 2

try:
	o = Observation("../objects/OB0011/", "../output/OB0011")
	o._clean_files()

	if o.name == "OB0011":
		NAME_ERR2 = None
	else:
		NAME_ERR2 = "%s is not the correct name"%o.name     
except Exception as e:
	tb = traceback.format_exc()
	NAME_ERR2 = tb   

error_dict['NAME_ERR2'] = NAME_ERR2


## regular initialization: check that all files are in the right place


# object files stored correctly
try:
	o = Observation("../objects/OB0011", "../output/OB0011")
	all_obj_files = []
	for key in o.object_files:
		all_obj_files += o.object_files[key]

	orig_obj_files = os.listdir("../objects/OB0011/object")
	if len(all_obj_files) == len(orig_obj_files):
		OBJ_STORE = None
		for f in all_obj_files:
			if f not in orig_obj_files:
				OBJ_STORE = "Wrong object files stored"
				break
	else:
		OBJ_STORE = "Number of object files is wrong"
except Exception as e:
	tb = traceback.format_exc()
	OBJ_STORE = tb
error_dict['OBJ_STORE'] = OBJ_STORE

# arc files stored correctly
try:
	all_arc_files = []
	for key in o.arc_files:
		all_arc_files += o.arc_files[key]

	orig_arc_files = os.listdir("../objects/OB0011/arc")
	if len(all_arc_files) == len(orig_arc_files):
		ARC_STORE = None
		for f in all_arc_files:
			if f not in orig_arc_files:
				ARC_STORE = "Wrong arc files stored"
	else:
		ARC_STORE = "Arc files not recorded correctly"
except Exception as e:
    tb = traceback.format_exc()
    ARC_STORE = tb
error_dict['ARC_STORE'] = ARC_STORE


# flat field objects stored correctly
try:

	all_flat_files = []
	for key in o.flat_files:
		all_flat_files += o.flat_files[key]

	orig_flat_files = os.listdir("../objects/OB0011/flat")
	if len(all_flat_files) == len(orig_flat_files):
		FLAT_STORE = None
		for f in all_flat_files:
			if f not in orig_flat_files:
				FLAT_STORE = "Wrong flat files stored"
	else:
		FLAT_STORE = "Flat field files not recorded correctly"
		o._clean_files()
except Exception as e:
	tb = traceback.format_exc()
	FLAT_STORE = tb
	
error_dict['FLAT_STORE'] = FLAT_STORE

# initialize 

try:
	o.initialize("object", "HK")
	data_to_analyze = os.path.join(o.emir_path, "data")

	if 10+len(o.object_files['HK']) == len(os.listdir(data_to_analyze)):
		OBJ_INIT = None
	else:
		OBJ_INIT = "Object analysis not initialized correctly"
	o._clean_emir_directory()
except Exception as e:
    tb = traceback.format_exc()
    OBJ_INIT = tb

error_dict['OBJ_INIT'] = OBJ_INIT

try:
	o.initialize("arc", "HK")
	data_to_analyze = os.path.join(o.emir_path, "data")

	if 10+len(o.arc_files['HK']) == len(os.listdir(data_to_analyze)):
		ARC_INIT = None
	else:
		ARC_INIT = "Arc analysis not initialized correctly"
	o._clean_emir_directory()
except Exception as e:
    tb = traceback.format_exc()
    ARC_INIT = tb

error_dict['ARC_INIT'] = ARC_INIT

try:
	o.initialize("flat", "HK")

	data_to_analyze = os.path.join(o.emir_path, "data")
	if 10+len(o.flat_files['HK']) == len(os.listdir(data_to_analyze)):
		FLAT_INIT = None
	else:
		FLAT_INIT = "Object analysis not initialized correctly"
	o._clean_emir_directory()
except Exception as e:
    tb = traceback.format_exc()
    FLAT_INIT = tb

error_dict['FLAT_INIT'] = FLAT_INIT


# control file is generated

try:
	o._generate_control_file('HK')
	if os.path.exists(os.path.join(o.emir_path, "control_HK.yaml")):
		CNTRL_FILE = None
	else:
		CNTRL_FILE = "No control.yaml file generated."
	o._clean_emir_directory()
except Exception as e:
    tb = traceback.format_exc()
    CNTRL_FILE = tb
error_dict['CNTRL_FILE'] = CNTRL_FILE


# observation result files are generated correctly

try:
	o.initialize("object", 'HK')
	o._generate_obs_res("object", 'HK')
	with open(os.path.join(o.emir_path, "object_obs_res_HK.yaml"), 'r') as f:
		name_line = f.readline()
		file_list = []
		for line in f:
			if line[:3] == " - ":
				file_list += [line[3:]]

		f.close()
	
	if name_line.split(':')[-1] != o.name:
		OBSRES_NAME_O = None
	else:
		OBSRES_NAME_O = "Object's Observation Result file does not have the right name"

	if len(file_list) == len(o.object_files['HK']):
		OBSRES_FILES_O = None
		for f in file_list:
			if f[:-1] not in o.object_files['HK']:
				OBSRES_FILES_O = "Wrong file in object's obsres file"
	else:
		OBSRES_FILES_O = "Wrong number of files in object's observation result file."
except FileNotFoundError:
	OBSRES_NAME_O = "No object_obs_res.yaml file generated"
	OBSRES_FILES_O = "No object_obs_res.yaml file generated"
except Exception as e:
	tb = traceback.format_exc()
	OBSRES_NAME_O = tb
	OBSRES_FILES_O = tb

error_dict['OBSRES_NAME_O'] = OBSRES_NAME_O
error_dict['OBSRES_FILES_O'] = OBSRES_FILES_O

try:
	o.initialize("arc", 'HK')
	o._generate_obs_res("arc", 'HK')
	with open(os.path.join(o.emir_path, "arc_obs_res_HK.yaml"), 'r') as f:
		name_line = f.readline()
		file_list = []
		for line in f:
			if line[:3] == " - ":
				file_list += [line[3:]]

		f.close()
	if name_line.split(':')[-1] != o.name:
		OBSRES_NAME_A = None
	else:
		OBSRES_NAME_A = "Arc's Observation Result file does not have the right name"

	num_files = len(os.listdir('../objects/OB0011/arc/'))
	if num_files == 10+len(o.arc_files):
		OBSRES_FILES_A = None
	else:
		OBSRES_FILES_A = "Wrong number of files in arc's observation result file."
	o._clean_files()
except FileNotFoundError:
	OBSRES_NAME_A = "No arc_obs_res.yaml file generated"
	OBSRES_FILES_A = "No arc_obs_res.yaml file generated"
except Exception as e:
	tb = traceback.format_exc()
	OBSRES_NAME_A = tb
	OBSRES_FILES_A = tb

error_dict['OBSRES_NAME_A'] = OBSRES_NAME_A
error_dict['OBSRES_FILES_A'] = OBSRES_FILES_A

# flat field files generated correctly

try:
	o = Observation("../objects/OB0011", "../output/OB0011")

	if os.path.exists("../output/OB0011/EMIR_files/data/master_flat_HK.fits"):
		FLAT_FILES = None
	else:
		FLAT_FILES = "No flat file generated"
	o._clean_files()
except Exception as e:
	tb = traceback.format_exc()
	FLAT_FILES = tb

error_dict['FLAT_FILES'] = FLAT_FILES

try:
	o = Observation("../objects/OB0011_no_flat", "../output/OB0011_no_flat")
	if os.path.exists("../output/OB0011_no_flat/EMIR_files/data/master_flat_spec.fits"):
		FLAT_FILES_N = None
	else:
		FLAT_FILES_N = "No flat file generated when there is no initial flat files"
	o._clean_files()
except Exception as e:
	tb = traceback.format_exc()
	FLAT_FILES_N = tb

error_dict['FLAT_FILES_N'] = FLAT_FILES_N


print_test_summary("FILE SYSTEM CHECKS", error_dict)




