import shutil
import os
import subprocess
import numpy as np
from scipy.interpolate import interp1d, splrep, splev

from astropy.io import fits

from EMIR_file_processing import *
from templates import *
from reduction_tools import *

DEFAULT_FLAT_FIELD_FILE = 'master_flat_spec.fits'
FLAT_FIELD_FILENAME = 'master_flat_%s.fits'

class Observation(object):

	def __init__(self, observation_dir, result_dir=None, conda_check=True):

		'''
		Constructor for Observation object. 
	
		Instance Variables:
			- name (str): the name of the observation directory	
			- obs_dir (str): the directory where the raw istrument data is
			- result_dir (str): the directory where all outputs will be stored 
					(defaults to [obs_dir]/results)
			- emir_path (str): the directory where the EMIR inputs/outputs are processed
			- qc (Quality_Control): stored information from the quality control file
			- object_files (dict): a dictionary of the raw object data files (key: GRISM, value: list of files)
			- arc_files (dict): a dictionary of the raw arc data files (key: GRISM, value: list of files)
			- flat_files (dict): a dictionary of the raw flat data files (key: GRISM, value: list of files)
			- must_generate_flat_field (bool): flag to determine if the flat field must be generated or copied from defaults
	
		'''

		# check that conda environment is active 
		if os.environ['CONDA_DEFAULT_ENV'] == 'base' and conda_check:
			raise EnvironmentError("EMIR environment has not been activated")

		self.name = observation_dir.split('/')[-1]
		if self.name == '':
			self.name = observation_dir.split('/')[-2]
	
		self.obs_dir = observation_dir
		self.result_dir = self._get_result_dir(result_dir)

		self.qc = self._get_quality_control_info()
		self.object_files, self.arc_files, self.flat_files = self._get_files()
		
		self.must_generate_flat_field = (len(self.flat_files) != 0 )

		if not os.path.exists(os.path.join(self.result_dir, "EMIR_directory")):
			self._build_EMIR_directory()
		if os.path.exists(os.path.join(self.result_dir, "ABBA")):
			self.abba_dir = os.path.join(self.result_dir, "ABBA")
		if os.path.exists(os.path.join(self.result_dir, "response_curve")):
			self.response_curve_dir = os.path.join(self.result_dir, "response_curve")

	def _get_result_dir(self, result_dir):
		'''
		Finds the results directory and creates it if it doesn't yet exist.
		'''

		result_dir = result_dir if result_dir is not None else os.path.join(self.obs_dir, "results")

		if not os.path.exists(self.obs_dir):
			raise FileNotFoundError("observation directory %s does not exist"%self.obs_dir)

		if not os.path.exists(result_dir):
			os.mkdir(result_dir)

		return result_dir

	def _get_quality_control_info(self):
		'''
		Finds and parses the quality control file. 
		'''		

		for f in os.listdir(self.obs_dir):
			if f[-4:] == ".txt" and f[:3] == "GTC":
				qc_filename = os.path.join(self.obs_dir,f)
				break
		try:	
			return Quality_Control(qc_filename)
		except UnboundLocalError:
			raise FileNotFoundError("No valid quality control file for observation %s"%self.name)
			

	def _get_files(self):
		'''
		Gets a tuple of the file names for each obervation type and filter. 
		Checks that there are enough files for analysis.
		'''

		tables = ['object_table', 'arc_table', 'flat_table']
		files = [{}, {}, {}]
		for i, tname in enumerate(tables):

			t = getattr(self.qc, tname)
			if "FILTER" not in t.columns or 'IMAGE' not in t.columns:
				raise ValueError('quality control file not parsed correctly for obs %s'%self.name)

			files[i] = {key: [] for key in list(set(t['FILTER']))}
			for row in t:
				files[i][row['FILTER']] += [row['IMAGE']]

			pathtype = tname[:-len('_table')]
			path = os.path.join(self.obs_dir,pathtype)
			setattr(self, "%s_path"%pathtype, path)
		
		return tuple(files)		
		


	def _build_EMIR_directory(self):
		'''
		Builds the necessary files to do EMIR runs.
		'''

		setattr(self, 'emir_path', os.path.join(self.result_dir, "EMIR_files"))
		if not os.path.exists(self.emir_path):
			os.mkdir(self.emir_path)
		if not os.path.exists(os.path.join(self.emir_path, "data")):
			os.mkdir(os.path.join(self.emir_path, "data"))
		
		for fil in self.object_files:
			self._generate_control_file(fil)
			for analysis_type in ['arc', 'object']:
				self._generate_obs_res(analysis_type, fil)

		emir_defaults_dir=os.path.join(os.environ['EMIR_PIPE'], 'default_EMIR_files')
		for f in os.listdir(emir_defaults_dir):

			if f == DEFAULT_FLAT_FIELD_FILE and self.must_generate_flat_field:
				self._generate_flat_field()
			else:
				source_file = os.path.join(emir_defaults_dir, f)
				target_file = os.path.join(self.emir_path, 'data', f)
				shutil.copyfile(source_file, target_file)

	def _generate_flat_field(self): 
		'''
		Takes flat field files and generates the flat field file for EMIR analysis.
		'''

		file_base = os.path.join(self.obs_dir, 'flat/')

		for fil in self.flat_files:
			for i, f in enumerate(self.flat_files[fil]):
				try:
					im = fits.open(file_base+f)[1].data
				except IndexError:
					im = fits.open(file_base+f)[0].data
				if i == 0:
					ff_coadd = im
				else:
					ff_coadd += im

			emir_defaults_dir=os.path.join(os.environ['EMIR_PIPE'], 'default_EMIR_files')
			default_fits = fits.open(os.path.join(emir_defaults_dir, DEFAULT_FLAT_FIELD_FILE))
			default_fits.data = ff_coadd
			default_fits.writeto(os.path.join(self.emir_path, 'data', FLAT_FIELD_FILENAME%fil), overwrite=True)


	def initialize(self, analysis_type, filter_type):
		'''
		Prepares the emir directory to run the analysis type.
		'''

		self._copy_files(analysis_type, filter_type)
		self._generate_control_file(filter_type)
		self._generate_obs_res(analysis_type, filter_type)

	def _copy_files(self, analysis_type, filter_type):
		try:
			file_dict = getattr(self, "%s_files"%analysis_type)
			data_file_list = file_dict[filter_type]
			data_path = os.path.join(self.obs_dir, analysis_type)
			for fname in data_file_list:
				source_path = os.path.join(data_path, fname)
				target_path = os.path.join(self.emir_path, "data", fname)
				shutil.copyfile(source_path, target_path)

		except FileNotFoundError:
			raise ValueError("%s is not a valid analysis type (options: 'object', 'arc', 'flat')"%(analysis_type))
		except KeyError:
			raise KeyError("%s is not a valid filter type (options: %s)"\
							%(filter_type, str(list(file_dict.keys()))))


	def _generate_control_file(self, fil):
		'''
		Generates the control file used in the EMIR analysis.
		'''

		with open(os.path.join(self.emir_path, 'control_%s.yaml'%fil), 'w') as f:
			flat_fname = "master_flat_%s.fits"%fil
			f.write(CONTROL_YAML_TEMPLATE%(flat_fname,flat_fname))
			f.close()

	def _generate_obs_res(self, analysis_type, filter_type):
		'''
		Generates all observation result files for all GRISMs and analysis types.
		'''
		file_list = '' 
		with open(os.path.join(self.emir_path, '%s_obs_res_%s.yaml'%(analysis_type,filter_type)), 'w') as f:
			try:
				for i, filename in enumerate(getattr(self, "%s_files"%analysis_type)[filter_type]):
					if analysis_type == 'arc':
						file_list += " - %s\n"%filename
					else:
						file_list = " - %s\n"%filename
						text = OBS_RES_TEMPLATE%(self.name+"_"+filename.split('-')[0], file_list)
						f.write(text)
						if i < len(getattr(self, "%s_files"%analysis_type)[filter_type])-1:
							f.write('---\n')
				
				if analysis_type == 'arc':
					text = OBS_RES_TEMPLATE%(self.name+'_arc', file_list)
					f.write(text)
				
			except KeyError:
				#raise ValueError("attempting to generate an observation result file without proper files")
				pass

			f.close()

	def run_analysis(self, filter_type):
	
		self.rectify_and_calibrate(filter_type)
		self.ABBA_subtract(filter_type)
		self.convert_flux(filter_type)
		# save final results

	def rectify_and_analyze(self, analysis_type, filter_type):
		'''
		Run the numina recipe for rectification and calibration.
		'''

		orig_cwd = os.getcwd()
		os.chdir(self.emir_path)

		result = subprocess.run('numina run %s_obs_res_%s.yaml -r control_%s.yaml'\
								%(analysis_type,filter_type,filter_type),
							shell=True, capture_output=True, text=True)

		with open("emir_out.txt", 'w') as f:
			f.write(result.stdout+'\n\n')
			f.write(result.stderr+'\n\n')
			f.close()		

		os.chdir(orig_cwd)

	def ABBA_subtract(self, filter_type):

		self._prepare_ABBA(filter_type)
		self._subtract_and_save_ABBA(filter_type)

	def _prepare_ABBA(self, filter_type):
		'''
		Prepares the directory to save ABBA information, checks that 
		there is rectified and calibrated data.
		'''
		# if observation result file not generated, must run initialization
		obs_res_file = os.path.join(self.emir_path, "object_obs_res_%s.yaml"%filter_type)
		if not os.path.exists(obs_res_file):
			raise FileNotFoundError("Must generate object result file from _initialize")
		# if 'obsid%s_%s_results/reduced_mos.fits' not generated, must run rectification and calibration
		has_run_rect_and_cal = False
		for f in os.listdir(self.emir_path):
			comps = f.split('_') 
			if len(f) < 3: continue
			if self.name in comps[0] and comps[1].isdigit() and "results" in comps[2]:
				has_run_rect_and_cal = True
		if not has_run_rect_and_cal:
			raise FileNotFoundError("Must run rectify and calibrate on object file")

		abba_dir = os.path.join(self.result_dir, 'ABBA/')
		setattr(self, "abba_dir", abba_dir)
		if not os.path.exists(abba_dir):
			os.mkdir(abba_dir)


	def _get_calibrated_data(self, filter_type): 
		
		A = [] ; B = []
		fname_base = os.path.join(self.emir_path,'obsid%s_%s_results/reduced_mos.fits')

		for row in self.qc.object_table:
			if row['FILTER'] != filter_type: continue

			fID = row['IMAGE'].split('-')[0]
			fname = fname_base%(self.name, fID)
			if not os.path.exists(fname): continue

			if row['TELPOS'] == 'NOD_A':
				A += [fits.open(fname)[0].data.astype(float)]
			else:
				fits_file = fits.open(fname)
				header = fits_file[0].header
				B += [fits.open(fname)[0].data.astype(float)]

		if len(A) == 0 or len(B) == 0:
			raise ValueError("quality control file has incorrect names (observation %s)"%self.name)

		return A, B, header

	def _subtract_and_save_ABBA(self, filter_type, only_AB=False):
		'''
		Performs the actual subtraction and saves the results.
		'''
		A, B, header = self._get_calibrated_data(filter_type)

		if len(A) != 2 and len(B) != 2:
			self._subtract_and_save_non_ABBA(filter_type)

		fname = os.path.join(self.abba_dir, "ABBA_subtracted_%s.fits"%filter_type)

		primary_hdu = fits.ImageHDU(data=A[0] + A[1] - B[0] - B[1], header=header)
		primary_hdu.writeto(fname, overwrite=True)


	def _subtract_and_save_non_ABBA(self, filter_type):

		A, B, header = self._get_calibrated_data(filter_type)
		fname = os.path.join(self.abba_dir, "ABBA_subtracted_%s.fits"%filter_type)
		
		final_image = np.zeros(A[0].shape)
		for a,b in zip(A,B):
			final_image += (a-b)

		primary_hdu = fits.ImageHDU(data=final_image, header=header)
		primary_hdu.writeto(fname, overwrite=True)

	def _get_tabulated_interpolation(self, res_curve_path):

		fits_file = fits.open(res_curve_path)
		f = fits_file[0]

		flux_tab = f.data
		wavelengths_tab = pixels_to_wavelength(f.header)

		fits_file.close()

		return interp1d(wavelengths_tab, flux_tab, kind="linear")

	def _make_prelim_resp_curve(self, tabulated_spectrum_path, counts, wave, header, filter_type, **kwargs):

		tabulated_interpolation = self._get_tabulated_interpolation(tabulated_spectrum_path)
		response_prelim = tabulated_interpolation(wave) / counts
		if 'prelim_smooth_radius' in kwargs:
			response_prelim = smooth(response_prelim, radius=kwargs['prelim_smooth_radius']) 

		hdu_sp = fits.PrimaryHDU(data=response_prelim, header=header)
		hdu_sp.writeto(self.response_curve_dir+'/%s_response_curve.fits'%filter_type, overwrite=True)

		return response_prelim

	def _add_telluric_correction(self, file_path, filter_type, **kwargs):

		with fits.open(file_path) as f: 
			header = f[0].header
			flux = f[0].data
			wave = pixels_to_wavelength(header)

		telluric = get_atmospheric_spectrum(wave)
		response_notell = flux / telluric
		if 'telluric_smooth_radius' in kwargs:
			response_notell = smooth(response_notell, radius=kwargs['telluric_smooth_radius'])

		hdu_sp = fits.PrimaryHDU(data=response_notell, header=header)
		hdu_sp.writeto(self.response_curve_dir+'/%s_response_curve_telluric_corr.fits'%filter_type, overwrite=True)

		return response_notell

	def _do_abba_fitting(self, counts, wave, filter_type, header, **kwargs):

		spacing = kwargs.get('sample_spacing', len(counts)//25)
		if 'feature_mask' in kwargs:
			available_pixels = np.arange(0,len(counts), 1)[kwargs['feature_mask']]
			sample = np.array([pix for i,pix in enumerate(available_pixels) if i%spacing==0])
		else:
			sample = np.array([i for i in range(len(counts)) if i%spacing==0])

		sampled_counts = counts[sample]
		sampled_wave = wave[sample]

		tck = splrep(sampled_wave, sampled_counts)
		counts = splev(wave, tck)
	
		hdu_sp = fits.PrimaryHDU(data=counts, header=header)
		hdu_sp.writeto(self.response_curve_dir+'/%s_splined_response.fits'%filter_type, overwrite=True)

		return counts


	def make_response_curve(self, tabulated_spectrum_path, filter_type, **kwargs):
		'''Makes a response curve of the data.'''

		setattr(self, 'response_curve_dir', os.path.join(self.result_dir, "response_curve"))
		if not os.path.exists(self.response_curve_dir): os.mkdir(self.response_curve_dir)

		ABBA_file = os.path.join(self.abba_dir,'ABBA_subtracted_%s.fits'%filter_type)
		counts, wave, header = get_spectrum(ABBA_file, filter_type, kwargs.get('SN_position', None))

		hdu_sp = fits.PrimaryHDU(data=counts, header=header)
		hdu_sp.writeto(self.response_curve_dir+'/%s_spectrum_raw.fits'%filter_type, overwrite=True)

		if 'raw_smooth_radius' in kwargs:
			counts = smooth(counts, kwargs['raw_smooth_radius'])

		if kwargs.get("best_fit", True):
			counts = self._do_abba_fitting(counts, wave, filter_type, header, **kwargs)

		self._make_prelim_resp_curve(tabulated_spectrum_path, counts, wave, header, filter_type, **kwargs)
		if kwargs.get('correct_telluric_lines', True):
			self._add_telluric_correction(self.response_curve_dir+'/%s_response_curve.fits'%filter_type, 
									  	  filter_type, **kwargs)

	def get_response_curve(self, filter_type, telluric_corr=False):
		
		if telluric_corr:
			fname = self.response_curve_dir+'/%s_response_curve_telluric_corr.fits'%filter_type
		else:
			fname = self.response_curve_dir+'/%s_response_curve.fits'%filter_type

		with fits.open(fname) as f:
			header = f[0].header
			flux = f[0].data
			wave = pixels_to_wavelength(header)

		return wave, flux

	def get_reduced_spectrum(self, calibration_obs, filter_type, **kwargs):
		
		wave, resp_curve = calibration_obs.get_response_curve(filter_type, kwargs.get("telluric_correction", True))
		wave, raw_data = self.get_raw_spectrum(filter_type, **kwargs)

		return wave, raw_data*resp_curve

	def get_raw_spectrum(self, filter_type, **kwargs):
		raw_data, wave, header = get_spectrum(self.abba_dir+'/ABBA_subtracted_%s.fits'%filter_type, filter_type, 
											  SN_position=kwargs.get('SN_position', None), debug=kwargs.get('raw_debug', False))
		return wave, raw_data

	def _clean_emir_directory(self):
		for f in os.listdir(os.path.join(self.emir_path, "data/")):
			if f[0].isalpha(): continue
			os.remove(os.path.join(self.emir_path, "data/", f))

	def _clean_files(self):
		try:
			shutil.rmtree(self.emir_path)
			shutil.rmtree(self.result_dir)
		except FileNotFoundError:
			pass
	


