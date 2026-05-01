import shutil
import os
import subprocess
import sys
import numpy as np
from scipy.interpolate import interp1d, splrep, splev

from astropy.io import fits

from EMIR_file_processing import *
from templates import *
from reduction_tools import *
from interactive_reduction_tools import *

DEFAULT_FLAT_FIELD_FILE = 'master_flat_spec.fits'
FLAT_FIELD_FILENAME = 'master_flat_%s.fits'

zp_fluxes_vega = {'J': 3.129e-10, 'H': 1.133e-10, 'K': 4.283e-11}

class Observation(object):

	def __init__(self, observation_dir, result_dir=None, conda_check=True, qc_fname=None):

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

		self.qc = self._get_quality_control_info(qc_fname)
		self.object_files, self.arc_files, self.flat_files = self._get_files()
		self.grisms = self._get_available_grisms()
		
		self.must_generate_flat_field = (len(self.flat_files) != 0 )

		self.emir_path = os.path.join(self.result_dir, 'EMIR_files')
		self.abba_dir = os.path.join(self.result_dir, "ABBA")
		self.sens_function_dir = os.path.join(self.result_dir, "sens_function")
		if not os.path.exists(os.path.join(self.result_dir, "EMIR_directory")):
			self._build_EMIR_directory()

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

	def _get_quality_control_info(self, qc_filename=None):
		'''
		Finds and parses the quality control file. 
		'''		

		if qc_filename is None:
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
		Gets a tuple of the file names for each obervation type and grism. 
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
		
	def _get_available_grisms(self):
		filters = list(set(self.qc.object_table['FILTER']))
		return filters

	def _build_EMIR_directory(self):
		'''
		Builds the necessary files to do EMIR runs.
		'''

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
			elif f == DEFAULT_FLAT_FIELD_FILE:
				source_file = os.path.join(emir_defaults_dir, f)
				target_file = os.path.join(self.emir_path, 'data', FLAT_FIELD_FILENAME)

				for grism in self.grisms:
					shutil.copyfile(source_file, target_file%grism)
			else:
				source_file = os.path.join(emir_defaults_dir, f)
				target_file = os.path.join(self.emir_path, 'data', f)
				shutil.copyfile(source_file, target_file)

		open(os.path.join(self.result_dir, "EMIR_directory"), 'w').close()

	def _generate_flat_field(self): 
		'''
		Takes flat field files and generates the flat field file for EMIR analysis.
		'''

		file_base = os.path.join(self.obs_dir, 'flat/')

		for fil in self.flat_files:
			for i, f in enumerate(self.flat_files[fil]):
				fi = fits.open(file_base+f)
				im = fi[-1].data
				if i == 0:
					ff_coadd = im
				else:
					ff_coadd += im
				fi.close()

			emir_defaults_dir=os.path.join(os.environ['EMIR_PIPE'], 'default_EMIR_files')
			default_fits = fits.open(os.path.join(emir_defaults_dir, DEFAULT_FLAT_FIELD_FILE))
			default_fits[0].data = ff_coadd
			default_fits.writeto(os.path.join(self.emir_path, 'data', FLAT_FIELD_FILENAME%fil), overwrite=True)
			default_fits.close()


	def initialize(self, analysis_type, grism_type):
		'''
		Prepares the emir directory to run the analysis type.
		'''

		self._copy_files(analysis_type, grism_type)
		self._generate_control_file(grism_type)
		self._generate_obs_res(analysis_type, grism_type)

	def _copy_files(self, analysis_type, grism_type):
		try:
			file_dict = getattr(self, "%s_files"%analysis_type)
			data_file_list = file_dict[grism_type]
			data_path = os.path.join(self.obs_dir, analysis_type)
			if not os.path.exists(os.path.join(self.obs_dir, analysis_type)):
				raise ValueError("'%s' is not a valid analysis type (options: 'object', 'arc', 'flat')"%(analysis_type))

			for fname in data_file_list:
				if not os.path.exists(os.path.join(self.obs_dir, analysis_type, fname)):
					print('WARNING: Qualty control file lists a file that does not exist (%s)'%fname)
					continue
				source_path = os.path.join(data_path, fname)
				target_path = os.path.join(self.emir_path, "data", fname)
				shutil.copyfile(source_path, target_path)

		except KeyError:
			raise KeyError("%s is not a valid grism type for %s analysis (options: %s)"\
							%(grism_type, analysis_type, str(list(file_dict.keys()))))


	def _generate_control_file(self, fil):
		'''
		Generates the control file used in the EMIR analysis.
		'''

		with open(os.path.join(self.emir_path, 'control_%s.yaml'%fil), 'w') as f:
			flat_fname = "master_flat_%s.fits"%fil
			f.write(CONTROL_YAML_TEMPLATE%(flat_fname,flat_fname,flat_fname,flat_fname))
			f.close()

	def _generate_obs_res(self, analysis_type, grism_type):
		'''
		Generates all observation result files for all GRISMs and analysis types.
		'''
		file_list = '' 
		with open(os.path.join(self.emir_path, '%s_obs_res_%s.yaml'%(analysis_type,grism_type)), 'w') as f:
			try:
				for i, filename in enumerate(getattr(self, "%s_files"%analysis_type)[grism_type]):
					fpath = os.path.join(self.obs_dir, analysis_type, filename)
					if not os.path.exists(fpath): 
						print("WARNING: skipping file %s"%filename)
						continue

					if analysis_type == 'arc':
						file_list += " - %s\n"%filename
					else:
						file_list = " - %s\n"%filename
						text = OBS_RES_TEMPLATE%(self.name+"_"+filename.split('-')[0], file_list)
						f.write(text)
						if i < len(getattr(self, "%s_files"%analysis_type)[grism_type])-1:
							f.write('---\n')
				
				if analysis_type == 'arc':
					text = OBS_RES_TEMPLATE%(self.name+'_arc', file_list)
					f.write(text)
				
			except KeyError:
				#raise ValueError("attempting to generate an observation result file without proper files")
				pass

			f.close()

	def rectify_and_analyze(self, analysis_type, grism_type):
		'''
		Run the numina recipe for rectification and calibration.
		'''

		numina_bin = os.path.join(os.path.dirname(sys.executable), 'numina')
		result = subprocess.run(
			[numina_bin, 'run',
			 '%s_obs_res_%s.yaml' % (analysis_type, grism_type),
			 '-r', 'control_%s.yaml' % grism_type],
			capture_output=True, text=True, cwd=self.emir_path)

		with open(os.path.join(self.emir_path, "emir_out.txt"), 'w') as f:
			f.write(result.stdout+'\n\n')
			f.write(result.stderr+'\n\n')

	def ABBA_subtract(self, grism_type):

		self._prepare_ABBA(grism_type)
		self._subtract_and_save_ABBA(grism_type)

	def _prepare_ABBA(self, grism_type):
		'''
		Prepares the directory to save ABBA information, checks that 
		there is rectified and calibrated data.
		'''
		# if observation result file not generated, must run initialization
		obs_res_file = os.path.join(self.emir_path, "object_obs_res_%s.yaml"%grism_type)
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

		if not os.path.exists(self.abba_dir):
			os.mkdir(self.abba_dir)


	def _get_calibrated_data(self, grism_type): 
		
		A = [] ; B = []
		fname_base = os.path.join(self.emir_path,'obsid%s_%s_results/reduced_mos.fits')

		for row in self.qc.object_table:
			if row['FILTER'] != grism_type: continue

			fID = row['IMAGE'].split('-')[0]
			fname = fname_base%(self.name, fID)
			if not os.path.exists(os.path.join(self.emir_path,'obsid%s_%s_results/'%(self.name, fID)))\
			  and not os.path.exists(fname) \
			  and os.path.exists(os.path.join(self.obs_dir, 'object', row['IMAGE'])):
				print(os.path.join(self.obs_dir, 'object', row['IMAGE']))
				raise ValueError('EMIR did not run correctly, look at log file.')
			elif not os.path.exists(fname): continue

			if row['TELPOS'] == 'NOD_A':
				fi = fits.open(fname)
				A += [fi[0].data.astype(float)]
			elif row['TELPOS'] == 'NOD_B':
				fi = fits.open(fname)
				header = fi[0].header
				B += [fi[0].data.astype(float)]
			else:
				raise ValueError("TELPOS should be NOD_A or NOD_B not %s (observation %s)"%(row['TELPOS'],self.name))
			fi.close()

		if len(A) == 0 or len(B) == 0:
			err_msg = "%i files in A position, %i in B, %i files in total (observation %s)"%(len(A), len(B), len(self.qc.object_table), self.name)
			err_msg2 = "\n if you're getting this message, the quality control file might have the wrong file names"
			raise ValueError(err_msg+err_msg2)

		return A, B, header

	def _subtract_and_save_ABBA(self, grism_type, only_AB=False):
		'''
		Performs the actual subtraction and saves the results.
		'''
		A, B, header = self._get_calibrated_data(grism_type)

		if len(A) != 2 or len(B) != 2:
			self._subtract_and_save_non_ABBA(grism_type)

		fname = os.path.join(self.abba_dir, "ABBA_subtracted_%s.fits"%grism_type)

		primary_hdu = fits.ImageHDU(data=A[0] + A[1] - B[0] - B[1], header=header)
		primary_hdu.writeto(fname, overwrite=True)


	def _subtract_and_save_non_ABBA(self, grism_type):

		A, B, header = self._get_calibrated_data(grism_type)
		fname = os.path.join(self.abba_dir, "ABBA_subtracted_%s.fits"%grism_type)
		
		final_image = np.zeros(A[0].shape)
		for a,b in zip(A,B):
			final_image += (a-b)

		primary_hdu = fits.ImageHDU(data=final_image, header=header)
		primary_hdu.writeto(fname, overwrite=True)

	###### RESPONSE CURVE FUNCTIONS #############

	def make_sens_function(self, grism_type, **kwargs):
		'''Makes a sensitivity curve for the CCD.'''

		tabulated_spectrum_path = os.path.join(os.environ['EMIR_PIPE'],
												"extra_files/uka0v.fits")
		if not os.path.exists(self.sens_function_dir): os.mkdir(self.sens_function_dir)

		ABBA_file = os.path.join(self.abba_dir,'ABBA_subtracted_%s.fits'%grism_type)
		counts, wave, header = get_spectrum(ABBA_file, grism_type, 
											kwargs.get('SN_position', None),
											kwargs.get('raw_debug', False))

		self._make_prelim_sens_fn(counts, wave, header, grism_type, **kwargs)
		
	def get_sens_function(self, grism_type):
		
		fname = self.sens_function_dir+'/%s_sens_function.fits'%grism_type
		wave, flux = fits_to_data(fname)

		return wave, flux


	def _make_prelim_sens_fn(self, counts, wave, header, grism_type, **kwargs):

		tabulated_spectrum_path = os.path.join(os.environ['EMIR_PIPE'],
                                                "extra_files/uka0v.fits")
		counts = self._do_abba_fitting(counts, wave, grism_type, header, **kwargs)
		tabulated_flux = self._get_tablulated_flux(wave, **kwargs)
		
		sens = tabulated_flux / counts

		telluric_corr = self._add_telluric_correction(counts, grism_type, **kwargs)
		sens = sens / telluric_corr

		if 'prelim_smooth_radius' in kwargs:
			sens = smooth(sens, radius=kwargs['prelim_smooth_radius']) 

		hdu_sp = fits.PrimaryHDU(data=sens, header=header)
		hdu_sp.writeto(self.sens_function_dir+'/%s_sens_function.fits'%grism_type, overwrite=True)

		return sens

	def _do_abba_fitting(self, counts, wave, grism_type, header, **kwargs):

		if 'sample_points' in kwargs:
			sample = []
			for point in kwargs['sample_points']:
				sample += [np.where(np.abs(wave-point) == min(np.abs(wave-point)))[0][0]]

			sample = np.array(sample)
		else:
			spacing = kwargs.get('sample_spacing', len(counts)//25)
			sample = np.array([i for i in range(len(counts)) if i%spacing==0])

		sampled_counts = counts[sample]
		sampled_wave = wave[sample]

		tck = splrep(sampled_wave, sampled_counts)
		counts = splev(wave, tck)
	
		hdu_sp = fits.PrimaryHDU(data=counts, header=header)
		hdu_sp.writeto(self.sens_function_dir+'/%s_splined_sens.fits'%grism_type, overwrite=True)

		return counts

	def _get_tablulated_flux(self, wave, **kwargs):
		
		tabulated_spectrum_path = os.path.join(os.environ['EMIR_PIPE'],
                                                "extra_files/uka0v.fits")
		tabulated_interpolation = self._get_tabulated_interpolation(tabulated_spectrum_path)
		tabulated_flux = tabulated_interpolation(wave)

		if 'sample_points' not in kwargs:
			spacing = kwargs.get('sample_spacing', len(wave)//25)
			sample = np.array([i for i in range(len(wave)) if i%spacing==0])
		else:
			sample = []
			for point in kwargs['sample_points']:
				sample += [np.where(np.abs(wave-point) == min(np.abs(wave-point)))[0][0]]

			sample = np.array(sample)

		sampled_flux = tabulated_interpolation(wave)[sample]
		sampled_wave = wave[sample]

		tck = splrep(sampled_wave, sampled_flux)
		tabulated_flux = splev(wave, tck)

		return tabulated_flux

	def _get_tabulated_interpolation(self, res_curve_path):

		fits_file = fits.open(res_curve_path)
		f = fits_file[0]

		flux_tab = f.data
		wavelengths_tab = pixels_to_wavelength(f.header)

		fits_file.close()

		return interp1d(wavelengths_tab, flux_tab, kind="linear")


	def _add_telluric_correction(self, abba_fit, grism_type, **kwargs):

		wave, flux = self.get_raw_spectrum(grism_type)
		wave, flux = self._remove_hydrogen(flux, wave, grism_type)

		telluric_correction = flux/abba_fit

		return telluric_correction

	def _remove_hydrogen(self, flux, wave, grism_type, widths=None, in_angstrom=False):
	
		all_lines, widths = get_hydrogen_lines(widths)
		if len(widths) != len(all_lines): raise ValueError('Widths length should be %i '%len(all_lines)+\
															'(current length: %i)'%(len(widths)))
		for l, w in zip(all_lines, widths):
			idx = np.where(np.abs(wave-l)==min(np.abs(wave-l)))[0][0]
			if idx != 0 and idx != len(wave)-1:
				x_vals = [wave[idx-w], wave[idx+w]]
				y_vals = [flux[idx-w], flux[idx+w]]
				interp = interp1d(x_vals, y_vals, kind="linear")

				x = wave[idx-w:idx+w]
				flux[idx-w:idx+w] = interp(x)

		return wave, flux

	##### FINAL REDUCTION FUNCTIONS #######

	def get_reduced_spectrum(self, tabulated_magnitude, calibration_obs, grism_type, **kwargs):
		
		wave, sens_fn = calibration_obs.get_sens_function(grism_type)
		wave, raw_calib_data = calibration_obs.get_raw_spectrum(grism_type, **get_parameters(calibration_obs, grism_type))
		wave, raw_data = self.get_raw_spectrum(grism_type, **kwargs)
		flux_scaling = self._calculate_flux_scaling(tabulated_magnitude, wave, raw_calib_data*sens_fn, grism_type, **kwargs)

		final_spectrum = flux_scaling*raw_data*sens_fn
		
		final_sps = getattr(self, 'final_sps', None)
		final_wvs = getattr(self, 'final_wvs', None)
		if final_sps is None: final_sps = {}
		if final_wvs is None: final_wvs = {}

		final_sps[grism_type] = final_spectrum
		final_wvs[grism_type] = wave
		setattr(self, 'final_sps', final_sps)
		setattr(self, 'final_wvs', final_wvs)

		return wave, final_spectrum

	def _calculate_flux_scaling(self, tabulated_magnitude, wave, SN_data, grism_type, **kwargs):

		overlap_mask, resp_flux = self._get_filter_response_values(wave, grism_type)
		return self._flux_scale(SN_data[overlap_mask], tabulated_magnitude, resp_flux, wave[overlap_mask], grism_type)
		
	def _get_filter_response_values(self, wave, grism_type):	
		wave_overlap = {} ; flux = {} 
		fresp_dir = os.getenv('EMIR_PIPE') + '/filter_response_functions/'
		for fil in list(grism_type):
			if fil == 'Y': continue
			fname = fil+'.txt'
			fr_wave = [] ; fr_flux = []
			with open(fresp_dir+fname, 'r') as fh:
				for line in fh:
					if not line[0].isdigit(): continue
					w, fr = line.split()
					fr_wave += [float(w)*10000]
					fr_flux += [float(fr)]
    
			interp = interp1d(fr_wave, fr_flux, kind="linear")
			data_overlap = (wave>min(fr_wave)) & (wave<max(fr_wave))

			flux[fil] = interp(wave[data_overlap])

		return data_overlap, flux

	def _flux_scale(self, data, tabulated_mag, flux, wave, filter_type):

		f = filter_type[-1] # either J or K

		area_f = np.sum(flux[f])*(wave[1]-wave[0]) # total area of the filter
		area_c = np.sum(data*flux[f])*(wave[1]-wave[0]) # data convolved with filter

		m_zp = (2.5*np.log10(np.abs(area_c/area_f)) 
				- 2.5*np.log10(zp_fluxes_vega[f]) 
				+ tabulated_mag[f] )
		f_zp = 10**(-m_zp/2.5) 

		return f_zp

	def get_raw_spectrum(self, grism_type, **kwargs):
		raw_data, wave, header = get_spectrum(self.abba_dir+'/ABBA_subtracted_%s.fits'%grism_type, grism_type, 
											  SN_position=kwargs.get('SN_position', None), debug=kwargs.get('raw_debug', False))
		return wave, raw_data

	def save_spectrum(self, fname):
	
		final_sps = getattr(self, 'final_sps', None)
		final_wvs = getattr(self, 'final_wvs', None)

		if final_sps is None: 
			raise ValueError("No final spectra calculated yet. Run get_reduced_spectra to calculate.")

		self._write_file(fname, final_wvs, final_sps)

	def _write_file(self, fname, final_wvs, final_sps):

		hh = get_header(self, 1)

		grism_types = list(final_sps.keys())
		file_info = "# TELESCOPE: GTC \n# INSTRUMENT: EMIR \n# GRISM: %s \n# DATE OBS: %s\n# N ABBA CYCLES: %s\n# EXP TIME PER EXPOSURE: %s\n# TOTAL EXP TIME: %s\n# AIRMASS: %.08F\n"
		fi = file_info%(str(grism_types), 
                        hh['DATE-OBS'],
						str([len(self.object_files[gr])//4 for gr in grism_types]),
                        str([hh['EXPTIME']]*len(grism_types)),
                        str([hh['EXPTIME']*len(self.object_files[gr]) for gr in grism_types]),
                        hh['AIRMASS']
                        )

		with open(fname, 'w') as f:
			f.write(fi)
			for wk,sk in zip(final_wvs, final_sps):
				wave, sp = final_wvs[wk], final_sps[sk]
				for w, s in zip(wave, sp):
					if np.isnan(s): continue
					f.write('%E\t%E\n'%(w,s))
	





