
from scipy.interpolate import interp1d, splrep, splev
import matplotlib.pyplot as plt
import numpy as np

import mpl_interactions.ipyplot as iplt
from matplotlib.widgets import Slider, Button

from astropy.io import fits
import numpy as np
import os
import pickle

from Observation import Observation
from reduction_tools import *

def get_parameters(o, filter_type):

	if not os.path.exists(o.result_dir+"/interactive_parameters.pkl"): return {}
	
	kwa = {}
	with open(o.result_dir+"/interactive_parameters.pkl", 'rb') as f:
		d = pickle.load(f)
		if 'SN_positions' in d and filter_type in d['SN_positions']:
			kwa['SN_position'] = d['SN_positions'][filter_type]
		if 'sample_points' in d and filter_type in d['sample_points']:
			kwa['sample_points'] = d['sample_points'][filter_type]
	
	return kwa

def interactive_SN_position_finder(o, filter_type):

	kwa = get_parameters(o, filter_type) 
	raw_wave, raw = o.get_raw_spectrum(filter_type)

	abba_fname = os.path.join(o.abba_dir, 'ABBA_subtracted_%s.fits'%filter_type)
	with fits.open(abba_fname) as f:
		abba = f[1].data
		abba_low = f[1].data[:,:300]
		abba_high = f[1].data[:,-300:]

	fig, axs = plt.subplots(2, 1)

	# plot histogram of pixel intensities
	y_flat = np.sum(abba, axis=1)
	y_flat_low = np.sum(abba_low, axis=1)
	y_flat_high = np.sum(abba_high, axis=1)

	low = 1000 ; high = 1500
	axs[1].plot(y_flat, label='full sum')
	axs[1].plot(y_flat_low, label='left side of image')
	axs[1].plot(y_flat_high, label='right side of image')

	axs[1].set_xlim((low, high))
	axs[1].set_ylim((min(y_flat[low:high])*0.75, max(y_flat[low:high])*1.25))

	# create interactive controls
	def f(x, low, high):
		kwargs = {'SN_position': (high, low)}
		w, raw = o.get_raw_spectrum(filter_type, **kwargs)
		return raw

	idx_range = range(low, high)
	
	plt.subplots_adjust(bottom=0.3, top=0.95)

	high_ax = plt.axes([0.45, 0.15, 0.4, 0.03])
	high_slider = Slider(high_ax, label="A (positive peak) position", valmin=low, valmax=high, valstep=idx_range, valinit=kwa.get('SN_position', [1])[0])
	low_ax = plt.axes([0.45, 0.10, 0.4, 0.03])
	low_slider = Slider(low_ax, label="B (negative peak) position", valmin=low, valmax=high, valstep=idx_range, valinit=kwa.get('SN_position', [1,1])[1])

	ax = plt.axes([0.01, 0.01, 0.1, 0.075])
	button = Button(ax, 'Finalize',color="yellow")
	def close_and_save(val):
		fname = o.result_dir+'/interactive_parameters.pkl'
		if os.path.exists(fname):
			with open(o.result_dir+'/interactive_parameters.pkl', 'rb') as f:
				d = pickle.load(f)
		else:
			d = {}

		if d.get('SN_positions') is None: d['SN_positions'] = {}
		d['SN_positions'][filter_type] = (ctrls.params['high'], ctrls.params['low'])

		with open(fname, 'wb') as f:
			pickle.dump(d, f)

		plt.close()
	button.on_clicked(close_and_save)
	
	ctrls = iplt.plot(raw_wave, f, low=low_slider, high=high_slider, ax=axs[0], label='raw data')
	iplt.axvline(ctrls["low"], ax=axs[1], c="k")
	iplt.axvline(ctrls["high"], ax=axs[1], c="k")
	[a.legend(loc='upper right') for a in axs]
	plt.show()
	return button, ctrls.params	

def interactive_response_curve_fit(o_calib, filter_type, tabulated_flux_path, num_points=7, **kwargs):

	with fits.open(tabulated_flux_path) as f:
		tab_data = f[0].data
		header = f[0].header
		tab_wave = pixels_to_wavelength(header)

	kwa = get_parameters(o_calib, filter_type)

	wave, raw = o_calib.get_raw_spectrum(filter_type, **(kwargs|kwa))
	
	interp = interp1d(wave, raw, kind="linear")
	tab_interp = interp1d(tab_wave, tab_data, kind="linear")

	fig, ax = plt.subplots(3+num_points, sharex=True, height_ratio=[num_points+2]*3+[1]*num_points)#, figsize=(20,15))
	#plt.subplots_adjust(bottom=0.35, top=0.9)
	
	button, sliders = init_widgets(wave, num_points, ax[3:])
	plot_raw_data(wave, raw, tab_interp, ax)

	if 'sample_points' in kwa:
		init_points = kwa['sample_points']
	else:
		init_points = [wave[0]+(wave[-1]-wave[0])*i/num_points for i in range(num_points)]
	pts, tab_pts, resp_pts = plot_sample_points(init_points, interp, tab_interp, ax)
	line, tab_line, resp_line = plot_spline_fits(wave, o_calib, tabulated_flux_path, filter_type, 
												 init_points, interp, tab_interp, ax, **(kwargs|kwa))
	[a.legend() for a in ax[:3]]

	def update(val):
		global x
		x = sorted([sliders[key].val for key in sliders])
		
		update_sample_points(x, interp, tab_interp, pts, tab_pts, resp_pts)
		flux, counts = update_splines(x, wave, interp, tab_interp, line, tab_line)

		kw = kwargs | {'sample_points': x}
		o_calib.make_response_curve(tabulated_flux_path, filter_type, **kw)
		resp_wave, response = o_calib.get_response_curve(filter_type)
		resp_line.set_xdata(resp_wave)
		resp_line.set_ydata(response)

	[sliders[key].on_changed(update) for key in sliders]

	def close_and_save(val):
		fname = o_calib.result_dir+'/interactive_parameters.pkl'
		if os.path.exists(fname):
			with open(o_calib.result_dir+'/interactive_parameters.pkl', 'rb') as f:
				d = pickle.load(f)
		else:
			d = {}

		if d.get('sample_points') is None: d['sample_points'] = {}
		d['sample_points'][filter_type] = sorted([sliders[key].val for key in sliders])

		with open(fname, 'wb') as f:
			pickle.dump(d, f)

		plt.close()
	button.on_clicked(close_and_save)
	
	plt.show()
	return button, sliders

def init_widgets(wave, num_points, slider_axs):
	#slider_axs = [plt.axes([0.15, i+0.05, 0.65, 0.03]) 
    # 	          for i in np.linspace(0,0.1, num_points)]
	sliders = {'point %i'%i: Slider(slider_axs[i], label="", 
    	              valmin=wave[0], valmax=wave[-1], 
        	          valinit=wave[0]+((wave[-1]-wave[0])*i/num_points)) 
           		for i in range(num_points)}

	ax = plt.axes([0.01, 0.01, 0.1, 0.075])
	cl_b = Button(ax, 'Finalize',color="yellow")

	return cl_b, sliders

def plot_raw_data(wave, raw, tab_interp, ax):
	ax[0].plot(wave, raw, color=plt.cm.Set1(1), label='Raw Spectral Data')
	ax[1].plot(wave, tab_interp(wave), color=plt.cm.Set1(1), label='Tabulated Data')

def plot_sample_points(init_points, interp, tab_interp, ax):
	pts, = ax[0].plot(init_points, interp(init_points), 'X', color=plt.cm.Set1(0))
	tab_pts, = ax[1].plot(init_points, tab_interp(init_points), 'X', color=plt.cm.Set1(0))
	resp_pts, = ax[2].plot(init_points, tab_interp(init_points)/interp(init_points), 'X', color=plt.cm.Set1(0))

	return pts, tab_pts, resp_pts

def update_sample_points(x, interp, tab_interp, pts, tab_pts, resp_pts):
	pts.set_xdata(x)
	tab_pts.set_xdata(x)
	resp_pts.set_xdata(x)

	pts.set_ydata(interp(x))
	tab_pts.set_ydata(tab_interp(x))
	resp_pts.set_ydata(tab_interp(x)/interp(x))

def plot_spline_fits(wave, o_calib, tabulated_flux_path, filter_type, init_points, dat_interp, tab_interp, ax, **kwargs):
	
	tck = splrep(init_points, dat_interp(init_points))
	counts = splev(wave, tck)
	tck = splrep(init_points, tab_interp(init_points))
	flux = splev(wave, tck)
	
	o_calib.make_response_curve(tabulated_flux_path, filter_type, **kwargs)
	resp_wave, response = o_calib.get_response_curve(filter_type)

	line, = ax[0].plot(wave, counts, color=plt.cm.Set1(0), label='spline fit')
	tab_line, = ax[1].plot(wave, flux, color=plt.cm.Set1(0), label='spline fit')
	resp_line, = ax[2].plot(resp_wave, response, label='response curve')
	
	return line, tab_line, resp_line

def update_splines(x, wave, interp, tab_interp, line, tab_line):
	flux = get_interpolation(x, wave, tab_interp)
	counts = get_interpolation(x, wave, interp)
	line.set_ydata(counts)
	tab_line.set_ydata(flux)
	return flux, counts	

def get_interpolation(x, wave, interpolation):
	tck = splrep(x, interpolation(x))
	y = splev(wave, tck)
	return y
