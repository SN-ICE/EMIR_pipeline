

CONTROL_YAML_TEMPLATE = "version: 1 \n\
products:\n\
  EMIR:\n\
  - {id: 2, type: 'MasterBadPixelMask', tags: {}, content: 'master_bpm.fits'}\n\
  - {id: 3, type: 'MasterDark', tags: {}, content: 'master_dark_zeros.fits'}\n\
  - {id: 4, type: 'MasterIntensityFlat', tags: {}, content: '%s'}\n\
  - {id: 5, type: 'MasterSpectralFlat', tags: {}, content: '%s'}\n\
  - {id: 11, type: 'MasterRectWave', tags: {grism: J, filter: J}, content: 'rect_wpoly_MOSlibrary_grism_J_filter_J.json'}\n\
  - {id: 12, type: 'MasterRectWave', tags: {grism: H, filter: H}, content: 'rect_wpoly_MOSlibrary_grism_H_filter_H.json'}\n\
  - {id: 13, type: 'MasterRectWave', tags: {grism: K, filter: Ksp}, content: 'rect_wpoly_MOSlibrary_grism_K_filter_Ksp.json'}\n\
  - {id: 14, type: 'MasterRectWave', tags: {grism: LR, filter: YJ}, content: 'rect_wpoly_MOSlibrary_grism_LR_filter_YJ.json'}\n\
  - {id: 15, type: 'MasterRectWave', tags: {grism: LR, filter: HK}, content: 'rect_wpoly_MOSlibrary_grism_LR_filter_HK.json'}\n\
  - {id: 21, type: 'RefinedBoundaryModelParam', tags: {grism: J, filter: J}, content: 'final_multislit_bound_param_grism_J_filter_J.json'}\n\
  - {id: 22, type: 'RefinedBoundaryModelParam', tags: {grism: H, filter: H}, content: 'final_multislit_bound_param_grism_H_filter_H.json'}\n\
  - {id: 23, type: 'RefinedBoundaryModelParam', tags: {grism: K, filter: Ksp}, content: 'final_multislit_bound_param_grism_K_filter_Ksp.json'}\n\
  - {id: 24, type: 'RefinedBoundaryModelParam', tags: {grism: LR, filter: YJ}, content: 'final_multislit_bound_param_grism_LR_filter_YJ.json'}\n\
  - {id: 25, type: 'RefinedBoundaryModelParam', tags: {grism: LR, filter: HK}, content: 'final_multislit_bound_param_grism_LR_filter_HK.json'}\n\
requirements:\n\
  EMIR:\n\
    default:\n\
     FULL_DITHERED_IMAGE:\n\
      - {name: 'x_offsets', tags: {}, content: 'ref_object_pos.txt'}\n"



OBS_RES_TEMPLATE = "id: %s\n\
instrument: EMIR\n\
mode: GENERATE_RECTWV_COEFF\n\
frames:\n\
%s\
enabled: True\n"


