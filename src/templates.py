
# numina 0.35.2 requires products nested under instrument profile UUIDs.
# Both known EMIR configs are listed so the same control file works for
# data taken before and after the 2023 detector upgrade.
_EMIR_PRODUCTS = "\
    - {id: 2, type: 'MasterBadPixelMask', tags: {}, content: 'master_bpm.fits'}\n\
    - {id: 3, type: 'MasterDark', tags: {}, content: 'master_dark_zeros.fits'}\n\
    - {id: 4, type: 'MasterIntensityFlat', tags: {}, content: '%s'}\n\
    - {id: 5, type: 'MasterSpectralFlat', tags: {}, content: '%s'}\n\
    - {id: 11, type: 'MasterRectWave', tags: {grism: J, filter: J}, content: 'rect_wpoly_MOSlibrary_grism_J_filter_J.json'}\n\
    - {id: 12, type: 'MasterRectWave', tags: {grism: H, filter: H}, content: 'rect_wpoly_MOSlibrary_grism_H_filter_H.json'}\n\
    - {id: 13, type: 'MasterRectWave', tags: {grism: K, filter: Ksp}, content: 'rect_wpoly_MOSlibrary_grism_K_filter_Ksp.json'}\n\
    - {id: 14, type: 'MasterRectWave', tags: {grism: LR, filter: YJ}, content: 'rect_wpoly_MOSlibrary_grism_LR_filter_YJ.json'}\n\
    - {id: 15, type: 'MasterRectWave', tags: {grism: LR, filter: HK}, content: 'rect_wpoly_MOSlibrary_grism_LR_filter_HK.json'}\n"

CONTROL_YAML_TEMPLATE = ("version: 1\n"
"products:\n"
"  EMIR:\n"
"    225fcaf2-7f6f-49cc-972a-70fd0aee8e96:\n"  # original config (2016–2023-03)
+ _EMIR_PRODUCTS
+ "    443fc0d1-e09a-48cc-a0fd-02be6f399da2:\n"  # updated config (2023-07+)
+ _EMIR_PRODUCTS
+ "requirements:\n"
"  EMIR:\n"
"    default:\n"
"     FULL_DITHERED_IMAGE:\n"
"      - {name: 'x_offsets', tags: {}, content: 'ref_object_pos.txt'}\n")


OBS_RES_TEMPLATE = "id: %s\n\
instrument: EMIR\n\
mode: GENERATE_RECTWV_COEFF\n\
frames:\n\
%s\
enabled: True\n"
