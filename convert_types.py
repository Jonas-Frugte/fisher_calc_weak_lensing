#!/usr/bin/env python3
"""
Script to convert remaining char* type parameters to int in data_importer_new.pyx
"""

import re

# Read the file
with open('code/data_importer_new.pyx', 'r') as f:
    content = f.read()

# List of all remaining function pairs that need conversion
functions_to_convert = [
    ('lbs_omch2_m', 'lps_omch2_m'),
    ('lbs_ns_p', 'lps_ns_p'),
    ('lbs_ns_m', 'lps_ns_m'),
    ('lbs_As_p', 'lps_As_p'),
    ('lbs_As_m', 'lps_As_m'),
    ('lbs_tau_p', 'lps_tau_p'),
    ('lbs_tau_m', 'lps_tau_m'),
    ('lbs_mnu_p', 'lps_mnu_p'),
    ('lbs_mnu_m', 'lps_mnu_m'),
    ('lbs_w0_p', 'lps_w0_p'),
    ('lbs_w0_m', 'lps_w0_m'),
    ('lbs_logT_AGN_p', 'lps_logT_AGN_p'),
    ('lbs_logT_AGN_m', 'lps_logT_AGN_m'),
]

# Replace all remaining char* type1, char* type2, char* type3 patterns in function signatures
for lbs_func, lps_func in functions_to_convert:
    # Pattern for lbs function with 3 type params
    pattern_lbs = f'cdef double {lbs_func}\\(int k1, int k2, int k3, char\\* type1, char\\* type2, char\\* type3, int num_samples, bint pb_correction\\)'
    replacement_lbs = f'cdef double {lbs_func}(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction)'
    content = re.sub(pattern_lbs, replacement_lbs, content)
    
    # Pattern for lps function with 2 type params
    pattern_lps = f'cdef double {lps_func}\\(int l, char\\* type1, char\\* type2\\)'
    replacement_lps = f'cdef double {lps_func}(int l, int type1, int type2)'
    content = re.sub(pattern_lps, replacement_lps, content)

# Also convert lps_der and lbs_der signatures
content = re.sub(
    r'cdef double lps_der\(int k, char\* type1, char\* type2, char\* par\)',
    'cdef double lps_der(int k, int type1, int type2, char* par)',
    content
)

content = re.sub(
    r'cdef double lbs_der\(int k1, int k2, int k3, char\* type1, char\* type2, char\* type3, int num_samples, bint pb_correction, char\* par\)',
    'cdef double lbs_der(int k1, int k2, int k3, int type1, int type2, int type3, int num_samples, bint pb_correction, char* par)',
    content
)

# Write back
with open('code/data_importer_new.pyx', 'w') as f:
    f.write(content)

print("Conversion complete!")
