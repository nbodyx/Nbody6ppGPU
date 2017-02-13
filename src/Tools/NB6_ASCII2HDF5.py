#
# Convert NBODY6++GPU ASCII output to HDF5 (Cai et al. 2015, ApJS, 219, 31).
#
# Author: Maxwell Xu CAI (NAOC/KIAA)
# Feedback/Bug report: maxwell.x.cai@gmail.com
#

import h5py
import os
import numpy as np
import glob


working_dir = os.path.dirname(os.path.realpath(__file__))
f_singles = sorted(glob.glob(os.path.join(working_dir, 'single.40_*')))

col_dtypes_s = {'names': ('NAM', 'M', 'X1', 'X2', 'X3', 'V1', 'V2', 'V3', 'RS', 'L', 'TE', 'MC', 'RC', 'KW'),
                'formats': ('int', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'int')}
col_dtypes_b = {'names': ('NAM1', 'NAM2', 'NAMC', 'M1', 'M2', 'XC1', 'XC2', 'XC3', 'VC1', 'VC2', 'VC3', 'XR1', 'XR2', 'XR3', 'VR1', 'VR2', 'VR3', 'A', 'ECC', 'P', 'RS1', 'RS2', 'L1', 'L2', 'TE1', 'TE2', 'MC1', 'MC2', 'RC1', 'RC2', 'KW1', 'KW2', 'KWC'),
                'formats': ('int', 'int', 'int', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'int', 'int', 'int')}
col_dtypes_m = {'names':('NAM1', 'NAM2', 'NAM3', 'NAMC', 'M1', 'M2', 'M3', 'XC1', 'XC2', 'XC3', 'VC1', 'VC2', 'VC3', 'XR01', 'XR02', 'XR03', 'VR01', 'VR02', 'VR03', 'XR11', 'XR12', 'XR13', 'VR11', 'VR12', 'VR13', 'A0', 'ECC0', 'P0', 'A1', 'ECC1', 'P1', 'KW1', 'KW2', 'KW3', 'KWC', 'RS1', 'RS2', 'RS3', 'L1', 'L2', 'L3', 'TE1', 'TE2', 'TE3', 'MC1', 'MC2', 'MC3', 'RC1', 'RC2', 'RC3'),
                'formats': ('int', 'int', 'int', 'int', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'int', 'int', 'int', 'int', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float')}

h5f = None
for f_single_id in range(len(f_singles)):
    f_single = f_singles[f_single_id]
    t_current_str = f_single.split('single.40_')[1]
    t_current = float(t_current_str)
    t_step = np.floor(t_current)
    t_substep = t_current - t_step
    f_binary = os.path.join(working_dir, ('binary.40_%s' % t_current_str))
    f_merger = os.path.join(working_dir, ('merger.40_%s' % t_current_str))
    if t_substep == 0.0:
        # create new file
        if h5f is not None:
            h5f.close()  # close previously opened file
        # create new h5 file
        h5f = h5py.File(os.path.join(working_dir, ('snap.40_%d.h5part' % t_step)), 'w')
        substep_id = 0
    else:
        substep_id += 1
    print 'Processing, t = %f, Step#%d, h5f = %s' % (t_current, substep_id, ('snap.40_%d.h5part' % t_step))
    h5g = h5f.create_group('Step#%d' % substep_id)
    h5g_bin = h5g.create_group('Binaries')
    h5g_mer = h5g.create_group('Mergers')
    data_s = np.loadtxt(f_single, dtype=col_dtypes_s, comments='#')
    for col_name in col_dtypes_s['names']:
        h5g.create_dataset(col_name, data=data_s[col_name])
    if os.path.isfile(f_binary):
        data_b = np.loadtxt(f_binary, dtype=col_dtypes_b, comments='#')
        for col_name in col_dtypes_b['names']:
            h5g_bin.create_dataset(col_name, data=data_b[col_name])
    if os.path.isfile(f_merger):
        data_m = np.loadtxt(f_merger, dtype=col_dtypes_m, comments='#')
        for col_name in col_dtypes_m['names']:
            h5g_mer.create_dataset(col_name, data=data_m[col_name])

if h5f is not None:
    h5f.close()  # close any previously opened h5 file


