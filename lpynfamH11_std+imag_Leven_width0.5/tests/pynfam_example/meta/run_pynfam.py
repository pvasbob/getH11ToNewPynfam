#!/usr/bin/env python
from __future__ import unicode_literals
from pynfam import pynfam_mpi_calc
import numpy as np


pynfam_inputs = {
 'directories': {
     'outputs' : 'pynfam_example',
     'exes'    : './exes',
     'scratch' : './tests'
     },

 'nr_parallel_calcs': None,

 'rerun_mode': 0,

 'hfb_mode': {
     'gs_def_scan'    : (+2, (0.225,)),
     'dripline_mode'  : 0,
     'ignore_nonconv' : 0,
     },

 'fam_mode': {
     'fam_contour': 'None',
     'beta_type'  : '-',
     'fam_ops'    : '1+', 

     }
}

override_settings = {
 'hfb' : {'proton_number'    : 60, 
          'neutron_number'   : 90,
          'number_of_shells' : 16,
          'number_gauss'     : 40,
          'number_laguerre'  : 40,
          'number_legendre'  : 80,
          'functional'       : 'SKM*',
          'force_parity'     : False,
          },


 'ctr' : {'energy_min' : 2.6,
          'energy_max' : 22.6,
          'half_width' : 0.5,
          'nr_points'  : 201,

     },

 'fam' : {'interaction_name' : 'SKM*',
          'vpair_t0'         : 0.0,
          'convergence_epsilon' : 1e-2,
          'max_iter'         : 100,
          'require_gauge_invariance': False,
          'override_cs0'        : '',
          'override_csr'        : '',
          'override_cds'        : '',
          'override_ct'         : '',
          'override_cgs'        : '',
          'override_cf'         : '',
 
     },

 'psi' : {}
}


if __name__ == '__main__':
    pynfam_mpi_calc(pynfam_inputs, override_settings, check=False)
