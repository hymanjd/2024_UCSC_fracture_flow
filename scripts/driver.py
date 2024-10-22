#"""
#   :synopsis: Driver run file for TPL example
#   :version: 2.0
#   :maintainer: Jeffrey Hyman
#.. moduleauthor:: Jeffrey Hyman <jhyman@lanl.gov>
#"""

from pydfnworks import *
import os

src_path = os.getcwd()
jobname = src_path + "/output"
dfnFlow_file = src_path + '/dfn_explicit_multi_material.in'

DFN = DFNWORKS(jobname,
               dfnFlow_file=dfnFlow_file,
               flow_solver="PFLOTRAN",
               ncpu=8)

DFN.params['domainSize']['value'] = [10000.0, 10000.0, 10000.0]
DFN.params['h']['value'] = 100
DFN.params['disableFram']['value'] = True
DFN.params['keepIsolatedFractures']['value'] = True


DFN.add_user_fract(shape='rect',
                   radii=5000,
                   translation=[0, 0, 1000],
                   normal_vector=[0, 0, 1],
                   permeability=1.0e-12)

DFN.add_user_fract(shape='rect',
                   radii=5000,
                   translation=[0, 0, 4000],
                   normal_vector=[0, 0, 1],
                   permeability=1.0e-12)

DFN.add_user_fract(shape='rect',
                   radii=2000,
                   translation=[0, 0, 3000],
                   normal_vector=[1, 0, 0],
                   permeability=1.0e-12)

DFN.add_user_fract(shape='rect',
                   radii=2000,
                   translation=[0, 0, 3000],
                   normal_vector=[0, 1, 0],
                   permeability=1.0e-12)

DFN.make_working_directory(delete=True)
DFN.check_input()
DFN.create_network()
DFN.mesh_network()

DFN.map_to_continuum(l=1000, orl=3)

# DFN.upscale(mat_perm=1e-15, mat_por=0.01)

# DFN.zone2ex(zone_file='all')

# DFN.pflotran()
# DFN.parse_pflotran_vtk_python()
# DFN.pflotran_cleanup()
