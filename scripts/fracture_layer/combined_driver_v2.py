#"""
#   :synopsis: Driver run file for TPL example
#   :version: 2.0
#   :maintainer: Jeffrey Hyman
#.. moduleauthor:: Jeffrey Hyman <jhyman@lanl.gov>
#"""

from pydfnworks import *
import os
import numpy as np

def translate_DFN(DFN, z0, offset = 0):
    for i in range(DFN.num_frac):
        DFN.centers[i] += [0,0,z0]

    for i in range(1,DFN.num_frac+1):
        key = f'fracture-{i}'
        for j in range(len(DFN.polygons[key])):
            DFN.polygons[key][j] += [0,0,z0]
        if offset > 0:
            new_key = f'fracture-{i + offset}'
            DFN.polygons[new_key] = DFN.polygons[key] 

    return DFN 


jobname = "combined_UDFM"
DFN = DFNWORKS(jobname,
               ncpu=12)
DFN.params['domainSize']['value'] = [10000.0, 10000.0, 10000.0]
DFN.params['h']['value'] = 100
DFN.params['disableFram']['value'] = True
DFN.params['keepIsolatedFractures']['value'] = True


# Add one family in layer #1
DFN.add_fracture_family(shape="ell",
                        distribution="tpl",
                        kappa=0.1,
                        theta=0.0,
                        phi=0.0,
                        alpha=1.2,
                        min_radius=100.0,
                        max_radius=1000.0,
                        p32=0.001,
                        hy_variable="aperture",
                        hy_function="constant",
                        hy_params={"mu": 1e-4})

DFN.h = 0.1
DFN.x_min = -5000
DFN.y_min = -5000
DFN.z_min = -5000
DFN.x_max = 5000
DFN.y_max = 5000
DFN.z_max = 5000

DFN.domain = {"x": 10000, "y": 10000, "z": 10000 }

src_dir = os.getcwd()



DFN.make_working_directory(delete=True)
DFN.check_input()

MIDDLE_DFN = DFNWORKS(pickle_file = f"{src_dir}/output_pd_layer_v2/middle_layer.pkl")

## combine DFN
DFN.num_frac = MIDDLE_DFN.num_frac  
DFN.centers = MIDDLE_DFN.centers
DFN.polygons = MIDDLE_DFN.polygons
DFN.normal_vectors = MIDDLE_DFN.normal_vectors

os.symlink(f"{src_dir}/output_pd_layer_v2/reduced_mesh.inp", "reduced_mesh.inp")


DFN.map_to_continuum(l = 1000, orl = 3)