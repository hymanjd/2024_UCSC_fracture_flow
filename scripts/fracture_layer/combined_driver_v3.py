#"""
#   :synopsis: Driver run file for TPL example
#   :version: 2.0
#   :maintainer: Jeffrey Hyman
#.. moduleauthor:: Jeffrey Hyman <jhyman@lanl.gov>
#"""

from pydfnworks import *
import os
import numpy as np
import subprocess
import pandas as pd 


jobname = "combined_UDFM"
DFN = DFNWORKS(jobname,
               ncpu=12)
DFN.params['domainSize']['value'] = [10000.0, 10000.0, 10000.0]
DFN.params['h']['value'] = 100
DFN.params['disableFram']['value'] = True
DFN.params['keepIsolatedFractures']['value'] = True


# DUMMY FRACTURE FAMILY, needed to trick code 
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

DFN.h = 1000 
DFN.x_min = -5000
DFN.y_min = -5000
DFN.z_min = -5000
DFN.x_max = 5000
DFN.y_max = 5000
DFN.z_max = 5000

DFN.domain = {"x": 10000, "y": 10000, "z": 10000 }

src_dir = os.getcwd()


# os.symlink('middle_layer/middle_layer.inp', 'middle_layer.inp')
# os.symlink('faults/faults.inp', 'faults.inp')


lagrit_script = """"
## prior to running you need to copy the reduced_mesh from the top & bottom DFN here (or symbolic link)
## also copy the *pkl files
# to run 
# lagrit < combine_mesh.lgi 

# read in mesh 1 
read / middle_layer.inp / mo_middle

# read in mesh 2 
read / faults.inp / mo_faults /

# combine mesh 1 and mesh 2 to make final mesh
addmesh / merge / mo_dfn / mo_middle / mo_faults

# write to file 
dump / combined_dfn.inp / mo_dfn 
dump / reduced_mesh.inp / mo_dfn 

finish 
"""

with open('combine_dfn.lgi', 'w') as fp:
    fp.write(lagrit_script)

import subprocess
subprocess.call('lagrit < combine_dfn.lgi', shell = True)


DFN.make_working_directory(delete=True)
DFN.check_input()

FAULT_DFN = DFNWORKS(pickle_file = f"{src_dir}/faults/faults.pkl")
MIDDLE_DFN = DFNWORKS(pickle_file = f"{src_dir}/middle_layer/middle_layer.pkl")

## combine DFN
## combine DFN
DFN.num_frac = FAULT_DFN.num_frac + MIDDLE_DFN.num_frac 
DFN.centers = np.concatenate((FAULT_DFN.centers, MIDDLE_DFN.centers))
DFN.polygons = FAULT_DFN.polygons.copy() 
DFN.polygons = DFN.polygons| MIDDLE_DFN.polygons
DFN.normal_vectors = np.concatenate((FAULT_DFN.normal_vectors, MIDDLE_DFN.normal_vectors))

os.symlink(f"{src_dir}/reduced_mesh.inp", "reduced_mesh.inp")

DFN.map_to_continuum(l = 1000, orl = 3)
DFN.upscale(mat_perm=1e-15, mat_por=0.01)

# load z values 
with open('octree_dfn.inp') as finp:
    header = finp.readline().split()
    num_nodes = int(header[0])
    print(num_nodes)
    x = np.zeros(num_nodes)
    y = np.zeros(num_nodes)
    z = np.zeros(num_nodes)
    for i in range(num_nodes):
        line = finp.readline().split()
        x[i] = float(line[1])
        y[i] = float(line[2])
        z[i] = float(line[3])
material_id = np.genfromtxt("tag_frac.dat").astype(int)

df = pd.DataFrame({'x': x, 'y': y, 'z': z, 'material': material_id, })
print(df)
df.to_pickle('octree_nodes.pkl')

lagrit_script = """
read / octree_dfn.inp / mo1
cmo / addatt / mo1 / frac_index / vdouble / scalar / nnodes
cmo / setatt / mo1 / frac_index / 1 0 0 / 1
cmo / readatt / mo1 / frac_index / 1, 0, 0 / tag_frac.dat 
dump / tmp.inp / mo1 
finish
"""
with open("color_mesh.lgi", "w") as fp:
    fp.write(lagrit_script)

subprocess.call('lagrit < color_mesh.lgi', shell = True)

