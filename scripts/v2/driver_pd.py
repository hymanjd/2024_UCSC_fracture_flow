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

src_path = os.getcwd()
jobname = src_path + "/output_pd"
dfnFlow_file = os.getcwd() + '/fehmn.files'

DFN = DFNWORKS(jobname,
               dfnFlow_file=dfnFlow_file,
               flow_solver="FEHM",
               ncpu=8)

DFN.params['domainSize']['value'] = [10000.0, 10000.0, 10000.0]
DFN.params['h']['value'] = 100
DFN.params['disableFram']['value'] = True
DFN.params['keepIsolatedFractures']['value'] = True


## Layer boundaries - set to matrix perm 
mat_perm = 1.0e-15
DFN.add_user_fract(shape='rect',
                   radii=5000,
                   translation=[0, 0, 1000],
                   normal_vector=[0, 0, 1],
                   permeability=mat_perm)

DFN.add_user_fract(shape='rect',
                   radii=5000,
                   translation=[0, 0, 4000],
                   normal_vector=[0, 0, 1],
                   permeability=mat_perm)

# Individual fractures
DFN.add_user_fract(shape='rect',
                   radii=2000,
                   translation=[0, 0, 3000],
                   normal_vector=[1, 0, 0],
                   aperture=1.0e-4)

DFN.add_user_fract(shape='rect',
                   radii=2000,
                   translation=[0, 0, 3000],
                   normal_vector=[0, 1, 0],
                   aperture=1.0e-4)

DFN.make_working_directory(delete=True)
DFN.check_input()
DFN.create_network()
DFN.mesh_network()

DFN.map_to_continuum(l=1000, orl=3)
DFN.upscale(mat_perm=mat_perm, mat_por=0.01)


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

exit()
####
lagrit_script = """
## read in dfn 
read / reduced_mesh.inp / mo_dfn
# make element set named elt_fracs s.t. itetclr (matid <= 3)
eltset / elt_fracs /  itetclr / lt / 3 
## remove the layer fractures
rmpoint / element / eltset, get, elt_fracs 
# read in the octree mesh 
read / octree_dfn.inp / mo_udfm


cmo / addatt / mo_udfm / tag / vdouble / scalar / nnodes
cmo / setatt / mo_udfm / tag / 1 0 0 / 1
cmo / readatt / mo_udfm / tag / 1, 0, 0 / tag_frac.dat


intersect_elements / mo_udfm / mo_dfn / fractures
eltset / eltfrac / fractures / ge 1 
pset / pfrac /  eltset / eltfrac 

pset/ p1 / geom / xyz / 1, 0, 0 /  & 
    -6000 -6000 -5000  /  6000 6000 1000 /  0,0,0

pset/ p2 / geom / xyz / 1, 0, 0 /  & 
    -6000 -6000 1000  /  6000 6000 4000 /  0,0,0

pset/ p3 / geom / xyz / 1, 0, 0 /  & 
    -6000 -6000 4000  /  6000 6000 6000 /  0,0,0

cmo / addatt / mo_udfm / layer / vdouble / scalar / nnodes
cmo / setatt / mo_udfm / layer / 1 0 0 / 0
cmo / setatt / mo_udfm / layer / pset, get, p1 / 1
cmo / setatt / mo_udfm / layer / pset, get, p2 / 2
cmo / setatt / mo_udfm / layer / pset, get, p3 / 3

dump / tmp.inp / mo_udfm 

finish 
"""






with open("paint_layers.lgi", "w") as fp:
    fp.write(lagrit_script)
    fp.flush() 


subprocess.call('lagrit < paint_layers.lgi' ,shell = True)

exit()



####

material_id = np.genfromtxt("tag_frac.dat").astype(int)
num_nodes = len(material_id)
perm = np.genfromtxt("perm_fehm.dat", skip_header = 1)[:,-1]
print(perm)
exit()

porosity = np.zeros(num_nodes)

# load z values 
with open('full_mesh.uge') as fuge:
    header = fuge.readline().split()
    num_cells = int(header[1])
    print(num_cells)
    z = np.zeros(num_cells)
    for i in range(num_cells):
        line = fuge.readline().split()
        z[i] = float(line[-2])

layer_1_matrix_perm = 1e-16 
layer_2_matrix_perm = 1e-12 
layer_3_matrix_perm = 1e-16

fracture_perm = 1e-10

for i in range(num_cells):
    if material_id[i] == 1: ## Matrix Cell
        if  z[i] < 1000:
            perm[i] = layer_1_matrix_perm
        elif 1000 < z[i] < 4000:
            perm[i] = layer_2_matrix_perm
        elif 4000 < z[i]:
            perm[i] = layer_3_matrix_perm

    # elif material_id[i] == 2:
    #     perm[i] = fracture_perm

with open("perm_layer.dat", "a") as fperm:
    fperm.write("perm\n")

for i in range(1, num_nodes + 1):
    with open("perm_layer.dat", "a") as fperm:
        fperm.write(
            str(i) + " " + str(i) + " " + "1" + " " +
            str(perm[i - 1]) + " " + str(perm[i - 1]) + " " +
            str(perm[i - 1]) + "\n")

with open("perm_layer.dat", "a") as fperm:
    fperm.write("")


with open("rock_fehm.dat", "a") as g:
    g.write(str(i) + " " + str(i) + " " + "1" + " " + "2165." + " " +
        "931." + " " + str(por_var[i - 1]) + "\n")



np.savetxt("my_perm.dat", perm)
lagrit_script = """
read / octree_dfn.inp / mo1
cmo / addatt / mo1 / perm / vdouble / scalar / nnodes
cmo / setatt / mo1 / perm / 1 0 0 / 1
cmo / readatt / mo1 / perm / 1, 0, 0 / my_perm.dat 
dump / perm.inp / mo1 
finish
"""
with open("color_mesh.lgi", "w") as fp:
    fp.write(lagrit_script)


subprocess.call('lagrit < color_mesh.lgi', shell = True)      
        
os.symlink(DFN.jobname + os.sep + "dfnGen_output", 'dfnGen_output')


# DFN.correct_stor_file()
# os.symlink(dfnFlow_file, 'fehmn.files')
# DFN.fehm()


# DFN.zone2ex(zone_file='all')

# DFN.pflotran()
# DFN.parse_pflotran_vtk_python()
# DFN.pflotran_cleanup()
