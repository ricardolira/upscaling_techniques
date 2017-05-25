# coding=utf-8
import numpy as np
from pymoab import core
from pymoab import types
from pymoab import topo_util


def get_block_by_ijk(i, j, k, n_i, n_j, n_k):
    """
    track down the block from its (i,j,k) position. Selects more than 1 (in z-direction) if 'k'
    is set >0

    example:
    --------
        get_block_by_ijk(0, 0, 0, 2, 2, 2)

        For a 2x2x2 mesh it will return vertical (with respect to 'z') blcks from i = 0, j = 0
        and all z from 0 to 2

        get_block_by_ijk(3, 2, 5, 3, 30, 30)

        For a 30x30x30 mesh it will return vertical (with respect to 'z') blcks from i = 3, j = 2
        and all z from 5 to 30

    It can simulate a vertical well that is completed in all 'z' selected

    """
    block = np.arange(k-1,n_k)*n_i*n_j+((i-1)+(j-1)*n_i)

    return block
# Parâmetros
#outfile_h5m = "multiscale.h5m"
outfile_vtk = "multiscale.vtk"

#Get multiscale parameters
coarse_ratio_x = 4
coarse_ratio_y = 6
coarse_ratio_z = 7

mesh_size_x = 18
mesh_size_y = 22
mesh_size_z = 15

mesh_size_x_coarse = mesh_size_x/coarse_ratio_x
mesh_size_y_coarse = mesh_size_y/coarse_ratio_y
mesh_size_z_coarse = mesh_size_z/coarse_ratio_z

block_size_x, block_size_y ,block_size_z = 1, 1, 1  # (ft)
volume_dims = (block_size_x, block_size_y ,block_size_z)

coarse_ids_x = [i // (coarse_ratio_x) for i in xrange(mesh_size_x)]
coarse_ids_y = [i // (coarse_ratio_y) for i in xrange(mesh_size_y)]
coarse_ids_z = [i // (coarse_ratio_z) for i in xrange(mesh_size_z)]
print coarse_ids_x, coarse_ids_y, coarse_ids_z

new_coarse_x = coarse_ids_x[mesh_size_x//coarse_ratio_x*coarse_ratio_x:]
if len(new_coarse_x) < mesh_size_x//coarse_ratio_x:
    new_coarse_x = np.repeat(max(coarse_ids_x)-1,len(new_coarse_x)).tolist()
    coarse_ids_x = coarse_ids_x[:mesh_size_x//coarse_ratio_x*coarse_ratio_x]+new_coarse_x
else:
    coarse_ids_x = coarse_ids_x

new_coarse_y = coarse_ids_y[mesh_size_y//coarse_ratio_y*coarse_ratio_y:]
if len(new_coarse_y) < mesh_size_y//coarse_ratio_y:
    new_coarse_y = np.repeat(max(coarse_ids_y)-1,len(new_coarse_y)).tolist()
    coarse_ids_y = coarse_ids_y[:mesh_size_y//coarse_ratio_y*coarse_ratio_y]+new_coarse_y
else:
    coarse_ids_y = coarse_ids_y

new_coarse_z = coarse_ids_z[mesh_size_z//coarse_ratio_z*coarse_ratio_z:]
if len(new_coarse_z) < mesh_size_z//coarse_ratio_y:
    new_coarse_z = np.repeat(max(coarse_ids_z)-1,len(new_coarse_z)).tolist()
    coarse_ids_z = coarse_ids_z[:mesh_size_z//coarse_ratio_z*coarse_ratio_z]+new_coarse_z
else:
    coarse_ids_y = coarse_ids_y

print coarse_ids_x, coarse_ids_y, coarse_ids_z

n_fine = mesh_size_x*mesh_size_y*mesh_size_z
n_coarse = (max(coarse_ids_x)+1)*(max(coarse_ids_y)+1)*(max(coarse_ids_z)+1)

max_mesh_size = max(mesh_size_z*block_size_z, mesh_size_y*block_size_y, mesh_size_x*block_size_x)

# MOAB Startup
mb = core.Core()
root_set = mb.get_root_set()
mesh_topo_util = topo_util.MeshTopoUtil(mb)


coords = np.array([
    (i, j, k) for k in (np.arange(mesh_size_z+1, dtype='float64')*block_size_z/max_mesh_size)
    for j in (np.arange(mesh_size_y+1, dtype='float64')*block_size_y/max_mesh_size)
    for i in (np.arange(mesh_size_x+1, dtype='float64')*block_size_x/max_mesh_size)
], dtype='float64')
verts = mb.create_vertices(coords.flatten())


# Tags
physical_tag = mb.tag_get_handle(
    "MATERIAL_SET", 1, types.MB_TYPE_INTEGER, True)

restriction_operator = mb.tag_get_handle(
    "RESTRICTION_OPERATOR", 1, types.MB_TYPE_INTEGER, True)

gid_tag = mb.tag_get_handle(
    "GLOBAL_ID", 1, types.MB_TYPE_INTEGER, True)

phi_tag = mb.tag_get_handle(
    "phi", 1, types.MB_TYPE_DOUBLE, True)

coarse_phi_tag = mb.tag_get_handle(
    "coarse_phi", 1, types.MB_TYPE_DOUBLE, True)

perm_tag = mb.tag_get_handle(
    "perm", 9, types.MB_TYPE_DOUBLE, True)

coarse_perm_tag = mb.tag_get_handle(
    "coarse_perm", 9, types.MB_TYPE_DOUBLE, True)

abs_perm_x_tag = mb.tag_get_handle(
    "abs_perm_x", 1, types.MB_TYPE_DOUBLE, True)
abs_coarse_perm_x_tag = mb.tag_get_handle(
    "abs_coarse_perm_x", 1, types.MB_TYPE_DOUBLE, True)

#create_wells_by_meshset
injection_tag = mb.tag_get_handle(
    "injection_well", 1, types.MB_TYPE_INTEGER, True, types.MB_TAG_SPARSE)

production_tag = mb.tag_get_handle(
    "production_well", 1, types.MB_TYPE_INTEGER, True, types.MB_TAG_SPARSE)

injection_wells = {}
production_wells = {}
injection_wells[1] = mb.create_meshset()

production_wells[1] = mb.create_meshset()
production_wells[2] = mb.create_meshset()
production_wells[3] = mb.create_meshset()
production_wells[4] = mb.create_meshset()

phi_values = []
with open('spe_phi.dat') as phi:
    for line in phi:
        phi_values.extend(line.rstrip().split('        	'))
phi_values = [float(val) for val in phi_values]

perm_values = []
i = 0
with open('spe_perm.dat') as perm:
    for line in perm:
        line_list = line.rstrip().split('        	')
        if len(line_list) > 1:
            perm_values.extend(line_list)
perm_values = [float(val) for val in perm_values]
abs_perm_x = perm_values[0:mesh_size_x*mesh_size_y*mesh_size_z]

#Primal grid generation

cur_id = 0
elems = []
primals = {}

#Create fine grid
for k, idz in zip(xrange(mesh_size_z), coarse_ids_z):
    print round(float(k+1)/float(mesh_size_z)*100,2), ' %'
    for j, idy in zip(xrange(mesh_size_y), coarse_ids_y):
        for i, idx in zip(xrange(mesh_size_x), coarse_ids_x):

            hexa = [verts[(i)+(j*(mesh_size_x+1))+(k*((mesh_size_x+1)*(mesh_size_y+1)))], # (i, j, k)
                    verts[(i+1)+(j*(mesh_size_x+1))+(k*((mesh_size_x+1)*(mesh_size_y+1)))], # (i+1, j, k)
                    verts[(i+1)+(j+1)*(mesh_size_x)+(j+1)+(k*((mesh_size_x+1)*(mesh_size_y+1)))], # (i+1, j+1, k)
                    verts[(i)+(j+1)*(mesh_size_x)+(j+1)+(k*((mesh_size_x+1)*(mesh_size_y+1)))], # (i, j+1, k)

                    verts[(i)+(j*(mesh_size_x+1))+((k+1)*((mesh_size_x+1)*(mesh_size_y+1)))], # (i, j, k+1)
                    verts[(i+1)+(j*(mesh_size_x+1))+((k+1)*((mesh_size_x+1)*(mesh_size_y+1)))], # (i+1, j, k+1)
                    verts[(i+1)+(j+1)*(mesh_size_x)+(j+1)+((k+1)*((mesh_size_x+1)*(mesh_size_y+1)))], # (i+1, j+1, k+1)
                    verts[(i)+(j+1)*(mesh_size_x)+(j+1)+((k+1)*((mesh_size_x+1)*(mesh_size_y+1)))]] # (i, j+1, k+1)

            el = mb.create_element(types.MBHEX, hexa)

#Show fine scale properties
            mb.tag_set_data(gid_tag, el, cur_id)
            mb.tag_set_data(phi_tag, el, phi_values[cur_id])
            mb.tag_set_data(perm_tag, el, [perm_values[cur_id], 0, 0,
                                           0, perm_values[cur_id+mesh_size_x*mesh_size_y*mesh_size_z], 0,
                                           0, 0, perm_values[cur_id+2*mesh_size_x*mesh_size_y*mesh_size_z]])
            mb.tag_set_data(abs_perm_x_tag,el,abs_perm_x[cur_id])
            cur_id += 1

            elems.append(el)

# Create primal coarse grid
            try:
                primal = primals[(idx, idy, idz)]
                mb.add_entities(primal, [el])
            except KeyError:
                primal = mb.create_meshset()
                primals[(idx, idy, idz)] = primal
                mb.add_entities(primal, [el])


#setting wells
elems = np.array(elems).reshape(-1, 1)

for el in elems[get_block_by_ijk(1, 1, 1, mesh_size_x, mesh_size_y, mesh_size_z)]:
    mb.add_entities(production_wells[1], el)
    mb.tag_set_data(production_tag, el, 1)

for el in elems[get_block_by_ijk(1, mesh_size_y, 1, mesh_size_x, mesh_size_y, mesh_size_z)]:
    mb.add_entities(production_wells[2], el)
    mb.tag_set_data(production_tag, el, 2)

for el in elems[get_block_by_ijk(mesh_size_x, mesh_size_y, 1, mesh_size_x, mesh_size_y, mesh_size_z)]:
    mb.add_entities(production_wells[3], el)
    mb.tag_set_data(production_tag, el, 3)

for el in elems[get_block_by_ijk(mesh_size_x, 1, 1, mesh_size_x, mesh_size_y, mesh_size_z)]:
    mb.add_entities(production_wells[4], el)
    mb.tag_set_data(production_tag, el, 4)

# Esse esta errado:
for el in elems[get_block_by_ijk(mesh_size_x//2, mesh_size_y//2, 1, mesh_size_x, mesh_size_y, mesh_size_z)]:
    mb.add_entities(injection_wells[1], el)
    mb.tag_set_data(injection_tag, el, 1)

coarse_perm={}
i, j, k = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
perm_xx = []
perm_yy = []
perm_zz = []

meshset_id = 0
restriction_op = np.zeros((n_coarse,n_fine))

min_coarse_ids = np.array([0, 0, 0])
max_coarse_ids = np.array([max(coarse_ids_x), max(coarse_ids_y), max(coarse_ids_z)])

#Get absolute permeability values for each direction
for primal_id, primal in primals.iteritems():
    #attribute physical properties from fine grid to meshset (or coarse grid)
    fine_elems_in_coarse = mb.get_entities_by_type(primal,types.MBHEX)
    elems_ids = mb.tag_get_data(gid_tag, fine_elems_in_coarse)
    fine_phi_values = mb.tag_get_data(phi_tag, fine_elems_in_coarse)
    coarse_mean_phi = fine_phi_values.mean()

    mb.tag_set_data(coarse_phi_tag, primal, coarse_mean_phi)
    mb.tag_set_data(coarse_phi_tag, fine_elems_in_coarse,
                    np.repeat(coarse_mean_phi, len(fine_elems_in_coarse)))

    fine_perm_values = mb.tag_get_data(perm_tag, fine_elems_in_coarse)
    coarse_perm.update({eid: tensor.reshape(3,3) for eid, tensor in
                        zip(fine_elems_in_coarse, fine_perm_values)
                        })

    perm_xx = [(np.dot(np.dot(coarse_perm[elem_id],i),i))
            for elem_id in fine_elems_in_coarse]
    #perm_xy = [(np.dot(np.dot(coarse_perm[elem_id],i),j))
    #        for elem_id in fine_elems_in_coarse]
    #perm_xz = [(np.dot(np.dot(coarse_perm[elem_id],i),k))
    #        for elem_id in fine_elems_in_coarse]

    #perm_yx = perm_xy
    perm_yy = [(np.dot(np.dot(coarse_perm[elem_id],j),j))
            for elem_id in fine_elems_in_coarse]
    #perm_yz = [(np.dot(np.dot(coarse_perm[elem_id],j),k))
    #        for elem_id in fine_elems_in_coarse]
    #perm_zx = perm_xz
    #perm_zy = perm_yz
    perm_zz = [(np.dot(np.dot(coarse_perm[elem_id],k),k))
            for elem_id in fine_elems_in_coarse]
#Take average value - upscaling with arithmetic mean
    coarse_perm_xx = np.mean(perm_xx)
    coarse_perm_yy = np.mean(perm_yy)
    coarse_perm_zz = np.mean(perm_zz)

    mb.tag_set_data(abs_coarse_perm_x_tag, primal, coarse_perm_xx)
    mb.tag_set_data(abs_coarse_perm_x_tag, fine_elems_in_coarse,
                    np.repeat(coarse_perm_xx, len(fine_elems_in_coarse)))
    #coarse_perm_tensor = [[coarse_perm_xx, coarse_perm_xy, coarse_perm_xz],
    #                      [coarse_perm_yx, coarse_perm_yy, coarse_perm_yz],
    #                      [coarse_perm_zx, coarse_perm_zy, coarse_perm_zz]
    #                      ]





# Escreve o arquivo
#mb.write_file(outfile_h5m)
mb.write_file(outfile_vtk)
