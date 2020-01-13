import numpy as np
import yaml

#================================================================

# Read in the options saved in the yaml file using the loader
with open('mesh_options.yaml') as fp:
    mesh_dict = yaml.load(fp, Loader = yaml.FullLoader)

#================================================================

# Include endpoints or not
if mesh_dict['include_endpoints']:
    # Create a mesh with the specified properties
    mesh = np.linspace(mesh_dict['x0'], mesh_dict['x1'], mesh_dict['nx'])
else:
    # Create a mesh with the specified properties
    mesh = np.linspace(mesh_dict['x0'], mesh_dict['x1'], mesh_dict['nx']+2)
    mesh = mesh[1:-1]

#================================================================

# Open the mesh file for writing
fp = open('mesh_files/xmesh%s' % (mesh_dict['filetype']), 'w+')

# Write all values to the file
fp.write('# Meshfile\n')
fp.write('Node ID, Value\n')
for k, x in enumerate(mesh):
    fp.write('%d, %.3f\n' % (k, x))

# Close the mesh file when finished
fp.close()

#================================================================
