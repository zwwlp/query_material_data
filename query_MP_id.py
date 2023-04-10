from pymatgen import MPRester
import json
import os
from pydash import py_
import numpy as np
import pandas as pd

API_KEY = "MWmp5Bws0oq6S9h4"
FILEPATH = "data/mateial_project_stable/"
chunk_size=1000
"""
query all the data_id, and return the set of cif files
"""
with MPRester(API_KEY) as m:
    docs=m.query({},["structure","task_id","formation_energy_per_atom", "band_gap", "nelements", "spacegroup", "structure", "nsites", "total_magnetization","volume"])
#save all data as npz
    for n, chunk in enumerate(py_.chunk(docs, chunk_size)):
        file_path = os.path.join(FILEPATH, "mp_data_{:03d}.npz".format(n))
        np.savez_compressed(file_path, materials=chunk)
    for d in docs:
        d["structure"].to(filename=d['task_id']+'.cif',fmt="CIF")

"""
query all the data_id, and return the task_id files
"""
task_id_list=[]
with MPRester(API_KEY) as m:
    docs=m.query({},["structure","task_id"])
    for d in docs:
        task_id_list.append(d['task_id'])
FILEPATH = "data/sample/"
task_id_list=json.dumps(task_id_list, sort_keys=True, indent=4, separators=(',', ': '))
with open(FILEPATH + '/task_id_list.json', 'w') as file:
    file.write(task_id_list)

with open(FILEPATH + '/task_id_list.json', 'r') as file:
    data = json.load(file)

