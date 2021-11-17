import h5py,time

from numpy.lib.financial import ipmt
from pathlib import Path
import os

p = os.path.dirname(str(Path(os.path.abspath(__file__))))
h = h5py.File(os.path.join(os.path.dirname(str(Path(os.path.abspath(__file__)))), './sift-128-euclidean.hdf5'),'r')

labels =  h['neighbors'][:]
query = h['test'][:]
embeddings = h['train'][:]

def calc_time(f, des="", *args,**kwargs):
    print(des," cost time:")
    s = time.time()
    v = f(*args,**kwargs)
    e = time.time()
    print(e-s)
    return v

def recall(key, value):
    a = set(labels[key])
    b = set(value)
    return len(a&b)/100
