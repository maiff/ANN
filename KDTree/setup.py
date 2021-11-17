from setuptools import setup, Extension
import os
from pybind11.setup_helpers import ParallelCompile, Pybind11Extension
cxx_std = int(os.environ.get("CMAKE_CXX_STANDARD", "14"))
functions_module = Pybind11Extension(
    name='kdtree',
    sources=['KDTree.cpp'],
    extra_compile_args=["-O3","-fPIC"],
     cxx_std=cxx_std,
#     include_dirs=['/Users/admin/opt/miniconda3/lib/python3.8/site-packages/pybind11/include/'],
)

setup(ext_modules=[functions_module])