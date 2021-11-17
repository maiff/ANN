# KDTree

Simple C++ static KD-Tree implementation with minimal functionality.

- points are given as STL vectors (and inserted in their own STL vector) so supports n-dimensional points for any n
- makes full trees, (i.e. does not cut-off the branching at some arbitrary level) giving the nearest neighbor query have (strong) logarithmic complexity.
- builds the tree in one go (does not support adding nodes, the tree is built from a list of points and cannot be altered afterwards)
- points are assumed to be STL vectors
- it provides the following queries:
	- nearest neighbor
	- neighbors within a given distance
	- k nearest neighbors (The project added with python wrapper)

### How to use
```
python setup.py build
mv build/*.so ./test_KDTree
cd test_KDTree
python test_build.py
```

The python wrapper as follow:
```
PYBIND11_MODULE(kdtree, m) {
    py::class_<KDTree>(m, "KDTree")
        .def(py::init<pointVec &>())
        .def("nearest_topk_point_array", &KDTree::nearest_topk_point_array)
        .def("nearest_topk_index", &KDTree::nearest_topk_index);
}
```

Meanwhile, the C++ examples which build witch cmake you can find in tests.



## License and copyright
© Maiff
© J. Frederico Carvalho

Licensed under the [BSD3 License](LICENSE)
