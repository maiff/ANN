## ANN Algorithm

### Data
You can use [ann-benchmarks](https://github.com/erikbern/ann-benchmarks) to test.
This project uses the [SIFT](http://ann-benchmarks.com/sift-128-euclidean.hdf5) dataset.

### Algorithm
Every Algorithm is implemented by C++ and having a Python wrapper.


✅ [KD-Tree](./KDTree)
- nearest neighbors algorithm
- k-nearest neighbors algorithm
- based [KDTree](https://githubmemory.com/repo/ertosns/KDTree)


✅ [Annoy](./test_annoy)

……

### Time Cost
All the algorithms are running in MacOS  with those settings:

![](http://test.maiff.cn/6bLox1LeGEimage.png)

Building the trees on SIFT dataset which has 1000000 embeddings with 128 dimension. Querying is on the same dataset which has 10000 embeddings with 128 dimension.

| Algorithm | params    | build time cost(s) | query knn time cost(s) | recall |
| --------- | --------- | ------------------ | ---------------------- | ------ |
| KD-Tree   | None      | 22±4              | 6800±10               | 1      |
| Annoy     | 10 trees  | 19±3              | 9±3                   | 0.42   |
| Annoy     | 100 trees | 72±3              | 49±4                  | 0.88   |
