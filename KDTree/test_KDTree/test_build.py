import sys
import kdtree

sys.path.append("../..")
from data import embeddings, query, calc_time, recall

f = 128
print(embeddings.shape)
def main():
    t = kdtree.KDTree(embeddings)
    return t

def testquery(t):
    r = 0
    for i,q in  enumerate(query[:100]):
        v = t.nearest_topk_index(q, 100)
        tt = recall(i, v)
        r+=tt
    return r/len(query[:100])

def testrange():
    t = calc_time(main, 'build kdtree trees')
    r = calc_time(testquery, 'test query', t)
    print("recall: ", r)

testrange()
