from annoy import AnnoyIndex
import sys

sys.path.append("..")
from data import embeddings, calc_time

f = 128
print(embeddings.shape)

def main(l):
    t = AnnoyIndex(f, 'angular')  # Length of item vector that will be indexed
    for i,v in enumerate(embeddings):
        t.add_item(i, v)

    t.build(l) # 100 trees
    return t

def testrange(r=[10,100]):
    for i in r:
        t = calc_time(main, 'build %s trees anno index'%i, i)
        t.save('test%s.ann'%i)
testrange()
