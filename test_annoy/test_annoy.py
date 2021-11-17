from annoy import AnnoyIndex
import sys
sys.path.append("..")
from data import query, calc_time,recall

f=128

def test(l):
    u = AnnoyIndex(f, 'angular')
    u.load('test%s.ann'%l) # super fast, will just mmap the file
    r = 0
    for i,q in  enumerate(query):
        v = u.get_nns_by_vector(q, 100)
        t = recall(i, v)
        r+=t
    return r/len(query)

def testrange(r=[10,100]):
    for i in r:
        t = calc_time(test, "test %s trees"%i, i)
        print("%s trees recall:"%i, t)
testrange()

