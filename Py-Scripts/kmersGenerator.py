import numpy as np

def kRecur(alphas, prfx, l, k):
    if k==0:
        print(prfx)
    else:
        for i in range(l):
            newPrfx = prfx + alphas[i]
            kRecur(alphas, newPrfx, l, k-1)

def combos(alphas, k ):
    l = len(alphas)
    kRecur(alphas, "", l, k)

alphabet = 'ACGT'

# combos(alphabet, 3)
results=[]

def KMerGenerator(kmer, k):
    if k==0:
        # print(kmer)
        results.append(kmer)
    else:
        for c in alphabet:
            newKmer = kmer + c
            KMerGenerator(newKmer, k-1)

for l in range(5):
    KMerGenerator("", l)
    print( "len={:n} -> {:,}".format(l, len(results)))
    results=[]

l1 = np.array(['AAA', 'ACC', 'CAA', 'CCT', 'GTT'])
l2 = np.array(['AAC', 'GTT', 'TTA', 'ACC', 'CAA', 'CAT', 'GGT'])

print(np.intersect1d(l1, l2))