""" extract.py, Chris Monico.
    Read a csv file of triples produced by the main program
    and produce a tabular summary. The idea is that one should
    concatenate all of the individually produced files into a
    single one and then run this program to generate a summary.
"""
import pandas as pd
import numpy as np
from tabulate import tabulate

FILE='results/results.csv'

df = pd.read_csv(FILE, header=None, names=['q', 'k', 'maxchainlen'])
A = df.to_numpy()
kvals = np.unique(A[:,1])
mvals = np.unique(A[:,2])

Results = {}
maxQrow = ['max q']

for k in kvals:
    theseK = A[A[:,1]==k]
    Q = max(theseK[:,0])
    maxQrow.append(str(Q))
    print(f"For k={k}, q up to and including {Q}:")
    chainlenvals = np.unique(theseK[:,2])
    for m in chainlenvals:
        theseresults = theseK[theseK[:,2]==m]
        M = max(theseresults[:,0])
        print(f"k={k}, m={m+1}, largest q = {M}")
        Results[(k,m)] = M


Output = [['m'] + [f'k={k}' for k in kvals]]
Output.append(maxQrow)

for m in mvals:
    row=[f'{m+1}']
    for k in kvals:
        if (k,m) in Results:
            row.append(str(Results[k,m]))
        else:
            row.append('')
    Output.append(row)

print("\n\n")
print(tabulate(Output, headers="firstrow"))
