#!/usr/bin/env python3
# Author : Diana Aguilar
'''
Only keep rows that contain elements from a column'''


import pandas as pd
import sys

intable=sys.argv[1]
column_number=int(sys.argv[2])
columnfile=sys.argv[3]
outfile=sys.argv[4]

annotation=pd.read_table(intable,header=None)
seqs2keep=pd.read_table(columnfile,header=None)
annotation2keep=annotation[annotation[column_number].isin(seqs2keep[0])]
annotation2keep.to_csv(outfile,sep="\t",index=False,header=False)
