#!/usr/bin/env python3
# Author : Diana Aguilar
'''
Use: filtertable.py table outfile column_number value action header
'''
import sys
import pandas as pd
import numpy as np

intable=sys.argv[1]
outfilename=sys.argv[2]
column_number=sys.argv[3]
value=int(sys.argv[4])
action=sys.argv[5]
header=sys.argv[6]

assert action in ["exclude","leq","geq"], "action must be :exclude, leq or geq"

if header=="True":
	table2filter=pd.read_table(intable)
else:
	table2filter=pd.read_table(intable,header=None)
	column_number=int(column_number)

if action=="exclude":
	table2keep=table2filter[float(table2filter[column_number])!=value]
elif action=="leq":
	table2keep=table2filter[float(table2filter[column_number])<=value]
elif action=="geq":
	table2keep=table2filter[float(table2filter[column_number])>=value]

table2keep.to_csv(outfilename,sep="\t",index=False,header=header)
print(len(table2filter),len(table2keep))
