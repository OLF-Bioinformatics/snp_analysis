#!/usr/bin/env python

import pandas as pd
from sys import argv

script, inputTable, cleanAlignment, sortedTable, organizedTable = argv

# There must be a top row header name to match for join to work, "reference_pos"
mytable = pd.read_csv(inputTable, delimiter="\t")

# mylist must be ordered with 2 top columns "reference_pos" and "reference_call"
mylist = pd.read_csv(cleanAlignment)

# axis drops the entire column if column is empty
mytable.dropna(axis='columns', how='all', inplace='True')

# mytable orders to mylist (alignment file)
# This sorts the rows vertically in the best position provided by RAxML
mytable = mylist.merge(mytable, on='reference_pos', how="outer")

mytable.to_csv(sortedTable, sep="\t", index=False)

col_num = mytable.shape[1]
#print ("There are %s columns in table." % col_num)

row_num = mytable.shape[0]
#print ("There are %s rows in table." % row_num)

# Counting the number of SNPs/position and group
# Iterate through each column of table
# If sample call is equal to reference call (cell [0,x]), count finding
# Place sum in new row for each column
count=0
# Get a column number
for each_column in range(1,col_num):
    # Iterate each cell in column
    for each_cell in mytable.ix[0:,each_column]:
        if each_cell == mytable.ix[0,each_column]:
            count += 1
    mytable.ix[row_num + 1,each_column] = count
    #mytable.ix[1, 1]="myvalue"
    count=0

# Iterate through each column of table
# Count distance SNP accures from reference call
count=0
# Get a column number
for each_column in range(1,col_num):
    # Iterate each cell in column
    for each_cell in mytable.ix[0:,each_column]:
        if each_cell == mytable.ix[0,each_column]:
            count += 1
        else:
            mytable.ix[row_num + 2,each_column] = count
            count=0
            break

mytrans = mytable.transpose()
col_num = mytrans.shape[1]
mytable = mytrans.sort_values([col_num, col_num - 1], ascending=[True, True]).transpose()

# Put the last column to the front
cols = mytable.columns.tolist()
cols = cols[-1:] + cols[:-1]
mytable = mytable[cols]

# Remove the last row with number counts
mytable = mytable[:-2]

mytable.to_csv(organizedTable, sep="\t", index=False)
