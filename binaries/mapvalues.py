#!/usr/bin/env python

import pandas as pd
import numpy as np
from sys import argv

# infile arg used to make compatible for both sorted and organized tables
script, infile, inquality, transposedTable, finisedTable = argv

mytable = pd.read_csv(infile, sep='\t')
quality = pd.read_csv(inquality, sep='\t')

# set index to "reference_pos" so generic index does not transpose
mytable = mytable.set_index('reference_pos')
mytable = mytable.transpose()

# write to csv to import back with generic index again
# seems like a hack that can be done better
mytable.to_csv(transposedTable, sep="\t", index_label='reference_pos')

# can't merge on index but this newly imported transpose is formated correctly
mytable = pd.read_csv(transposedTable, sep='\t')
mytable = mytable.merge(quality, on='reference_pos', how='inner')

# set index to "reference_pos" so generic index does not transpose 
mytable = mytable.set_index('reference_pos')
mytable = mytable.transpose()

# since "reference_pos" was set as index it needs to be explicitly written into csv
mytable.to_csv(finisedTable, sep="\t", index_label='reference_pos')
