"""Christa Caggiano 
Jan 2017 
Takes in coverage plot data (signal summed across a given number of genes)
and manipulates from a .result file to a csv for use in R. Sums the total signal
for each bp analyzed (default 1000), and averages any replicates together."""


import pandas as pd
from builtins import input
import glob
import csv 

# functions that remove extraneous columns and sum all values in one columns
def remove_columns(df):

    for i in range(0, 5):
        df = df.drop(i, 1) # 
    return df

def sum_columns(df):
    list_column_sums = [sum([row[column] for _, row in df.iterrows()]) for column in df]
    return list_column_sums
def append_columns(d, df, i):
    print(len(d))
    df[i] = d
    return df

# accepts user input for files they want modified - best if not running code en masse 
# file_number = int(input('How many files do you have?: '))
# for i in range(file_number):
#     file_names.append(input('Please type the name of file {}: ' .format(i+1)))

# automatically pulls files with desired extension from current directory 
file_names = []
file_names += (glob.glob("*.sort_BothStrand.result"))


# bp_number = int((input('How many base pairs around the peak did you analyze? Hit 0 to use default (500bp) ')))
bp_number = 0 

if not (bp_number == 500 or bp_number == 0):
    axis = list(range(-int((bp_number / 2)), int((bp_number / 2)), 1))
    column_sums_for_graphing = pd.DataFrame(axis, index=range(bp_number)) 
else:
    axis = list(range(-500, 502, 1))
    column_sums_for_graphing = pd.DataFrame(axis, index=range(1002))

lists = []

# performs the manipulation of each file the user inputted
for i in file_names:
    i_df = pd.read_table(i, sep='\t', header=None)
    i_df = remove_columns(i_df)
    i_sum = sum_columns(i_df)
    column_sums_for_graphing = append_columns(i_sum, column_sums_for_graphing, i)
    lists.append(i_sum)

#saves as a CSV
column_sums_for_graphing.to_csv('graphing_output.csv', sep=',')

