import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Imports the TRILEGAL data file as a table


def load_info():
	dataFile = os.path.join(os.getcwd(), os.path.normpath(
		"data/TRILEGAL.table"))  # File location
	info = pd.read_table(dataFile, delim_whitespace=True)
	return info

# returns a specific comlumn (string)


def info_col(column):
	info = load_info()
	return info[[column]]

# returns number of rows in dataset


def info_len():
	info = load_info()
	return len(info.index)
# if run as script, do something
if __name__ == '__main__':
	df = load_info()
	print df
