import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

dataFile = os.path.join(os.getcwd(), os.path.normpath("data/TRILEGAL.table"))	# File location


# Imports the TRILEGAL data file as a table
def load():
	info = pd.read_table(dataFile, delim_whitespace=True)
	return info
