import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

dataFile = os.path.join(os.getcwd(), os.path.normpath("data/TRILEGAL.dat"))

def trImport():
    trInfo = pd.read_table(dataFile, delim_whitespace=True)
    return trInfo