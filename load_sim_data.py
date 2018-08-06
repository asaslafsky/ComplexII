# import the pysb module and all its methods and functions
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np
import pandas as pd
import seaborn as sns
from pysb.simulator import StochKitSimulator
from pysb.simulator import SimulationResult
from pysb import *

df = SimulationResult.load('filename').dataframe

print(df)