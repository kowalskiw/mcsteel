from os import listdir, getcwd, chdir, popen, mkdir, path, replace
import json as js
from pynput.keyboard import Key, Controller
import time
from sys import argv
import numpy as np
from numpy import sqrt, log, random, pi
from datetime import datetime as dt
from export import Export
from fires import Fires
import ezdxf
import shapely.geometry as sh


'''Read geometry and map it to the fire (choose the most exposed sections)'''


class Single:
    def __init__(self):
        pass

    # import geometry data
    def read_dxf(self):
        pass

    # map fire to geometry - choose fire scenario
    def map(self):
        pass



# import fire coordinates, HRR-time and Diameter-time tables


# generates set of n scenarios
def generate_set(n):
    pass

# generates initial files (elem.in, prof.in, locafi.txt)

def generate_in(sim_data):
   pass

