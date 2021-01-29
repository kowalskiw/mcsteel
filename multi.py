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

'''Run single simulation'''


class RunSim:
    def __init__(self):
        # import FEA path
        pass

    # run SAFIR T2D
    def t2d(self):
        pass

    # calculate mean temperature of the profile
    def mean_temp(self):
        pass

    # create temperature-time table of the section and return it
    def main(self):
        pass


'''Run node queue'''


class Queue:
    def __init__(self):
        # import simulation configuration (including stop conditions)
        # import queue
        pass

    # run simulation queue
    def run(self):
        # loop
        #   check for the tasks in queue
        #   run first task if queue is not empty
        #   edit queue when finished

        pass

    # choose theta_a,max and t_theta,a,cr and add those values to case.res
    def save_res(self):
        pass


'''Run multisimulation on cluster'''


class Cluster:
    def __init__(self):
        # import cluster configuration
        pass

    # send task to the nodes' queues
    def assign(self):
        pass

    # initialize calculations on nodes
    def wake_nodes(self):
        pass

    # check stop condition (sim number or RMSE), stop analysis or generate more scenarios if needed
    def check_stop(self):
        pass
