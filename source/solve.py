import time
import os
import csv
from multiprocessing import Process, Manager
from pysat.solvers import Glucose42, Cadical195
from ver1_glu import solve_inc, compute_upper_bound

TIMEOUT = 600

def main():
