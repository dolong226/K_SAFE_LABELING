from typing import List, Tuple
from pysat.formula import CNF
from pysat.solvers import Glucose3
import time
from multiprocessing import Process, Manager
import pandas as pd
import numpy as np

# n: vertex
# edges
# k
# UB: upper bound

TIMEOUT = 600
def var_id(v: int,l: int,UB: int) -> int:
    return v * UB + l
def control_var_id(s: int, n:int, UB:int) -> int:
    return n * UB + s
def build_base_cnf(n: int, edges: List[Tuple[int, int]], k:int, UB:int, root_fix: bool = True) -> CNF:
    cnf = CNF()

    # EXO per vertex:
    for v in range(n):
        #ALO 
        cnf.append([var_id(v,i,UB) for i in range(1,UB+1)])
        #AMO binomial
        for a in range (1,UB + 1):
            for b in range (a+1, UB + 1):
                cnf.append([-var_id(v,a,UB), -var_id(v,b,UB)])

    # All Different per label
    for l in range(1, UB+1):
        for u in range(n):
            for v in range(u+1, n):
                cnf.append([-var_id(u,l,UB), -var_id(v,l,UB)])
    
    # k-safe
    for (u,v) in edges:
        if u > v: u,v = v,u
        for l1 in range(1, UB + 1):
            low = max (1, l1 - (k - 1))
            high = min (UB, l1 + (k - 1))
            for l2 in range(low, high + 1):
                cnf.append([-var_id(u,l1,UB), - var_id(v,l2, UB)])
    
    # symmetry breaking
    if root_fix:
        cnf.append([var_id(0,1,UB)])

    # inc
    for S in range (1, UB + 1):
        Cs = control_var_id(S,n,UB)
        for v in range(n):
            for l in range(S+1, UB + 1):
                cnf.append([-var_id(v,l,UB), Cs])
    
    return cnf

def solve_inc(n: int, edges:List[Tuple[int,int]], k:int, UB:int, shared):
    cnf = build_base_cnf(n,edges, k, UB)
    solver = Glucose3(bootstrap_with=cnf.clauses)

    best_S, best_model = None, None
    high = UB
    while high > 0:
        C = control_var_id(high, n, UB)
        sat = solver.solve(assumptions=[-C])

        if sat: 
            best_S = high
            best_model = solver.get_model()
            shared["best_S"] = best_S
            shared["best_model"] = best_model
            high = best_S - 1
        else: 
            break
    # decode labels
    labels = [-1] * n
    modelset = set(l for l in best_model if l > 0)
    for v in range(n):
        for l in range(1,UB + 1):
            if var_id(v,l,UB) in modelset:
                labels[v] = l
                break
    shared["labels"] = labels
    print ("Optimal UB =", best_S)
    for v in range(n):
        print(f"Vertex {v} -> label {labels[v]}")
    
    
    return best_S, labels


def main():
    n = int(input())
    m = int(input())
    edges = [tuple(map(int,input().split())) for _ in range (m)]
    k = int(input())
    UB = k * (n - 1)+1

    start = time.perf_counter()

    with Manager() as manager:
        shared = manager.dict()
        p = Process(target=solve_inc, args = (n,edges,k,UB,shared))
        p.start()
        p.join(timeout = TIMEOUT)

        if p.is_alive():
            print(f"Timeout sau {TIMEOUT} s")
            p.terminate()
            p.join()

            if "best_S" in shared and shared["best_S"] is not None:
                print("Last UB before timeout =", shared["best_S"])
            else:
                print("No solution was found before the timeout")
        
        end = time.perf_counter()
        print(f"Time: {end - start:.2f} s")

if __name__ == "__main__":
    main()


