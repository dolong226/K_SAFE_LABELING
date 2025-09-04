from typing import List, Tuple
from pysat.formula import CNF
from pysat.solvers import Glucose42, Cadical195
import time
from multiprocessing import Process, Manager
import pandas as pd
import numpy as np
import networkx as nx
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

def compute_upper_bound(n, edges, k):

    # tạo đồ thị gốc
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(edges)

    # tạo complete graph
    G2 = nx.complete_graph(n)
    E_c = list(nx.non_edges(G))
    edge_checked = set()
    clique_unchanged = 0
    clique_changed_flag = 0
    lst2 = []

    while E_c:
        q_old = nx.graph_clique_number(G2)

        e_incident = None
        for e in E_c:
            if e[0] not in edge_checked and e[1] not in edge_checked:
                e_incident = e
                edge_checked.add(e[0])
                edge_checked.add(e[1])
                break

        if e_incident is None:
            if len(edge_checked) != n:
                for e in E_c:
                    if e[0] in edge_checked and e[1] in edge_checked:
                        e_incident = e
                        break
            else:
                edge_checked = set()

        if e_incident is None:
            edge_checked = set()
            continue

        E_c.remove(e_incident)
        if G2.has_edge(*e_incident):
            G2.remove_edge(*e_incident)
        elif G2.has_edge(e_incident[1], e_incident[0]):
            G2.remove_edge(e_incident[1], e_incident[0])

        q_now = nx.graph_clique_number(G2)

        if q_now == q_old and q_now >= int(n / 2) and clique_changed_flag == 0:
            clique_unchanged = q_old
            clique_changed_flag = 1

        if q_now >= clique_unchanged:
            upper_bound = (k - 1) * q_now + n + 1 - k
        else:
            upper_bound = (k - 1) * clique_unchanged + n + 1 - k

        lst2.append(upper_bound)

    return lst2[-1] if lst2 else None
  
def solve_inc(n: int, edges:List[Tuple[int,int]], k:int, UB:int, shared, solver_cls = Glucose42):
    cnf = build_base_cnf(n,edges, k, UB)
    solver = solver_cls(bootstrap_with=cnf.clauses)

    # print clauses, vars
    print("Variables: ", solver.nof_vars())
    print("Clauses: ", solver.nof_clauses())
    
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
    UB = compute_upper_bound(n,edges, k)

    start = time.perf_counter()

    with Manager() as manager:
        shared = manager.dict()
        p = Process(target=solve_inc, args = (n,edges,k,UB,shared,Glucose42))
        # p = Process(target=solve_inc, args = (n,edges,k,UB,shared,Cadical195))
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


