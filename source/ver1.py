from pysat.formula import CNF
from typing import List, Tuple
from pysat.solvers import Glucose42, Cadical195
import time
from multiprocessing import Process, Manager
import networkx as nx

TIMEOUT = 600

def X_id(i: int, j: int, UB: int) -> int:
    return i * UB + j

def control_var_id(S: int, n: int, UB: int) -> int:
    return n * UB + S

def symmetry_find(n: int, edges:List[Tuple[int,int]]) -> int:
    fmax = 0
    sym = 0
    arr = [0] * n
    for (u,v) in edges:
        arr[u] += 1
        arr[v] += 1
        if arr[u] > fmax:
            fmax = arr[u]
            sym = u
        if arr[v] > fmax:
            fmax = arr[v]
            sym = v
    return sym

def build_base_cnf(n: int, edges:List[Tuple[int,int]], k:int, UB: int, root_fix:bool = True) -> CNF:
    cnf = CNF()

    # X_i,UB = 1
    for i in range(n):
        cnf.append([X_id(i,UB,UB)])

    # X_i,j => X_i,j+1
    for i in range (n):
        for j in range (1, UB):
            cnf.append([-X_id(i,j,UB), X_id(i,j+1,UB)])

    # -X_i1,1 V - X_i2,1
    for i1 in range(n):
        for i2 in range(i1 + 1, n):
            cnf.append([-X_id(i1,  1, UB), -X_id(i2, 1, UB)])

    # -(K_i1,j ∧ K_i2,j)
    # -X_i1,j V X_i1,j-1 V -X_i2,j V X_i2,j-1
    for i1 in range(n):
        for i2 in range(i1 + 1, n):
            for j in range(2, UB + 1):
                cnf.append([-X_id(i1,j,UB), X_id(i1,j-1,UB), -X_id(i2,j,UB), X_id(i2,j-1,UB)])

    # K_u,val => X_v,val-k V -X_v,val+k-1
    # -X_u,val V X_u,val-1 V X_v,val-k V -X_v,val+k-1
    directed_edges = edges + [(v,u) for u,v in edges]
    for u,v in directed_edges:
        for val in range(1,UB + 1): 
            condition_lits = []
            if val - k >= 1:
                condition_lits.append(X_id(v,val-k,UB))
            if val + k - 1 <= UB:
                condition_lits.append(-X_id(v,val + k - 1, UB))
            
            if not condition_lits:
                if val == 1:
                    clause = [-X_id(u,1,UB)]
                else:
                    clause = [-X_id(u,val,UB), X_id(u, val-1, UB)]
            else:
                if val == 1:
                    clause = [-X_id(u,1,UB)] + condition_lits
                else:
                    clause = [-X_id(u,val,UB), X_id(u,val-1,UB)] + condition_lits
            
            if clause:
                cnf.append(clause)
    
    # symmetry breaking ()
    if root_fix:
        s = symmetry_find(n,edges)
        cnf.append([X_id(s,1,UB)])

    # incremental
    for S in range(1, UB + 1):
        Cs = control_var_id(S,n,UB)
        for i in range(n):
            cnf.append([Cs, X_id(i,S,UB)])

    return cnf

def solve_inc(n: int, edges: List[Tuple[int, int]], k: int, UB: int, shared, solver_cls = Glucose42):
    cnf = build_base_cnf(n, edges, k,UB)
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
    if best_S is None:
        print("UNSAT")
        return None, None
    # decode labels
    labels = [-1] * n
    modelset = set(l for l in best_model if l > 0)
    for v in range(n):
        for l in range(1,UB + 1):
            if X_id(v,l,UB) in modelset:
                labels[v] = l
                break
    shared["labels"] = labels
    print("Optimal UB =", best_S)
    for v in range(n):
        print(f"Vertex {v} -> label {labels[v]}")

    return best_S, labels

def compute_upper_bound(n, edges, k):

    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(edges)

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


def main():
    n,m,k = map(int, input().split())

    edges = [tuple(map(int, input().split())) for _ in range(m)]
    UB = compute_upper_bound(n,edges, k)
    

    start = time.perf_counter()

    with Manager() as manager:
        shared = manager.dict()
        p = Process(target=solve_inc, args= (n,edges,k,UB, shared, Glucose42))
        # Process(target=solve_inc, args= (n,edges,k,UB, shared, Cadical195))

        p.start()
        p.join(timeout = TIMEOUT)

        if p.is_alive():
            print(f"TIMEOUT")
            p.terminate()
            p.join()
        
            if "best_S" in shared and shared["best_S"] is not None:
                print("Last UB before timeout =", shared["best_S"])
            else:
                print("No solution was found before timeout")
        end = time.perf_counter()
        print(f"Time: {end - start:.2f} s")

if __name__ == "__main__":
    main()