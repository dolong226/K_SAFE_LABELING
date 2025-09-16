import os
import time
import argparse
import pandas as pd
from multiprocessing import Process, Manager
import importlib

from pysat.solvers import Glucose42, Cadical195
# from ver1_glu import compute_upper_bound, solve_inc, build_base_cnf
# from direct import compute_upper_bound,solve_inc,build_base_cnf
from ver2 import compute_upper_bound, solve_inc, build_base_cnf
def parse_instance(path: str):
    with open(path, "r", encoding="utf-8") as f:
        first = ""
        while first == "":
            line = f.readline()
            if not line:
                raise ValueError(f"{path}: file trống hoặc thiếu dòng đầu tiên (n m k)")
            first = line.strip()
            if first.startswith("#"):
                first = ""
        parts = first.split()
        if len(parts) < 3:
            raise ValueError(f"{path}: dòng đầu phải có 3 số n m k")
        n, m, k = map(int, parts[:3])

        edges = []
        while len(edges) < m:
            line = f.readline()
            if not line:
                break
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            u, v = map(int, s.split())
            if u == v:
                continue  
            edges.append((u, v))
        if len(edges) != m:
            raise ValueError(f"{path}: mong {m} cạnh, đọc được {len(edges)}")
    return n, m, k, edges


def cnf_stats(n: int, edges, k: int, UB: int, mod):
    cnf = mod.build_base_cnf(n, edges, k, UB)
    clauses = len(cnf.clauses)
    max_var = 0
    for cl in cnf.clauses:
        for lit in cl:
            v = abs(lit)
            if v > max_var:
                max_var = v
    return clauses, max_var


def choose_solver(name: str):
    name = name.lower()
    if name in ("glucose", "g3", "glucose42"):
        return Glucose42
    if name in ("cadical", "cadical195", "c195"):
        return Cadical195
    raise ValueError(f"Invalid solver: {name}")


def run_one(path: str, solver_name: str, timeout_sec: float, method: str):
    n, m, k, edges = parse_instance(path)

    # import module theo method
    if method == "direct":
        mod = importlib.import_module("direct")
    elif method == "ver1":
        mod = importlib.import_module("ver1")
    elif method == "ver2":
        mod = importlib.import_module("ver2")
        raise ValueError(f"Method chưa hỗ trợ: {method}")

    UB = mod.compute_upper_bound(n, edges, k)
    if UB is None:
        UB = k * (n - 1) + 1

    clauses, vars_ = cnf_stats(n, edges, k, UB, mod)

    start = time.perf_counter()
    with Manager() as manager:
        shared = manager.dict()
        SolverCls = choose_solver(solver_name)

        p = Process(target=mod.solve_inc, args=(n, edges, k, UB, shared, SolverCls))
        p.start()
        p.join(timeout_sec)

        status = None
        span = None
        if p.is_alive():
            # Timeout
            p.terminate()
            p.join()
            status = "TO"
            if "best_S" in shared and shared["best_S"] is not None:
                span = int(shared["best_S"])
            elapsed = timeout_sec
        else:
            if "best_S" in shared and shared["best_S"] is not None:
                span = int(shared["best_S"])
                status = "SAT"
            else:
                status = "UNSAT"
            elapsed = time.perf_counter() - start

    return {
        "Instance": f"{os.path.basename(path)}[{method}]",  # in kèm tên method
        "n": n,
        "m": m,
        "k": k,
        "UB": int(UB),
        "Span": (None if span is None else int(span)),
        "Clauses": int(clauses),
        "Vars": int(vars_),
        "Status": status,
        "Time(s)": round(elapsed, 2),
    }


def collect_paths(input_path: str, ext: str):
    paths = []
    if os.path.isdir(input_path):
        for name in sorted(os.listdir(input_path)):
            if name.lower().endswith(ext.lower()):
                paths.append(os.path.join(input_path, name))
    else:
        if input_path.lower().endswith(ext.lower()) or ext == "*":
            paths = [input_path]
    if not paths:
        raise ValueError(f"Không tìm thấy instance nào trong '{input_path}' với ext '{ext}'")
    return paths


def main():
    parser = argparse.ArgumentParser(description="Solve k-safe labeling instances và xuất CSV")
    parser.add_argument("--input", required=True, help="Thư mục chứa các instance hoặc đường dẫn 1 file")

    # file name
    parser.add_argument("--out", default="results_ver2_cadical.csv", help="Đường dẫn file CSV output")
    
    parser.add_argument("--solver", default="glucose", help="glucose | cadical")
    parser.add_argument("--timeout", type=float, default=None, help="Timeout mỗi instance (giây)")
    parser.add_argument("--ext", default=".txt", help="Phần mở rộng file instance khi --input là thư mục (vd .txt)")
    parser.add_argument("--method", default="default", help="Phiên bản thuật toán: default | ver1_glu | ...")

    args = parser.parse_args()

    # nếu không truyền timeout thì lấy TIMEOUT từ module method
    if args.timeout is None:
        if args.method == "default":
            args.timeout = importlib.import_module("direct").TIMEOUT
        elif args.method == "ver1":
            args.timeout = importlib.import_module("ver1_glu").TIMEOUT
        elif args.method == "ver2":
            args.timeout == importlib.import_module("ver2").TIMEOUT

    paths = collect_paths(args.input, args.ext)

    rows = []
    for p in paths:
        try:
            row = run_one(p, args.solver, args.timeout, args.method)
            rows.append(row)
            print(f"{row['Instance']}: {row['Status']} span={row['Span']} time={row['Time(s)']}s")
        except Exception as e:
            print(f"Lỗi với {os.path.basename(p)}: {e}")

    df = pd.DataFrame(rows, columns=[
        "Instance", "n", "m", "k", "UB", "Span", "Clauses", "Vars", "Status", "Time(s)"
    ])
    df.to_csv(args.out, index=False)
    print(f"\nĐã ghi {args.out} với {len(df)} hàng.")


if __name__ == "__main__":
    main()

# ví dụ chạy:
# python solve.py --input ./instances --solver cadical --method default
# python solve.py --input ./instances --solver glucose --method ver1_glu

# python source/solve2.py --input ./instance --solver cadical --method ver2