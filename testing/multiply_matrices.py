import numpy as np
import time

def read_matrix(filename):
    rows=[]
    with open(filename) as f:
        for line in f:
            line=line.strip()
            if not line:
                continue
            if line.endswith(";"):
                line=line[:-1]
            rows.append([int(x.strip()) for x in line.split(",")])
    return np.array(rows,dtype=np.int64)

A=read_matrix("matrix.csv")
B=read_matrix("matrix.csv")

start = time.perf_counter()
C=A @ B
end = time.perf_counter()

print(f"Time: {end - start}s")

print(C)
