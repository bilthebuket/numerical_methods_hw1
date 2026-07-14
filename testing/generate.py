import random

N = 2048
random.seed(1)

with open("matrix.csv", "w") as f:
    for i in range(N):
        row = ", ".join(str(random.randint(0, 99)) for _ in range(N))
        if i != N - 1:
            f.write(row + ";\n")
        else:
            f.write(row + "\n")
