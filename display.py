import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("output.csv")

x = df.iloc[:, 0]

plt.figure(figsize=(8, 6))

for column in df.columns[1:]:
    plt.semilogy(x, df[column], marker='o', linewidth=2, label=column)

plt.xlabel("Iteration #")
plt.ylabel("Residual Norm")
plt.title("Residual Norm vs Iteration Count for Jacobi, Gauss-Seidel, and SOR (omega = .8) Methods")
plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.legend()

plt.tight_layout()
plt.show()
