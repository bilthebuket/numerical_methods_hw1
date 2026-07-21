import matplotlib.pyplot as plt

# File containing comma-separated values
filename = "plot2.txt"

# Read the values
with open(filename, "r") as f:
    contents = f.read()

arr = [float(x.strip()) for x in contents.split(",") if x.strip()]

# Generate x-values (indices)
x = range(len(arr))

# Create semi-log plot (logarithmic y-axis)
plt.semilogy(x, arr)

plt.xlabel("Iteration Number")
plt.ylabel("inf_norm(v_k - v_1)")
plt.title("Semilogy plot of inf_norm(v_k - v_1) vs Iteration Number")
plt.grid(True, which="both", linestyle="--")

plt.show()
