import numpy as np
import matplotlib.pyplot as plt

# Load data from text files
n_values = np.loadtxt('n.txt')
p_values = np.loadtxt('p.txt')
v_values = np.loadtxt('v.txt')
x_values = np.loadtxt('x.txt')

# Check that all lists have the same length
if not (len(n_values) == len(p_values) == len(v_values) == len(x_values)):
    raise ValueError("All files must have the same number of entries.")

# Calculate start and end indices range edit
def get_subset_indices(data):
    total_len = len(data)
    start_idx = int(total_len*0)
    end_idx = int(total_len)
    return start_idx, end_idx

# Select subset of values
n_start, n_end = get_subset_indices(n_values)
p_start, p_end = get_subset_indices(p_values)
v_start, v_end = get_subset_indices(v_values)
x_start, x_end = get_subset_indices(x_values)

n_subset = n_values[n_start:n_end]
p_subset = p_values[p_start:p_end]
v_subset = v_values[v_start:v_end]
x_subset = x_values[x_start:x_end]

# Plot n vs x
plt.figure(figsize=(10, 6))
plt.plot(x_subset, n_subset, label='n vs x', color='blue', marker='o')
plt.xlabel('x')
plt.ylabel('n(x)')
plt.title('n vs x')
plt.legend()
plt.grid(True)
plt.savefig("nplot.png")
plt.show()

# Plot p vs x
plt.figure(figsize=(10, 6))
plt.plot(x_subset, p_subset, label='p vs x', color='red', marker='s')
plt.xlabel('x')
plt.ylabel('p(x)')
plt.title('p vs x')
plt.legend()
plt.grid(True)
plt.savefig("pplot.png")
plt.show()

# Plot v vs x
plt.figure(figsize=(10, 6))
plt.plot(x_subset, v_subset, label='v vs x', color='violet', marker='^')
plt.xlabel('x')
plt.ylabel('v(x)')
plt.title('v vs x')
plt.legend()
plt.grid(True)
plt.savefig("vplot.png")
plt.show()

# Plot n and p vs x
plt.figure(figsize=(10, 6))
plt.plot(x_subset, n_subset, label='n vs x', color='blue', marker='o')
plt.plot(x_subset, p_subset, label='p vs x', color='red', marker='s')
plt.xlabel('x')
plt.ylabel('n(x), p(x)')
plt.title('n and p vs x')
plt.legend()
plt.grid(True)
plt.savefig("npplot.png")
plt.show()
