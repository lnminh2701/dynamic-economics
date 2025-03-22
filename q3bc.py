import numpy as np
import matplotlib.pyplot as plt

try:
    import pandas as pd
except ModuleNotFoundError:
    print("Warning: 'pandas' module not found. Install it using: pip install pandas")

def rouwenhorst(N, rho, sigma, mu):
    """
    Implements Rouwenhorst's method to approximate an AR(1) process.
    """
    y_mean = mu / (1 - rho)
    y_sd = np.sqrt(sigma**2 / (1 - rho**2))
    state_range = y_sd * np.sqrt(N - 1)
    state_values = np.linspace(y_mean - state_range, y_mean + state_range, N)
    
    def build_transition_matrix(n, p, q):
        if n == 2:
            return np.array([[p, 1 - p], [1 - q, q]])
        prev_matrix = build_transition_matrix(n - 1, p, q)
        pmat = np.zeros((n, n))
        pmat[:n-1, :n-1] += p * prev_matrix
        pmat[:n-1, 1:] += (1 - p) * prev_matrix
        pmat[1:, :-1] += (1 - q) * prev_matrix
        pmat[1:, 1:] += q * prev_matrix
        pmat[1:n-1, :] /= 2  # Normalize rows
        return pmat
    
    p = (1 + rho) / 2
    transition_matrix = build_transition_matrix(N, p, p)
    return state_values, transition_matrix


def simulate(grid, pmat, T):
    """
    Simulates a discrete Markov chain given a transition matrix.
    """
    N = len(grid)
    state0 = np.random.choice(N, p=np.ones(N) / N)
    cmat = np.cumsum(pmat, axis=1)
    y = np.zeros(T)
    for i in range(T):
        y[i] = grid[state0]
        state0 = np.searchsorted(cmat[state0, :], np.random.uniform())
    return y

# Part (b): Discretizing AR(1) Process Using Rouwenhorst's Method
gamma_1 = 0.85
N = 7
mu = 0.5
sigma = 1
T = 50
seed = 2025

# Set seed
np.random.seed(seed)

# Compute state values and transition matrix
state_values_b, transition_matrix_b = rouwenhorst(N, gamma_1, sigma, mu)

# Print the discretized state values
print("Discretized State Values (Grid):")
print(np.round(state_values_b, 4))

# Print the transition matrix
print("\nTransition Matrix:")
if 'pd' in globals():
    print(pd.DataFrame(np.round(transition_matrix_b, 4), 
                       index=np.round(state_values_b, 2), 
                       columns=np.round(state_values_b, 2)))
else:
    print(np.round(transition_matrix_b, 4))


# Part (c): Simulating the Markov Chain and Plotting
# Simulate the Markov chain
sim_b = simulate(state_values_b, transition_matrix_b, T)

# Plot results
plt.figure(figsize=(8, 5))
plt.plot(range(T), sim_b, marker='', linestyle='-', label="Simulated AR(1)")
plt.axhline(y=np.mean(state_values_b), color='red', linestyle='--', label="Grid Mean")
plt.xlabel('Time')
plt.ylabel('Y')
plt.title("AR(1) Process with Grid (Rouwenhorst)")
plt.legend()
plt.grid(True, linestyle=':', linewidth=0.7)
plt.show()
