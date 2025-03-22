from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import uniform

def rouwenhorst(N, rho, sigma, mu):
    """
    Implements the Rouwenhorst method to approximate an AR(1) process.
    
    Parameters:
    N : int - Number of states
    rho : float - Autoregressive coefficient
    sigma : float - Standard deviation of the shock
    mu : float - Mean of the process
    
    Returns:
    state_values : numpy array - Discretized state space
    transition_matrix : numpy array - Transition probability matrix
    """
    # Ensure rho is a single float value
    if isinstance(rho, list):
        raise ValueError("rho should be a single float, not a list.")

    # Compute the mean and standard deviation of the stationary distribution
    y_mean = mu / (1 - rho)
    y_sd = np.sqrt(sigma**2 / (1 - rho**2))

    # Compute the state range
    state_range = y_sd * np.sqrt(N - 1)

    # Define state bounds
    upper_bound = y_mean + state_range
    lower_bound = y_mean - state_range

    # Generate state values
    state_values = np.linspace(lower_bound, upper_bound, N)

    def build_transition_matrix(n, p, q):
        """Recursively constructs the transition matrix for the Rouwenhorst method."""
        if n == 2:
            return np.array([[p, 1 - p], [1 - q, q]])
        
        prev_matrix = build_transition_matrix(n - 1, p, q)
        
        p1 = np.zeros((n, n))
        p2 = np.zeros((n, n))
        p3 = np.zeros((n, n))
        p4 = np.zeros((n, n))
        
        p1[:n - 1, :n - 1] = p * prev_matrix
        p2[:n - 1, 1:] = (1 - p) * prev_matrix
        p3[1:, :-1] = (1 - q) * prev_matrix
        p4[1:, 1:] = q * prev_matrix
        
        pmat = p1 + p2 + p3 + p4
        pmat[1:n - 1, :] /= 2  # Normalize rows
        
        return pmat

    # Compute transition probabilities
    p = (1 + rho) / 2
    q = p

    # Build transition matrix
    transition_matrix = build_transition_matrix(N, p, q)

    return state_values, transition_matrix


def simulate(grid, pmat, T):
    """
    Simulates a discrete Markov chain given a transition matrix.
    
    Parameters:
    grid : numpy array - Discretized state values
    pmat : numpy array - Transition probability matrix
    T : int - Number of time periods
    
    Returns:
    y : numpy array - Simulated AR(1) process
    """
    N = len(grid)
    state0 = np.random.choice(N, p=np.ones(N) / N)  # Random initial state

    # Compute cumulative transition probabilities
    cmat = np.cumsum(pmat, axis=1)

    # Initialize container for simulation
    y = np.zeros(T)

    # Simulation loop
    for i in range(T):
        y[i] = grid[state0]  # Store current state value
        r = np.random.uniform()
        state0 = np.searchsorted(cmat[state0, :], r)  # Determine next state

    return y


# Input parameters
N = 7
mu = 0.5
rho_list = [0.75, 0.85, 0.95, 0.99]  # Multiple rho values
sigma = 1
T = 50

# Seed for reproducibility
seed = 2025
np.random.seed(seed)

# Dictionary to store simulations
sim = {}

# Loop through different rho values and simulate AR(1) processes
plt.figure(figsize=(10, 6))

for rho in rho_list:
    state_values, transition_matrix = rouwenhorst(N, rho, sigma, mu)
    sim[f"rho_{rho}"] = simulate(state_values, transition_matrix, T)
    
    # Plot the simulated time series
    plt.plot(range(T), sim[f"rho_{rho}"], label=f"Rho = {rho}")

# Plot state grid as reference
for state in state_values:
    plt.axhline(y=state, color='gray', linestyle=':', alpha=0.5)

plt.xlabel('Time')
plt.ylabel('Y')
plt.title("AR(1) Process with Rouwenhorst for Different Rho Values")
plt.legend()
plt.grid()
plt.show()
