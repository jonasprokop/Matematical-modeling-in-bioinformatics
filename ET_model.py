from plotly import graph_objects as go
import numpy as np
import matplotlib.pyplot as plt

class ET_model():
    # Class that implements the ET model computations,
    # adds a function for inserting either effector or target cells
    # and plots the results both in a scatter and trajectory plot

    def __init__(self, E, T, p, m, n, r, k, c, u, v, s, d, delta_t = 0.01):
        self._E0 = E  # Initial number of effector cells (immune cells)
        self._T0 = T  # Initial number of target cells (disease cells)
        self._p = p  # Growth rate of effector cells based on target cell count
        self._m = m  # Half-saturation constant for target cell stimulation
        self._n = n  # Hill coefficient for effector cell proliferation 
        self._r = r  # Natural growth rate of target cells
        self._k = k  # Killing rate of target cells by effector cells
        self._c = c  # Half-saturation constant for effector cell self-regulation 
        self._u = u  # Rate at which the immune system responds to changes in T
        self._v = v  # Hill coefficient for target cell influence on effector cell growth
        self._s = s  # Rate of effector cell self-renewal
        self._d = d  # Death rate of effector cells
        self._delta_t = delta_t # Timestep for one iteration in the model
        self._num_of_iterations = 0 # Number of iterations
        self._T = np.array([]) # Array for T results colection
        self._E = np.array([]) # Array for E results colection

    def _update(self, T_prev, E_prev):
        # Computes rates of change for T and E
        dT_dt = self._r * T_prev - self._k * T_prev * E_prev
        dE_dt = (
            self._p * ((T_prev**self._u) / (self._m**self._v + T_prev**self._v)) +
            self._s * (E_prev**self._n / (self._c**self._n + E_prev**self._n)) -
            self._d * E_prev
        )

        # Applys updates based on timestep 
        T_next = max(0, T_prev + dT_dt * self._delta_t)
        E_next = max(0, E_prev + dE_dt * self._delta_t)

        return T_next, E_next
        
    def _prepare_modulation_schedule(self, iterations, modulation):
        # Creates modulation arrays for E and T
        T_modulation = np.zeros(iterations)
        E_modulation = np.zeros(iterations)

        for target, point, strength in modulation:
            if target == "T":
                T_modulation[point] = strength
            elif target == "E":
                E_modulation[point] = strength

        return T_modulation, E_modulation
    

    def observe(self, iterations, modulation=[]):
        # Observing n states of the system
        self._iterations = iterations
        self._iterations_aray = np.arange(self._iterations)
        self._T = np.zeros(iterations)
        self._E = np.zeros(iterations)
        self._T[0] = self._T0
        self._E[0] = self._E0
        T_modulation, E_modulation = self._prepare_modulation_schedule(iterations, modulation)
        
        for iteration in range(1, iterations):
            # Planned modulation is applied based on preprepared modulation arrays
            T_prev = max(0, self._T[iteration - 1] + T_modulation[iteration])
            E_prev = max(0, self._E[iteration - 1] + E_modulation[iteration])
            # And the system updates to the next state
            self._T[iteration], self._E[iteration] = self._update(T_prev, E_prev)
        
        self._iterations_aray = np.arange(self._iterations)
        
        return self._T, self._E
    
    def plot_scatter(self):
        # Scatter plot 

        fig_xy = go.Figure()

        fig_xy.add_trace(go.Scatter(x=self._iterations_aray, y=self._E, mode='lines+markers', name='Effector cells (E)', line=dict(color='blue')))
        fig_xy.add_trace(go.Scatter(x=self._iterations_aray, y=self._T, mode='lines+markers', name='Target cells (T)', line=dict(color='green')))

        fig_xy.update_layout(
            title="Effector cells (E) and Target cells (T) over Iterations (Time)",
            xaxis_title="Iterations (Time)",
            yaxis_title="Cell counts",
            legend_title="Cell type"
        )
        fig_xy.show()

    def plot_trajectory(self):
        # Trajectory plot

        plt.figure(figsize=(10, 8))

        plt.plot(self._T, self._E, label="Trajectory", color='blue')
        
        plt.scatter(self._T[0], self._E[0], color='green', label='Start (T0, E0)')
        plt.scatter(self._T[-1], self._E[-1], color='red', label='End (T, E)')
        
        plt.xlabel("Target cells (T)")
        plt.ylabel("Effector cells (E)")
        plt.title("Phase plane trajectory")
        plt.legend()
        plt.grid()
        plt.show()