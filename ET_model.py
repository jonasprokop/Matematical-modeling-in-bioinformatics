from plotly import graph_objects as go
import numpy as np

class ET_model():
    # Class that implements the ET model computations,
    # adds a function for inserting either effector or target cells
    # and plots the results both in a scatter and phase plot

    def __init__(self, E, T, p, m, n, r, k, c, u, v, s, d, delta_t = 0.01):
        self._E = E  # Initial number of effector cells (immune cells)
        self._T = T  # Initial number of target cells (disease cells)
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
        self._num_of_iterations=0 # Number of iterations
        self._results_array=np.array([]) # Array for results colection

    def _update(self, T_prev, E_prev):
        # Compute rates of change for T and E
        dT_dt = self._r * T_prev - self._k * T_prev * E_prev
        dE_dt = (
            self._p * ((T_prev**self._u) / (self._m**self._v + T_prev**self._v)) +
            self._s * (E_prev**self._n / (self._c**self._n + E_prev**self._n)) -
            self._d * E_prev
        )

        # Apply updates based on timestep 
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
        self._results_array = np.zeros((iterations, 2))
        T = np.zeros(iterations)
        E = np.zeros(iterations)
        T[0] = self._T
        E[0] = self._E
        T_modulation, E_modulation = self._prepare_modulation_schedule(iterations, modulation)
        
        for iteration in range(1, iterations):
            # First the initial state of the system is noted
            self._results_array[iteration] = [self._E, self._T]
            # Then planned modulation is applied based on preprepared modulation arrays
            T_prev = max(0, T[iteration - 1] + T_modulation[iteration])
            E_prev = max(0, E[iteration - 1] + E_modulation[iteration])
            # At last the system updates to the next state
            T[iteration], E[iteration] = self._update(T_prev, E_prev)
        
        
        self._iterations_aray = np.arange(self._iterations)
        self._results_array = np.vstack((E, T)).T
            
        return self._results_array

    
    def plot_scatter(self):
        # Scatter plot 

        E_values, T_values = self._results_array[:, 0], self._results_array[:, 1]

        fig_xy = go.Figure()
        fig_xy.add_trace(go.Scatter(x=self._iterations_aray, y=E_values, mode='lines+markers', name='Effector cells (E)', line=dict(color='blue')))
        fig_xy.add_trace(go.Scatter(x=self._iterations_aray, y=T_values, mode='lines+markers', name='Target cells (T)', line=dict(color='green')))

        fig_xy.update_layout(
            title="Effector cells (E) and target cells (T) over time",
            xaxis_title="Iterations (Time)",
            yaxis_title="Cell counts",
            legend_title="Cell type"
        )
        fig_xy.show()

        return fig_xy


    def plot_phase_diagram(self):
        # Phase diagram plot
        
        E_values, T_values = self._results_array[:, 0], self._results_array[:, 1]

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=T_values, y=E_values, mode='lines', name='Trajectory',
            line=dict(color='royalblue', width=2)
        ))

        fig.update_layout(
            title="Phase diagram: Counts of effector cells (E) vs target cells (T)",
            xaxis_title="Count of target cells (T)",
            yaxis_title="Count of effector cells (E)",
            showlegend=False
        )

        fig.show()

        return fig


simple_model = ET_model(E=1, T=5*10**4, p=0.7, m=1, n=3, r=0.15, k=0.1, c=1, u = 1, v=2, s=2, d=1, delta_t = 0.01)
simple_model.observe(10**3, modulation=[["E", 5*10**2, 5*10**1]])
simple_model.plot_scatter()
simple_model.plot_phase_diagram()