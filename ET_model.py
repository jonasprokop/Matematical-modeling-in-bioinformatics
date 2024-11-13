from plotly import graph_objects as go

import numpy as np
class ET_model():
    # Class that implements the ET model computations and adds a function for inserting either effector or target cells

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
        self._results_array=np.array # Array for results colection
        

    def _update(self):
        # Compute rates of change for T and E
        dT_dt = self._r * self._T - self._k * self._T * self._E
        dE_dt = (
            self._p * ((self._T**self._u) / (self._m**self._v + self._T**self._v)) +
            self._s * (self._E**self._n / (self._c**self._n + self._E**self._n)) -
            self._d * self._E
        )
        self._T += dT_dt * self._delta_t
        self._E += dE_dt * self._delta_t
        
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
        self._results_array = np.zeros((iterations, 2))
        T_modulation, E_modulation = self._prepare_modulation_schedule(iterations, modulation)
        
        for iteration in range(iterations):
            # First the initial state of the system is noted
            self._results_array[iteration] = [self._E, self._T]
            # Then planned modulation is applied
            self._T = max(0, self._T + T_modulation[iteration])
            self._E = max(0, self._E + E_modulation[iteration])
            # The system updates to the next state
            self._update()
            
            
        return self._results_array

    
    def plot_scatter(self):
        # Scatter plot 
        E_values, T_values = self._results_array[:, 0], self._results_array[:, 1]
        fig_xy = go.Figure()
        fig_xy.add_trace(go.Scatter(x=self._iterations_aray, y=E_values, mode='lines+markers', name='Effector Cells (E)', line=dict(color='blue')))
        fig_xy.add_trace(go.Scatter(x=self._iterations_aray, y=T_values, mode='lines+markers', name='Target Cells (T)', line=dict(color='green')))
        fig_xy.update_layout(
            title="Effector Cells (E) and Target Cells (T) Over Time",
            xaxis_title="Iterations (Time)",
            yaxis_title="Cell Counts",
            legend_title="Cell Type"
        )
        fig_xy.show()


    def plot_phase_diagram(self):
        
        E_values = self._results_array[:, 0]
        T_values = self._results_array[:, 1]
        

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=T_values, y=E_values, mode='lines', name='Trajectory',
            line=dict(color='royalblue', width=2)
        ))

        fig.update_layout(
            title="Phase Diagram: Counts of Effector Cells (E) vs Target Cells (T), both axes are reversed",
            xaxis_title="Count of Target Cells (T), reversed",
            yaxis_title="Count of Effector Cells (E), reversed",
            showlegend=False
        )

        fig.show()



simple_model = ET_model(E=0.1, T=1000, p=0.7, m=1, n=3, r=0.15, k=0.1, c=1, u = 1, v=2, s=2, d=1, delta_t = 0.01)
simple_model.observe(10**4, modulation=[])
simple_model.plot_scatter()
simple_model.plot_phase_diagram()

