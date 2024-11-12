from plotly import graph_objects as go

import numpy as np
class ET_model():
    # Class that implements the ET model computations and adds a function for inserting either effector or target cells

    def __init__(self, E, T, p, m, n, r, k, c, v, s, d, delta_t = 0.01):
        self.E = E  # Initial number of effector cells (immune cells)
        self.T = T  # Initial number of target cells (disease cells)
        self.p = p  # Growth rate of effector cells based on target cell count
        self.m = m  # Half-saturation constant for target cell stimulation
        self.n = n  # Hill coefficient for effector cell proliferation 
        self.r = r  # Natural growth rate of target cells
        self.k = k  # Killing rate of target cells by effector cells
        self.c = c  # Half-saturation constant for effector cell self-regulation 
        self.v = v  # Hill coefficient for target cell influence on effector cell growth
        self.s = s  # Rate of effector cell self-renewal
        self.d = d  # Death rate of effector cells
        self.delta_t = delta_t # Timestep for one iteration in the model
        self.num_of_iterations=0 # Number iterations
        self.results=[] # Dict for results colection

    def _update(self):
        # Compute rates of change for T and E
        dT_dt = self.r * self.T - self.k * self.T * self.E
        dE_dt = (
            self.p * ((self.T**self.v) / (self.m**self.v + self.T**self.v)) +
            self.s * (self.E**self.n / (self.c**self.n + self.E**self.n)) -
            self.d * self.E
        )
        self.T += dT_dt * self.delta_t
        self.E += dE_dt * self.delta_t

    def observe(self, iterations, modulation):
        # Observing n states of the system
        modulation_points, modulation_schedulge = self._plan_modulation(modulation)
        
        for iteration in range(iterations):
            self._update()
            if iteration in modulation_points:
                self._insert_injection(modulation_schedulge[iteration])
            self.results.append([self.E, self.T])
            self.num_of_iterations += 1
        return self.results
    
    
    def _plan_modulation(modulation):
        # Adjust the modulation triplets into apropriate form
        modulation_schedulge = {}
        for triplet in modulation:
            modulation_schedulge[triplet[1]] = [triplet[0], triplet[2]]
        modulation_points = modulation_schedulge.keys()
        return modulation_points, modulation_schedulge
    
    def _insert_injection(self, modulation):
        # Inserts an injection of either effector or target cells,
        # While keeping the null condition for cell count for negative injections (eg. killing each cell type)

        if modulation[0] == "E":
            self.E += modulation[1]
            if self.E < 0:
                self.E = 0

        if modulation[0] == "T":
            self.T += modulation[1] 
            if self.E < 0:
                self.E = 0
    
    def plot_scatter(self):
        # Scatter plot výsledků
        E_values = [res[0] for res in self.results]
        T_values = [res[1] for res in self.results]
        iterations = list(range(1, len(self.results)))
        fig_xy = go.Figure()
        fig_xy.add_trace(go.Scatter(x=iterations, y=E_values, mode='lines+markers', name='Effector Cells (E)', line=dict(color='blue')))
        fig_xy.add_trace(go.Scatter(x=iterations, y=T_values, mode='lines+markers', name='Target Cells (T)', line=dict(color='green')))
        fig_xy.update_layout(
            title="Effector Cells (E) and Target Cells (T) Over Time",
            xaxis_title="Iterations (Time)",
            yaxis_title="Cell Counts",
            legend_title="Cell Type"
        )
        fig_xy.show()


# přidat fázový diagram