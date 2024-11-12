from plotly import graph_objects as go

import numpy as np
class ET_model():
    def __init__(self, E, T, p, m, n, r, k, c, u, v, s, d, time_points, strength, delta_t = 0.01, treatment_present = False):
        self.E = E # num efektor
        self.T = T # num target
        self.p = p 
        self.m = m
        self.n = n
        self.r = r
        self.k = k
        self.c = c
        self.u = u
        self.v = v
        self.s = s
        self.d = d
        self.delta_t = delta_t
        self._treatment_present = treatment_present
        self._time_points = time_points
        self.strength = strength

        self.num_of_iterations=0
        self.results=[]

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

    def _observe(self, iterations):
        # observing n states of the system
        for count in range(iterations):
            self._update()
            if self._treatment_present:
                self._insert_injection()
            self.results.append([self.E, self.T])
            self.num_of_iterations += 1
        return self.results

    def _observe_repated_infection(self, iterations, num_of_infections, dose):
        infection_step = int(iterations * (1/(num_of_infections+1)))
        # observing n states of the system
        for count in range(iterations):
            if count % infection_step == 0 and count not in [0, iterations]:
                self.T += dose
            if self._treatment_present:
                self._insert_injection()
            self._update()
            self.results.append([self.E, self.T])
            self.num_of_iterations += 1
        return self.results

    def _insert_injection(self):
        # dosing the treatment
        time = self.delta_t * self.num_of_iterations
        for time_point in self._time_points:
            if time_point == time:
                self.T -= self.strength     
                if self.T < 0:
                    self.T = 0       
    
    def _plot_scatter(self):
        #scatter plot výsledků
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
# bude to potřeba překopat tak, že spojíš oba dva modely, a nastavíš tam variable time injectcí
# ideál bude překopat tu funkci tak, aby se tam dali volit 2x2 mody operandi - eg. volit mezi injekci T a E, a volit mezi regularizovaným dosingem a point time dosingem
#