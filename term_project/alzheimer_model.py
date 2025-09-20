import numpy as np
import pandas as pd
import yaml
from matplotlib import pyplot as plt 

class AlzheimerModel:
    def __init__(self, path_params: str, path_initial_state: str, end_time: float = 20,  start_time: float = 0.0, time_step: float = 0.25):
        
        params = self.load_config(path_params)
        initial_state = self.load_config(path_initial_state)
        
        self.params = params
        self.initial_state = initial_state
        self.state = initial_state.copy()
        self.time = start_time 
        self.end_time = end_time
        self.time_step = time_step
        self.steps = (self.end_time - start_time) / time_step

        # initialize the constasnt govering the systems dynamics
        self.initialize_params(params)

        # initialize the state of the system at time 0 
        self.initialize_state(initial_state)

    def initialize_params(self, params: list):
        """
        Initialize the model parameters based on the provided list.
        """
        params = params["values"]
        self.alfa1 = params[0]
        self.alfa2 = params[1]
        self.alfa3 = params[2]
        self.alfa4 = params[3]
        self.alfa5 = params[4]
        self.alfa6 = params[5]
        self.alfa7 = params[6]
        self.alfa8 = params[7]
        self.alfa9 = params[8]
        self.alfa10 = params[9]
        self.alfa11 = params[10]
        self.alfa12 = params[11]
        self.alfa13 = params[12]
        self.alfa14 = params[13]
        self.alfa15 = params[14]
        self.alfa16 = params[15]
        self.alfaR = params[16]
    
    def initialize_state(self, initial_state: list):
        """
        Initialize the model's state variables based on the provided initial conditions.
        The initial state of the system is expected to be in the following order:
        IC0, IC1, IC3, IC4, IC5, IC6, IC7
        where:
        IC0 = Neuronal survival (NS)    
        IC1 = Neuron death (ND)
        IC3 = Astrocyte quiescent (AQ)
        IC4 = Astrocyte rapidly proliferating (AP)
        IC5 = Microglia type 1 (M1)
        IC6 = Microglia type 2 (M2)
        IC7 = Amyloid beta (AB)

        And then saves them into state transition dicts
        """
        initial_state = initial_state["values"]

        # initialize the system states
        self.NS = initial_state[0]
        self.ND = initial_state[1]
        self.AQ = initial_state[2]
        self.AP = initial_state[3]
        self.M1 = initial_state[4]
        self.M2 = initial_state[5]
        self.AB = initial_state[6]

        # initialise the state transition arrays
        self.states_NS = {}
        self.states_ND = {}
        self.states_AQ = {}
        self.states_AP = {}
        self.states_M1 = {}
        self.states_M2 = {}
        self.states_AB = {}

    def save_state(self, iteration):
        self.states_NS[iteration] = self.NS
        self.states_ND[iteration] = self.ND
        self.states_AQ[iteration] = self.AQ
        self.states_AP[iteration] = self.AP
        self.states_M1[iteration] = self.M1
        self.states_M2[iteration] = self.M2
        self.states_AB[iteration] = self.AB

    def update(self):
        """
        Update the model state based on the current parameters and conditions.
        """
        dNS = (self.alfa1 * self.AQ - self.alfa2 * self.AP - self.alfa3 * self.M1) 

        dND = - dNS

        dAQ = (self.alfa4 * self.M2 - self.alfa5 * self.M1) 

        dAP = -dAQ

        dM2 = ((self.alfa6 + self.alfa11) * self.NS
            - self.alfa10 * self.ND
            + (self.alfa7 + self.alfa12) * self.AQ
            - self.alfa9 * self.M1
            + self.alfa14 * self.M2
            - (self.alfa8 + self.alfa13) * self.AB)
        
        dM1 = - dM2 

        dAB = (self.alfa15 * self.NS
            - self.alfa16 * self.M2
            - self.alfaR * self.AB)
        

        # Update each population (Euler method)
        self.NS = max(0, self.NS + dNS * self.time_step)
        self.ND = max(0, self.ND + dND * self.time_step)
        self.AQ = max(0, self.AQ + dAQ * self.time_step)
        self.AP = max(0, self.AP + dAP * self.time_step)
        self.M2 = max(0, self.M2 + dM2 * self.time_step)
        self.M1 = max(0, self.M1 + dM1 * self.time_step)
        self.AB = max(0, self.AB + dAB * self.time_step)


    def simulate(self):
        """
        Simulates the time evolution of the system, each time the t-1 state is saved 
        """
        for iteration in range(int(self.steps)):
            self.save_state(iteration)
            self.update()

    def load_config(self, path):
        with open(path, 'r') as file:
            data = yaml.safe_load(file)     
            return data
        
    def show_results(self):
        return {
        "Neuron survival": self.states_NS,
        "Neuronal death": self.states_ND,
        "Astrocytes quiescent": self.states_AQ,
        "Astrocytes proliferating": self.states_AP,
        "Microglia M1": self.states_M1,
        "Microglia M2": self.states_M2,
        "Amyloid Beta": self.states_AB
        }


    def plot_scatter(self):
        population_NS = self.states_NS
        population_M1 = self.states_M1
        population_AB = self.states_AB

        plt.scatter(population_NS.keys(), population_NS.values(), color='red', label='NS')
        
        plt.scatter(population_M1.keys(), population_M1.values(), color='blue', label='M1')
        
        plt.scatter(population_AB.keys(), population_AB.values(), color='black', label='AB')

        plt.xlabel('Time Steps')
        plt.ylabel('Population')
        plt.title('Scatter Plot of Populations Over Time')
        plt.legend()
        plt.tight_layout()
        plt.show()


    def perturb_param(self, param_name: str, factor: float = 2.0):
        """
        Perturb a single parameter (e.g., 'alfa3') by multiplying it with a factor.
        Returns a new AlzheimerModel instance with the perturbed parameter.
        """
        perturbed = AlzheimerModel.__new__(AlzheimerModel)
        perturbed.__dict__ = self.__dict__.copy() 

        if hasattr(perturbed, param_name):
            setattr(perturbed, param_name, getattr(perturbed, param_name) * factor)
        else:
            raise ValueError(f"Parameter {param_name} not found in model.")

        return perturbed


    def perturb_initial(self, species: str, factor: float = 10.0):
        """
        Perturb initial condition of a species (e.g., 'NS', 'M1') by a factor.
        Returns a new AlzheimerModel instance with perturbed initial state.
        """
        perturbed = AlzheimerModel.__new__(AlzheimerModel)
        perturbed.__dict__ = self.__dict__.copy()

        if hasattr(perturbed, species):
            setattr(perturbed, species, getattr(perturbed, species) * factor)
        else:
            raise ValueError(f"Species {species} not found in model.")

        return perturbed


    def compare_runs(self, models: dict, variable: str):
        """
        Compare trajectories of a given variable across multiple model runs.
        'models' should be a dict: {"label": model_instance, ...}
        """
        for label, model in models.items():
            model.simulate()
            data = getattr(model, f"states_{variable}")
            plt.plot(data.keys(), data.values(), label=label)

        plt.xlabel('Time Steps')
        plt.ylabel(variable)
        plt.title(f'Comparison of {variable} under Perturbations')
        plt.legend()
        plt.tight_layout()
        plt.show()