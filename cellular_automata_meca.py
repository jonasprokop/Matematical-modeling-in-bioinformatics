import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt


# Class that simulates and visualises the Meca automaton
class Meca_Cellular_Automata:
    def __init__(self, number_of_epochs, len_intial_grid
                ):
        
        self._number_of_epochs = number_of_epochs
        self._intial_grid = self._create_intial_grid(len_intial_grid)
        self._all_results = {}

    # Creates random intial grid to evolve of the selected length
    def _create_intial_grid(self, len_intial_grid):
        return np.array([random.randint(0, 1) for _ in range(len_intial_grid)])
    
    # Generate the value grid for rule
    def _generate_value_grid(self, rule_number):
        binary_repr = f"{rule_number:08b}" 
        return {i: int(binary_repr[7 - i]) for i in range(8)}

    # Compare the binary inputs and return result based on inputed rule set
    def _convert_to_binary_compare(self, input_list, value_grid):
        binary_int = (input_list[0] << 2) | (input_list[1] << 1) | input_list[2]
        return value_grid[binary_int]


    # Evolve grid from state t to state t+1
    def _evolve_grid_meca(self, value_grid, prev_grid, curr_grid):
        left_neighbors = np.roll(prev_grid, -1)
        right_neighbors = np.roll(prev_grid, 1)

        new_grid = np.array([
            self._convert_to_binary_compare([left_neighbors[i], prev_grid[i], right_neighbors[i]], value_grid)
            for i in range(len(curr_grid))
        ])

        return new_grid

    # Iterate through the grid evolution n-times
    def _evolve_grid_n_times_meca(self, epochs, intial_grid, value_grid):
        prev_grid = intial_grid.copy()
        curr_grid = intial_grid.copy()
        results = [intial_grid.copy()]

        for _ in range(epochs):
            new_grid = self._evolve_grid_meca(value_grid, prev_grid, curr_grid)
            prev_grid = curr_grid.copy()
            curr_grid = new_grid.copy()
            results.append(curr_grid)

        return results

    # Visualize the evolution of the grid
    def _plot_binary_map(self, results, rule_number):
        plt.figure(figsize=(15, 10))
        plt.imshow(results, cmap='binary', interpolation='nearest', aspect='auto')
        plt.colorbar(label="State (0 or 1)")
        plt.xlabel("Cell Index")
        plt.ylabel("Epoch")
        value_grid = self._generate_value_grid(rule_number)
        plt.title(f"Evolution of Memory-Based Cellular Automaton based on rule number: {rule_number} \n {value_grid}")
        plt.show()

    # Simulate all of the defined rules  
    def _simulate_all_rules(self):
        for rule_number in range(256):
            value_grid = self._generate_value_grid(rule_number)
            results = self._evolve_grid_n_times_meca(self._number_of_epochs, self._intial_grid, value_grid)
            self._all_results[rule_number] = results
            if rule_number % 10 == 0 and rule_number > 0:
                print(f"I have solved grids for rules num {rule_number-10} - {rule_number}")
            elif rule_number == 255:
                print("I have finished solving the whole space of grids")


    # Visualises the selected rule based on precomputed results
    def view_rule_number(self, rule_number):
        results = self._all_results[rule_number]
        if not results:
            print("Grids must be solved prior to visualisation")
        self._plot_binary_map(results, rule_number)

    # Visualises all rules based on precomputed term (use with caution due to memory requirements)
    def iterate_and_view_all_rules(self):
        for rule_num in range(256):
            self.view_rule_number(rule_num)