import random
import matplotlib.pyplot as plt
import numpy as np



# Rule 90
value_grid = {0:0, 1:1, 2:0, 3:1, 4:1, 5:0, 6:1, 7:0}


# dict implementation

# entry_grid = {i: random.randint(0, 1) for i in range(500)}

# def convert_to_binary_compare(input_list, value_grid):
#     binary_str = ''.join(['1' if bit else '0' for bit in input_list])

#     binary_int = int(binary_str, 2)

#     if binary_int in value_grid:
#         result = value_grid[binary_int]

#     return result


# # evolve grid from state t to state t+1
# def evolve_grid(value_grid, entry_grid):
#     new_grid = {}
#     val_1 = False
#     val_2 = False
#     n = len(entry_grid)
#     for key, value in entry_grid.items():
#         val_1 = entry_grid[(key - 2) % n]  
#         val_2 = entry_grid[(key + 2) % n]  

#         new_val = convert_to_binary_compare([val_1, value, val_2], value_grid)  

#         new_grid[key] = new_val

#     return new_grid

# # iterate evolution n times
# def evolve_grid_n_times(epochs, value_grid, entry_grid):
#     results = []
#     for epoch in range(epochs):
#         entry_grid = evolve_grid(value_grid, entry_grid)
#         results.append(list(entry_grid.values()))
#     return results


# # Running the evolution for 500 iterations
# results = evolve_grid_n_times(500, value_grid, entry_grid)


# # Visualize the evolution of the grid
# plt.figure(figsize=(15, 10))
# plt.imshow(results, cmap='binary', interpolation='nearest', aspect='auto')
# plt.colorbar(label="State (0 or 1)")
# plt.xlabel("Cell Index")
# plt.ylabel("Epoch")
# plt.title("Evolution of Cellular Automaton")
# plt.show()



# numpy implemntation

entry_grid = np.array([random.randint(0, 1) for _ in range(500)])

def convert_to_binary_compare(input_list, value_grid):
    binary_str = ''.join(['1' if bit else '0' for bit in input_list])

    binary_int = int(binary_str, 2)

    if binary_int in value_grid:
        result = value_grid[binary_int]

    return result


# evolve grid from state t to state t+1
def evolve_grid_numpy(value_grid, entry_grid):
    left_neighbors = np.roll(entry_grid, -2)
    right_neighbors = np.roll(entry_grid, 2)

    new_grid = np.array([
        convert_to_binary_compare([left_neighbors[i], entry_grid[i], right_neighbors[i]], value_grid)
        for i in range(len(entry_grid))
    ])

    return new_grid

# iterate evolution n times
def evolve_grid_n_times_numpy(epochs, value_grid, entry_grid):
    results = [entry_grid]  
    for _ in range(epochs):
        entry_grid = evolve_grid_numpy(value_grid, entry_grid)
        results.append(entry_grid)  
    return np.array(results) 


results = evolve_grid_n_times_numpy(500, value_grid, entry_grid)


# Visualize the evolution of the grid
plt.figure(figsize=(15, 10))
plt.imshow(results, cmap='binary', interpolation='nearest', aspect='auto')
plt.colorbar(label="State (0 or 1)")
plt.xlabel("Cell Index")
plt.ylabel("Epoch")
plt.title("Evolution of Cellular Automaton")
plt.show()
