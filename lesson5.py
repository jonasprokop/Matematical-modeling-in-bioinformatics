import random
import matplotlib.pyplot as plt

value_grid = {0:0, 1:1, 2:0, 3:1, 4:1, 5:0, 6:1, 7:0}
entry_grid = {i: random.randint(0, 1) for i in range(500)}

def convert_to_binary_compare(input_list, value_grid):
    binary_str = ''.join(['1' if bit else '0' for bit in input_list])

    binary_int = int(binary_str, 2)

    if binary_int in value_grid:
        result = value_grid[binary_int]

    return result

def evolve_grid(n, value_grid, entry_grid):
    new_grid = {}
    val_1 = False
    val_2 = False
    for key, value in entry_grid.items():
        val_1 = entry_grid[(key - 2) % n]  
        val_2 = entry_grid[(key + 2) % n]  

        new_val = convert_to_binary_compare([val_1, value, val_2], value_grid)  

        new_grid[key] = new_val

    return new_grid

def evolve_grid_n_times(epochs, value_grid, entry_grid):
    results = []
    for epoch in range(epochs):
        entry_grid = evolve_grid(len(entry_grid), value_grid, entry_grid)
        results.append(list(entry_grid.values()))
    return results

# Running the evolution for 500 iterations
results = evolve_grid_n_times(500, value_grid, entry_grid)

# Visualize the evolution of the grid
plt.figure(figsize=(15, 10))
plt.imshow(results, cmap='binary', interpolation='nearest', aspect='auto')
plt.colorbar(label="State (0 or 1)")
plt.xlabel("Cell Index")
plt.ylabel("Epoch")
plt.title("Evolution of Cellular Automaton")
plt.show()