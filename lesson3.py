iterations = 20
n0 = 10

def observe(n0):

    n1 = n0/2
    n2 = n1/3
    n3 = n2*6
    return n1, n2, n3

n1s = []
n2s = []
n3s = []

for _ in range(iterations):
    n1, n2, n3 = observe(n0)
    n3 = n0
    n1s.append(n1)
    n2s.append(n2)
    n3s.append(n3)


print(n1s, n2s, n3s)


values = []
def observe(n0):

    n1 = -0.3 * n0
    n2 = 0.2 * n1
    n3 = 0.02 * n2
    n4 = 2 * n3
    return n1, n2, n3, n4

for _ in range(iterations):
    n1, n2, n3, n4 = observe(n0)
    
    values.append([n1, n2, n3, n4])

    n0 = n4  
 
print(values)


from plotly import graph_objects as go
def _plot_scatter(results):
    #scatter plot výsledků
    E_values = [res[0] for res in results]
    T_values = [res[1] for res in results]
    g_values = [res[2] for res in results]
    h_values = [res[3] for res in results]
    iterations = list(range(1, len(results)))
    fig_xy = go.Figure()
    fig_xy.add_trace(go.Scatter(x=iterations, y=E_values, mode='lines+markers', name='Effector Cells (E)', line=dict(color='blue')))
    fig_xy.add_trace(go.Scatter(x=iterations, y=T_values, mode='lines+markers', name='Target Cells (T)', line=dict(color='green')))
    fig_xy.add_trace(go.Scatter(x=iterations, y=g_values, mode='lines+markers', name='Target Cells (T)', line=dict(color='yellow')))
    fig_xy.add_trace(go.Scatter(x=iterations, y=h_values, mode='lines+markers', name='Target Cells (T)', line=dict(color='red')))
    fig_xy.update_layout(
        title="Effector Cells (E) and Target Cells (T) Over Time",
        xaxis_title="Iterations (Time)",
        yaxis_title="Cell Counts",
        legend_title="Cell Type"
    )
    fig_xy.show()








_plot_scatter(values)