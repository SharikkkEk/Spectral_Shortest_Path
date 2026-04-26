import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from SpectralPath import rand_sparse_graph, spectral_history, get_adjacency_list, dijkstra_history, eigenvector_2D, eigenvector_2D_history, rand_tree, rand_dense_graph, spectral, dijkstra

from visualize_eigenvector import Eig_visualizer
matrix = rand_dense_graph(50)

class Path_visualizer:
    def __init__(self, matrix, path, fig, ax):
        self.matrix = matrix
        self.path = path
        self.fig = fig
        self.ax = ax
        self.ani = None
        
    def visualize_path(self):        
        G = nx.Graph()
        for i, neighbors in enumerate(get_adjacency_list(self.matrix)):
            for neighbor, weight in neighbors:
                G.add_edge(i, neighbor, weight=weight)
        
        pos = nx.spring_layout(G, k = 0.7, seed = 42)
        path_vertices = [v for v, _ in self.path]
        path_times = [t for _, t in self.path]

        self.ax.set_xlim(-1.2, 1.2)
        self.ax.set_ylim(-1.2, 1.2)
        self.ax.axis('off')
        
        base_edges = nx.draw_networkx_edges(G, pos, ax=self.ax, edge_color='gray', 
                                            alpha=0.3, width=1)
        base_nodes = nx.draw_networkx_nodes(G, pos, ax=self.ax, node_color='lightgray', 
                                            node_size=200)
        
        labels = {i: str(i) for i in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels, ax=self.ax, font_size=8)
        edge_labels = nx.get_edge_attributes(G, 'weight')
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels,
                                    font_size=8, font_color='red')
        
        info_text = self.ax.text(-1.1, 1.1, '', fontsize=10)
        
        def update(frame):
            if frame > 0 and frame < len(path_vertices):
                prev_v = path_vertices[frame - 1]
                curr_v = path_vertices[frame]
                
                x1, y1 = pos[prev_v]
                x2, y2 = pos[curr_v]
                self.ax.plot([x1, x2], [y1, y2], 'b-', lw=3, alpha=0.7)
            elif frame == len(path_vertices):
                if self.ani:
                    self.ani.event_source.stop()                        
                return
            
            current_node = path_vertices[frame]
            nx.draw_networkx_nodes(G, pos, ax=self.ax, nodelist=[current_node], 
                                   node_color='blue', node_size=300)
            
            info_text.set_text(f'Vertex: {current_node}\n'
                             f'Step: {frame + 1}/{len(path_vertices)}\n'
                             f'Sum Dijkstra: {int(dijkstra(self.matrix))}\n'
                             f'Sum Spectral: {int(spectral(self.matrix))}')
            
            return info_text

        plt.title('Path Visualization')
        plt.tight_layout()

        interval_ms = max(50, int(path_times[-1] / len(path_vertices) * 5))
        self.ani = animation.FuncAnimation(
            self.fig, update, frames=len(path_vertices)+1,
            interval=interval_ms, blit=False
        )

if __name__ == '__main__':
    fig_d, ax_d = plt.subplots(figsize=(10, 10))
    mgr = plt.get_current_fig_manager()
    mgr.window.wm_geometry("800x900+0+0")
    
    dijkstra_vis = Path_visualizer(matrix, dijkstra_history(matrix), fig_d, ax_d)
    dijkstra_vis.visualize_path()
    
    fig_s, ax_s = plt.subplots(figsize=(10, 10))
    mgr = plt.get_current_fig_manager()
    mgr.window.wm_geometry("800x900+800+0")

    spectral_vis = Path_visualizer(matrix, spectral_history(matrix), fig_s, ax_s)
    eig_vis = Eig_visualizer(matrix, eigenvector_2D_history(matrix), eigenvector_2D(matrix), fig_s, ax_s, spectral_vis)
    eig_vis.visualize_eigenvector()

    plt.show()
    
    
