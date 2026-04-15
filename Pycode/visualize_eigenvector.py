import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class Eig_visualizer:
    def __init__(self, matrix, history, eigvector, fig, ax, path_vis):
        self.matrix = matrix
        self.history = history
        self.eigvector = eigvector
        self.fig = fig
        self.ax = ax
        self.anim = None
        self.path_vis = path_vis
        
    def visualize_eigenvector(self):
        iterations = [np.array(vec) for vec, _ in self.history]
        times = [t for _, t in self.history]
        
        self.ax.set_xlim(-1.5, 1.5)
        self.ax.set_ylim(-1.5, 1.5)
        self.ax.set_aspect('equal')
        self.ax.grid(True, alpha=0.3)
        
        circle = plt.Circle((0, 0), 1, fill=False, color='gray', linestyle='-', linewidth=1.5)
        self.ax.add_patch(circle)
        
        self.ax.annotate('', xy=(self.eigvector[0], self.eigvector[1]), xytext=(0, 0),
                    arrowprops=dict(arrowstyle='->', color='red', linestyle='--', lw=2))

        current_arrow = self.ax.annotate('', xy=(1, 1), xytext=(0, 0),
                                   arrowprops=dict(arrowstyle='->', color='blue', lw=2.5,
                                                   mutation_scale=20))
        
        info_text = self.ax.text(-1.4, 1.3, '', fontsize=10)

        self.ax.plot([], [], 'r--', lw=2, label='Настоящий собственный вектор')
        self.ax.plot([], [], 'b--', lw=4, label='Текущая итерация')
        self.ax.legend(loc='upper right')
        
        def update(frame):
            if frame == len(iterations):
                if self.anim:
                    self.anim.event_source.stop()
                    plt.cla()
                    self.path_vis.visualize_path()
                return
            
            nonlocal current_arrow
            current_arrow.remove()
            
            current_arrow = self.ax.annotate('', xy=(iterations[frame][0], iterations[frame][1]), xytext=(0, 0),
                                       arrowprops=dict(arrowstyle='->', color='blue', lw=2.5,
                                                       mutation_scale=20))
            
            info_text.set_text(f'Iteration: {frame + 1}/{len(iterations)}\nTime: {times[frame]:.2f} ms')
            
            return current_arrow, info_text

        interval_ms = 50
        self.anim = animation.FuncAnimation(
            self.fig, update, frames=len(iterations)+1,
            interval=interval_ms, blit=False
        )
        
        plt.title('Eigenvector Finding Visualization')
