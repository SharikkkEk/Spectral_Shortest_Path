import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from SpectralPath import rand_sparse_graph, eigenvector, eigenvector_history, test



def visualize_eigenvector(matrix):
    history = eigenvector_history(matrix)
    final_eigenvector = eigenvector(matrix)
    plt.figure()
    plt.plot([0, 0], [final_eigenvector[0], final_eigenvector[1]])
    plt.show()

    import matplotlib; import matplotlib.pyplot as plt; plt.figure(); plt.plot([1,2],[1,2]); plt.savefig('test.png'); print('OK')
    '''iterations = [np.array(vec) for vec, _ in history]
    times = [t for _, t in history]
    
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    
    circle = plt.Circle((0, 0), 1, fill=False, color='gray', linestyle='-', linewidth=1.5)
    ax.add_patch(circle)
    
    ax.annotate('', xy=(final_eigenvector[0], final_eigenvector[1]), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color='red', linestyle='--', lw=2))
    
    current_arrow = ax.annotate('', xy=(iterations[0][0], iterations[0][1]), xytext=(0, 0),
                               arrowprops=dict(arrowstyle='->', color='blue', lw=2.5,
                                               mutation_scale=20))
    
    info_text = ax.text(-1.4, 1.3, '', fontsize=10)
    
    ax.plot([], [], 'r--', lw=2, label='True eigenvector')
    ax.legend(loc='upper right')
    
    def update(frame):
        current_arrow.remove()
        
        current_arrow = ax.annotate('', xy=(iterations[frame][0], iterations[frame][1]), xytext=(0, 0),
                                   arrowprops=dict(arrowstyle='->', color='blue', lw=2.5,
                                                   mutation_scale=20))
        
        info_text.set_text(f'Iteration: {frame + 1}/{len(iterations)}\nTime: {times[frame]:.2f} ms')
        
        return current_arrow, info_text
    
    anim = animation.FuncAnimation(
        fig, update, frames=len(iterations),
        interval=max(50, 5000 // len(iterations)), blit=False
    )
    
    plt.title('Eigenvector Finding Visualization')
    plt.show()'''

m = eigenvector([[1, 2], [2, 3]])
print(m)
