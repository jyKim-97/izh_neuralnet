import matplotlib.pyplot as plt
import numpy as np


def draw_raster_plot(tspk, xlim=None, ylim=None, colors=None, cell_types=None, s=1):
    if colors is None:
        colors = np.array([[203, 67, 53], [36, 113, 163]])/255
    
    num_cells = len(tspk)
    if cell_types is None:
        cell_types = np.ones(num_cells, dtype=np.int)
    
    # get scatter array
    x, y, c = [], [], []    
    for n in range(num_cells):
        num_fire = len(tspk[n])
        
        x.extend(tspk[n])
        y.extend(n * np.ones(num_fire))
        c.extend(np.tile(colors[cell_types[n]], [num_fire, 1]))
            
    c = np.array(c)
        
    plt.scatter(x, y, s=s, c=c)
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is None:
        ylim = [0, num_cells]
    plt.ylim(ylim)
    plt.xlabel('time (ms)', fontsize=12)
    yt = [0, ylim[1]//2, ylim[1]]
    plt.yticks(yt, labels=['#%d'%(i) for i in yt])


