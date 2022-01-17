import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append('./include/')
import nnpy.izh_tools as it
import nnpy.visu_tools as visu


if __name__ == "__main__":
    tag = "./data/single_ntk"
    nid = 10
    
    obj = it.IzhReader(tag)

    kwargs = {"lw":1, "color":"k"}
    
    def set_ax():
        plt.xlabel("time (ms)", fontsize=10)
        plt.xticks(fontsize=8)
        plt.yticks(fontsize=8)
    
    plt.figure(dpi=150, figsize=(9, 5))
    plt.subplot(231)
    plt.plot(obj.ts, obj.vs[nid], **kwargs)
    set_ax()
    plt.ylabel("v (mV)", fontsize=10)
    
    plt.subplot(232)
    plt.plot(obj.ts, obj.us[nid], **kwargs)
    set_ax()
    plt.ylabel("u", fontsize=10)
    
    plt.subplot(233)
    plt.plot(obj.ts, obj.ics[nid], **kwargs)
    set_ax()
    plt.ylabel("$i_{cell}$", fontsize=10)
    
    plt.subplot(212)
    visu.draw_raster_plot(obj.t_spks, xlim=[0, obj.tmax], ylim=[-1, obj.num_cells], cell_types=obj.cell_types)
    set_ax()
    plt.ylabel("cell id", fontsize=10)
    
    plt.tight_layout()
    plt.show()
