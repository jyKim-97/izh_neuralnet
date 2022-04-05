import numpy as np
import matplotlib.pyplot as plt
import nnpy.izh_tools as it


if __name__=="__main__":
    import sys

    if len(sys.argv) == 2:
        tag = sys.argv[1]
        print(f"Print result of the {tag}")
    else:
        print("Please input correct tag argument")
    
    obj = it.IzhReader(tag)

    fig = plt.figure(dpi=150, figsize=(8, 5))
    it.draw_single_summary(obj, xlim=(obj.tmax-1000, obj.tmax), xplim=(1000, obj.tmax-1000))
    fname = "./figs/"+tag[tag.rfind("/"):]+".png"
    fig.savefig(fname)
    