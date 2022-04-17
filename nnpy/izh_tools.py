from struct import unpack
import numpy as np
import matplotlib.pyplot as plt
import json
import os


class IzhReader:
    def __init__(self, tag):
        self.tag = tag
        self.vs = None
        self.us = None
        self.rs = None
        self.ics = None
        self.t_spks = None
        
        if os.path.isfile(f'{self.tag}_env.txt'):
            self._read_info_txt()
        else:
            self._read_info()
        self.ts = np.arange(0, self.tmax, self.dt)+self.dt
        self.read_all_data()
        self.read_summary()

    def _read_info_txt(self):
        with open(f'{self.tag}_env.txt', 'r') as fid:
            self.num_cells = int(fid.readline().split('=')[1][:-1])
            self.num_bck = int(fid.readline().split('=')[1][:-1])
            self.tmax = float(fid.readline().split('=')[1][:-1])
            self.dt = float(fid.readline().split('=')[1][:-1])
            self.cell_types = [int(i) for i in fid.readline().split(",")[:-1]]

    def _read_info(self):
        with open(f'{self.tag}_env.json', 'r') as fid:
            self.info = json.load(fid)
            self.num_cells = self.info["num_cells"]
            self.num_bck = self.info["num_bck"]
            self.dt = self.info["dt"]
            self.tmax = self.info["tmax"]
            self._gen_types()
    
    def _gen_types(self):
        ratio = self.info["cell_type_ratio"]
        self.cell_types = np.zeros(self.num_cells, dtype=int)
        tp, n_end = 0, int(self.num_cells * ratio[0])
        for n in range(self.num_cells):
            if n > n_end:
                tp += 1
                n_end += int(self.num_cells * ratio[tp])
            self.cell_types[n] = tp

    def read_all_data(self):
        var_names = ["vs", "us", "ics"]
        fnames = ["fv", "fu", "fi"]
        for n in range(3):
            if os.path.isfile(self.tag+"_%s.dat"%(fnames[n])):
                self.read_dat(var_names[n], fnames[n])
        
        if os.path.isfile(self.tag+"_ft_spk.dat"):
            self.read_tspk_dat()
        elif os.path.isfile(self.tag+"_ft_spk.txt"):
            self.read_tspk()

    def read_dat(self, varname, fname):
        x = read_byte_data(f"{self.tag}_%s.dat"%(fname), self.num_cells)      
        exec("self.%s = x"%(varname))
        
    def read_tspk(self):
        self.t_spks = [[] for i in range(N)]
        with open(f"{self.tag}_ft_spk", "r") as fid:
            line = fid.readline() 
        data = line.split(",")[:-1]
        for pair in data:
            nstep, nid = pair.split("-")
            self.t_spks[int(nid)].append(ts[int(nstep)])
        self.num_spks = [len(t) for t in self.t_spks]

    def read_tspk_dat(self):
        with open(self.tag+"_ft_spk.info", "r") as fp:
            self.num_spks = [int(n) for n in fp.readline().split(",")[:-1]]
        
        with open(self.tag+"_ft_spk.dat", "rb") as fp:
            t_spks_flat = np.fromfile(fp, dtype=np.int32)*self.dt
        # align t_spks
        n0 = 0
        self.t_spks = []
        for n in range(self.num_cells):
            self.t_spks.append(t_spks_flat[n0:n0+self.num_spks[n]])
            n0 += self.num_spks[n]

    def read_summary(self):
        fname = self.tag+"_summary"
        if not os.path.isfile(fname+".info"):
            print("There is no summary file")
            return

        with open(fname+".info", "r") as f:
            vals = f.readline().split(",")
        nt = int(vals[0].split("=")[1])
        nf  = int(vals[1].split("=")[1])

        with open(fname+".dat", "rb") as f:
            data = np.fromfile(f, dtype=np.double)
        self.ts_s = data[:nt]
        self.vm = data[nt:2*nt]
        self.rk = data[2*nt:3*nt]
        self.freq = data[3*nt:3*nt+nf]
        self.yf = data[3*nt+nf:3*nt+2*nf]


def read_byte_data(fname, N):
    with open(fname, "rb") as fid:
        data = np.fromfile(fid, dtype=np.dtype(np.double))
    nline = len(data) // N
    # return data.reshape([N, nline], order='C')
    return data.reshape([nline, N], order='C').T


# def read_n


def get_raster_plot(t_spks, cell_types):
    cell_colors = ['#C62828', '#1565C0', '#D68910', '#2ECC71']
    ts = []
    id_spks = []
    cs = []
    
    for i, ctp in enumerate(cell_types):
        ts.extend(t_spks[i])
        for j in range(len(t_spks[i])):
            id_spks.append(i)
            cs.append(cell_colors[ctp])
    
    return ts, id_spks, cs


def get_spike_hist(obj, tbin=5):
    num_exc = int(obj.num_cells * 0.8)
    edges = np.arange(0, obj.ts[-1]+tbin/2, tbin)
    t_spks_vec = np.zeros(len(edges)-1)
    for n in range(num_exc):
        bins, _ = np.histogram(obj.t_spks[n], edges)
        bins[bins > 0] = 1
        t_spks_vec += bins
    t_hist = (edges[1:] + edges[:-1])/2
    return t_spks_vec/num_exc, t_hist


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

    x = np.array(x)
    y = np.array(y)
    c = np.array(c)

    if xlim is not None:
        idx = (x>=xlim[0]) & (x<xlim[1])
    else:
        idx = np.ones(len(x), dtype=bool)
        
    plt.scatter(x[idx], y[idx], s=s, c=c[idx])
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is None:
        ylim = [0, num_cells]
    plt.ylim(ylim)
    plt.xlabel('time (ms)', fontsize=12)
    yt = [0, ylim[1]//2, ylim[1]]
    plt.yticks(yt, labels=['#%d'%(i) for i in yt])
    
    return x, y, c


def draw_single_summary(obj, xlim=None, flim=(5,105), clim=None, vlim=None, title=None, xplim=None, ha='center', markersize=1):
    
    from matplotlib.patches import Rectangle
    import nnpy.pyeeg as pyeeg

    def set_axes(xt=()):
        plt.xticks(xt, fontsize=8)
        plt.yticks(fontsize=8)
    
    spk_hist, t_hist = get_spike_hist(obj, tbin=5)
    psd, fpsd, tpsd = pyeeg.get_stfft(obj.vm, obj.ts_s, 2000, wbin_t=0.5, mbin_t=0.01, f_range=(4, 200))

    ax1 = plt.axes((0.1, 0.65, 0.6, 0.2))
    if xlim is None: xlim = (0, obj.tmax)

    draw_raster_plot(obj.t_spks, cell_types=obj.cell_types, xlim=xlim, s=markersize)
    plt.xlabel("")
    plt.xlim(xlim)
    plt.ylim([-30, obj.num_cells+30])
    set_axes(np.arange(xlim[0], xlim[1]+1, 100))
    plt.ylabel("cell id", fontsize=10)
    ax = plt.twinx()
    plt.plot(t_hist, spk_hist, 'k', lw=1)
    set_axes(np.arange(xlim[0], xlim[1]+1, 100))
    plt.ylim([-0.1, 1.1])
    plt.ylabel(r"$p_{e,fire}$")
    
    if vlim is None: vlim = (min(obj.vm-10), max(obj.vm+10))
    if xplim is None: xplim = (0, obj.tmax)
    rect_args = [[xlim[0], vlim[0]], xlim[1]-xlim[0], vlim[1]-vlim[0]]
    rect_kwargs = {"edgecolor": "k", "facecolor": "k", "fill":True, "alpha": 0.2, "linewidth":1}
    
    xt = np.arange(xplim[0], xplim[1]+1, 2000)
    ax2 = plt.axes((0.1, 0.35, 0.6, 0.22))
    plt.plot(obj.ts_s, np.array(obj.vm), 'k', lw=1)
    plt.xlim(xplim)
    plt.ylim(vlim)
    set_axes()
    plt.ylabel(r"$\langle v \rangle$", fontsize=10)
    ax = plt.twinx()
    plt.plot(obj.ts_s, obj.rk, 'r', lw=1)
    plt.ylabel(r"$r_{k}$")
    plt.ylim([-0.1, 1.1])
    ax.spines["right"].set_color("r")
    ax.yaxis.label.set_color("r")
    ax.tick_params(axis="y", color="r")
    plt.yticks(color="r")
    plt.xticks(xt, labels=[])
    plt.yticks(fontsize=8, color="r")
    ax2.add_patch(Rectangle(*rect_args, **rect_kwargs))
    
    rect_args[0][1] = flim[0]
    rect_args[2] = flim[1]-flim[0]
    
    ax3 = plt.axes((0.1, 0.1, 0.6, 0.25))
    yt = np.arange((flim[0]//10)*10, flim[-1]+10, 20)
    if clim is None:
        clim = [None, None]
    plt.imshow(psd, origin='lower', extent=(tpsd[0]*1e3, tpsd[-1]*1e3, fpsd[0], fpsd[-1]), aspect='auto', cmap='jet', vmax=clim[1], vmin=clim[0])
    plt.xticks(xt, fontsize=8)
    plt.yticks(yt, fontsize=8)
    plt.xlim(xplim)
    plt.ylim(flim)
    plt.xlabel("time (ms)", fontsize=10)
    plt.ylabel("frequency (Hz)", fontsize=10)
    ax3.add_patch(Rectangle(*rect_args, **rect_kwargs))

    plt.axes((0.78, 0.1, 0.15, 0.75))
    plt.yticks([])
    plt.xlabel('frequency (Hz)', fontsize=10)
    plt.xticks(fontsize=8)
    set_axes(yt)
    
    plt.twinx()
    nf = np.argmax(obj.yf)
    yl = obj.yf[nf]
    yl = [-0.05*yl, 1.2*yl]
    
    plt.plot(obj.freq, obj.yf, 'k', lw=1)
    plt.plot(obj.freq[nf], 1.02*obj.yf[nf], 'rv', markersize=5)
    plt.text(obj.freq[nf], 1.06*obj.yf[nf], r"$%d$ Hz"%(obj.freq[nf]), fontsize=8, ha=ha)
    plt.xticks(yt, fontsize=8)
    plt.xlim(flim)
    plt.ylim(yl)
    plt.yticks(fontsize=8)
    tw = ((obj.tmax-5000)/1000, obj.tmax/1000)
    plt.ylabel("fft result from %d~%ds"%(tw), fontsize=10)
    
    if title is None:
        title = obj.tag
    plt.suptitle(title)
    

### Network tools ###
def read_graph(graph_name):
    graph = dict()
    graph["weight"] = []
    graph["adj_list"] = []
    
    with open(graph_name, "r") as fid:
        line= fid.readline()
        while line:
            _, n_pre, n_post, w, _ = line[:-1].split(",")
            n_pre, n_post, w = int(n_pre), int(n_post), float(w)
            if len(graph["adj_list"]) < n_pre+1:
                graph["adj_list"].append([n_post])
                graph["weight"].append([w])
            else:
                graph["adj_list"][n_pre].append(n_post)
                graph["weight"][n_pre].append(w)
            line = fid.readline()
    return graph


def convertGraph2DL(graph, fdl):
    num_cells = len(graph["adj_list"])
    
    with open(fdl, "w") as fid:
        # write node
        fid.write("nodedef>name VARCHAR\n")
        for n in range(num_cells):
            fid.write("%d\n"%(n))
        fid.write("edgedef>node1,node2,weight DOUBLE,directed BOOLEAN\n")
        
        for n in range(num_cells):
            for i in range(len(graph["adj_list"][n])):
                n_post = graph["adj_list"][n][i]
                w = graph["weight"][n][i]
                fid.write("%s,%s,%s,true\n"%(n,n_post,w))
    print("Done\n")


def readDL(fdl):
    graph_viz = {"id":[], "x":[], "y":[], "edge_line":[[],[]], "line_id": []}
    with open(fdl, "r") as fid:
        info = fid.readline()
        line = fid.readline()
        while "edge" not in line:
            vals =  line.split(",")
            graph_viz["id"].append(int(vals[0]))
            graph_viz["x"].append(float(vals[4]))
            graph_viz["y"].append(float(vals[5]))
            line = fid.readline()
        # read edge info
        graph_viz["line_id"] = [[-1, 0] for n in range(len(graph_viz["id"]))]
        
        line = fid.readline()
        while line:
            vals = line.split(",")
            n1, n2 = int(vals[0]), int(vals[1])
            ex = [graph_viz["x"][n1], graph_viz["x"][n2], np.nan]
            ey = [graph_viz["y"][n1], graph_viz["y"][n2], np.nan]
            
            n_line = len(graph_viz["edge_line"][0])
            if graph_viz["line_id"][n1][0] == -1:
                graph_viz["line_id"][n1][0] = n_line
            graph_viz["line_id"][n1][1] = n_line+2
            graph_viz["edge_line"][0].extend(ex)
            graph_viz["edge_line"][1].extend(ey)
            line = fid.readline()
    return graph_viz


class SpikeVisu:
    def __init__(self, obj, tbin=0.1):
        self.tbin = tbin # ms
        self.num_cells = obj.num_cells
        self.id_t = np.zeros(obj.num_cells, dtype=int)
        self.t = 0
        self.t_spks = obj.t_spks
        self.num_spks = [len(t) for t in obj.t_spks]
        self._set_ctype_cs()
        
    def _set_ctype_cs(self):
        self.colors = [[0.8, 0.2, 0.2], [0.2, 0.2, 0.8]]
        self.ctypes = []
        for n in range(self.num_cells):
            if n < 0.8*self.num_cells:
                self.ctypes.append(0)
            else:
                self.ctypes.append(1)
    
    def reset(self):
        for n in range(self.num_cells):
            self.id_t[n] = 0
    
    def set_time(self, t):
        self.t = t
        for n in range(self.num_cells):
            self.id_t[n] = np.where(np.array(obj.t_spks[n]) < self.t)[0][-1]
            
    def __next__(self):
        self.t += self.tbin
        cs = np.ones([self.num_cells, 3])
        act_node = []
        for n in range(self.num_cells):
            flag = False
            while (self.id_t[n]+1 < self.num_spks[n]) and (self.t_spks[n][self.id_t[n]+1] < self.t):
                self.id_t[n] += 1
                flag = True
                
            if flag:
                cs[n,:] = self.colors[self.ctypes[n]]
                act_node.append(n)
                
        return cs, act_node