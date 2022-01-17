from struct import unpack
import numpy as np
import matplotlib.pyplot as plt


class IzhReader:
    def __init__(self, tag):
        self.tag = tag
        self.vs = None
        self.us = None
        self.rs = None
        self.ics = None
        self.t_spks = None
        self._read_info()
        
    def _read_info(self):
        with open(f'{self.tag}_env.txt', 'r') as fid:
            self.num_cells = int(fid.readline().split('=')[1][:-1])
            self.num_bck = int(fid.readline().split('=')[1][:-1])
            self.tmax = float(fid.readline().split('=')[1][:-1])
            self.dt = float(fid.readline().split('=')[1][:-1])
            self.cell_types = [int(i) for i in fid.readline().split(",")[:-1]]
            self.file_id = fid.readline()
        self.ts = np.arange(0, self.tmax, self.dt) + self.dt
        self.read_file()
        
    def read_file(self):
        var_names = ["vs", "us", "ics"]
        fnames = ["fv", "fu", "fi"]
        for n in range(3):
            if int(self.file_id[n]) == 1:
                self.read_dat(var_names[n], fnames[n])
        if int(self.file_id[3]) == 1:
            self.read_tspk()
        
    def read_dat(self, varname, fname):
        x = read_byte_data(f"{self.tag}_%s.dat"%(fname), self.num_cells)      
        exec("self.%s = x"%(varname))
        
    def read_tspk(self):
        self.t_spks = read_tspk(f"{self.tag}_ft_spk.txt", self.num_cells, self.ts)


def read_tspk(fname, N, ts):
    t_spks = [[] for i in range(N)]
    with open(fname, "r") as fid:
        line = fid.readline()
    data = line.split(",")[:-1]
    for pair in data:
        nstep, nid = pair.split("-")
        t_spks[int(nid)].append(ts[int(nstep)])
    return t_spks


def read_byte_data(fname, N):
    with open(fname, "rb") as fid:
        data = np.fromfile(fid, dtype=np.dtype(np.double))
    nline = len(data) // N
    return data.reshape([nline, N], order='C').T


def get_phase(t_spks, N, ts):
    phase = np.zeros([N, len(ts)])
    phase[:] = np.nan 

    dt = ts[1]-ts[0]

    for n in range(N):
        for t in ts:
            i = int(t / dt)


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

