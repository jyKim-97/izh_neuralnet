from struct import unpack
import numpy as np
import matplotlib.pyplot as plt


class Reader:
    def __init__(self, tag):
        self.tag = tag
        self._read_info()
        self.vs = None
        self.us = None
        self.rs = None
        self.ics = None
        self.t_spks = None

    def _read_info(self):
        with open(f'{self.tag}_env.txt', 'r') as fid:
            line = fid.readline().split(',')
            self.N, self.num_syns, self.tmax, self.dt = int(line[0]), int(line[1]), float(line[2]), float(line[3])
            self.cell_types = [int(i) for i in fid.readline().split(',')]
            self.id_presyn  = [int(i) for i in fid.readline().split(',')]
        self.ts = np.arange(0, self.tmax+self.dt/2, self.dt)

    def read_uv(self):
        self.vs = self._read_byte_data(f'{self.tag}_v.dat')
        self.us = self._read_byte_data(f'{self.tag}_u.dat')

    def read_ic(self):
        self.ics = self._read_byte_data(f'{self.tag}_i.dat')

    def read_r(self):
        self.rs = self._read_byte_data(f'{self.tag}_r.dat')
    
    def read_tspk(self):
        self.t_spks = [[] for i in range(self.N)]
        with open(f'{self.tag}_t.txt', "r") as fid:
            for i in range(len(self.ts)):
                ids = [int(i) for i in fid.readline().split(',')[:-1]]
                for j in ids:
                    self.t_spks[j].append(self.ts[i])
        for i in range(len(self.t_spks)):
            self.t_spks[i] = np.array(self.t_spks[i])

    def _read_byte_data(self, fname):
        with open(fname, "rb") as fid:
            # line = fid.read()
            data = np.fromfile(fid, dtype=np.dtype(np.double))
        nline = len(data) // self.N            
        # nline = len(line) // self.N // 8  # double
        # data = np.array(unpack(f"{self.N*nline}d", line))
        return data.reshape([nline, self.N])


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
