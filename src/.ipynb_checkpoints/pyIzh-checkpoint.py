import numpy as np
from copy import copy
import json

# Contents
# 1. FILE IO from/to C
# 2. Izhikevich single neuron model to analyze dynamics

# define default params
default_cell_params = {
    'fs': [0.1, 0.2, -65, 2],
    'rs': [0.02, 0.2, -65, 8],
    'ib': [0.02, 0.2, -55, 4],
    'ch': [0.02, 0.2, -50, 2]}

class InfoIzh:
    """
    Description:
        - NumCells: the # of the cells, N
        - NumTypes: the # of the cell types, n
        - CntProbs: Connection probability betweeen cell types, n x n
        - gSyns: initial connection strength between cells, n x n
        - tauSyns: time constant of the synapses, n
        - veqSyns: equilibrium potential of the synapses, n
        - NumIn: # of the input source
        - fin: input frequency
    """

    def __init__(self):
        # set default keys
        self.default_keys = ['NumCells', 'NumTypes', 'CellParams', # single cell info
            'CntProbs', 'gSyns', 'tauSyns', 'veqSyns', # network info
            'NumIn', 'fIn', 'gIn', 'tauIn', 'pIn'] # network input
        self.info = dict()

    def add_params(self, **kwargs):
        for key in kwargs.keys():
            # check is the key correpsonds to default keys
            if not key in self.default_keys:
                raise ValueError(
                    f'"{key}" is not included in default key values {self.default_keys}'
                )
            self.info[key] = kwargs[key]
        
    def export_info(self, fname):
        # check the variable
        self._inspect_info()
        with open(fname, 'w', encoding='utf-8') as f:
            json.dump(self.info, f, ensure_ascii=False, indent='\t')

    def find_missed_key(self):
        missed_keys = copy(self.default_keys)
        for key in self.info.keys():
            missed_keys.remove(key)
        return missed_keys

    def _inspect_info(self):
        # check all the default_keys are in info
        missed_keys = self.find_missed_key()
        if len(missed_keys) > 0:
            raise ValueError(f'Some keys are missed, {missed_keys}')

        # check type # of params are corresponds to "NumTypes"
        n = self.info['NumTypes']
        for key in ['CellParams', 'gIn', 'pIn', 'tauSyns', 'veqSyns']:
            if len(self.info[key]) != n:
                raise ValueError(
                    f'The number of {key} is not matched to the # of cell types'
                )

        for key in ['CntProbs', 'gSyns']:
            if len(self.info[key]) != n*n:
                raise ValueError(
                    f'The number of {key} is not matched to (the # of cell types)^2'
                )

    def print_info(self):
        pass

class IzhNeuron:
    def __init__(self, a, b, c, d, dt=0.01):
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.dt = dt
    
    def run(self, tmax, ic, t_on=0, t_off=None):
        nitr = int(tmax/self.dt)
        if t_off is None: t_off = tmax
        self.ts = np.arange(0, tmax+self.dt/2, self.dt)
        self.vs = np.zeros(nitr+1)
        self.us = np.zeros(nitr+1)
        # set init params
        self.vs[0] = self.c
        self.us[0] = self.b*self.c
        
        for n in range(nitr):
            if (n*self.dt>=t_on) and (n*self.dt<=t_off):
                ic_input = ic
            else:
                ic_input = 0
            self.vs[n+1], self.us[n+1] = self._update(self.vs[n], self.us[n], ic_input)
        print("Calculate done")
    
    def _update(self, v, u, ic):
        # update v
        dv = self._rk4(self._fv, v, u, ic)
        du = self._rk4(self._fu, u, v)
        v += dv
        u += du
        if (v >= 30):
            v = self.c
            u += self.d
        return v, u
    
    def _fv(self, v, u, ic):
        return 0.04*v**2 + 5*v + 140 - u + ic
    
    def _fu(self, u, v):
        return self.a*(self.b*v - u)
    
    def _rk4(self, f, x, *args):
        dx1 = f(x, *args)*self.dt
        dx2 = f(x+dx1/2, *args)*self.dt
        dx3 = f(x+dx2/2, *args)*self.dt
        dx4 = f(x+dx3, *args)*self.dt
        return (dx1 + 2*dx2 + 2*dx3 + dx4)/6

if __name__ == "__main__":
    # test single neuron model 
    import matplotlib.pyplot as plt

    dt = 0.005
    rs_izh = IzhNeuron(0.02, 0.2, -65, 8, dt)
    ib_izh = IzhNeuron(0.02, 0.2, -55, 4, dt)
    ch_izh = IzhNeuron(0.02, 0.2, -50, 2, dt)
    fs_izh = IzhNeuron(0.1,  0.2, -65, 2, dt)

    tmax = 500
    ic = 10
    t_on = 20
    t_off = 60

    obj_izh = [rs_izh, ib_izh, ch_izh, fs_izh]
    for obj in obj_izh:
        obj.run(tmax, ic, t_on=t_on, t_off=t_off)

    xl = [0, 100]
    obj_name = ["RS", "IB", "CH", "FS"]
    colors = np.array([[217,45,45], [34,49,223]])/255

    plt.figure(dpi=300, figsize=(6,6))
    for i in range(len(obj_izh)):
        ax1 = plt.subplot(4,1,i+1)
        plt.plot(obj_izh[i].ts, obj_izh[i].vs, lw=1.5, c=colors[0])
        plt.grid(True, lw=0.5)
        plt.ylabel('v (mV)', fontsize=12)
        ax2 = ax1.twinx()
        plt.plot(obj_izh[i].ts, obj_izh[i].us, lw=1, c=colors[1])
        plt.xlim(xl)
        plt.ylabel('u', fontsize=12)
        plt.title(obj_name[i], fontsize=13, fontweight='bold')

    ax1.set_xlabel('time (ms)', fontsize=10)

    plt.tight_layout()
    plt.show()