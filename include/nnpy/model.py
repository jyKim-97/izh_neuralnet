import numpy as np
from abc import *
from time import time
from functools import wraps

def wrap_init_net(func):
    @wraps(func)
    # initialize network parameters
    def init_net(self, cell_params, cell_types, weight_mat):
        # save parameters
        self.cell_params = cell_params
        self.cell_types = cell_types
        self.weight_mat  = weight_mat
        # get # of the cells
        self.N = len(cell_types)
        self.num_types = len(cell_params)
        self.v = np.zeros(self.N)
        self.u = np.zeros(self.N)
        # set params
        func(self, cell_params, cell_types, weight_mat)
    return init_net


def update_x_rk4(func, dt, x, *args):
    k1 = func(x, *args) * dt
    k2 = func(x+k1/2, *args) * dt
    k3 = func(x+k2/2, *args) * dt
    k4 = func(x+k3, *args) * dt
    return x + (k1+2*k2+2*k3+k4)/6 # return x_{n+1}


# test current
def const_current(N, n_on, amp, cell_types=None):
    if cell_types is None:
        n0 = N
    else:
        n0 = 0
        while cell_types[n0] == 0: n0 += 1
    
    n = 0
    while True:
        n += 1
        if n > n_on:
            ic = np.zeros(N)
            ic[:n0] = amp
            yield ic
        else:
            yield 0


class SingleCell():
    def __init__(self, a, b, c, d, dt=0.01):
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.ic_t = None
        self.ic_amp = 1
        self.dt = dt

    def set_square_current(self, input_time, amp):
        self.ic_t = input_time
        self.ic_amp = amp

    def run(self, tmax, ic=None):
        nitr = int(tmax/self.dt)+1
        self.vs = np.zeros(nitr)
        self.us = np.zeros(nitr)
        self.ts = np.arange(nitr) * self.dt

    def update(self, n, ic):
        vnew = update_x_rk4(self._fv, self.dt, self.v, self.u, ic+self.isyn)
        unew = update_x_rk4(self._fu, self.dt, self.u, self.v, self.a, self.b)
        

    def _fv(self, v, u, I):
        return 0.04*v**2 + 5*v + 140 - u + I

    def _fu(self, u, v, a, b):
        return a * (b*v - u)

    def find_spk(self, vs):
        ids = vs >= 30
        return ids



# Network models
class NeuralNet(metaclass=ABCMeta):
    def __init__(self, dt=0.01):
        self.dt = dt

    @abstractmethod
    def compile_nn(self, cell_params, adj_mat, weight_mat):
        pass

    @abstractmethod
    def update(self):
        pass
    
    @abstractmethod
    def save_history(self, n):
        # n: current simulation step number
        pass

    @abstractmethod
    def init_run(self, nt):
        # nt: the # of the times
        pass

    @abstractmethod
    def get_syn_current(self):
        # get synaptic current from the presynaptic neuron
        pass

    def _check_current_type(self, Is):
        # Is (1, nt) or (self.N, nt)
        self.is_gen = False
        if Is == 0:
            return
        elif Is.__class__.__name__ == 'generator':
            self.is_gen = True
        elif len(Is.shape) == 1:
            Is.resize(1, len(Is))

    def _get_ext_current(self, Is, n):
        if self.is_gen:
            return next(Is)
        else:
            return Is[:, n]
    
    # run simulation
    def run(self, Is=0, tmax=1000):

        self._check_current_type(Is)

        print("Start run newtork")
        tic = time()
        nt = int(tmax/self.dt)
        self.init_run(nt)
        for n in range(nt):
            iext = self._get_ext_current(Is, n)
            self.update(iext)
            self.save_history(n+1)            
        
        print("Execution time = %.2fs"%(time()-tic))

    def find_spk(self, vs):
        ids = vs >= 30
        return ids


class IzhSimpleNet(NeuralNet):
    @wrap_init_net
    def compile_nn(self, cell_params, cell_types, weight_mat):
        # weight_mat [post_neuron, pre_neuron]
        # set params - a, b, c, d
        self._alloc_cell_params(cell_params, cell_types)
        
    def _alloc_cell_params(self, cell_params, cell_types):
        for str_var in ['a', 'b', 'c', 'd']:
            self.__dict__[str_var] = np.array(
                [cell_params[ctp][str_var] for ctp in cell_types]
                )

    def _fv(self, v, u, I):
        return 0.04*v**2 + 5*v + 140 - u + I

    def _fu(self, u, v, a, b):
        return a * (b*v - u)

    def _init_cell(self, nt):
        # intializing state variables
        self.v = np.ones(self.N) * (-50)
        self.u = np.zeros(self.N)
        self.is_fire = np.zeros(self.N, dtype=bool)
        self.ts = np.arange(nt+1)*self.dt # time array (ms)
        
        # initializing save matrices
        self.vs = np.zeros([self.N, nt+1])
        self.us = np.zeros([self.N, nt+1])
        self.isyns = np.zeros([self.N, nt+1])
        self.vs[:,0], self.us[:,0] = self.v, self.u
        self.spk_mat = np.zeros([self.N, nt+1])

    def _update_cell(self, iext):
        self.get_syn_current()
        vnew = update_x_rk4(self._fv, self.dt, self.v, self.u, iext+self.isyn)
        unew = update_x_rk4(self._fu, self.dt, self.u, self.v, self.a, self.b)
        self.is_fire = self.find_spk(vnew)
        vnew[self.is_fire] = self.c[self.is_fire]
        unew[self.is_fire] += self.d[self.is_fire]
        self.v, self.u = vnew, unew

    def _save_cell_history(self, n):
        self.vs[:, n] = self.v
        self.us[:, n] = self.u
        self.isyns[:, n] = self.isyn
        self.spk_mat[:, n] = self.is_fire

    def get_syn_current(self):
        if not any(self.is_fire):
            return np.zeros(self.N)
        self.isyn = np.sum(self.weight_mat[:, self.is_fire], axis=1)

    def save_history(self, n):
        self._save_cell_history(n)

    def update(self, iext):
        self._update_cell(iext)

    def init_run(self, nt):
        self._init_cell(nt)
        self.save_history(0)


class IzhExpNet(IzhSimpleNet):
    @wrap_init_net
    def compile_nn(self, cell_params, cell_types, weight_mat):
        # weight_mat [post_neuron, pre_neuron]
        # set params - a, b, c, d
        self._alloc_cell_params(cell_params, cell_types)
        # set synaptic parameters - tau, ev (equilibrium potential)
        self._alloc_syn_params(cell_params, cell_types)
    
    def _alloc_syn_params(self, cell_params, cell_types):
        # allocate synaptir parameters
        for str_var in ['tau', 'ev']:
            self.__dict__[str_var] = np.array(
                [cell_params[ctp][str_var] for ctp in cell_types]
                )
        self.ev = np.tile(self.ev, [self.N, 1])
        self.D = 0.05

    def _fsyn(self, r, tau, is_fired, D):
        return (-r + D * is_fired)/tau

    def _update_syn(self):
        self.r = update_x_rk4(self._fsyn, self.dt, self.r, self.tau, self.is_fire, self.D)

    def _init_syn(self, nt):
        self.r = np.zeros(self.N)
        self.rs = np.zeros([self.N, nt+1])

    def _save_syn_history(self, n):
        self.rs[:,n] = self.r

    def init_run(self, nt):
        self._init_cell(nt)
        self._init_syn(nt)

    def save_history(self, n):
        self._save_cell_history(n)
        self._save_syn_history(n)        

    def get_syn_current(self):
        vpost = np.tile(self.v, [self.N, 1]).T
        isyn_tmp = self.weight_mat * (vpost - self.ev)
        self.isyn = -np.dot(isyn_tmp, self.r)

    def update(self, iext):
        self._update_cell(iext)
        self._update_syn()

    
class IzhExpNet_poisson_input(IzhExpNet):
    def set_Poisson_input(self, Next, gext, tauext, pext, pf, n_on=0, targets=None):
        # gext (list, 1 x num_types): maximal conductance
        # pext (list, 1 x num_types): connection prob
        # pf (1, ,firing prob)
        self.Next = Next
        self.tauext = tauext # constant
        self.pext_spk = pf
        if targets is None:
            self.id_target = np.ones(self.N, dtype=bool)
        else:
            self.id_target = np.zeros(self.N, dtype=bool)
            self.id_target[targets] = True
        
        self._gen_poisson_weight(gext, pext)
        self.rext = np.zeros(self.Next)
        self.n_on = n_on
    
    def _gen_poisson_weight(self, gext, pext):
        self.weight_pos = np.zeros([self.N, self.Next])
        for i in range(self.N):
            ctp = self.cell_types[i]
            for j in range(self.Next):
                if not self.id_target[i]:
                    continue

                if np.random.rand() < pext[ctp]:
                    self.weight_pos[i, j] = gext[ctp]

    def _update_ext_syn(self):
        is_fire = np.random.rand(self.Next) < self.pext_spk
        self.rext = update_x_rk4(self._fsyn, self.dt, self.rext, self.tauext, is_fire, self.D)

    def _get_ext_current(self, Iext, n):
        if n > self.n_on:
            self._update_ext_syn()
            gr = np.dot(self.weight_pos, self.rext)
            return Iext - self.v * gr
        else:
            return Iext
