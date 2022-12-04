from lzma import CHECK_SHA256
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
    # k1 = func(x, *args) * dt
    # k2 = func(x+k1/3, *args) * dt
    # k3 = func(x-k1/3+k2, *args) * dt
    # k4 = func(x+k1-k2+k3, *args) * dt
    # return x + (k1+3*k2+3*k3+k4)/8 # return x_{n+1}


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


class SingleCell:
    def __init__(self, a, b, c, d, ic_time=None, ic_amp=None, dt=0.01):
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        self.dt = dt

    def init_cell(self, tmax):
        nitr = int(tmax/self.dt)+1

        self.vs = np.zeros(nitr)
        self.us = np.zeros(nitr)
        self.ics = np.zeros(nitr-1)
        self.nid_spk = []
        self.ts = np.arange(nitr) * self.dt
        self.vs[0] = self.c


    def run(self, tmax, ic_custom=None, ic_sq_t=None, ic_sq_amp=None):
        nitr = int(tmax/self.dt)+1

        ic_sq_n = ic_sq_t
        if ic_sq_t is not None:
            ic_sq_n[0] /= self.dt
            ic_sq_n[1] /= self.dt
            ic_square = True
        else:
            ic_square = False

        self.init_cell(tmax)
        for n in range(nitr-1):
            ic = 0
            if ic_custom is not None:
                ic += ic_custom[n]
            if ic_square and (ic_sq_n[0] <= n < ic_sq_n[1]):
                ic += ic_sq_amp

            self.update(n, ic=ic)
            self.ics[n] = ic


    def update(self, n, ic=0):
        self.vs[n+1] = update_x_rk4(self._fv, self.dt, self.vs[n], self.us[n], ic)
        self.us[n+1] = update_x_rk4(self._fu, self.dt, self.us[n], self.vs[n], self.a, self.b)
        self.check_spk(n+1)

    def _fv(self, v, u, I):
        return 0.04*v**2 + 5*v + 140 - u + I

    def _fu(self, u, v, a, b):
        return a * (b*v - u)

    def check_spk(self, n):
        is_spk = self.vs[n] >= 30
        if is_spk:
            self.vs[n] = self.c
            self.us[n] += self.d
            self.nid_spk.append(n)


def draw_single_result(obj, idt_show=None):

    import matplotlib.pyplot as plt

    if idt_show is None:
        idt_show = np.ones(len(obj.ts), dtype=bool)

    def set_axes():
        plt.xticks(fontsize=9)
        plt.yticks(fontsize=9)
    
    def draw_current(obj):
        plt.twinx()
        plt.plot(obj.ts[1:], obj.ics, 'r', lw=2)
        yl = [min(obj.ics)-1, max(obj.ics)*10]
        plt.ylim(yl)
        plt.ylabel("input current", fontsize=11)
        set_axes()

    plt.figure(dpi=150, figsize=(8,8))

    plt.subplot(211)
    plt.plot(obj.vs[idt_show], obj.us[idt_show], 'k', lw=1)
    plt.xlabel("v", fontsize=12)
    plt.ylabel("u", fontsize=12)
    plt.title("a=%.3f, b=%.3f, c=%.3f, d=%.3f"%(obj.a, obj.b, obj.c, obj.d), fontsize=14)

    plt.subplot(223)
    plt.plot(obj.ts[idt_show], obj.vs[idt_show], 'k', lw=1)
    plt.xlabel("time (ms)", fontsize=12)
    plt.ylabel("voltage (v) (mV)", fontsize=12)
    set_axes()
    draw_current(obj)
    plt.ylabel("")

    plt.subplot(224)
    plt.plot(obj.ts[idt_show], obj.us[idt_show], 'k', lw=1)
    plt.xlabel("time (ms)", fontsize=12)
    plt.ylabel("u", fontsize=12)
    set_axes()
    draw_current(obj)
    
    plt.tight_layout()


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


if __name__=="__main__":

    import matplotlib.pyplot as plt

    # a = 0.1
    # b = 0.2
    # c = -65
    # d = 2
    # ic_sq_amp = 3.87
    a = 0.02
    b = 0.25
    c = -65
    d = 2
    ic_sq_amp = 0.6

    dt = 0.05
    tmax = 1200
    ic_sq_t = [500, 1000]

    ic_custom = np.zeros(int(tmax/dt))
    for n in range(int(tmax/dt)):
        if (520 <= n*dt < 550):
            ic_custom[n] = 1e-1
    
    obj = SingleCell(a, b, c, d, dt=dt)
    obj.run(tmax, ic_sq_t=ic_sq_t, ic_custom=ic_custom, ic_sq_amp=ic_sq_amp)

    idt_show = 10 <= obj.ts
    draw_single_result(obj, idt_show)
    plt.savefig("./test_single_cell.png")
    plt.close()
    
    # plt.figure(dpi=150, figsize=(8,8))

    # plt.subplot(211)
    # plt.plot(obj.vs[idt_show], obj.us[idt_show], 'k', lw=1)
    # plt.xlabel("v", fontsize=12)
    # plt.ylabel("u", fontsize=12)
    # plt.title("a=%.3f, b=%.3f, c=%.3f, d=%.3f"%(obj.a, obj.b, obj.c, obj.d), fontsize=14)

    # plt.subplot(223)
    # plt.plot(obj.ts[idt_show], obj.vs[idt_show], 'k', lw=1)
    # plt.xlabel("time (ms)", fontsize=12)
    # plt.ylabel("voltage (v) (mV)", fontsize=12)
    # set_axes()
    # draw_current(obj)
    # plt.ylabel("")

    # plt.subplot(224)
    # plt.plot(obj.ts[idt_show], obj.us[idt_show], 'k', lw=1)
    # plt.xlabel("time (ms)", fontsize=12)
    # plt.ylabel("u", fontsize=12)
    # set_axes()
    # draw_current(obj)
    
    # plt.tight_layout()
    # plt.savefig("./test_single_cell.png")
    # plt.close()