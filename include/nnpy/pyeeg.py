import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, sosfilt, hilbert, sosfiltfilt, resample

try:
    from spectral_connectivity import Multitaper
    from spectral_connectivity import Connectivity
    # sc_load_fail = False
except:
    # sc_load_fail = True
    print("Spectral connectivity toolbox is not completely loaded, check conda env")
    

def bandpass_filt(x, fs, frange, order=5):
    # bandpass signal
    sos = butter(order, frange, fs=fs, btype='bandpass', output='sos', analog=False)
    return sosfiltfilt(sos, x)


class PAC:
    # Phase Amplitude coupling
    def __init__(self, xl, xh, fs, order=5):
        self.xl = xl
        self.xh = xh
        self.fs = fs
        self.order = order
        
    def run(self, flows, fhighs, idt=None, nbin=50):
        if idt is None:
            idt = np.ones_like(self.xl, dtype=bool)
        
        self.nbin = nbin
        self.flows, self.fhighs = flows, fhighs
        nl, nh = len(flows), len(fhighs)
        self.xphs = np.zeros([nl, sum(idt)])
        self.yamp = np.zeros([nh, sum(idt)])

        # bandpass filtering
        for i in range(nl):
            tmp, _ = self._get_filtered(self.xl, flows[i])
            self.xphs[i,:] = tmp[idt]
        for i in range(nh):
            _, tmp = self._get_filtered(self.xh, fhighs[i])
            self.yamp[i,:] = tmp[idt]
        
        self.amp_align = np.zeros([nl, nh, self.nbin])
        for i in range(nl):
            self.amp_align[i,:,:] = self._align_amp(self.xphs[i], self.yamp)
            
        # normalize
        self.prob = self.amp_align / np.nansum(self.amp_align, 2).reshape([nl, nh, 1])
        
        # get MI
        self.mi_max = np.log(self.nbin)
        self.prob[self.prob==0] = 1
        self.mi = -np.nansum(self.prob * np.log(self.prob), 2)
        self.mi = (self.mi_max - self.mi) / self.mi_max
        self.mi = self.mi.T
        
    def show(self, **kwargs):
        plt.figure(dpi=150, figsize=(4,4))
        plt.imshow(self.mi, cmap='jet', extent=(self.flows[0][0], self.flows[-1][1], self.fhighs[0][0], self.fhighs[-1][1]), aspect='auto', **kwargs)
        plt.xlabel("low frequency (phase, Hz)", fontsize=11)
        plt.ylabel("high frequency (amplitude, Hz)", fontsize=11)
        plt.title("Modulation Index", fontsize=13)
        plt.colorbar()
        plt.show()
    
    def _align_amp(self, xl, yhs):
        edges = np.linspace(-np.pi, np.pi, self.nbin+1)
        ids = np.digitize(xl, edges)-1
        amp_align = np.zeros([len(yhs), self.nbin])
        
        for i in range(self.nbin):
            tmp = yhs[:, ids==i]
            if np.shape(tmp)[1] == 0:
                continue
            amp_align[:,i] = np.nanmean(tmp, 1)
        
        return amp_align
    
    def _get_filtered(self, x, frange):
        filt = bandpass_filt(x, self.fs, frange, order=self.order)
        h = hilbert(filt)
        amp = np.abs(h)
        phs = np.angle(h)
        return phs, amp
    

def get_fft_with_t(x, t, dt=None, t_range=None):
    # use Fourier transform without 
    if t_range is None:
        t_range = [t[0], t[-1]]
    
    idt = (t >= t_range[0]) & (t < t_range[1])
    N = sum(idt)
    
    freq = np.linspace(0, 1/(2*dt), N//2)
    yf = np.fft.fft(x[idt])
    yf = 2/N * np.abs(yf[:N//2])
    
    return freq, yf


def downsample_signal(x, t, fs_org, fs_new, window=None):
    x_new, t_new = resample(x, int(len(x)*fs_new/fs_org), t=t, window=window, domain='time')
    return x_new, t_new


def measure_granger_causality(data, srate, method="pairwise_spectral_granger"):
    # measure granger causality
    # data (time series, samples, data class)
    m = Multitaper(data,
                   sampling_frequency=srate,
                   time_halfbandwidth_product=1,
                   start_time=0)
    
    c = Connectivity(fourier_coefficients=m.fft(),
                     frequencies=m.frequencies,
                     time=m.time)
    
    gc = c.pairwise_spectral_granger_prediction()
    return gc, c.frequencies
    
    
if __name__ == "__main__":
    # test bandpass filtering
    import matplotlib.pyplot as plt
    
    fs = 1000

    t_test = np.arange(0, 100, 1/fs)
    y_test = 5*np.sin(2*np.pi*t_test) + 4*np.cos(10*2*np.pi*t_test)*(np.cos(2*np.pi*t_test)+0)

    frange = [8, 12]
    order = 10

    y_filt = bandpass_filt(y_test, fs, frange, order)
    y_hilt = hilbert(y_filt, len(y_filt))
    
    plt.figure(dpi=150)
    plt.plot(t_test, y_test, 'k')
    plt.plot(t_test, y_filt, 'r')
    plt.plot(t_test, np.abs(y_hilt), 'b')
    plt.xlim([20, 22])
    plt.ylim([-8, 8])
    plt.yticks(np.arange(-8, 8.1, 1))
    plt.grid(True)
    plt.show()