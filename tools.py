import numpy as np

class Tools():
    """docstring for Utils"""
    def __init__(self):
        super(Tools, self).__init__()

        self.broad_bib = {
        'Lorentzian':self.lorentz,
        'Gaussian':self.gaussian
        }
        self.procedures_bib = {
        'procedure_01':self.proc_01,
        'procedure_02':self.proc_02,
        'procedure_03':self.proc_03,
        'procedure_04':self.proc_04,
        'procedure_05':self.proc_05,
        }
        
    def read_file(self):
        f = open(self.file, 'r')
        lines = f.readlines()
        f.close()
        return lines

    def get_lines_interval(self, interval):
        f = open(self.file, 'r')
        lines = f.readlines()
        f.close()
        return lines[interval[0]:interval[1]]

    def get_lines(self, string, index=False, maxlines=None):
        '''
        Read lines and retur an array with 
        the lines and the index
        '''
        lines = self.read_file()
        if maxlines is None:
            maxlines = len(lines)
        storeLines = []
        storeIndex = []

        if len(lines) < maxlines:
            maxlines = len(lines)
            
        for n in range(maxlines):
            if string.lower() in lines[n].lower():
                storeLines.append(lines[n])
                storeIndex.append(n)

        if index:
            return storeLines, storeIndex
        else:
            return storeLines

    def which_calculation(self, maxlines=150):
        for key in self.calc_type.keys():
            if self.get_lines(key, maxlines=maxlines) != []:
                self.calc_type[key] = True
        self.calc = [ key for key in self.calc_type.keys() if self.calc_type[key]]
   
    def get_all_data(self):
        for calc_ in self.calc_type:
            if self.calc_type[calc_]:
                self.calc_bib[calc_]()

    def set_emax(self, emax):
        if len(emax) == 1:
            for n in range(self.nstates):
                if self.energies[n] > emax[0]:
                    self.set_groups('1-{}'.format(n))
                    return
        elif len(emax) == 2:
            for n in range(self.nstates):
                if self.energies[n] > self.energies[emax[0]-1]*emax[1]:
                    self.set_groups('1-{}'.format(n))
                    return

    def lorentz(self, x, x0, y0, g):
        g /= 2
        g2 = g ** 2
        dx2 = (x - x0) ** 2
        return (y0 / np.pi) * g / (dx2 + g2)

    def gaussian(self, x, x0, y0, g):
        g /= 2
        g2 = g ** 2
        dx2 = (x - x0) ** 2
        log2 = np.log(2)
        c1 = y0 * np.sqrt(log2 / np.pi) / g
        return c1 * np.exp(-log2 * dx2 / g2)

    def convolution_averaged(self):
        self.spectra_y = np.zeros((self.ngroups, self.steps))
        for ngroups in range(self.ngroups):
            for states in self.groups[ngroups]:
                x0max = np.average(self.energies[states], weights=self.intensities[states])
                y0max = np.average(self.intensities[states])
                self.spectra_y[ngroups] += self.broad_func(self.spectra_x, x0max, y0max, self.fwhm[ngroups])
        self.y_max = np.max(self.spectra_y)

    def convolution_configurations(self):
        self.spectra_y = np.zeros((self.ngroups, self.steps))
        for ngroups in range(self.ngroups):
            for states in self.groups[ngroups]:
                nmax = len(self.energies[states])
                for n in range(nmax):
                    x0 = self.energies[states][n]
                    y0 = self.intensities[states][n]
                    self.spectra_y[ngroups] += self.broad_func(self.spectra_x, x0, y0, self.fwhm[ngroups])
            self.spectra_y[ngroups] /= nmax
        self.y_max = np.max(self.spectra_y)

    def split_in_groups(self):
        new_states = []
        for ngroups in range(self.ngroups):
            for states in self.groups[ngroups]:
                new_states += [str(states)]
        new_group = ';'.join(new_states)
        self.set_groups(new_group)

    def proc_01(self):
        self.split_in_groups()
        if len(self.fwhm) < (self.ngroups):
            nmax = self.ngroups - len(self.fwhm)
            for n in range(nmax):
                self.fwhm += [0.4]
        self.convolution_averaged()

    def proc_02(self):
        self.split_in_groups()
        self.def_fhwm_std()
        self.convolution_averaged()

    def proc_03(self):
        self.convolution_configurations()

    def proc_04(self):
        self.split_in_groups()
        self.def_fhwm_std()
        self.convolution_configurations()

    def proc_05(self):
        from scipy.optimize import brentq
        self.spectra_y = np.zeros((self.ngroups, self.steps))
        self.fwhm_optimal = []
        for ngroups in range(self.ngroups):
            fwhm = self.fwhm[ngroups]
            self.fwhm_optimal += [brentq(self.min_proc_05, 0.001, 3.0, args=[fwhm, ngroups])]
        self.y_max = np.max(self.spectra_y)

    def min_proc_05(self, x, args=None):
        fwhm = args[0]
        ngroups = args[1]
        self.spectra_y[ngroups] = np.zeros(self.steps)
        for states in self.groups[ngroups]:
            nmax = len(self.energies[states])
            for n in range(nmax):
                x0 = self.energies[states][n]
                y0 = self.intensities[states][n]
                self.spectra_y[ngroups] += self.broad_func(self.spectra_x, x0, y0, x)
        self.spectra_y[ngroups] /= nmax
        f = self.measure_fwhm(self.spectra_x, self.spectra_y[ngroups])
        return fwhm - f
        
    def measure_fwhm(self, x, y):
        ymax = np.max(y)
        yr = ymax / 2
        d = np.sign(yr - y[0:-1]) - np.sign(yr - y[1:])
        left_idx = np.argmax(d)
        right_idx = np.argmin(d)
        f = x[right_idx] - x[left_idx]
        return f

    def def_fhwm_std(self):
        std = []
        for ngroups in range(self.ngroups):
            for states in self.groups[ngroups]:
                t = self.weighted_std(self.energies[states], self.intensities[states]) * 2 * np.sqrt( 2 * np.log(2))
                std += ['{0:4.2f}'.format(t)]
        new_fwhm = ';'.join(std)
        self.set_fwhm(new_fwhm)

    def weighted_std(self, values, weights):
        # For simplicity, assume len(values) == len(weights)
        # assume all weights > 0
        sum_of_weights = np.sum(weights)
        values = np.array(values)
        weights = np.array(weights)
        weighted_average = np.sum(values * weights) / sum_of_weights
        n = len(weights)
        numerator = np.sum(n * weights * (values - weighted_average) ** 2.0)
        denominator = (n - 1) * sum_of_weights
        weighted_std = np.sqrt(numerator / denominator)
        return weighted_std      

    def set_key(self, dictionary, key, value):
        if key not in dictionary:
            dictionary[key] = [value]
        else:
            dictionary[key] = [*dictionary[key], value]