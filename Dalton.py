import numpy as np
import tools

AU2EV = 27.2114

class Initialize(tools.Tools):
    """
    Class to deal with Dalton output file
    1) To get the energies and intensities from TD-DFT calculations 
       to calculate OPA or TPA UV/VIS spectrum
    """
    def __init__(self, file):
        super(Initialize, self).__init__()
        tools.Tools().__init__()
        self.name = file
        self.file = file
        self.prog = 'dal'
        self.calc_type = {'.SINGLE RESIDUE':False,
                          '.TWO-PHOTON'    :False,}
        self.calc_bib  = {'.SINGLE RESIDUE':self.get_exitations,
                          '.TWO-PHOTON'    :self.get_a2f_exitations}

        self.which_calculation(maxlines=1000)

    def get_exitations(self):
        self.get_energies()
        self.get_intensities()

    def get_a2f_exitations(self, polar='  Linear  '):
        lines, _ = self.get_lines(polar)
        nmax = len(lines)
        self.energies = {}
        self.intensities = {}
        for n in range(nmax):
            lsplit = lines[n].split()
            state = int(lsplit[1])
            energy = float(lsplit[2])
            intensity = float(lsplit[7]) * 0.05 * np.pi
            self.set_key(self.energies, state, energy)
            self.set_key(self.intensities, state, intensity)

    def get_energies(self):
        lines, _ = self.get_lines('@ Excitation energy : ')
        state_lines, _ = self.get_lines('@ Excited state no')
        nmax = len(lines)
        self.energies = {}
        for n in range(nmax):
            state = int(state_lines[n].split()[4])
            energy = float(lines[n].split()[4]) * AU2EV
            self.set_key(self.energies, state, energy)
        

    def get_intensities(self):
        lines, _ = self.get_lines('@ Oscillator strength (LENGTH) ')
        self.intensities = {}
        if len(lines) % 3:
            print('The length of the Oscillator strength is not divisible by 3')
            return()
        state_lines, _ = self.get_lines('@ Excited state no')
        nmax = len(lines)
        nstates = int(nmax / 3)
        intensities = np.zeros(nmax)
        self.intensities = {}
        for n in range(nmax):
            intensities[n] = lines[n].split()[5]
        intensities = [ np.linalg.norm(x) for x in np.reshape(intensities, (nstates,3)) ]
        for n in range(nstates):
            state = int(state_lines[n].split()[4])
            self.set_key(self.intensities, state, intensities[n])