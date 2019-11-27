import numpy as np
import tools

class Initialize(tools.Tools):
    """
    Class to deal with Gaussian output file
    1) To get the energies and intensities from TD-DFT caluclations 
       to calculate the UV/VIS spectrum
    """
    def __init__(self, file):
        super(Initialize, self).__init__()
        tools.Tools().__init__()
        self.name = file
        self.file = file
        self.prog = 'gau'
        self.calc_type = {'TD':False,
                          'CIS':False,
                          #'CIS(D)':False,
                          #'IR':False,
                          #'RAMAN':False
                          }
        self.calc_bib  = {'TD':self.get_exitations,
                          'CIS':self.get_exitations,
                          #'CIS(D)':False,
                          #'IR':False,
                          #'RAMAN':False
                          }

        self.which_calculation()
            
    def get_exitations(self):
        lines = self.get_lines('Excited State  ')
        nmax = len(lines)
        self.energies = {}
        self.intensities = {}
        for n in range(nmax):
            lsplit = lines[n].split()
            state = int(lsplit[2][:-1])
            energy = float(lsplit[4])
            intensity = float(lsplit[8][2:])
            self.set_key(self.energies, state, energy)
            self.set_key(self.intensities, state, intensity)

    def get_transition_dipoles(self):
        string = 'Ground to excited state transition electric dipole moments (Au)'
        _, line1 = self.get_lines(string, index=True)
        line1 = line1[0] + 2
        if not self.nstates:
            string = 'Ground to excited state transition velocity dipole moments (Au)'
            _, line2 = self.get_lines(string, index=True)
            line1 = line2[0]
        else:
            line2 = line1 + self.nstates
        lines = self.get_lines_interval((line1, line2))
        self.tr_dipoles = np.array([np.array(x.split()[1:4], float) for x in lines])

    def get_permanent_dipoles(self):
        string = 'Dipole moment (field-independent basis, Debye)'
        _, index = self.get_lines(string, index=True)
        lines = []
        for n in index:
            lines += self.get_lines_interval((n + 1, n + 2))
        #   X=             -9.8637    Y=              3.6963    Z=             -1.1599  Tot=             10.5972
        lines = [x.replace('X=','').replace('Y=','').replace('Z=','') for x in lines]
        self.pe_dipoles = np.array([np.array(x.split()[:3], float) for x in lines])
