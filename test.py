import Gaussian

file = Gaussian.Initialize('a.log')

file.get_orbitals()

print(file.orbitals_energies)
print(file.orbitals_intensities)