'''
au : atomic unit (properly written a.u.)
kcal ~ kcal/mol is a handy unit of energy
cm   ~ cm^{-1} or "wavenumber" is a handy unit of energy
ang ~ Ångström is unit of length
eV  ~ electron-volt, unit of energy

femtoSec2au, kcal2cm, ang2au, au2cm : conversion factor to change units as it is written.

mass_n_au : mass of the neutron in a.u.
mass_p_au : mass of the neutron in a.u.
mass_e_au : mass of the neutron in a.u.
mX : mass of the atom X in a.u (X= H, Mu, etc...) Mu is for Muonium and H for Hydrogen.
mass_n_au : mass of the neutron in a.u.
'''
femtoSec2au = 41.341374578
kcal2cm = 349.75
au2cm = 27.211*8065.5
ang2au = 1.88973
eV2cm = 8065.6
mass_n_au = 1838.6854702982287
mass_p_au = 1836.1528047787378
mass_mu_au = (1.0/9.0)*mass_p_au
mass_e_au = 1.0
mH = mass_p_au+1
mMu = 0.113*(mass_p_au+1)
