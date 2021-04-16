#! /usr/bin/env python
""" """

import numpy as np
import matplotlib.pyplot as plt

from toolbox import *

def ifile(ekey, pola, deckey, tupkey):
    return f'tuples/tuple_{deckey}_{tupkey}_{ekey}{pola}.root'

def lepton_theta_phi(cand):
    """ """
    pass

def plot_mpipi(data):
    """ """
    data['costh'] = data.apply(costh, axis=1)
    data_posi = data[data.pdgid_mc > 0]
    data_nega = data[data.pdgid_mc < 0]

    edge = np.cos(-np.pi/18)

    plt.hist(data_posi.costh, histtype='step', bins=100, density=True, range=[-edge, edge])
    plt.hist(data_nega.costh, histtype='step', bins=100, density=True, range=[-edge, edge])

    plt.xlim([-edge, edge])
    plt.show()

def plot_pi0_mass(data):
    plt.hist(data.tau_pi0_M, histtype='step', bins=100)
    plt.xlim((0.05, 0.30))
    plt.show()

def plot_pipi0_mass(data):
    plt.hist(data.M, histtype='step', bins=100)
    plt.show()

def make_tau_mc_momentum(data):
    data['tau_px_mc'] = data.tau_pi_px_mc + data.tau_pi0_gamma_px_mc + data.tau_pi0_gamma0_px_mc
    data['tau_py_mc'] = data.tau_pi_py_mc + data.tau_pi0_gamma_py_mc + data.tau_pi0_gamma0_py_mc
    data['tau_pz_mc'] = data.tau_pi_pz_mc + data.tau_pi0_gamma_pz_mc + data.tau_pi0_gamma0_pz_mc
    data['tau_ptot'] = np.sqrt(data.tau_px_mc**2 + data.tau_py_mc**2 + data.tau_pz_mc**2)
    data['tau_costh'] = data.tau_pz_mc / data.tau_ptot
    data['tau_phi'] = np.arctan2(data.tau_py_mc, data.tau_px_mc)

def plot_pipi_costh(data):
    data_posi = data[data.tau_pi_pdgid_mc > 0]
    data_nega = data[data.tau_pi_pdgid_mc < 0]
    plt.hist(data_posi.tau_costh, histtype='step', bins=100)
    plt.hist(data_nega.tau_costh, histtype='step', bins=100)
    plt.show()

def plot_pipi_phi(data):
    data_posi = data[data.tau_pi_pdgid_mc > 0]
    data_nega = data[data.tau_pi_pdgid_mc < 0]
    plt.hist(data_posi.tau_phi, histtype='step', bins=40)
    plt.hist(data_nega.tau_phi, histtype='step', bins=40)
    plt.show()

def main():
    """ """
    # for pola in ['posi', 'nega', 'zero']:
    #     muons = open_data(ifile('thr', pola, 'taulnu_taupipinu', 'muon'), 'muons')
    #     plot_mpipi(muons)

    pola = 'zero'
    taus = open_data(ifile('thr', pola, 'taulnu_taupipinu', 'tau'), 'taus')
    taus = taus[np.abs(taus.tau_pi_pdgid_mc) == 211].groupby('evtn').first()
    make_tau_mc_momentum(taus)
    
    # plot_pi0_mass(taus)
    # plot_pipi0_mass(taus)
    # plot_pipi_costh(taus)
    plot_pipi_phi(taus)

if __name__ == '__main__':
    main()
