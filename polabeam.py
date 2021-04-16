#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from toolbox import *

def cos_delta_th(pcl1, pcl2):
    mom1, mom2 = [three_momentum(item) for item in [pcl1, pcl2]]
    return mom1 @ mom2 / np.sqrt(mom1 @ mom1 * mom2 @ mom2)

def plot_costh(data):
    data['costh'] = data.apply(costh, axis=1)
    data_posi = data[data.pdgid_mc > 0]
    data_nega = data[data.pdgid_mc < 0]

    plt.hist(data_posi.costh, histtype='step', bins=40, density=True)
    plt.hist(data_nega.costh, histtype='step', bins=40, density=True)

    plt.xlim((-1, 1))
    plt.show()

def plot_delta_costh(data):
    dcosth = []
    for _, idata in data.groupby('evtn'):
        if idata.shape[0] != 2:
            continue
        dcosth.append(cos_delta_th(idata.iloc[0], idata.iloc[1]))

    print(len(dcosth))
    plt.hist(dcosth, histtype='step', bins=100, density=True)
    plt.xlim((-1, 1))
    plt.show()

if __name__ == '__main__':
    ifiles = [
        ['tupletaukaon_nega.root', 'kaons'],
        ['tupletaukaon_posi.root', 'kaons'],
        ['tupletaukaon_zero.root', 'kaons'],
        # ['tupletaupion_nega.root', 'pions'],
        # ['tupletaupion_posi.root', 'pions'],
        # ['tupletaupion_zero.root', 'pions'],
    ]

    for fname, tname in ifiles:
        data = open_data('/'.join(['tuples', fname]), tname)
        plot_costh(data)
        # plot_delta_costh(data)
