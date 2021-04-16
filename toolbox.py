import numpy as np
import uproot

mass_dict = {321 : 0.493677, 211 : 0.13957061}

def open_data(ifile, treename):
    df = uproot.open(ifile)[treename].pandas.df()
    print(df.head())
    print(df.shape)
    return df

def mass(pcl):
    return mass_dict[np.abs(pcl.pigid_mc)]

def three_momentum(pcl):
    return np.array([pcl.px_mc, pcl.py_mc, pcl.pz_mc])

def costh(pcl):
    momentum = three_momentum(pcl)
    return pcl.pz_mc / np.sqrt(momentum @ momentum)
