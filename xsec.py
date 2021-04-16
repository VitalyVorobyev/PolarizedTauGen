""" """

mtau = 1.77682


def xsec_max():
    """ maximal sigma(e+e- -> tau+ tau-) (nb) """
    return 3.562


def total_xsec(beta, w1, w2):
    """ sigma(e+e- -> tau+ tau-) in units of maximal xsec """
    return beta * (1 - beta**2) * (3 - beta**2) * (1 + w1 * w2)


def main():
    import matplotlib.pyplot as plt
    import numpy as np

    sqrts = np.linspace(2*mtau, 6)
    gamma = 0.5 * sqrts / mtau
    beta = np.sqrt(1 - 1./gamma**2)
    plt.plot(sqrts, total_xsec(beta, 0.8, 0.0))
    plt.show()

if __name__ == '__main__':
    main()
