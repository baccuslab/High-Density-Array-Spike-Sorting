# -*- coding:utf-8 -*-

##---IMPORTS

import scipy as sp
from scipy import linalg as sp_la
from common import TimeSeriesCovE
from plot import P, mcdata
from nsim.scene.noise import ArNoiseGen
from database import MunkSession


##---CONSTANTS

DB = MunkSession()


##---FUNCTIONS

def load_cmx():
    """load the a the covariance matrix for L011 t7 trial20"""

    return sp.load('/home/phil/Data/Munk/cmx_L011_t7_trl20.npz')['C']


def load_ce(ana=1216, trl=20):

    return DB.get_covariance_estimator(ana, trl)


def AR_from_C(C, nc, max_order=10):
    """build covariance sequence for lags [0,max_order+1] from block toeplitz matrix
    
    Here the input matrix C is assumed to exhibit a block structure with respect to the
    multichanneled data and the temporal lag, where the temporal per channel forms the
    Toeplitz structure, and the multichanneled part the block structure.
    
    The output is a sequence of one nc x nc covariance matrix per lag in the toeplitz
    structure of C. Also the AR coefficients for the order specified
    """

    # inits
    assert C.shape[0] == C.shape[1], 'C is not square'
    tf = C.shape[0] / float(nc)
    assert tf == sp.round_(tf), 'nc does not match tf'
    tf, nc = int(tf), int(nc)
    order = min(tf - 1, max_order + 1)
    print 'doing lags upto:', order
    gamma_y = sp.zeros((order, nc, nc))
    GammaY0 = sp.empty(((order - 1) * nc, (order - 1) * nc))
    # build gamma_y_seq from C
    for o in xrange(order):
        for c1 in xrange(nc):
            for c2 in xrange(nc):
                gamma_y[o, c1, c2] = C[c1 * tf + o, c2 * tf]
    # build GammaY0 from gamma_y_seq, blockwise
    for i in xrange(order - 1):
        for j in xrange(i, order - 1):
            GammaY0[i * nc:(i + 1) * nc, j * nc:(j + 1) * nc] = gamma_y[j - i]
            if j > i:
                GammaY0[j * nc:(j + 1) * nc, i * nc:(i + 1) * nc] = gamma_y[j - i].T
    # build calculation matrices
    GammaY0inv = sp_la.inv(GammaY0)
    gamma_y_seq = sp.hstack(gamma_y[1:])

    # return
    return sp.dot(gamma_y_seq, GammaY0inv), gamma_y


def Sigma_u_from_C0(CE, chan_set=(0, 1, 2, 3)):
    """produce innovation process covariance matrix from C0"""

    Gamma_y0 = CE.get_cmx(tf=1, chan_set=chan_set)
    U, s, Vh = sp_la.svd(Gamma_y0)
    Sigma_u_hat = sp.dot(sp.dot(sp.diag(sp.sqrt(s)), U), Vh)
    return Sigma_u_hat


def matshow_stuff(*args):
    """matshow all matrices passed"""

    do_names = False
    if args[0] is True:
        do_names = True
    for i in xrange(int(do_names), len(args), 1 + int(do_names)):
        try:
            item = args[i]
            P.matshow(item)
            P.colorbar()
            if do_names:
                print args[i + 1]
                P.title(str(args[i + 1]))
            print item.shape
        except:
            pass


def test():

    CS = (0, 1, 2, 3)
    TF = 65
    ORDER = 16
    noise_ce = load_ce()
    noise_cmx = noise_ce.get_cmx(tf=TF, chan_set=CS)
    A, gamma = AR_from_C(noise_cmx, 4, ORDER)
    Sigma_u_hat = Sigma_u_from_C0(noise_ce, chan_set=CS)
    NG = ArNoiseGen((A, Sigma_u_hat), check=True)
    x = NG.query(100000)
    mcdata(x, show=False)
    CE = TimeSeriesCovE(tf_max=TF, nc=4)
    CE.new_chan_set(CS)
    CE.update(x)
    matshow_stuff(True,
                  noise_ce.get_cmx(tf=TF, chan_set=CS), 'Original Noise-Covariance',
                  CE.get_cmx(tf=TF, chan_set=CS), 'AR Noise-Covariance',
                  NG.coeffs, 'AR Coefficients',
                  NG.sigma, 'AR Innovation-Covariance')
    P.show()


##---MAIN

if __name__ == '__main__':

    test()
