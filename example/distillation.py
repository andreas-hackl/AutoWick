import numpy as np
import autowick
from autowick.contractions import combined_contraction, cross_contraction, sequential, loop, PerambBaryon

#### EXAMPLE OF DISTILLATION CODE ####
"""
TYPICAL OVERHEAD

class NAME(PerambBaryon):
    
    def __init__(self, t0, ALL ARGUMENTS FOR DISTILLATION OBJECTS AND GAMMA MATRICES, tval=None):
        # tval defines the time slices that are computed
        self.t0 = t0 # time slice of the source

        HANDLING OF ALL ARGUMENTS

        self.pions = {
            NAME OF PION: lambda t: (MOMENTUM INSERTION, GAMMA MATRIX),
        }
    
    # AUTOGENERATED CODE FROM CONTRACTOR

    def corr(self, t):
        return SIGN * PI0_PREFACTOR * sum(self.diagrams(t))

"""



class PerambNucleonPion2NucleonPion(PerambBaryon):

    def __init__(self, eps_src, eps_snk, perambs, mins_src, mins_snk, t0, Gamma_src=None, Gamma_snk=None, 
                 Gamma_pi_src=None, Gamma_pi_snk=None, P_src=None, P_snk=None, tval=None):
        self.t0 = t0
        self.eps_src = np.conjugate(eps_src[self.t0])
        self.eps_snk = eps_snk
        self.perambs = perambs
        self.Nt = perambs[0].shape[0]
        self.tval = range(self.Nt) if tval is None else tval

        self.mins_src = mins_src
        self.mins_snk = mins_snk
        
        self.P_src = np.eye(4) if P_src is None else P_src
        self.P_snk = np.eye(4) if P_snk is None else P_snk

        self.Gamma_src = autowick.mat.nucleon['Cg5'] if Gamma_src is None else Gamma_src
        self.Gamma_snk = autowick.mat.nucleon['Cg5'] if Gamma_snk is None else Gamma_snk

        self.Gamma_pi_src = autowick.mat.pion['g5'] if Gamma_pi_src is None else Gamma_pi_src
        self.Gamma_pi_snk = autowick.mat.pion['g5'] if Gamma_pi_snk is None else Gamma_pi_snk

        self.pions = {
            "pi+/src": lambda t: (self.mins_src[t], self.Gamma_pi_src),
            "pi+/snk": lambda t: (self.mins_snk[t], self.Gamma_pi_snk)
        }

    #######################################################################
    ### FROM HERE ON THE CODE IS AUTOGENERATED WITH examples/npi_2pt.py ###
    ### OUTPUT OF PRINT(dist)                                           ###
    #######################################################################

    def seq0(self, t): 
        return sequential([ self.perambs[t][t], self.perambs[self.t0][t] ], [self.pions["pi+/snk"](t) ])

    def seq1(self, t): 
        return sequential([ self.perambs[self.t0][t], self.perambs[self.t0][self.t0] ], [self.pions["pi+/src"](self.t0) ])

    def seq2(self, t): 
        return sequential([ self.perambs[t][t], self.perambs[self.t0][t], self.perambs[self.t0][self.t0] ], [self.pions["pi+/snk"](t), self.pions["pi+/src"](self.t0) ])

    def seq3(self, t): 
        return sequential([ self.perambs[self.t0][t], self.perambs[t][self.t0], self.perambs[self.t0][t] ], [self.pions["pi+/src"](self.t0), self.pions["pi+/snk"](t) ])

    def seq4(self, t): 
        return loop([ self.perambs[t][self.t0], self.perambs[self.t0][t] ], [self.pions["pi+/snk"](t), self.pions["pi+/src"](self.t0) ])

    def diagrams(self, t):
        S0 = self.seq0(t)
        S1 = self.seq1(t)
        S2 = self.seq2(t)
        S3 = self.seq3(t)
        S4 = self.seq4(t)


        def contraction(func, ops, factor):
            d0, d1 = func(self.eps_src, self.eps_snk[t], ops, self.P_snk, self.P_src)
            d0 *= factor
            d1 *= factor
            return d0, d1


        # shortcuts

        L = self.perambs[self.t0][t]
        # GX
        GS3 = self.Gamma_snk @ S3
        GS1 = self.Gamma_snk @ S1
        GL = self.Gamma_snk @ L

        # XG 
        LG = L @ self.Gamma_src
        S2G = S2 @ self.Gamma_src

        # GXG 
        GS1G = self.Gamma_snk @ S1 @ self.Gamma_src


        diags = []

        # combined contractions
        for pre, ops in [(-1, [LG, GL, S2 ]), (-1, [S0, GS1G, L ]), (-1, [S2G, GL, L ]), (-1, [LG, GS3, L ])]:
            di, dj = contraction(combined_contraction, ops, factor=pre)
            diags.append(di)
            diags.append(dj)


        # cross contractions
        for pre, ops in [(-1, [S0, LG, GS1 ])]:
            di, dj = contraction(cross_contraction, ops, factor=pre)
            diags.append(di)
            diags.append(dj)


        # LOOP CONTRACTIONS
        # combined contractions
        for pre, ops, loops in [(1, [LG, GL, L ], [S4])]:
            factor = pre * np.prod(loops)
            di, dj = contraction(combined_contraction, ops, factor=factor)
            diags.append(di)
            diags.append(dj)


        return diags
    
    ##########################################################################################
    #### END OF AUTOGENERATED CODE FROM examples/npi_2pt.py                                ###
    ##########################################################################################

    def corr(self, t):
        # IMPORTANT HERE: IF ONE USES PI0 ONE HAS TO MULTIPLY THE SUM WITH A 1/SQRT(2) PREFACTOR FOR EACH PI0 
        return - sum(self.diagrams(t))  # Minus sign comes from O_{pi+}^\dagger(p) = - O_{pi-}(-p)
    

if __name__ == '__main__':
    eps = np.load('elemental.npy')
    perambs = {t: np.load(f'peramb_t{t}.npy') for t in range(4)}
    mins = np.load("mins.npy")

    cont = PerambNucleonPion2NucleonPion(eps, eps, perambs, mins, mins, 0)
    print(cont.corr(0))