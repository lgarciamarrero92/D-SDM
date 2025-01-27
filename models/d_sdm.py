import numpy as np
from scipy.special import lambertw
from core.params import Iph, Isat, eta, Rs, Rsh, Ns, nNsVth
from core.iv import IVModel
from utils.constants import BOLTZMANN_CONSTANT, ELECTRON_VOLT
from utils.helpers import celsius_to_kelvin
from models.sdm import SDM
from models.sdm_4param import SDM4PARAM
import pandas as pd

class D_SDM(IVModel):

    def __init__(self,tm=None,G=None,alpha=None,ns=None,Voc=None,Isc=None,method="lambertw_v_from_i",estimate_Iph2_fromG=False,no_estimate_eta=False, bounds={}):
        self.name = "D-SDM"
        self.method = method
        self.bounds = bounds
        self.tm = tm
        self.Voc = Voc
        self.G = G
        self.Isc = Isc
        self.alpha = alpha
        self.ns = ns
        self.estimate_Iph2_fromG = estimate_Iph2_fromG
        self.no_estimate_eta = no_estimate_eta

        if self.tm is not None and self.no_estimate_eta == False:
            self.params = [
                Iph(instance=1),
                Isat(instance=1),
                eta(instance=1,bounds=(1,2)),
                Rs(instance=1,bounds = (0.001,0.1)),
                Rsh(instance=1,bounds=(1,200)),
                Ns(instance=1,bounds=(1,ns)),
                Iph(instance=2) if self.estimate_Iph2_fromG == False else None
            ]
        else:
            self.params = [
                Iph(instance=1),
                Isat(instance=1),
                nNsVth(),
                Rs(instance=1,bounds = (0.001,0.1)),
                Rsh(instance=1,bounds=(1,200)),
                Ns(instance=1,bounds=(1,ns)),
                Iph(instance=2) if self.estimate_Iph2_fromG == False else None,
            ]

        # Filter out the None values if the conditional is False
        self.params = [param for param in self.params if param is not None]
        
        super().__init__(params=self.params)

    def get_all_params(self,vars):
        #This params refer to one cell
        if self.estimate_Iph2_fromG == False:
            Iph1, Isat, etaORnNsVth, Rs, Rsh, Ns1, Iph2 = vars
        else:
            Iph1, Isat, etaORnNsVth, Rs, Rsh, Ns1 = vars
            Iph2 = (self.G/1000)*(self.Isc*((Rsh+Rs)/Rsh)+self.alpha*(self.tm-25))
        
        ns1 = Ns1
        ns2 = self.ns - ns1

        Rs1 = Rs
        Rs2 = Rs

        Isat1 = Isat
        Isat2 = Isat1
        
        Rsh1 = Rsh

        return (Iph1, Isat1, etaORnNsVth, Rs1, Rsh1, ns1, Iph2, Isat2, etaORnNsVth, Rs2, ns2)

    def _calc(self,vars,v_points=None,i_points=None):
        
        Iph1, Isat1, etaORnNsVth, Rs1, Rsh1, ns1, Iph2, Isat2, etaORnNsVth, Rs2, ns2 = self.get_all_params(vars)
        
        if self.no_estimate_eta == False:
            sdm1 = SDM(self.tm,1,method='lambertw_v_from_i')
            sdm2 = SDM4PARAM(self.tm,1,method='v_from_i')
        else:
            sdm1 = SDM(None,1,method='lambertw_v_from_i')
            sdm2 = SDM4PARAM(None,1,method='v_from_i')

        veps = 0.1 #V
        
        if self.method == 'lambertw_v_from_i':
            v1 = sdm1.calc((Iph1, Isat1, etaORnNsVth, Rs1, Rsh1), v_points, i_points)[0]
            v2 = sdm2.calc((Iph2, Isat2, etaORnNsVth, Rs2), v_points, i_points)[0]
            condition = (v1 + veps >= 0) & (v2 + veps >= 0)
            v_out = np.full(v_points.shape, np.nan, dtype=np.float64)
            v_out[condition] = ns1*v1[condition] + ns2*v2[condition]

            return (v_out,i_points)
        
        elif self.method == 'lambertw_bisection_i_from_v':
            max_i = np.max(i_points)

            i_initial = np.full_like(i_points, 0.0)
            i_final = np.full_like(i_points, 2.0 * max_i)

            valid_i_final = np.full_like(i_points, False, dtype=bool)

            eps = 1e-4

            while np.any((i_final - i_initial) > eps):
                i_midd = (i_initial + i_final) / 2.0
                v1 = sdm1.calc((Iph1, Isat1, etaORnNsVth, Rs1, Rsh1), v_points, i_midd)[0]
                v2 = sdm2.calc((Iph2, Isat2, etaORnNsVth, Rs2), v_points, i_midd)[0]
                condition = (v1 + veps >= 0) & (v2 + veps >= 0)
                v_out = np.full(v_points.shape, -np.inf, dtype=np.float64)
                v_out[condition] = ns1*v1[condition] + ns2*v2[condition]

                large_v = (v_out > v_points)
                short_v = (v_out <= v_points)

                i_initial[large_v] = i_midd[large_v]
                i_final[short_v] = i_midd[short_v]
                valid_i_final[~condition & short_v] = False
                valid_i_final[condition & short_v] = True

            I = (i_initial + i_final) / 2.0
            I[valid_i_final == False] = np.nan
            return (np.copy(v_points), I)
        else:
            raise NotImplementedError("Method not implemented")
        
    def formatSol(self, x, error):
        Iph1, Isat1, etaORnNsVth, Rs1, Rsh1, ns1, Iph2, Isat2, etaORnNsVth, Rs2, ns2 = self.get_all_params(x)

        dict = {
            "$I_{ph{\text -}a}[A]$": [Iph1],
            "$I_s[A]$": [Isat1],
            "$eta$" if self.tm is not None and self.no_estimate_eta == False else "$nNsVth$": [etaORnNsVth],
            "$R_s[\Omega]$": [Rs1],
            "$R_{sh}[\Omega]$": [Rsh1],
            "$N_{s{\text -}a}$": [ns1],
            "$I_{ph{\text -}b}[A]$": [Iph2],
            "$Error$": [error]
        }
        
        return pd.DataFrame.from_dict(dict)
        