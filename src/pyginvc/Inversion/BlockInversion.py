import logging
import numpy as np
from scipy import optimize
from pyginvc.Inversion.BaseInversion import BaseInversion

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
    datefmt="%d-%M-%Y %H:%M:%S")

class BlockInversion(BaseInversion):
    def __init__(self, data, dict_weight={}, dict_bound={}):
        super().__init__(data, dict_weight, dict_bound)
            
    def assemble_design_matrix(self, nblock, G_gps_block):
        """
        Assemble desgin matrix for inversion.
        """

        nparam     = nblock*3
        W, _       = self.get_weight()
        G2I        = W @ G_gps_block
        G2R        = G_gps_block
        return G2I, G2R

    def run_inversion(self, nblock, G_gps_block):

        d2I, d2R = self.assemble_data_vector(0)
        
        # dimensions
        len_gps        = len(self.data.d_gps)

        # init parameters
        r              = np.zeros(len_gps)
        misfit         = 0.0
        euler          = np.zeros(3*nblock)

        nparam         = nblock*3
        bu             = np.full(nparam, np.inf)
        bl             = np.full(nparam,-np.inf)
        
        W, WSAR        = self.get_weight()

        G2I, G2R = self.assemble_design_matrix(nblock, G_gps_block)
        x, *_    = np.linalg.lstsq(G2I, d2I, rcond=None)
        logging.info('Unconstrained linear least square finished.')
            
        dhat     = G2R @ x
        r        = d2R -  W @ dhat
        misfit   = r.dot(r)
        logging.info('GPS Weighted Residual Sum of Squares (WRSS) %10.3f' %(misfit))
        logging.info('GPS: WRSS/(N) %f' %(misfit/len_gps))
        
        self.misfit     = misfit
        self.r          = r
        self.dhat       = dhat
        self.euler      = x
