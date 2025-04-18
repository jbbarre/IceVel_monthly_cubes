import numpy as np
from progress_bar import progress_bar

def cube_calib_using_low_speed(self,inplace=False):

    vx_med = np.nanmedian(self.vx_().filled(np.nan),axis=0)
    vy_med = np.nanmedian(self.vy_().filled(np.nan),axis=0)

    v_med = np.sqrt(vx_med**2+vy_med**2)

    if inplace:
        c2 = self
    else:
        c2 = self.deepcopy()

    for i, d in enumerate(c2.data):
        progress_bar(i/c2.nz*100.)
        if d.vx[v_med<100].count() > 100:
            difx = np.ma.median(d.vx[v_med<100]-vx_med[v_med<100])
            if difx > 10:
                d.vx = d.vx - difx
            
            dify = np.ma.median(d.vy[v_med<100]-vy_med[v_med<100])
            if dify > 10:
                d.vy = d.vy - dify
        else:
            if d.vx[v_med<200].count() > 100:
                difx = np.ma.median(d.vx[v_med<200]-vx_med[v_med<200])
                if difx > 30:
                    d.vx = d.vx - difx
            
                dify = np.ma.median(d.vy[v_med<100]-vy_med[v_med<100])
                if dify > 30:
                    d.vy = d.vy - dify
            #else:
            #    difx = np.ma.median(d.vx-vx_med)
            #    if difx > 250:
            #        d.vx = d.vx - difx
            
            #    dify = np.ma.median(d.vy-vy_med)
            #    if dify > 250:
            #        d.vy = d.vy - dify

    return c2