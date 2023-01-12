import numpy as np

def smoothstep(in_vec, lo_edge, hi_edge, lo_val, hi_val):
    out_vec = np.zeros(np.shape(in_vec))
    
    for k, x in enumerate(in_vec):
        if x < lo_edge:
            out_vec[k] = lo_val
        elif x >= hi_edge:
            out_vec[k] = hi_val
        else:
            xs = (x - lo_edge) / (hi_edge - lo_edge)
            
            out_vec[k] = (hi_val-lo_val) * xs * xs * xs * (xs * (xs * 6.0 - 15.0) + 10.0) + lo_val
            
    return out_vec

def linstep(in_vec, lo_edge, hi_edge, lo_val, hi_val):
    out_vec = np.zeros(np.shape(in_vec))
    
    for k, x in enumerate(in_vec):
        if x <= lo_edge:
            out_vec[k] = lo_val
        elif x > hi_edge:
            out_vec[k] = hi_val
        else:
            xs = (x - lo_edge) / (hi_edge - lo_edge)
            out_vec[k] = lo_val + (hi_val-lo_val) * xs
            
    return out_vec
