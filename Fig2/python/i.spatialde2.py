#!/usr/bin/env python3

import sys
import os
import pandas as pd
import numpy as np
import time
import scipy

# scipy version compatibility patch
minimal_patch = [
    'zeros_like', 'ones_like', 'argsort', 'array', 'arange',
    'sum', 'mean', 'var', 'std', 'min', 'max',
    'where', 'unique', 'sort', 'clip', 'dot',
    'exp', 'log', 'sqrt', 'square', 'abs', 'sign',
    'sin', 'cos', 'pi', 'inf', 'linspace', 'logspace',
    'log10', 'ones', 'zeros', 'eye',
    'transpose', 'reshape', 'flatten', 'ravel',
    'concatenate', 'hstack', 'vstack',
    'argmin', 'argmax', 'nonzero', 'isfinite', 'isnan', 'isinf'
]
for func_name in minimal_patch:
    if not hasattr(scipy, func_name) and hasattr(np, func_name):
        setattr(scipy, func_name, getattr(np, func_name))

import NaiveDE, SpatialDE


def main():
    if len(sys.argv) < 2:
        print("Error: SI value not provided")
        print("Usage: python i.spatialde2.py <SI_value>")
        sys.exit(1)
    
    si = float(sys.argv[1])
    sit = int(si) if si == int(si) else round(si, 2)
    
    print(f"Starting SpatialDE analysis for SI: {sit}")
    
    dpath = "data/sim.generated"
    rpath = "data/bhmk.res/i/SPATIALDE2"
    
    os.makedirs(rpath, exist_ok=True)
    
    try:
        count_path = os.path.join(dpath, "sim.isl", f"sim_isl_SI{sit}.csv")
        perm_path = os.path.join(dpath, "sim.isl", f"perm_isl_SI{sit}.csv")
        gt_path = os.path.join(dpath, "sim.isl", f"GT_SI{sit}.csv")
        coord_path = os.path.join(dpath, "coord.csv")
        
        if not os.path.exists(count_path):
            print(f"Error: File not found - {count_path}")
            sys.exit(1)
            
        counts = pd.read_csv(count_path, index_col=0)
        perm = pd.read_csv(perm_path, index_col=0)
        g = pd.read_csv(gt_path)
        coord = pd.read_csv(coord_path, index_col=0)
        coord_perm = coord.copy()
        

        t0 = time.time()
        # Calculate total_counts from ORIGINAL raw counts BEFORE any processing
        coord['total_counts'] = counts.sum(axis=0)
        
        # NaiveDE normalization
        counts = NaiveDE.stabilize(counts)
        
        # SpatialDE analysis
        resid_expr = NaiveDE.regress_out(coord, counts, 'np.log(total_counts)').T
        X = coord[['X1', 'X2']].values
        ressig = SpatialDE.run(X, resid_expr)
        
        t1 = time.time() - t0
        
        # Run on permuted data
        coord_perm['total_counts'] = perm.sum(axis=0)
        perm = NaiveDE.stabilize(perm)
        resid_expr_perm = NaiveDE.regress_out(coord_perm, perm, 'np.log(total_counts)').T
        resperm = SpatialDE.run(X, resid_expr_perm)
        
        gs = g['gs'].values
        gc = g['gc'].values
        thr = 0.05
        svgsim = ressig[ressig["qval"] < thr]["g"].values
        svgprm = resperm[resperm["qval"] < thr]["g"].values
        
        res = pd.DataFrame({
            'ID': ["Power_sim", "FPR_sim", "Power_prm", "FPR_prm"],
            'SPATIALDE2': [
                len(np.intersect1d(gs, svgsim)) / len(gs) if len(gs) > 0 else 0,
                len(np.intersect1d(gc, svgsim)) / len(gc) if len(gc) > 0 else 0,
                len(np.intersect1d(gs, svgprm)) / len(gs) if len(gs) > 0 else 0,
                len(np.intersect1d(gc, svgprm)) / len(gc) if len(gc) > 0 else 0
            ],
            'SI': [sit] * 4,
            'total_time': [t1] * 4,
            'mapping_time': [np.nan] * 4
        })
        
        res.to_csv(os.path.join(rpath, f"res_spatialde_SI{sit}.csv"), index=False)
        
    except Exception as e:
        print(f"Error during SpatialDE analysis: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()