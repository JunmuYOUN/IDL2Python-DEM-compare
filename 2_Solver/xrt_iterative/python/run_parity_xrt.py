import os, sys
import numpy as np
from scipy.io import readsav
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)

def load(oracle_dir, pid):
    rec=np.atleast_1d(readsav(os.path.join(oracle_dir,'probe_%s.sav'%pid),verbose=False)['chk'])[0]
    return {n.lower():np.asarray(rec[n]) for n in rec.dtype.names if n.lower()!='probe_id'}

def main():
    od, pd = sys.argv[1], sys.argv[2]
    os.makedirs(pd, exist_ok=True); os.environ['CHK_DIR']=pd
    from xrt_demstat import xrt_iter_demstat
    a=load(od,'00_input'); b=load(od,'00b_lines')
    dem=a['dem_in']; t=a['t']
    emis_arr=np.transpose(b['emis_arr'])         # (nt,n_line)->(n_line,nt)
    i_obs=b['i_obs']; i_err=b['i_err']
    n_line=emis_arr.shape[0]; nt=t.shape[0]
    from chk_dump import chk_dump
    chk_dump('00_input', dem_in=dem, t=t)
    chk_dump('00b_lines', emis_arr=emis_arr, i_obs=i_obs, i_err=i_err)
    line_t=np.broadcast_to(t,(n_line,nt)).copy() # native grid == common grid
    line_nt=np.full(n_line, nt)
    i_mod,di,chisq=xrt_iter_demstat(emis_arr,line_t,line_nt,i_obs,i_err,dem,t)
    print('n_line',n_line,'chisq',float(chisq),'i_mod',i_mod)

if __name__=='__main__': main()
