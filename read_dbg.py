import numpy as np
import matplotlib.pyplot as plt

def read_debug(fname) :
    f=open(fname,"rb")
    fsize,ipix,nside,nside_spec,n_sub,n_nu,n_pol,n_comp,n_spec_vary,n_samples=np.fromfile(f,dtype=np.int32,count=10)
    if fsize==4 :
        print "is float"
        dt=np.float32
    else :
        print "is double"
        dt=np.float64
    chain=(np.fromfile(f,dtype=dt,count=n_spec_vary*n_samples)).reshape(n_samples,n_spec_vary)

    for i in np.arange(n_spec_vary) :
        plt.plot(chain[:,i])
        plt.show()

    for i in np.arange(n_spec_vary) :
        for j in np.arange(n_spec_vary-i-1)+i+1 :
            plt.plot(chain[:,i],chain[:,j],'.',markersize=0.1); plt.show()

    f.close()


read_debug("test.dbg")
