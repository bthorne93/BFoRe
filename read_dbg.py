import numpy as np
import matplotlib.pyplot as plt

def read_debug(fname) :
    f=open(fname,"rb")
    fsize,ipix,nside,nside_spec,n_sub,n_nu,n_pol,n_comp,n_spec_vary,n_samples,nb_samples=np.fromfile(f,dtype=np.int32,count=11)
    if fsize==4 :
        print "is float"
        dt=np.float32
    else :
        print "is double"
        dt=np.float64
    chain=(np.fromfile(f,dtype=dt,count=n_spec_vary*(n_samples+nb_samples))).reshape(n_samples+nb_samples,n_spec_vary)

    for i in np.arange(n_spec_vary) :
        plt.hist(chain[nb_samples:,i],bins=100); plt.show()
        plt.plot(chain[nb_samples:,i]); plt.show()
        for j in np.arange(n_spec_vary-i-1)+i+1 :
            plt.plot(chain[nb_samples:,i],chain[nb_samples:,j],'.',markersize=0.1); plt.show()

    f.close()


read_debug("test.dbg")
