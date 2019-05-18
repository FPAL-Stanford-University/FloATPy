import numpy as np

def ddx1D(f,x):
    N = np.size(f)
    L = max(x)-min(x)
    dk = 2.0*np.pi/L
    k = [dk*wavenum for wavenum in range(-N/2,N/2-1)];
    fhat = np.fft.fft(f[0:-1])
    dfhat = fhat*1j*k;

    # zero the oddball
    dfhat[0] = 0
    df = np.fft.ifft(dfhat)
   
    # make periodic again
    df = np.hstack((df,df[0]))
    return df

# some tests
if __name__ == '__main__':

    Lx = 10
    Nx = 100
    x = np.linspace(0,Lx,Nx)

    f = np.sin(x)# + np.cos(3*x) + np.sin(x)**2
    dfdx_exact = np.cos(x) #- 3*np.sin(3*x) + 2*np.sin(x)*np.cos(x)

    # test ddx:
    dfdx = ddx1D(f,x)

    err = sum(np.sqrt(((dfdx-dfdx_exact)**2)))
    print(err)

