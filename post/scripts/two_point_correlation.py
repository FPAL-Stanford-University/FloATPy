import h5py
import numpy as np
import pyfftw
import sys

import floatpy.readers.padeops_reader as por

def two_point_correlation(Y_CO2, f, fhat, fft, ifft):
    # Remove mean and scale by variance
    f[:,:,:] = ( Y_CO2[:,:,:] - np.mean(Y_CO2,axis=2,keepdims=True) )

    # Get Fourier transform
    fft(input_array=f, output_array=fhat)

    # Get power spectral density
    psd = fhat * fhat.conjugate()

    # Get autocorrelation
    ifft(input_array=psd, output_array=f)

    # Define mask function for averaging in the x-y plane
    phi = 4.*Y_CO2*(1.-Y_CO2)
    phi[phi > 0.1] = 1.0
    phi[phi < 0.9] = 0.0

    # Compute averaged autocorrelation
    Y_corr = np.sum( f / (f[:,:,0:1] + 1.e-10) * phi, axis=(0,1) ) / np.sum( phi, axis=(0,1) )
    r_z = np.linspace(0., L_z, num=Y_corr.shape[0]+1)[:-1]

    return r_z, Y_corr


if len(sys.argv) < 3:
    print "Usage: "
    print "    python %s <resolution> <large span (True or False)> <step(s)>" %sys.argv[0]
    sys.exit()

resolution = int(sys.argv[1])
steps = [ int(sys.argv[i]) for i in range(3,len(sys.argv)) ]
# yslice = 90

large_span = False
if sys.argv[2] == 'True':
    large_span = True

L_y = 0.1143

if large_span:
    L_z = 1.0*L_y
    directory = '../LS_N%04d' %resolution
    output_prefix = 'two_point_correlation_data/LS_N%04d_' %resolution
else:
    L_z = 0.25*L_y
    directory = '../N%04d' %resolution
    output_prefix = 'two_point_correlation_data/N%04d_' %resolution

filename_prefix = directory + '/inclinedRM_'
reader = por.PadeopsReader(filename_prefix, periodic_dimensions=(False,False,True))

nx, ny, nz = reader.domain_size

f    = pyfftw.empty_aligned((nx, ny, nz), dtype='float64', order='F')
fhat = pyfftw.empty_aligned((nx, ny, nz/2+1), dtype='complex128', order='F')

fft  = pyfftw.FFTW(f, fhat, axes=(2,),direction='FFTW_FORWARD',  flags=('FFTW_MEASURE', ),                    threads=4, planning_timelimit=None)
ifft = pyfftw.FFTW(fhat, f, axes=(2,),direction='FFTW_BACKWARD', flags=('FFTW_MEASURE', ),                    threads=4, planning_timelimit=None)


# Only get a single x-z slice
# reader.sub_domain    = (0,yslice,0), (reader.domain_size[0]-1,    yslice, reader.domain_size[2]-1)
# x, y, z = reader.readCoordinates()
# steps = sorted(reader.steps)

for step in steps:
    outputfile = output_prefix + 'two_point_correlation_%04d.dat' %step
    reader.step = step
    u, v, w, Y_CO2 = reader.readData(('u', 'v', 'w', 'Massfraction_02'))

    # nx, ny, nz = Y_CO2.shape

    r_z, Y_corr = two_point_correlation(Y_CO2, f, fhat, fft, ifft)
    r_z, u_corr = two_point_correlation(u, f, fhat, fft, ifft)
    r_z, v_corr = two_point_correlation(v, f, fhat, fft, ifft)
    r_z, w_corr = two_point_correlation(w, f, fhat, fft, ifft)

    # Write to file
    with open(outputfile, 'w') as fout:
        fout.write('# %24s %26s %26s %26s %26s\n' %('r_z', 'Y_CO2', 'u', 'v', 'w'))
        for i in range(r_z.shape[0]):
            fout.write('%26.16e %26.16e %26.16e %26.16e %26.16e\n' %(r_z[i], Y_corr[i], u_corr[i], v_corr[i], w_corr[i]))

    print "Finished saving data to %s" %outputfile
