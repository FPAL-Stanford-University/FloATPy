from mpi4py import MPI
import numpy

from floatpy.parallel import transpose_wrapper
from floatpy.readers import moving_grid_data_reader, parallel_reader
import floatpy.readers.padeops_reader as por

## NEED TO CHANGE THE PATH OF DATA!
grid = 'grid_B'
directory_name = '/net/scratch1/manlong1/convert_data/3D_Poggi_RMI_RD/case_1_1/' + grid

output_filenames_k = [ \
    'spectrum_TKE_k_case_1_1_t_m0_05_' + grid, \
     'spectrum_TKE_k_case_1_1_t_0_05_' + grid, \
     'spectrum_TKE_k_case_1_1_t_0_40_' + grid, \
     'spectrum_TKE_k_case_1_1_t_0_75_' + grid, \
     'spectrum_TKE_k_case_1_1_t_1_10_' + grid, \
     'spectrum_TKE_k_case_1_1_t_1_20_' + grid, \
     'spectrum_TKE_k_case_1_1_t_1_30_' + grid, \
     'spectrum_TKE_k_case_1_1_t_1_40_' + grid, \
     'spectrum_TKE_k_case_1_1_t_1_50_' + grid, \
     'spectrum_TKE_k_case_1_1_t_1_60_' + grid, \
     'spectrum_TKE_k_case_1_1_t_1_75_' + grid]

output_filenames_E_x = [ \
    'spectrum_TKE_E_x_case_1_1_t_m0_05_' + grid, \
     'spectrum_TKE_E_x_case_1_1_t_0_05_' + grid, \
     'spectrum_TKE_E_x_case_1_1_t_0_40_' + grid, \
     'spectrum_TKE_E_x_case_1_1_t_0_75_' + grid, \
     'spectrum_TKE_E_x_case_1_1_t_1_10_' + grid, \
     'spectrum_TKE_E_x_case_1_1_t_1_20_' + grid, \
     'spectrum_TKE_E_x_case_1_1_t_1_30_' + grid, \
     'spectrum_TKE_E_x_case_1_1_t_1_40_' + grid, \
     'spectrum_TKE_E_x_case_1_1_t_1_50_' + grid, \
     'spectrum_TKE_E_x_case_1_1_t_1_60_' + grid, \
     'spectrum_TKE_E_x_case_1_1_t_1_75_' + grid]

output_filenames_E_y = [ \
    'spectrum_TKE_E_y_case_1_1_t_m0_05_' + grid, \
     'spectrum_TKE_E_y_case_1_1_t_0_05_' + grid, \
     'spectrum_TKE_E_y_case_1_1_t_0_40_' + grid, \
     'spectrum_TKE_E_y_case_1_1_t_0_75_' + grid, \
     'spectrum_TKE_E_y_case_1_1_t_1_10_' + grid, \
     'spectrum_TKE_E_y_case_1_1_t_1_20_' + grid, \
     'spectrum_TKE_E_y_case_1_1_t_1_30_' + grid, \
     'spectrum_TKE_E_y_case_1_1_t_1_40_' + grid, \
     'spectrum_TKE_E_y_case_1_1_t_1_50_' + grid, \
     'spectrum_TKE_E_y_case_1_1_t_1_60_' + grid, \
     'spectrum_TKE_E_y_case_1_1_t_1_75_' + grid]

output_filenames_E_z = [ \
    'spectrum_TKE_E_z_case_1_1_t_m0_05_' + grid, \
     'spectrum_TKE_E_z_case_1_1_t_0_05_' + grid, \
     'spectrum_TKE_E_z_case_1_1_t_0_40_' + grid, \
     'spectrum_TKE_E_z_case_1_1_t_0_75_' + grid, \
     'spectrum_TKE_E_z_case_1_1_t_1_10_' + grid, \
     'spectrum_TKE_E_z_case_1_1_t_1_20_' + grid, \
     'spectrum_TKE_E_z_case_1_1_t_1_30_' + grid, \
     'spectrum_TKE_E_z_case_1_1_t_1_40_' + grid, \
     'spectrum_TKE_E_z_case_1_1_t_1_50_' + grid, \
     'spectrum_TKE_E_z_case_1_1_t_1_60_' + grid, \
     'spectrum_TKE_E_z_case_1_1_t_1_75_' + grid]

time_steps = range(11)

# Create the parallel data reader.

comm = MPI.COMM_WORLD

## NEED TO CHANGE THE DATA READER!
serial_reader = moving_grid_data_reader.MovingGridDataReader(directory_name, \
    data_order = 'F')

reader = parallel_reader.ParallelDataReader(comm, serial_reader)
steps = reader.steps
reader.setStep(steps[time_steps[0]])

grid_partition = reader.grid_partition

domain_size = reader.domain_size
Nx_full = domain_size[0]
Ny_full = domain_size[1]
Nz_full = domain_size[2]

lo_c, hi_c = reader.full_chunk

x_coords, y_coords, z_coords = reader.readCoordinates()

dx = x_coords[1,0,0] - x_coords[0,0,0]
dy = y_coords[0,1,0] - y_coords[0,0,0]
dz = z_coords[0,0,1] - z_coords[0,0,0]

Lx_full = Nx_full*dx
Ly_full = Ny_full*dy
Lz_full = Nz_full*dz

del x_coords, y_coords, z_coords

# Set up the transpose wrappers.

tw_y = transpose_wrapper.TransposeWrapper(reader.grid_partition, direction=1, dimension=3)
lo_y, hi_y = tw_y.full_pencil
Nx_y = hi_y[0] - lo_y[0] + 1
Ny_y = hi_y[1] - lo_y[1] + 1
Nz_y = hi_y[2] - lo_y[2] + 1
assert Ny_y == Ny_full

tw_z = transpose_wrapper.TransposeWrapper(reader.grid_partition, direction=2, dimension=3)
lo_z, hi_z = tw_z.full_pencil
Nx_z = hi_z[0] - lo_z[0] + 1
Ny_z = hi_z[1] - lo_z[1] + 1
Nz_z = hi_z[2] - lo_z[2] + 1
assert Nz_z == Nz_full

M_SF6 = 146.055
M_air = 28.0135

for n in range(len(time_steps)):
    reader.setStep(steps[time_steps[n]])
    
    if comm.rank == 0:
        print("Time: %g" % reader.time)
    
    # Determine the location of mixing layer.
    
    Y0, = reader.readData('mass fraction 0')
    
    M = Y0/M_SF6 + (1.0 - Y0)/M_air
    M = 1.0/M
    
    X0 = Y0*M/M_SF6
    del Y0, M
    
    fcomm_yz = grid_partition.commyz()
    comm_yz = MPI.Comm.f2py(fcomm_yz)
    
    X0_sum_yz_local = numpy.sum(numpy.sum(X0, axis=2), axis=1)
    X0_sum_yz_global = numpy.empty(X0_sum_yz_local.shape)
    
    comm_yz.Allreduce([X0_sum_yz_local, MPI.DOUBLE], [X0_sum_yz_global, MPI.DOUBLE], op=MPI.SUM)
    phi = X0_sum_yz_global/(Ny_full*Nz_full)
    weight = 4.0*phi*(1.0 - phi)
    mixing_layer_indices, = numpy.where( numpy.array([e > 0.9 if ~numpy.isnan(e) else False for e in weight], dtype=bool) )
    
    within_mixing_layer_local = numpy.full(Nx_full, False, dtype=numpy.bool)
    within_mixing_layer_local[mixing_layer_indices + lo_c[0]] = True
    
    del X0
    
    # Compute the TKE spectra in the mixing layer and sum them up.
    
    rho, vel = reader.readData(('density', 'velocity'))
    
    rho_sum_yz_local = numpy.sum(numpy.sum(rho, axis=2), axis=1)
    rho_sum_yz_global = numpy.empty(rho_sum_yz_local.shape)
    comm_yz.Allreduce([rho_sum_yz_local, MPI.DOUBLE], [rho_sum_yz_global, MPI.DOUBLE], op=MPI.SUM)
    rho_mean_yz = rho_sum_yz_global/(Ny_full*Nz_full)
    
    rho_u = rho*vel[:,:,:,0]
    rho_u_sum_yz_local = numpy.sum(numpy.sum(rho_u, axis=2), axis=1)
    rho_u_sum_yz_global = numpy.empty(rho_u_sum_yz_local.shape)
    comm_yz.Allreduce([rho_u_sum_yz_local, MPI.DOUBLE], [rho_u_sum_yz_global, MPI.DOUBLE], op=MPI.SUM)
    rho_u_mean_yz = rho_u_sum_yz_global/(Ny_full*Nz_full)
    u_mean_yz = rho_u_mean_yz/rho_mean_yz
    del rho_u
    
    rho_v = rho*vel[:,:,:,1]
    rho_v_sum_yz_local = numpy.sum(numpy.sum(rho_v, axis=2), axis=1)
    rho_v_sum_yz_global = numpy.empty(rho_v_sum_yz_local.shape)
    comm_yz.Allreduce([rho_v_sum_yz_local, MPI.DOUBLE], [rho_v_sum_yz_global, MPI.DOUBLE], op=MPI.SUM)
    rho_v_mean_yz = rho_v_sum_yz_global/(Ny_full*Nz_full)
    v_mean_yz = rho_v_mean_yz/rho_mean_yz
    del rho_v
    
    rho_w = rho*vel[:,:,:,2]
    rho_w_sum_yz_local = numpy.sum(numpy.sum(rho_w, axis=2), axis=1)
    rho_w_sum_yz_global = numpy.empty(rho_w_sum_yz_local.shape)
    comm_yz.Allreduce([rho_w_sum_yz_local, MPI.DOUBLE], [rho_w_sum_yz_global, MPI.DOUBLE], op=MPI.SUM)
    rho_w_mean_yz = rho_w_sum_yz_global/(Ny_full*Nz_full)
    w_mean_yz = rho_w_mean_yz/rho_mean_yz
    del rho_w
    
    u = vel[:,:,:,0]
    v = vel[:,:,:,1]
    w = vel[:,:,:,2]
    
    del vel
    
    TKE_in_x = numpy.empty(rho.shape)
    for idx in range(rho.shape[0]):
        u_prime = u[idx,:,:] - u_mean_yz[idx]
        TKE_in_x[idx,:,:] = numpy.sqrt(rho[idx,:,:])*u_prime
    del u
    
    TKE_in_y = numpy.empty(rho.shape)
    for idx in range(rho.shape[0]):
        v_prime = v[idx,:,:] - v_mean_yz[idx]
        TKE_in_y[idx,:,:] = numpy.sqrt(rho[idx,:,:])*v_prime
    del v
    
    TKE_in_z = numpy.empty(rho.shape)
    for idx in range(rho.shape[0]):
        w_prime = w[idx,:,:] - w_mean_yz[idx]
        TKE_in_z[idx,:,:] = numpy.sqrt(rho[idx,:,:])*w_prime
    del w
    
    del rho
    
    # Perform FFT.
    
    N_max = numpy.ceil(numpy.sqrt((Ny_full/2)**2 + (Nz_full/2)**2))
    N_max = N_max.astype(numpy.int32)
    
    TKE_x_spectrum_sum_local = numpy.zeros([N_max], dtype=numpy.float64)
    TKE_y_spectrum_sum_local = numpy.zeros([N_max], dtype=numpy.float64)
    TKE_z_spectrum_sum_local = numpy.zeros([N_max], dtype=numpy.float64)
    
    # FFT for TKE in x-direction.
    
    TKE_z = tw_z.transposeToPencil(TKE_in_x)
    TKE_fft_z = numpy.fft.fft(TKE_z, axis=2)
    TKE_fft_re_z = numpy.asfortranarray(numpy.real(TKE_fft_z))
    TKE_fft_im_z = numpy.asfortranarray(numpy.imag(TKE_fft_z))
    
    del TKE_fft_z
    
    TKE_fft_re = tw_z.transposeFromPencil(TKE_fft_re_z)
    TKE_fft_im = tw_z.transposeFromPencil(TKE_fft_im_z)
    
    del TKE_fft_re_z, TKE_fft_im_z
    
    TKE_fft_re_y = tw_y.transposeToPencil(TKE_fft_re)
    TKE_fft_im_y = tw_y.transposeToPencil(TKE_fft_im)
    
    TKE_fft_y = TKE_fft_re_y + 1j*TKE_fft_im_y
    
    del TKE_fft_re_y, TKE_fft_im_y
    
    TKE_fft_y = numpy.fft.fft(TKE_fft_y, axis=1)
    TKE_fft_re_y = numpy.asfortranarray(numpy.real(TKE_fft_y))
    TKE_fft_im_y = numpy.asfortranarray(numpy.imag(TKE_fft_y))
    
    del TKE_fft_y
    
    TKE_fft_re = tw_y.transposeFromPencil(TKE_fft_re_y)
    TKE_fft_im = tw_y.transposeFromPencil(TKE_fft_im_y)
    
    del TKE_fft_re_y, TKE_fft_im_y
    
    TKE_in_x_fft = TKE_fft_re + 1j*TKE_fft_im
    
    del TKE_fft_re, TKE_fft_im
    
    # FFT for TKE in y-direction.
    
    TKE_z = tw_z.transposeToPencil(TKE_in_y)
    TKE_fft_z = numpy.fft.fft(TKE_z, axis=2)
    TKE_fft_re_z = numpy.asfortranarray(numpy.real(TKE_fft_z))
    TKE_fft_im_z = numpy.asfortranarray(numpy.imag(TKE_fft_z))
    
    del TKE_fft_z
    
    TKE_fft_re = tw_z.transposeFromPencil(TKE_fft_re_z)
    TKE_fft_im = tw_z.transposeFromPencil(TKE_fft_im_z)
    
    del TKE_fft_re_z, TKE_fft_im_z
    
    TKE_fft_re_y = tw_y.transposeToPencil(TKE_fft_re)
    TKE_fft_im_y = tw_y.transposeToPencil(TKE_fft_im)
    
    TKE_fft_y = TKE_fft_re_y + 1j*TKE_fft_im_y
    
    del TKE_fft_re_y, TKE_fft_im_y
    
    TKE_fft_y = numpy.fft.fft(TKE_fft_y, axis=1)
    TKE_fft_re_y = numpy.asfortranarray(numpy.real(TKE_fft_y))
    TKE_fft_im_y = numpy.asfortranarray(numpy.imag(TKE_fft_y))
    
    del TKE_fft_y
    
    TKE_fft_re = tw_y.transposeFromPencil(TKE_fft_re_y)
    TKE_fft_im = tw_y.transposeFromPencil(TKE_fft_im_y)
    
    del TKE_fft_re_y, TKE_fft_im_y
    
    TKE_in_y_fft = TKE_fft_re + 1j*TKE_fft_im
    
    del TKE_fft_re, TKE_fft_im
    
    # FFT for TKE in z-direction.
    
    TKE_z = tw_z.transposeToPencil(TKE_in_z)
    TKE_fft_z = numpy.fft.fft(TKE_z, axis=2)
    TKE_fft_re_z = numpy.asfortranarray(numpy.real(TKE_fft_z))
    TKE_fft_im_z = numpy.asfortranarray(numpy.imag(TKE_fft_z))
    
    del TKE_fft_z
    
    TKE_fft_re = tw_z.transposeFromPencil(TKE_fft_re_z)
    TKE_fft_im = tw_z.transposeFromPencil(TKE_fft_im_z)
    
    del TKE_fft_re_z, TKE_fft_im_z
    
    TKE_fft_re_y = tw_y.transposeToPencil(TKE_fft_re)
    TKE_fft_im_y = tw_y.transposeToPencil(TKE_fft_im)
    
    TKE_fft_y = TKE_fft_re_y + 1j*TKE_fft_im_y
    
    del TKE_fft_re_y, TKE_fft_im_y
    
    TKE_fft_y = numpy.fft.fft(TKE_fft_y, axis=1)
    TKE_fft_re_y = numpy.asfortranarray(numpy.real(TKE_fft_y))
    TKE_fft_im_y = numpy.asfortranarray(numpy.imag(TKE_fft_y))
    
    del TKE_fft_y
    
    TKE_fft_re = tw_y.transposeFromPencil(TKE_fft_re_y)
    TKE_fft_im = tw_y.transposeFromPencil(TKE_fft_im_y)
    
    del TKE_fft_re_y, TKE_fft_im_y
    
    TKE_in_z_fft = TKE_fft_re + 1j*TKE_fft_im
    
    del TKE_fft_re, TKE_fft_im
    
    for indices in mixing_layer_indices:
        fft_in_x_abs_sq = numpy.absolute(TKE_in_x_fft[indices, :, :])
        fft_in_x_abs_sq = numpy.square(fft_in_x_abs_sq)
        
        fft_in_y_abs_sq = numpy.absolute(TKE_in_y_fft[indices, :, :])
        fft_in_y_abs_sq = numpy.square(fft_in_y_abs_sq)
        
        fft_in_z_abs_sq = numpy.absolute(TKE_in_z_fft[indices, :, :])
        fft_in_z_abs_sq = numpy.square(fft_in_z_abs_sq)
        
        for k in range(max(0, lo_c[2]), min(Nz_full/2, hi_c[2] + 1)):
            for j in range(max(0, lo_c[1]), min(Ny_full/2, hi_c[1] + 1)):
                idx = numpy.round(numpy.sqrt(j**2 + k**2))
                idx = idx.astype(numpy.int32)
                TKE_x_spectrum_sum_local[idx] = TKE_x_spectrum_sum_local[idx] + fft_in_x_abs_sq[j - lo_c[1], k - lo_c[2]]
                TKE_y_spectrum_sum_local[idx] = TKE_y_spectrum_sum_local[idx] + fft_in_y_abs_sq[j - lo_c[1], k - lo_c[2]]
                TKE_z_spectrum_sum_local[idx] = TKE_z_spectrum_sum_local[idx] + fft_in_z_abs_sq[j - lo_c[1], k - lo_c[2]]
    
    del TKE_in_x_fft
    del TKE_in_y_fft
    del TKE_in_z_fft
    
    # Communicate.
    
    TKE_x_spectrum_sum_global = numpy.empty(TKE_x_spectrum_sum_local.shape, dtype=numpy.float64)
    comm.Reduce([TKE_x_spectrum_sum_local, MPI.DOUBLE], [TKE_x_spectrum_sum_global, MPI.DOUBLE], \
                op=MPI.SUM, root=0)
    
    TKE_y_spectrum_sum_global = numpy.empty(TKE_y_spectrum_sum_local.shape, dtype=numpy.float64)
    comm.Reduce([TKE_y_spectrum_sum_local, MPI.DOUBLE], [TKE_y_spectrum_sum_global, MPI.DOUBLE], \
                op=MPI.SUM, root=0)
    
    TKE_z_spectrum_sum_global = numpy.empty(TKE_z_spectrum_sum_local.shape, dtype=numpy.float64)
    comm.Reduce([TKE_z_spectrum_sum_local, MPI.DOUBLE], [TKE_z_spectrum_sum_global, MPI.DOUBLE], \
                op=MPI.SUM, root=0)
    
    within_mixing_layer_global = numpy.empty(within_mixing_layer_local.shape, dtype=numpy.bool)
    comm.Reduce([within_mixing_layer_local, MPI.BOOL], [within_mixing_layer_global, MPI.BOOL], \
                op=MPI.LOR, root=0)
    
    if comm.rank == 0:
        num_mixing_layer_indices = numpy.count_nonzero(within_mixing_layer_global == True)
        
        print(num_mixing_layer_indices)
        
        dk_y = 2.0*numpy.pi/Ly_full
        dk_z = 2.0*numpy.pi/Lz_full
        dk = 2.0*numpy.pi/Ly_full
        k = dk*numpy.arange(N_max)
        
        TKE_x_spectrum_global = (dk_y*dk_z/dk)*(dy*dy*dz*dz/((2.0*numpy.pi*Ly_full)*(2.0*numpy.pi*Lz_full))/num_mixing_layer_indices)* \
            TKE_x_spectrum_sum_global
        
        TKE_y_spectrum_global = (dk_y*dk_z/dk)*(dy*dy*dz*dz/((2.0*numpy.pi*Ly_full)*(2.0*numpy.pi*Lz_full))/num_mixing_layer_indices)* \
            TKE_y_spectrum_sum_global
        
        TKE_z_spectrum_global = (dk_y*dk_z/dk)*(dy*dy*dz*dz/((2.0*numpy.pi*Ly_full)*(2.0*numpy.pi*Lz_full))/num_mixing_layer_indices)* \
            TKE_z_spectrum_sum_global
        
        numpy.save(output_filenames_k[n], k)
        numpy.save(output_filenames_E_x[n], TKE_x_spectrum_global)
        numpy.save(output_filenames_E_y[n], TKE_y_spectrum_global)
        numpy.save(output_filenames_E_z[n], TKE_z_spectrum_global)

