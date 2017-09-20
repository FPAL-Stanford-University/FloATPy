import os
import subprocess


def BuildFortranObjects(sources, compiler='gfortran', compiler_id='GNU'):

    if compiler_id == 'GNU':
        compiler_flags = '-Wall -Wconversion -Wextra -Waliasing -ffree-form -ffree-line-length-none -ffast-math -march=native -funroll-loops -fno-protect-parens'
    elif compiler_id == 'Intel':
        compiler_flags = '-mkl -warn all -xhost -dynamic -qopt-report=2 -qopt-report-phase=vec -qoverride-limits'
    else:
        print "Unknown compiler_id %s. Only 'GNU' and 'Intel' supported." %compiler_id
    
    objects = []
    
    for source in sources:
        path_dir, name = source.rsplit(os.path.sep, 1)
        
        path_object = os.path.join(
            path_dir,
            os.path.splitext(name)[0] + '.o')
        
        objects.append(os.path.relpath(path_object))
        
        command_compile_fortran_mod = (
            compiler + ' -O3 -fPIC ' + compiler_flags + ' -J' + path_dir + ' '
            + source + ' -c -o ' + path_object)
        
        print(command_compile_fortran_mod)
        
        code = subprocess.check_output(command_compile_fortran_mod, shell=True)
    
    return objects
