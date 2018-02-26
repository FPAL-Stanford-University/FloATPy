import os
import subprocess

compiler_flags = '-Wall -Wconversion -Wextra -Waliasing -ffree-form -ffree-line-length-none -ffast-math -march=native -funroll-loops -fno-protect-parens'

def BuildFortranObjects(sources, compiler='gfortran'):
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
