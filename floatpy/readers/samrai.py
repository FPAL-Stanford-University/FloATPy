"""
Module for reading and handling samrai data.
"""

import copy
import h5py
import numpy

class samraiDataReader:
    """
    Class to read samrai data.
    """
    
    def __init__(self):
        """
        Constructor of the class.
        """
        
        self.__data_directory_path = ""
        
        self.__basic_info = {}
        self.__summary_loaded = False
        
        self.__x_coords = []
        self.__y_coords = []
        self.__z_coords = []
        self.__data = {}
        self.__data_loaded = False
    
    
    def setDataDirectoryPath(self, data_directory_path):
        """
        Set the absolute path of the data directory.
        """
        
        self.__data_directory_path = data_directory_path
    
    
    def getDataDirectoryPath(self):
        """
        Get the stored absolute path of the data directory.
        """
        
        return self.__data_directory_path
    
    
    def readSummary(self):
        """
        Get the basic information, patch extents and path map from the summary file.
        """
        
        # Check whether the data directory path is already set.
        
        if not self.__data_directory_path:
            raise RuntimeError('The data directory path is not set yet!')
        
        # Open the summary file.
        
        summary_file_path = self.__data_directory_path + '/summary.samrai'
        f_summary = h5py.File(summary_file_path)
        
        # Get the basic information.
        
        basic_info = f_summary['BASIC_INFO']
        
        # Get the number of file clusters.
        self.__basic_info['num_file_clusters'] = basic_info['number_file_clusters'].value[0]
        
        # Get the time and dimension of the data.
        
        self.__basic_info['t'] = basic_info['time'].value[0]
        self.__basic_info['t_dump'] = basic_info['time_of_dump'].value[0]
        self.__basic_info['n'] = basic_info['time_step_number'].value[0]
        self.__basic_info['dim'] = basic_info['number_dimensions_of_problem'].value[0]
        
        # Get and check the grid type.
        
        self.__basic_info['grid_type'] = basic_info['grid_type'].value[0]
        
        if numpy.char.strip(self.__basic_info['grid_type']) != 'CARTESIAN':
            raise RuntimeError("Grid type other than 'CARTESIAN' not supported!")
        
        # Get the number of levels and number of patches at different levels.
        
        self.__basic_info['num_levels'] = basic_info['number_levels'].value[0]
        self.__basic_info['num_patches'] = basic_info['number_patches_at_level'].value
        self.__basic_info['num_global_patches'] = basic_info['number_global_patches'].value[0]
        
        # Get the ratios to coarser levels at different levels.
        
        self.__basic_info['ratios_to_coarser_levels'] = basic_info['ratios_to_coarser_levels'].value
        
        # Get the variable names, number of variables and number of components in each variable.
        
        self.__basic_info['var_names'] = basic_info['var_names'].value
        self.__basic_info['num_variables'] = basic_info['number_visit_variables'].value[0]
        self.__basic_info['num_var_components'] = basic_info['var_number_components'].value
        
        # Get the geometry and check the dimension.
        
        self.__basic_info['x_lo'] = basic_info['XLO'].value
        self.__basic_info['dx'] = basic_info['dx'].value
        
        extents = f_summary['extents']
        
        # Get the patch extents.
        
        self.__patch_extents = extents['patch_extents'].value
        
        # Geth the patch map.
        
        self.__patch_map =  extents['patch_map'].value
        
        # Set the flag for loading summary file to be true.
        
        self.__summary_loaded = True
    
    
    def getBasicInfo(self):
        """
        Return the loaded basic information.
        """
        
        if not self.__summary_loaded:
            raise RuntimeError('The summary file is not read yet!')
        
        return self.__basic_info
    
    
    def getPatchExtents(self):
        """
        Return the loaded patch extents.
        """
        
        if not self.__summary_loaded:
            raise RuntimeError('The summary file is not read yet!')
        
        return self.__patch_extents
    
    
    def getPatchMap(self):
        """
        Return the loaded patch map.
        """
        
        if not self.__summary_loaded:
            raise RuntimeError('The summary file is not read yet!')
        
        return self.__patch_map
    
    
    def getDataCoordinates(self):
        """
        Return the coordinates of the loaded data.
        """
        
        if not self.__data_loaded:
            raise RuntimeError('No data is read yet!')
        
        return self.__x_coords, self.__y_coords, self.__z_coords
    
    
    def getData(self, var_name):
        """
        Return the loaded data.
        """
        
        if not self.__data_loaded:
            raise RuntimeError('No data is read yet!')
        
        return self.__data[var_name]
    
    
    def clearData(self):
        """
        Clear any loaded data.
        """
        
        self.__x_coords = []
        self.__y_coords = []
        self.__z_coords = []
        self.__data.clear()
        self.__data_loaded = False
    
    
    def clear(self):
        """
        Clear all data in the class.
        """
        
        self.__data_directory_path = ""
        
        self.__basic_info.clear()
        self.__summary_loaded = False
        
        self.__x_coords = []
        self.__y_coords = []
        self.__z_coords = []
        self.__data.clear()
        self.__data_loaded = False
    
    
    def getDomainSizeAtOneLevel(self, \
            level_num):
        """
        Get the domain size at one particular level.
        """
        
        # Read the summary file if it is not yet read.
        
        if not self.__summary_loaded:
            self.readSummary()
        
        dim = self.__basic_info['dim']
        num_levels = self.__basic_info['num_levels']
        num_patches = self.__basic_info['num_patches']
        
        # Check whether the required level is valid.
        
        if level_num >= num_levels:
            raise RuntimeError('Level number is greater than number of levels!')
        elif level_num < 0:
            raise RuntimeError('Level number is negative!')
        
        # Get the domain shape at the level.
        
        num_patches_level = num_patches[level_num]
        
        patch_level_start_idx = 0
        for level_idx in range(0, level_num):
            patch_level_start_idx = patch_level_start_idx + num_patches[level_idx]
        
        lo_level = self.__patch_extents[patch_level_start_idx][0]
        up_level = self.__patch_extents[patch_level_start_idx][1]
        
        for global_patch_idx in range(patch_level_start_idx + 1, patch_level_start_idx + num_patches_level):
            lo_level = numpy.minimum(lo_level, self.__patch_extents[global_patch_idx][0])
            up_level = numpy.maximum(up_level, self.__patch_extents[global_patch_idx][1])
        
        domain_shape = up_level[0:dim] - lo_level[0:dim] + numpy.ones(dim, dtype = numpy.int)
        
        return domain_shape
    
    
    def getRefinedDomainSize(self):
        """
        Get the full domain size refined to the highest level.
        """
        
        # Read the summary file if it is not yet read.
        
        if not self.__summary_loaded:
            self.readSummary()
        
        dim = self.__basic_info['dim']
        num_levels = self.__basic_info['num_levels']
        num_patches = self.__basic_info['num_patches']
        
        # Get the ratio from the coarest level to the finest level.
        ratio_of_coarest_to_finest = 1
        ratios_to_coarser_levels = self.__basic_info['ratios_to_coarser_levels']
        
        for level_idx in range(1, num_levels):
            ratio_of_coarest_to_finest = numpy.multiply(ratios_to_coarser_levels[level_idx], ratio_of_coarest_to_finest)
        
        # Get the refined domain shape at the root level.
        
        num_patches_root_level = num_patches[0]
        
        lo_root_level = self.__patch_extents[0][0]
        up_root_level = self.__patch_extents[0][1]
        
        for patch_idx in range(1, num_patches_root_level):
            lo_root_level = numpy.minimum(lo_root_level, self.__patch_extents[patch_idx][0])
            up_root_level = numpy.maximum(up_root_level, self.__patch_extents[patch_idx][1])
        
        domain_shape = up_root_level[0:dim] - lo_root_level[0:dim] + numpy.ones(dim, dtype = numpy.int)
        domain_shape = numpy.multiply(domain_shape, ratio_of_coarest_to_finest[0:dim])
        
        return domain_shape
    
    
    def readDataAtOneLevel(self, \
            var_names,
            level_num):
        """
        Read data at one particular level.
        """
        
        # Read the summary file if it is not yet read.
        
        if not self.__summary_loaded:
            self.readSummary()
        
        # Get the number of file clusters.
        
        num_file_clusters = self.__basic_info['num_file_clusters']
        
        # Get the dimension of the problem, number of levels and number of patches.
        
        dim = self.__basic_info['dim']
        num_levels = self.__basic_info['num_levels']
        num_patches = self.__basic_info['num_patches']
        
        # Check whether the required level is valid.
        
        if level_num >= num_levels:
            raise RuntimeError('Level number is greater than number of levels!')
        elif level_num < 0:
            raise RuntimeError('Level number is negative!')
        
        # Get the variable names.
        
        var_num_components = {}
        var_component_names = {}
        
        for var_name in var_names:
            var_idx = numpy.where(self.__basic_info['var_names'] == var_name)[0][0]
            var_num_components[var_name] = self.__basic_info['num_var_components'][var_idx]
            
            var_component_names[var_name] = [None]*var_num_components[var_name]
            if var_num_components[var_name] == 1:
                var_component_names[var_name] = [var_name]
            else:
                for component_idx in range(0, var_num_components[var_name]):
                    var_component_names[var_name][component_idx] = var_name + '.' + str(component_idx).zfill(2)
        
        # Get the domain shape.
        
        domain_shape = self.getDomainSizeAtOneLevel(level_num)
        
        # Get the index of the lower corner and the coordinates of the level.
        
        num_patches_level = num_patches[level_num]
        dx = self.__basic_info['dx']
        
        patch_level_start_idx = 0
        for level_idx in range(0, level_num):
            patch_level_start_idx = patch_level_start_idx + num_patches[level_idx]
        
        lo_level = self.__patch_extents[patch_level_start_idx][0]
        x_lo_level = self.__patch_extents[patch_level_start_idx][2]
        x_up_level = self.__patch_extents[patch_level_start_idx][3]
        
        for global_patch_idx in range(patch_level_start_idx + 1, patch_level_start_idx + num_patches_level):
            lo_level = numpy.minimum(lo_level, self.__patch_extents[global_patch_idx][0])
            x_lo_level = numpy.minimum(x_lo_level, self.__patch_extents[global_patch_idx][2])
            x_up_level = numpy.maximum(x_up_level, self.__patch_extents[global_patch_idx][3])
        
        if dim == 1:
            self.__x_coords = numpy.linspace(x_lo_level[0] + 0.5*dx[level_num][0], \
                x_up_level[0] - 0.5*dx[level_num][0], \
                num = domain_shape[0])
        
        elif dim == 2:
            self.__x_coords = numpy.linspace(x_lo_level[0] + 0.5*dx[level_num][0], \
                x_up_level[0] - 0.5*dx[level_num][0], \
                num = domain_shape[0])
            self.__y_coords = numpy.linspace(x_lo_level[1] + 0.5*dx[level_num][1], \
                x_up_level[1] - 0.5*dx[level_num][1], \
                num = domain_shape[1])
        
        elif dim == 3:
            self.__x_coords = numpy.linspace(x_lo_level[0] + 0.5*dx[level_num][0], \
                x_up_level[0] - 0.5*dx[level_num][0], \
                num = domain_shape[0])
            self.__y_coords = numpy.linspace(x_lo_level[1] + 0.5*dx[level_num][1], \
                x_up_level[1] - 0.5*dx[level_num][1], \
                num = domain_shape[1])
            self.__z_coords = numpy.linspace(x_lo_level[2] + 0.5*dx[level_num][2], \
                x_up_level[2] - 0.5*dx[level_num][2], \
                num = domain_shape[2])
        
        elif dim < 1 or dim > 3:
            raise RuntimeError('Problem dimension < 1 or > 3 not supported!')
               
        # Initialize container to store the data. The elements in the container are initialized as NAN values.
        
        for var_name in var_names:
            data_shape = numpy.insert(domain_shape, 0, var_num_components[var_name])
            self.__data[var_name] = numpy.empty(data_shape)
            self.__data[var_name][:] = numpy.NAN
        
        # Get the data from all patches at the specified level.
        
        if dim == 1:
            for process_idx in range(0, num_file_clusters):
                file_name = 'processor_cluster.' + str(process_idx).zfill(5) + '.samrai'
                full_path = self.__data_directory_path + '/' + file_name
                f_input = h5py.File(full_path, 'r')
                
                file_cluster = f_input['processor.' + str(process_idx).zfill(5)]
                file_cluster_level = file_cluster['level.' + str(level_num).zfill(5)]
                
                for var_name in var_names:
                    for patch_key in file_cluster_level:
                        patch_idx = int(patch_key.replace('patch.', ''))
                        global_patch_idx = patch_level_start_idx + patch_idx
                        file_cluster_patch = file_cluster_level[patch_key]
                        
                        lo_patch = patch_extents[global_patch_idx][0] - lo_level
                        up_patch = patch_extents[global_patch_idx][1] - lo_level
                        patch_shape = up_patch[0] - lo_patch[0] + 1
                        
                        x_start_idx = lo_patch[0]
                        x_end_idx = up_patch[0] + 1
                        
                        for component_idx in range(0, var_num_components[var_name]):
                            self.__data[var_name][component_idx, x_start_idx:x_end_idx] = \
                                file_cluster_patch[var_component_names[var_name][component_idx]].value.reshape(patch_shape, order = 'F')
                
                f_input.close()
        
        elif dim == 2:
            for process_idx in range(0, num_file_clusters):
                file_name = 'processor_cluster.' + str(process_idx).zfill(5) + '.samrai'
                full_path = self.__data_directory_path + '/' + file_name
                f_input = h5py.File(full_path, 'r')
                
                file_cluster = f_input['processor.' + str(process_idx).zfill(5)]
                file_cluster_level = file_cluster['level.' + str(level_num).zfill(5)]
                
                for var_name in var_names:
                    for patch_key in file_cluster_level:
                        patch_idx = int(patch_key.replace('patch.', ''))
                        global_patch_idx = patch_level_start_idx + patch_idx
                        file_cluster_patch = file_cluster_level[patch_key]
                        
                        lo_patch = self.__patch_extents[global_patch_idx][0] - lo_level
                        up_patch = self.__patch_extents[global_patch_idx][1] - lo_level
                        patch_shape = up_patch[0:2] - lo_patch[0:2] + numpy.ones(2, dtype = numpy.int)
                        
                        x_start_idx = lo_patch[0]
                        x_end_idx = up_patch[0] + 1
                        
                        y_start_idx = lo_patch[1]
                        y_end_idx = up_patch[1] + 1
                        
                        for component_idx in range(0, var_num_components[var_name]):
                            self.__data[var_name][component_idx, x_start_idx:x_end_idx, y_start_idx:y_end_idx] = \
                                file_cluster_patch[var_component_names[var_name][component_idx]].value.reshape(patch_shape, order = 'F')
                
                f_input.close()
        
        elif dim == 3:
            for process_idx in range(0, num_file_clusters):
                file_name = 'processor_cluster.' + str(process_idx).zfill(5) + '.samrai'
                full_path = self.__data_directory_path + '/' + file_name
                f_input = h5py.File(full_path, 'r')
                
                file_cluster = f_input['processor.' + str(process_idx).zfill(5)]
                file_cluster_level = file_cluster['level.' + str(level_num).zfill(5)]
                
                for var_name in var_names:
                    for patch_key in file_cluster_level:
                        patch_idx = int(patch_key.replace('patch.', ''))
                        global_patch_idx = patch_level_start_idx + patch_idx
                        file_cluster_patch = file_cluster_level[patch_key]
                        
                        lo_patch = self.__patch_extents[global_patch_idx][0] - lo_level
                        up_patch = self.__patch_extents[global_patch_idx][1] - lo_level
                        patch_shape = up_patch[0:3] - lo_patch[0:3] + numpy.ones(3, dtype = numpy.int)
                        
                        x_start_idx = lo_patch[0]
                        x_end_idx = up_patch[0] + 1
                        
                        y_start_idx = lo_patch[1]
                        y_end_idx = up_patch[1] + 1
                        
                        z_start_idx = lo_patch[2]
                        z_end_idx = up_patch[2] + 1
                        
                        for component_idx in range(0, var_num_components[var_name]):
                            self.__data[var_name][component_idx, x_start_idx:x_end_idx, y_start_idx:y_end_idx, z_start_idx:z_end_idx] = \
                                file_cluster_patch[var_component_names[var_name][component_idx]].value.reshape(patch_shape, order = 'F')
                
                f_input.close()
        
        else:
            raise RuntimeError('Problem dimension < 1 or > 3 not supported!')
        
        self.__data_loaded = True
    
    
    def readCombinedDataFromAllLevels(self, \
            var_names, \
            num_ghosts, \
            periodic_dimension):
        """
        Read data from all levels.
        """
        
        # Read the summary file if it is not yet read.
        
        if not self.__summary_loaded:
            self.readSummary()
        
        # Get the number of file clusters.
        
        num_file_clusters = self.__basic_info['num_file_clusters']
        
        # Get the dimension of the problem, number of levels and number of patches.
        
        dim = self.__basic_info['dim']
        num_levels = self.__basic_info['num_levels']
        num_patches = self.__basic_info['num_patches']
        
        # Check that cooresponding element of num_ghosts is greater than zero if the dimension is periodic.
        
        if (periodic_dimension[0] == True and num_ghosts[0] < 1):
            raise RuntimeError('The periodic dimension is not consistent with the number of ghost cells in the first direction!')
        
        if dim > 1:
            if (periodic_dimension[1] == True and num_ghosts[1] < 1):
                raise RuntimeError('The periodic dimension is not consistent with the number of ghost cells in the second direction!')
        
        if dim > 2:
            if (periodic_dimension[2] == True and num_ghosts[2] < 1):
                raise RuntimeError('The periodic dimension is not consistent with the number of ghost cells in the third direction!')
        
        # Get the variable names.
        
        var_num_components = {}
        var_component_names = {}
        
        for var_name in var_names:
            var_idx = numpy.where(self.__basic_info['var_names'] == var_name)[0][0]
            var_num_components[var_name] = self.__basic_info['num_var_components'][var_idx]

            var_component_names[var_name] = [None]*var_num_components[var_name]
            if var_num_components[var_name] == 1:
                var_component_names[var_name] = [var_name]
            else:
                for component_idx in range(0, var_num_components[var_name]):
                    var_component_names[var_name][component_idx] = var_name + '.' + str(component_idx).zfill(2)
                
        # Get the ratios of different levels to the finest level.
        
        ratios_to_coarser_levels = self.__basic_info['ratios_to_coarser_levels']
        ratios_to_finest_level = numpy.empty(ratios_to_coarser_levels.shape, dtype = ratios_to_coarser_levels.dtype)
        ratios_to_finest_level[num_levels - 1] = -1*ratios_to_coarser_levels[0]
        for level_idx in range(num_levels - 2, -1, -1):
            ratios_to_finest_level[level_idx] = numpy.multiply(ratios_to_coarser_levels[level_idx + 1], ratios_to_finest_level[level_idx + 1])
        
        # Get the domain shape refined to the finest level.
        
        domain_shape = self.getRefinedDomainSize()
        domain_shape = domain_shape + 2*num_ghosts[0:dim]
        
        # Get the coordinates of the domain refiend to the finest level.
        
        num_patches_root_level = num_patches[0]
        dx = self.__basic_info['dx']
        
        x_lo_root_level = self.__patch_extents[0][2]
        x_up_root_level = self.__patch_extents[0][3]
        
        for patch_idx in range(1, num_patches_root_level):
            x_lo_root_level = numpy.minimum(x_lo_root_level, self.__patch_extents[patch_idx][2])
            x_up_root_level = numpy.maximum(x_up_root_level, self.__patch_extents[patch_idx][3])
        
        if dim == 1:
            self.__x_coords = numpy.linspace(x_lo_root_level[0] + (0.5 - num_ghosts[0])*dx[-1][0], \
                x_up_root_level[0] + (num_ghosts[0] - 0.5)*dx[-1][0], \
                num = domain_shape[0])
        
        elif dim == 2:
            self.__x_coords = numpy.linspace(x_lo_root_level[0] + (0.5 - num_ghosts[0])*dx[-1][0], \
                x_up_root_level[0] + (num_ghosts[0] - 0.5)*dx[-1][0], \
                num = domain_shape[0])
            self.__y_coords = numpy.linspace(x_lo_root_level[1] + (0.5 - num_ghosts[1])*dx[-1][1], \
                x_up_root_level[1] + (num_ghosts[1] - 0.5)*dx[-1][1], \
                num = domain_shape[1])
        
        elif dim == 3:
            self.__x_coords = numpy.linspace(x_lo_root_level[0] + (0.5 - num_ghosts[0])*dx[-1][0], \
                x_up_root_level[0] + (num_ghosts[0] - 0.5)*dx[-1][0], \
                num = domain_shape[0])
            self.__y_coords = numpy.linspace(x_lo_root_level[1] + (0.5 - num_ghosts[1])*dx[-1][1], \
                x_up_root_level[1] + (num_ghosts[1] - 0.5)*dx[-1][1], \
                num = domain_shape[1])
            self.__z_coords = numpy.linspace(x_lo_root_level[2] + (0.5 - num_ghosts[2])*dx[-1][2], \
                x_up_root_level[2] + (num_ghosts[2] - 0.5)*dx[-1][2], \
                num = domain_shape[2])
        
        elif dim < 1 or dim > 3:
            raise RuntimeError('Problem dimension < 1 or > 3 not supported!')
        
        # Initialize containers to store the data at different levels. The elements in the containers 
        # are initialized as NAN values.
        
        level_data = {}
        
        for var_name in var_names:
            data_shape = numpy.insert(domain_shape, 0, var_num_components[var_name])
            data = numpy.empty(data_shape)
            data[:] = numpy.NAN
            
            level_data[var_name] = []
            for level_idx in range(num_levels):
                level_data[var_name].append(copy.deepcopy(data))
        
        # Get the data from all levels and upsampling them to the finest level.
        
        if dim == 1:
            for process_idx in range(0, num_file_clusters):
                file_name = 'processor_cluster.' + str(process_idx).zfill(5) + '.samrai'
                full_path = self.__data_directory_path + '/' + file_name
                f_input = h5py.File(full_path, 'r')
                
                for level_num in range(num_levels):
                    patch_level_start_idx = 0
                    for level_idx in range(0, level_num):
                        patch_level_start_idx = patch_level_start_idx + num_patches[level_idx]
                    
                    file_cluster = f_input['processor.' + str(process_idx).zfill(5)]
                    file_cluster_level = file_cluster['level.' + str(level_num).zfill(5)]
                    
                    for var_name in var_names:
                        for patch_key in file_cluster_level:
                            patch_idx = int(patch_key.replace('patch.', ''))
                            global_patch_idx = patch_level_start_idx + patch_idx
                            file_cluster_patch = file_cluster_level[patch_key]
                            
                            lo_patch = self.__patch_extents[global_patch_idx][0]
                            up_patch = self.__patch_extents[global_patch_idx][1]
                            patch_shape = up_patch[0] - lo_patch[0] + 1
                            
                            x_start_idx = lo_patch[0]*ratios_to_finest_level[level_num][0] + num_ghosts[0]
                            x_end_idx = (up_patch[0] + 1)*ratios_to_finest_level[level_num][0] + num_ghosts[0]
                            
                            for component_idx in range(0, var_num_components[var_name]):
                                patch_data = file_cluster_patch[var_component_names[var_name][component_idx]].value.reshape( \
                                    patch_shape, order = 'F')
                                patch_data = numpy.repeat(patch_data, ratios_to_finest_level[level_num][0], axis = 0)
                                
                                level_data[var_name][level_num][component_idx, x_start_idx: x_end_idx] = patch_data
                
                f_input.close()
        
        elif dim == 2:
            for process_idx in range(0, num_file_clusters):
                file_name = 'processor_cluster.' + str(process_idx).zfill(5) + '.samrai'
                full_path = self.__data_directory_path + '/' + file_name
                f_input = h5py.File(full_path, 'r')
                
                for level_num in range(num_levels):
                    patch_level_start_idx = 0
                    for level_idx in range(0, level_num):
                        patch_level_start_idx = patch_level_start_idx + num_patches[level_idx]
                    
                    file_cluster = f_input['processor.' + str(process_idx).zfill(5)]
                    file_cluster_level = file_cluster['level.' + str(level_num).zfill(5)]
                    
                    for var_name in var_names:
                        for patch_key in file_cluster_level:
                            patch_idx = int(patch_key.replace('patch.', ''))
                            global_patch_idx = patch_level_start_idx + patch_idx
                            file_cluster_patch = file_cluster_level[patch_key]
                            
                            lo_patch = self.__patch_extents[global_patch_idx][0]
                            up_patch = self.__patch_extents[global_patch_idx][1]
                            patch_shape = up_patch[0:2] - lo_patch[0:2] + numpy.ones(2, dtype = numpy.int)
                            
                            x_start_idx = lo_patch[0]*ratios_to_finest_level[level_num][0] + num_ghosts[0]
                            x_end_idx = (up_patch[0] + 1)*ratios_to_finest_level[level_num][0] + num_ghosts[0]
                            
                            y_start_idx = lo_patch[1]*ratios_to_finest_level[level_num][1] + num_ghosts[1]
                            y_end_idx = (up_patch[1] + 1)*ratios_to_finest_level[level_num][1] + num_ghosts[1]
                            
                            for component_idx in range(0, var_num_components[var_name]):
                                patch_data = file_cluster_patch[var_component_names[var_name][component_idx]].value.reshape( \
                                    patch_shape, order = 'F')
                                patch_data = numpy.repeat(patch_data, ratios_to_finest_level[level_num][0], axis = 0)
                                patch_data = numpy.repeat(patch_data, ratios_to_finest_level[level_num][1], axis = 1)
                                
                                level_data[var_name][level_num][component_idx, \
                                    x_start_idx:x_end_idx, y_start_idx:y_end_idx] = \
                                        patch_data
                
                f_input.close()
        
        elif dim == 3:
            for process_idx in range(0, num_file_clusters):
                file_name = 'processor_cluster.' + str(process_idx).zfill(5) + '.samrai'
                full_path = self.__data_directory_path + '/' + file_name
                f_input = h5py.File(full_path, 'r')
                
                for level_num in range(num_levels):
                    patch_level_start_idx = 0
                    for level_idx in range(0, level_num):
                        patch_level_start_idx = patch_level_start_idx + num_patches[level_idx]
                    
                    file_cluster = f_input['processor.' + str(process_idx).zfill(5)]
                    file_cluster_level = file_cluster['level.' + str(level_num).zfill(5)]
                    
                    for var_name in var_names:
                        for patch_key in file_cluster_level:
                            patch_idx = int(patch_key.replace('patch.', ''))
                            global_patch_idx = patch_level_start_idx + patch_idx
                            file_cluster_patch = file_cluster_level[patch_key]
                            
                            lo_patch = self.__patch_extents[global_patch_idx][0]
                            up_patch = self.__patch_extents[global_patch_idx][1]
                            patch_shape = up_patch[0:3] - lo_patch[0:3] + numpy.ones(3, dtype = numpy.int)
                            
                            x_start_idx = lo_patch[0]*ratios_to_finest_level[level_num][0] + num_ghosts[0]
                            x_end_idx = (up_patch[0] + 1)*ratios_to_finest_level[level_num][0] + num_ghosts[0]
                            
                            y_start_idx = lo_patch[1]*ratios_to_finest_level[level_num][1] + num_ghosts[1]
                            y_end_idx = (up_patch[1] + 1)*ratios_to_finest_level[level_num][1] + num_ghosts[1]
                            
                            z_start_idx = lo_patch[2]*ratios_to_finest_level[level_num][2] + num_ghosts[2]
                            z_end_idx = (up_patch[2] + 1)*ratios_to_finest_level[level_num][2] + num_ghosts[2]
                            
                            for component_idx in range(0, var_num_components[var_name]):
                                patch_data = file_cluster_patch[var_component_names[var_name][component_idx]].value.reshape( \
                                    patch_shape, order = 'F')
                                patch_data = numpy.repeat(patch_data, ratios_to_finest_level[level_num][0], axis = 0)
                                patch_data = numpy.repeat(patch_data, ratios_to_finest_level[level_num][1], axis = 1)
                                patch_data = numpy.repeat(patch_data, ratios_to_finest_level[level_num][2], axis = 2)
                                
                                level_data[var_name][level_num][component_idx, \
                                    x_start_idx:x_end_idx, y_start_idx:y_end_idx, z_start_idx:z_end_idx] = \
                                        patch_data
            
            f_input.close()
        
        else:
            raise RuntimeError('Problem dimension < 1 or > 3 not supported!')
        
        # Combine data at all levels.
        
        for var_name in var_names:
            self.__data[var_name] = level_data[var_name][0]
            
            for level_idx in range(1, num_levels):
                is_finite_idx = numpy.isfinite(level_data[var_name][level_idx])
                self.__data[var_name][is_finite_idx] = level_data[var_name][level_idx][is_finite_idx]
            
            # Apply the periodic boundary conditions on the ghost cells.
            
            if dim == 1:
                for component_idx in range(0, var_num_components[var_name]):
                    if periodic_dimension[0] == True:
                        self.__data[var_name][component_idx, 0:num_ghosts[0]] = \
                            self.__data[var_name][component_idx, -2*num_ghosts[0]:-num_ghosts[0]]
                        
                        self.__data[var_name][component_idx, -num_ghosts[0]:] = \
                            self.__data[var_name][component_idx, num_ghosts[0]:2*num_ghosts[0]]
            
            elif dim == 2:
                for component_idx in range(0, var_num_components[var_name]):
                    if periodic_dimension[0] == True:
                        self.__data[var_name][component_idx, 0:num_ghosts[0], :] = \
                            self.__data[var_name][component_idx, -2*num_ghosts[0]:-num_ghosts[0], :]
                        
                        self.__data[var_name][component_idx, -num_ghosts[0]:, :] = \
                            self.__data[var_name][component_idx, num_ghosts[0]:2*num_ghosts[0], :]
                    
                    if periodic_dimension[1] == True:
                        self.__data[var_name][component_idx, :, 0:num_ghosts[1]] = \
                            self.__data[var_name][component_idx, :, -2*num_ghosts[1]:-num_ghosts[1]]
                        
                        self.__data[var_name][component_idx, :, -num_ghosts[1]:] = \
                            self.__data[var_name][component_idx, :, num_ghosts[1]:2*num_ghosts[1]]
            
            elif dim == 3:
                for component_idx in range(0, var_num_components[var_name]):
                    if periodic_dimension[0] == True:
                        self.__data[var_name][component_idx, 0:num_ghosts[0], :, :] = \
                            self.__data[var_name][component_idx, -2*num_ghosts[0]:-num_ghosts[0], :, :]
                        
                        self.__data[var_name][component_idx, -num_ghosts[0]:, :, :] = \
                            self.__data[var_name][component_idx, num_ghosts[0]:2*num_ghosts[0], :, :]
                    
                    if periodic_dimension[1] == True:
                        self.__data[var_name][component_idx, :, 0:num_ghosts[1], :] = \
                            self.__data[var_name][component_idx, :, -2*num_ghosts[1]:-num_ghosts[1], :]
                        
                        self.__data[var_name][component_idx, :, -num_ghosts[1]:, :] = \
                            self.__data[var_name][component_idx, :, num_ghosts[1]:2*num_ghosts[1], :]
                    
                    if periodic_dimension[2] == True:
                        self.__data[var_name][component_idx, :, :, 0:num_ghosts[2]] = \
                            self.__data[var_name][component_idx, :, :, -2*num_ghosts[2]:-num_ghosts[2]]
                        
                        self.__data[var_name][component_idx, :, :, -num_ghosts[2]:] = \
                            self.__data[var_name][component_idx, :, :, num_ghosts[2]:2*num_ghosts[2]]
            
            else:
                raise RuntimeError('Problem dimension < 1 or > 3 not supported!')
        
        self.__data_loaded = True
    
    
    def readCombinedDataInSubdomainFromAllLevels(self, \
            var_names, \
            lo_subdomain, \
            up_subdomain, \
            num_ghosts, \
            periodic_dimension):
        """
        Read data in a sub-domain from all levels, refine the data to the finest level
        and combine the data from different levels.
        """
        
        # Read the summary file if it is not yet read.
        
        if not self.__summary_loaded:
            self.readSummary()
        
        # Check that lo_subdomain is smaller than or equal to up_subdomain.
        
        if numpy.all(numpy.less_equal(lo_subdomain, up_subdomain)) == False:
            raise RuntimeError('lo_subdomain is greater than up_subdomain!')
        
        # Get the number of file clusters.
        
        num_file_clusters = self.__basic_info['num_file_clusters']
        
        # Get the dimension of the problem, number of levels and number of patches.
        
        dim = self.__basic_info['dim']
        num_levels = self.__basic_info['num_levels']
        num_patches = self.__basic_info['num_patches']
        num_patches_root_level = num_patches[0]
        
        # Get the grid spacings at different levels.
        
        dx = self.__basic_info['dx']
        
        # Get the variable names.
        
        var_num_components = {}
        var_component_names = {}

        for var_name in var_names:
            var_idx = numpy.where(self.__basic_info['var_names'] == var_name)[0][0]
            var_num_components[var_name] = self.__basic_info['num_var_components'][var_idx]

            var_component_names[var_name] = [None]*var_num_components[var_name]
            if var_num_components[var_name] == 1:
                var_component_names[var_name] = [var_name]
            else:
                for component_idx in range(0, var_num_components[var_name]):
                    var_component_names[var_name][component_idx] = var_name + '.' + str(component_idx).zfill(2)
        
        ratios_to_coarser_levels = self.__basic_info['ratios_to_coarser_levels']
        ratios_to_finest_level = numpy.empty(ratios_to_coarser_levels.shape, dtype = ratios_to_coarser_levels.dtype)
        ratios_to_finest_level[num_levels - 1] = -1*ratios_to_coarser_levels[0]
        for level_idx in range(num_levels - 2, -1, -1):
            ratios_to_finest_level[level_idx] = numpy.multiply(ratios_to_coarser_levels[level_idx + 1], \
                                                ratios_to_finest_level[level_idx + 1])
        
        # Get the lower and upper indices of the domain.
        
        lo_root_level = self.__patch_extents[0][0]
        up_root_level = self.__patch_extents[0][1]
        x_lo_root_level = self.__patch_extents[0][2]
        x_up_root_level = self.__patch_extents[0][3]
        
        for patch_idx in range(1, num_patches_root_level):
            lo_root_level = numpy.minimum(lo_root_level, self.__patch_extents[patch_idx][0])
            up_root_level = numpy.maximum(up_root_level, self.__patch_extents[patch_idx][1])
            x_lo_root_level = numpy.minimum(x_lo_root_level, self.__patch_extents[patch_idx][2])
            x_up_root_level = numpy.maximum(x_up_root_level, self.__patch_extents[patch_idx][3])
        
        # Refine the the lower and upper indices of the domain to the highest level.
        
        lo_root_level_refined = numpy.multiply(lo_root_level[0:dim], ratios_to_finest_level[0][0:dim])
        up_root_level_refined = numpy.multiply(up_root_level[0:dim] + numpy.ones(dim, dtype = numpy.int), \
            ratios_to_finest_level[0][0:dim]) \
            - numpy.ones(dim, dtype = numpy.int)
        
        # Compute the shape of the domain refined to the highest level.
        
        domain_shape = up_root_level_refined[0:dim] - lo_root_level_refined[0:dim] + numpy.ones(dim, dtype = numpy.int)
        
        # Check whether the requested sub-domain is inside the computational domain.
        
        if numpy.all(numpy.greater_equal(lo_subdomain[0:dim], lo_root_level_refined)) == False:
            raise RuntimeError('Input sub-domain not inside the computational domain!')
        if numpy.all(numpy.less_equal(up_subdomain[0:dim], up_root_level_refined)) == False:
            raise RuntimeError('Input sub-domain not inside the computational domain!')
        
        # Include the ghost cells in the sub-domain.
        
        lo_subdomain = lo_subdomain[0:dim] - num_ghosts[0:dim]
        up_subdomain = up_subdomain[0:dim] + num_ghosts[0:dim]
        
        
        # Determine which file clusters to load.
        
        load_file_clusters = numpy.empty(num_file_clusters, dtype = numpy.bool)
        load_file_clusters.fill(False)
        
        for level_num in range(num_levels):
            patch_level_start_idx = 0
            for level_idx in range(0, level_num):
                patch_level_start_idx = patch_level_start_idx + num_patches[level_idx]
            
            for patch_idx in range(num_patches[level_num]):
                global_patch_idx = patch_level_start_idx + patch_idx
                lo_patch_refined = self.__patch_extents[global_patch_idx][0]*ratios_to_finest_level[level_num][0]
                up_patch_refined = (self.__patch_extents[global_patch_idx][1] + 1)*ratios_to_finest_level[level_num][0] - 1
                
                file_cluster_num = self.__patch_map[global_patch_idx][1]
                
                # Determine whether the sub-domain overlaps with the patches.
                
                if dim == 1:
                    if (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]):
                        load_file_clusters[file_cluster_num] = True
                    
                    # Check whether the sub-domain touches periodic boundaries.
                    # If the sub-domain touches the periodic boundaries, check whether ghost cell region
                    # overlaps with the patches.
                    
                    # Check the left boundary.
                    
                    if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True):
                        if (up_patch_refined[0] + num_ghosts[0] > up_root_level_refined[0]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the right boundary.
                    
                    if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True):
                        if (lo_patch_refined[0] - num_ghosts[0] < lo_root_level_refined[0]):
                            load_file_clusters[file_cluster_num] = True
                
                elif dim == 2:
                    if (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]) and \
                       (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]):
                        load_file_clusters[file_cluster_num] = True
                    
                    # Check whether the sub-domain touches periodic boundaries.
                    # If the sub-domain touches the periodic boundaries, check whether ghost cell region
                    # overlaps with the patches.
                    
                    # Check the left edge.
                    
                    if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True):
                        if (up_patch_refined[0] + num_ghosts[0] > up_root_level_refined[0]) and \
                           (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the right edge.
                    
                    if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True):
                        if (lo_patch_refined[0] - num_ghosts[0] < lo_root_level_refined[0]) and \
                           (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the bottom edge.
                    
                    if (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True):
                        if (up_patch_refined[1] + num_ghosts[1] > up_root_level_refined[1]) and \
                           (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the top edge.
                    
                    if (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True):
                        if (lo_patch_refined[1] - num_ghosts[1] < lo_root_level_refined[1]) and \
                           (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the left-bottom corner.
                    
                    if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True):
                        if (up_patch_refined[0] + num_ghosts[0] > up_root_level_refined[0]) and \
                           (up_patch_refined[1] + num_ghosts[1] > up_root_level_refined[1]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the left-top corner.
                    
                    if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True):
                        if (up_patch_refined[0] + num_ghosts[0] > up_root_level_refined[0]) and \
                           (lo_patch_refined[1] - num_ghosts[1] < lo_root_level_refined[1]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the right-bottom corner.
                    
                    if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True):
                        if (lo_patch_refined[0] - num_ghosts[0] < lo_root_level_refined[0]) and \
                           (up_patch_refined[1] + num_ghosts[1] > up_root_level_refined[1]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the right-top corner.
                    
                    if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True):
                        if (lo_patch_refined[0] - num_ghosts[0] < lo_root_level_refined[0]) and \
                           (lo_patch_refined[1] - num_ghosts[1] < lo_root_level_refined[1]):
                            load_file_clusters[file_cluster_num] = True
                
                elif dim == 3:
                    if (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]) and \
                       (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]) and \
                       (lo_patch_refined[2] <= up_subdomain[2]) and (up_patch_refined[2] >= lo_subdomain[2]):
                        load_file_clusters[file_cluster_num] = True
                    
                    # Check whether the sub-domain touches periodic boundaries.
                    # If the sub-domain touches the periodic boundaries, check whether ghost cell region
                    # overlaps with the patches.
                    
                    # Check the left face.
                    
                    if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]) and \
                       (lo_patch_refined[2] <= up_subdomain[2]) and (up_patch_refined[2] >= lo_subdomain[2]):
                        if (up_patch_refined[0] + num_ghosts[0] > up_root_level_refined[0]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the right face.
                    
                    if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]) and \
                       (lo_patch_refined[2] <= up_subdomain[2]) and (up_patch_refined[2] >= lo_subdomain[2]):
                        if (lo_patch_refined[0] - num_ghosts[0] < lo_root_level_refined[0]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the bottom face.
                    
                    if (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True) and \
                       (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]) and \
                       (lo_patch_refined[2] <= up_subdomain[2]) and (up_patch_refined[2] >= lo_subdomain[2]):
                        if (up_patch_refined[1] + num_ghosts[1] > up_root_level_refined[1]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the top face.
                    
                    if (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True) and \
                       (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]) and \
                       (lo_patch_refined[2] <= up_subdomain[2]) and (up_patch_refined[2] >= lo_subdomain[2]):
                        if (lo_patch_refined[1] - num_ghosts[1] < lo_root_level_refined[1]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the back face.
                    
                    if (lo_subdomain[2] == lo_root_level_refined[2] - num_ghosts[2]) and (periodic_dimension[2] == True) and \
                       (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]) and \
                       (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]):
                        if (up_patch_refined[2] + num_ghosts[2] > up_root_level_refined[1]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the front face.
                    
                    if (up_subdomain[2] == up_root_level_refined[2] + num_ghosts[2]) and (periodic_dimension[2] == True) and \
                       (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]) and \
                       (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]):
                        if (lo_patch_refined[2] - num_ghosts[2] < lo_root_level_refined[2]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the left-bottom edge.
                    
                    if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True) and \
                       (lo_patch_refined[2] <= up_subdomain[2]) and (up_patch_refined[2] >= lo_subdomain[2]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the left-top edge.
                    
                    if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True) and \
                       (lo_patch_refined[2] <= up_subdomain[2]) and (up_patch_refined[2] >= lo_subdomain[2]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the left-back edge.
                    
                    if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (lo_subdomain[2] == lo_root_level_refined[2] - num_ghosts[2]) and (periodic_dimension[2] == True) and \
                       (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the left-front edge.
                    
                    if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (up_subdomain[2] == up_root_level_refined[2] + num_ghosts[2]) and (periodic_dimension[2] == True) and \
                       (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the right-bottom edge.
                    
                    if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True) and \
                       (lo_patch_refined[2] <= up_subdomain[2]) and (up_patch_refined[2] >= lo_subdomain[2]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the right-top edge.
                    
                    if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True) and \
                       (lo_patch_refined[2] <= up_subdomain[2]) and (up_patch_refined[2] >= lo_subdomain[2]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the right-back edge.
                    
                    if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (lo_subdomain[2] == lo_root_level_refined[2] - num_ghosts[2]) and (periodic_dimension[2] == True) and \
                       (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the right-front edge.
                    
                    if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (up_subdomain[2] == up_root_level_refined[2] + num_ghosts[2]) and (periodic_dimension[2] == True) and \
                       (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the bottom-back edge.
                    
                    if (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True) and \
                       (lo_subdomain[2] == lo_root_level_refined[2] - num_ghosts[2]) and (periodic_dimension[2] == True) and \
                       (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the bottom-front edge.
                    
                    if (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True) and \
                       (up_subdomain[2] == up_root_level_refined[2] + num_ghosts[2]) and (periodic_dimension[2] == True) and \
                       (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the top-back edge.
                    
                    if (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True) and \
                       (lo_subdomain[2] == lo_root_level_refined[2] - num_ghosts[2]) and (periodic_dimension[2] == True) and \
                       (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the top-front edge.
                    
                    if (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True) and \
                       (up_subdomain[2] == up_root_level_refined[2] + num_ghosts[2]) and (periodic_dimension[2] == True) and \
                       (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the left-bottom-back corner.
                    
                    if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True) and \
                       (lo_subdomain[2] == lo_root_level_refined[2] - num_ghosts[2]) and (periodic_dimension[2] == True):
                        if (up_patch_refined[0] + num_ghosts[0] > up_root_level_refined[0]) and \
                           (up_patch_refined[1] + num_ghosts[1] > up_root_level_refined[1]) and \
                           (up_patch_refined[2] + num_ghosts[2] > up_root_level_refined[2]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the left-top-back corner.
                    
                    if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True) and \
                       (lo_subdomain[2] == lo_root_level_refined[2] - num_ghosts[2]) and (periodic_dimension[2] == True):
                        if (up_patch_refined[0] + num_ghosts[0] > up_root_level_refined[0]) and \
                           (lo_patch_refined[1] - num_ghosts[1] < lo_root_level_refined[1]) and \
                           (up_patch_refined[2] + num_ghosts[2] > up_root_level_refined[2]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the right-bottom-back corner.
                    
                    if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True) and \
                       (lo_subdomain[2] == lo_root_level_refined[2] - num_ghosts[2]) and (periodic_dimension[2] == True):
                        if (lo_patch_refined[0] - num_ghosts[0] < lo_root_level_refined[0]) and \
                           (up_patch_refined[1] + num_ghosts[1] > up_root_level_refined[1]) and \
                           (up_patch_refined[2] + num_ghosts[2] > up_root_level_refined[2]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the right-top-back corner.
                    
                    if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True) and \
                       (lo_subdomain[2] == lo_root_level_refined[2] - num_ghosts[2]) and (periodic_dimension[2] == True):
                        if (lo_patch_refined[0] - num_ghosts[0] < lo_root_level_refined[0]) and \
                           (lo_patch_refined[1] - num_ghosts[1] < lo_root_level_refined[1]) and \
                           (up_patch_refined[2] + num_ghosts[2] > up_root_level_refined[2]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the left-bottom-front corner.
                    
                    if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True) and \
                       (up_subdomain[2] == up_root_level_refined[2] + num_ghosts[2]) and (periodic_dimension[2] == True):
                        if (up_patch_refined[0] + num_ghosts[0] > up_root_level_refined[0]) and \
                           (up_patch_refined[1] + num_ghosts[1] > up_root_level_refined[1]) and \
                           (lo_patch_refined[2] - num_ghosts[2] < lo_root_level_refined[2]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the left-top-front corner.
                    
                    if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True) and \
                       (up_subdomain[2] == up_root_level_refined[2] + num_ghosts[2]) and (periodic_dimension[2] == True):
                        if (up_patch_refined[0] + num_ghosts[0] > up_root_level_refined[0]) and \
                           (lo_patch_refined[1] - num_ghosts[1] < lo_root_level_refined[1]) and \
                           (lo_patch_refined[2] - num_ghosts[2] < lo_root_level_refined[2]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the right-bottom-front corner.
                    
                    if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True) and \
                       (up_subdomain[2] == up_root_level_refined[2] + num_ghosts[2]) and (periodic_dimension[2] == True):
                        if (lo_patch_refined[0] - num_ghosts[0] < lo_root_level_refined[0]) and \
                           (up_patch_refined[1] + num_ghosts[1] > up_root_level_refined[1]) and \
                           (lo_patch_refined[2] - num_ghosts[2] < lo_root_level_refined[2]):
                            load_file_clusters[file_cluster_num] = True
                    
                    # Check the right-top-front corner.
                    
                    if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                       (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True) and \
                       (up_subdomain[2] == up_root_level_refined[2] + num_ghosts[2]) and (periodic_dimension[2] == True):
                        if (lo_patch_refined[0] - num_ghosts[0] < lo_root_level_refined[0]) and \
                           (lo_patch_refined[1] - num_ghosts[1] < lo_root_level_refined[1]) and \
                           (lo_patch_refined[2] - num_ghosts[2] < lo_root_level_refined[2]):
                            load_file_clusters[file_cluster_num] = True
                
                else:
                    raise RuntimeError('Problem dimension < 1 or > 3 not supported!')
        
        
        # Get the coordinates of the sub-domain refined to the finest level.
        
        full_domain_shape = up_root_level[0:dim] - lo_root_level[0:dim] + numpy.ones(dim, dtype = numpy.int)
        full_domain_shape = numpy.multiply(full_domain_shape, ratios_to_finest_level[0][0:dim]) + 2*num_ghosts
        
        if dim == 1:
            self.__x_coords = numpy.linspace(x_lo_root_level[0] + (0.5 - num_ghosts[0])*dx[-1][0], \
                x_up_root_level[0] + (num_ghosts[0] - 0.5)*dx[-1][0], \
                num = full_domain_shape[0])
            
            self.__x_coords = x_coords[lo_subdomain[0] + num_ghosts[0]:up_subdomain[0] + 1 + num_ghosts[0]]
        
        elif dim == 2:
            self.__x_coords = numpy.linspace(x_lo_root_level[0] + (0.5 - num_ghosts[0])*dx[-1][0], \
                x_up_root_level[0] + (num_ghosts[0] - 0.5)*dx[-1][0], \
                num = full_domain_shape[0])
            self.__y_coords = numpy.linspace(x_lo_root_level[1] + (0.5 - num_ghosts[1])*dx[-1][1],
                x_up_root_level[1] + (num_ghosts[1] - 0.5)*dx[-1][1], \
                num = full_domain_shape[1])
            
            self.__x_coords = self.__x_coords[lo_subdomain[0] + num_ghosts[0]:up_subdomain[0] + 1 + num_ghosts[0]]
            self.__y_coords = self.__y_coords[lo_subdomain[1] + num_ghosts[1]:up_subdomain[1] + 1 + num_ghosts[1]]
        
        elif dim == 3:
            self.__x_coords = numpy.linspace(x_lo_root_level[0] + (0.5 - num_ghosts[0])*dx[-1][0], \
                x_up_root_level[0] + (num_ghosts[0] - 0.5)*dx[-1][0], \
                num = full_domain_shape[0])
            self.__y_coords = numpy.linspace(x_lo_root_level[1] + (0.5 - num_ghosts[1])*dx[-1][1], \
                x_up_root_level[1] + (num_ghosts[1] - 0.5)*dx[-1][1], \
                num = full_domain_shape[1])
            self.__z_coords = numpy.linspace(x_lo_root_level[2] + (0.5 - num_ghosts[2])*dx[-1][2], \
                x_up_root_level[2] + (num_ghosts[2] - 0.5)*dx[-1][2], \
                num = full_domain_shape[2])
            
            self.__x_coords = self.__x_coords[lo_subdomain[0] + num_ghosts[0]:up_subdomain[0] + 1 + num_ghosts[0]]
            self.__y_coords = self.__y_coords[lo_subdomain[1] + num_ghosts[1]:up_subdomain[1] + 1 + num_ghosts[1]]
            self.__z_coords = self.__z_coords[lo_subdomain[2] + num_ghosts[2]:up_subdomain[2] + 1 + num_ghosts[2]]
        
        elif dim < 1 or dim > 3:
            raise RuntimeError('Problem dimension < 1 or > 3 not supported!')
        
        # Initialize containers to store the data at different levels. The elements in the containers 
        # are initialized as NAN values.
        
        level_data = {}
        
        for var_name in var_names:
            data_shape = up_subdomain[0:dim] - lo_subdomain[0:dim] + numpy.ones(dim, dtype = numpy.int)
            data_shape = numpy.insert(data_shape, 0, var_num_components[var_name])
            data = numpy.empty(data_shape)
            data[:] = numpy.NAN
            
            level_data[var_name] = []
            for level_idx in range(num_levels):
                level_data[var_name].append(copy.deepcopy(data))
        
        # Get the data from all levels and upsampling them to the finest level.
        
        file_clusters_to_load = numpy.arange(num_file_clusters)
        file_clusters_to_load = file_clusters_to_load[load_file_clusters]
        
        if dim == 1:
            for process_idx in file_clusters_to_load:
                file_name = 'processor_cluster.' + str(process_idx).zfill(5) + '.samrai'
                full_path = self.__data_directory_path + '/' + file_name
                f_input = h5py.File(full_path, 'r')
                
                for level_num in range(num_levels):
                    patch_level_start_idx = 0
                    for level_idx in range(0, level_num):
                        patch_level_start_idx = patch_level_start_idx + num_patches[level_idx]
                    
                    file_cluster = f_input['processor.' + str(process_idx).zfill(5)]
                    file_cluster_level = file_cluster['level.' + str(level_num).zfill(5)]
                    
                    for var_name in var_names:
                        for patch_key in file_cluster_level:
                            patch_idx = int(patch_key.replace('patch.', ''))
                            global_patch_idx = patch_level_start_idx + patch_idx
                            file_cluster_patch = file_cluster_level[patch_key]
                            
                            # Get the lower and upper indices of the current patch.
                            
                            lo_patch = patch_extents[global_patch_idx][0]
                            up_patch = patch_extents[global_patch_idx][1]
                            
                            # Get the shape of the patch.
                            
                            patch_shape = up_patch[0] - lo_patch[0] + 1
                            
                            for component_idx in range(0, var_num_components[var_name]):
                                # Get the patch data and upsample the data to the finest resolution.
                                
                                patch_data = file_cluster_patch[var_component_names[var_name][component_idx]].value.reshape( \
                                    patch_shape, order = 'F')
                                patch_data = numpy.repeat(patch_data, ratios_to_finest_level[level_num][0], axis = 0)
                                
                                # Get the indices after the refinement of the current patch.
                                
                                lo_patch_refined = numpy.empty(1, dtype = numpy.int)
                                up_patch_refined = numpy.empty(1, dtype = numpy.int)
                                
                                lo_patch_refined[0] = lo_patch[0]*ratios_to_finest_level[level_num][0]
                                
                                up_patch_refined[0] = (up_patch[0] + 1)*ratios_to_finest_level[level_num][0] - 1
                                
                                self.__loadDataFromPatchToSubdomain(lo_subdomain, up_subdomain, \
                                    lo_patch_refined, up_patch_refined, \
                                    level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check whether the sub-domain touches periodic boundaries.
                                # If the sub-domain touches the periodic boundaries, check whether ghost cell region
                                # overlaps with the patches.
                                
                                lo_subdomain_shifted = numpy.empty(1, dtype = numpy.int)
                                up_subdomain_shifted = numpy.empty(1, dtype = numpy.int)
                                
                                # Check the left boundary.
                                
                                if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True):
                                    if (up_patch_refined[0] + num_ghosts[0] > up_root_level_refined[0]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] + domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] + domain_shape[0]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the right boundary.
                                
                                if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True):
                                    if (lo_patch_refined[0] - num_ghosts[0]  < lo_root_level_refined[0]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] - domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] - domain_shape[0]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                
                f_input.close()
        
        elif dim == 2:
            for process_idx in range(0, num_file_clusters):
                file_name = 'processor_cluster.' + str(process_idx).zfill(5) + '.samrai'
                full_path = self.__data_directory_path + '/' + file_name
                f_input = h5py.File(full_path, 'r')
                
                for level_num in range(num_levels):
                    patch_level_start_idx = 0
                    for level_idx in range(0, level_num):
                        patch_level_start_idx = patch_level_start_idx + num_patches[level_idx]
                    
                    file_cluster = f_input['processor.' + str(process_idx).zfill(5)]
                    file_cluster_level = file_cluster['level.' + str(level_num).zfill(5)]
                    
                    for var_name in var_names:
                        for patch_key in file_cluster_level:
                            patch_idx = int(patch_key.replace('patch.', ''))
                            global_patch_idx = patch_level_start_idx + patch_idx
                            file_cluster_patch = file_cluster_level[patch_key]
                            
                            # Get the lower and upper indices of the current patch.
                            
                            lo_patch = self.__patch_extents[global_patch_idx][0]
                            up_patch = self.__patch_extents[global_patch_idx][1]
                            
                            # Get the shape of the patch.
                            
                            patch_shape = up_patch[0:2] - lo_patch[0:2] + numpy.ones(2, dtype = numpy.int)
                            
                            for component_idx in range(0, var_num_components[var_name]):
                                # Get the patch data and upsample the data to the finest resolution.
                                
                                patch_data = file_cluster_patch[var_component_names[var_name][component_idx]].value.reshape( \
                                    patch_shape, order = 'F')
                                patch_data = numpy.repeat(patch_data, ratios_to_finest_level[level_num][0], axis = 0)
                                patch_data = numpy.repeat(patch_data, ratios_to_finest_level[level_num][1], axis = 1)
                                
                                # Get the indices after the refinement of the current patch.
                                
                                lo_patch_refined = numpy.empty(2, dtype = numpy.int)
                                up_patch_refined = numpy.empty(2, dtype = numpy.int)
                                
                                lo_patch_refined[0] = lo_patch[0]*ratios_to_finest_level[level_num][0]
                                lo_patch_refined[1] = lo_patch[1]*ratios_to_finest_level[level_num][1]
                                
                                up_patch_refined[0] = (up_patch[0] + 1)*ratios_to_finest_level[level_num][0] - 1
                                up_patch_refined[1] = (up_patch[1] + 1)*ratios_to_finest_level[level_num][1] - 1
                                
                                self.__loadDataFromPatchToSubdomain(lo_subdomain, up_subdomain, \
                                    lo_patch_refined, up_patch_refined, \
                                    level_data[var_name][level_num][component_idx, :, :], patch_data)
                                
                                # Check whether the sub-domain touches periodic boundaries.
                                # If the sub-domain touches the periodic boundaries, check whether ghost cell region
                                # overlaps with the patches.
                                
                                lo_subdomain_shifted = numpy.empty(2, dtype = numpy.int)
                                up_subdomain_shifted = numpy.empty(2, dtype = numpy.int)
                                
                                # Check the left edge.
                                
                                if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True):
                                    if (up_patch_refined[0] + num_ghosts[0] > up_root_level_refined[0]) and \
                                       (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] + domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] + domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1]
                                        up_subdomain_shifted[1] = up_subdomain[1]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the right edge.
                                
                                if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True):
                                    if (lo_patch_refined[0] - num_ghosts[0] < lo_root_level_refined[0]) and \
                                       (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] - domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] - domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1]
                                        up_subdomain_shifted[1] = up_subdomain[1]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the bottom edge.
                                
                                if (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True):
                                    if (up_patch_refined[1] + num_ghosts[1] > up_root_level_refined[1]) and \
                                       (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0]
                                        up_subdomain_shifted[0] = up_subdomain[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] + domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] + domain_shape[1]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the top edge.
                                
                                if (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True):
                                    if (lo_patch_refined[1] - num_ghosts[1] < lo_root_level_refined[1]) and \
                                       (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0]
                                        up_subdomain_shifted[0] = up_subdomain[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] - domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] - domain_shape[1]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the left-bottom corner.
                                
                                if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True):
                                    if (up_patch_refined[0] + num_ghosts[0] > up_root_level_refined[0]) and \
                                       (up_patch_refined[1] + num_ghosts[1] > up_root_level_refined[1]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] + domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] + domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] + domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] + domain_shape[1]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the left-top corner.
                                
                                if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True):
                                    if (up_patch_refined[0] + num_ghosts[0] > up_root_level_refined[0]) and \
                                       (lo_patch_refined[1] - num_ghosts[1] < lo_root_level_refined[1]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] + domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] + domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] - domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] - domain_shape[1]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the right-bottom corner.
                                
                                if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True):
                                    if (lo_patch_refined[0] - num_ghosts[0] < lo_root_level_refined[0]) and \
                                       (up_patch_refined[1] + num_ghosts[1] > up_root_level_refined[1]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] - domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] - domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] + domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] + domain_shape[1]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the right-top corner.
                                
                                if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True):
                                    if (lo_patch_refined[0] - num_ghosts[0] < lo_root_level_refined[0]) and \
                                       (lo_patch_refined[1] - num_ghosts[1] < lo_root_level_refined[1]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] - domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] - domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] - domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] - domain_shape[1]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                
                f_input.close()
        
        elif dim == 3:
            for process_idx in range(0, num_file_clusters):
                file_name = 'processor_cluster.' + str(process_idx).zfill(5) + '.samrai'
                full_path = self.__data_directory_path + '/' + file_name
                f_input = h5py.File(full_path, 'r')
                
                for level_num in range(num_levels):
                    patch_level_start_idx = 0
                    for level_idx in range(0, level_num):
                        patch_level_start_idx = patch_level_start_idx + num_patches[level_idx]
                    
                    file_cluster = f_input['processor.' + str(process_idx).zfill(5)]
                    file_cluster_level = file_cluster['level.' + str(level_num).zfill(5)]
                    
                    for var_name in var_names:
                        for patch_key in file_cluster_level:
                            patch_idx = int(patch_key.replace('patch.', ''))
                            global_patch_idx = patch_level_start_idx + patch_idx
                            file_cluster_patch = file_cluster_level[patch_key]
                            
                            # Get the lower and upper indices of the current patch.
                            
                            lo_patch = patch_extents[global_patch_idx][0]
                            up_patch = patch_extents[global_patch_idx][1]
                            
                            # Get the shape of the patch.
                            
                            patch_shape = up_patch - lo_patch + numpy.ones(3, dtype = numpy.int)
                            
                            for component_idx in range(0, var_num_components[var_name]):
                                # Get the patch data and upsample the data to the finest resolution.
                                
                                patch_data = file_cluster_patch[var_component_names[var_name][component_idx]].value.reshape( \
                                    patch_shape, order = 'F')
                                patch_data = numpy.repeat(patch_data, ratios_to_finest_level[level_num][0], axis = 0)
                                patch_data = numpy.repeat(patch_data, ratios_to_finest_level[level_num][1], axis = 1)
                                patch_data = numpy.repeat(patch_data, ratios_to_finest_level[level_num][2], axis = 2)
                                
                                # Get the indices after the refinement of the current patch.
                                
                                lo_patch_refined = numpy.empty(3, dtype = numpy.int)
                                up_patch_refined = numpy.empty(3, dtype = numpy.int)
                                
                                lo_patch_refined[0] = lo_patch[0]*ratios_to_finest_level[level_num][0]
                                lo_patch_refined[1] = lo_patch[1]*ratios_to_finest_level[level_num][1]
                                lo_patch_refined[2] = lo_patch[2]*ratios_to_finest_level[level_num][2]
                                
                                up_patch_refined[0] = (up_patch[0] + 1)*ratios_to_finest_level[level_num][0] - 1
                                up_patch_refined[1] = (up_patch[1] + 1)*ratios_to_finest_level[level_num][1] - 1
                                up_patch_refined[2] = (up_patch[2] + 1)*ratios_to_finest_level[level_num][1] - 1
                                
                                self.__loadDataFromPatchToSubdomain(lo_subdomain, up_subdomain, \
                                    lo_patch_refined, up_patch_refined, \
                                    level_data[var_name][level_num][component_idx, :, :, :], patch_data)
                                
                                # Check whether the sub-domain touches periodic boundaries.
                                # If the sub-domain touches the periodic boundaries, check whether ghost cell region
                                # overlaps with the patches.
                                
                                lo_subdomain_shifted = numpy.empty(3, dtype = numpy.int)
                                up_subdomain_shifted = numpy.empty(3, dtype = numpy.int)
                                
                                # Check the left face.
                                
                                if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]) and \
                                   (lo_patch_refined[2] <= up_subdomain[2]) and (up_patch_refined[2] >= lo_subdomain[2]):
                                    if (up_patch_refined[0] + num_ghosts[0] > up_root_level_refined[0]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] + domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] + domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1]
                                        up_subdomain_shifted[1] = up_subdomain[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2]
                                        up_subdomain_shifted[2] = up_subdomain[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the right face.
                                
                                if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]) and \
                                   (lo_patch_refined[2] <= up_subdomain[2]) and (up_patch_refined[2] >= lo_subdomain[2]):
                                    if (lo_patch_refined[0] - num_ghosts[0] < lo_root_level_refined[0]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] - domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] - domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1]
                                        up_subdomain_shifted[1] = up_subdomain[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2]
                                        up_subdomain_shifted[2] = up_subdomain[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the bottom face.
                                
                                if (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True) and \
                                   (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]) and \
                                   (lo_patch_refined[2] <= up_subdomain[2]) and (up_patch_refined[2] >= lo_subdomain[2]):
                                    if (up_patch_refined[1] + num_ghosts[1] > up_root_level_refined[1]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0]
                                        up_subdomain_shifted[0] = up_subdomain[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] + domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] + domain_shape[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2]
                                        up_subdomain_shifted[2] = up_subdomain[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the top face.
                                
                                if (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True) and \
                                   (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]) and \
                                   (lo_patch_refined[2] <= up_subdomain[2]) and (up_patch_refined[2] >= lo_subdomain[2]):
                                    if (lo_patch_refined[1] - num_ghosts[1] < lo_root_level_refined[1]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0]
                                        up_subdomain_shifted[0] = up_subdomain[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] - domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] - domain_shape[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2]
                                        up_subdomain_shifted[2] = up_subdomain[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the back face.
                                
                                if (lo_subdomain[2] == lo_root_level_refined[2] - num_ghosts[2]) and (periodic_dimension[2] == True) and \
                                   (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]) and \
                                   (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]):
                                    if (up_patch_refined[2] + num_ghosts[2] > up_root_level_refined[1]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0]
                                        up_subdomain_shifted[0] = up_subdomain[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1]
                                        up_subdomain_shifted[1] = up_subdomain[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2] + domain_shape[2]
                                        up_subdomain_shifted[2] = up_subdomain[2] + domain_shape[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the front face.
                                
                                if (up_subdomain[2] == up_root_level_refined[2] + num_ghosts[2]) and (periodic_dimension[2] == True) and \
                                   (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]) and \
                                   (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]):
                                    if (lo_patch_refined[2] - num_ghosts[2] < lo_root_level_refined[2]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0]
                                        up_subdomain_shifted[0] = up_subdomain[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1]
                                        up_subdomain_shifted[1] = up_subdomain[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2] - domain_shape[2]
                                        up_subdomain_shifted[2] = up_subdomain[2] - domain_shape[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the left-bottom edge.
                                
                                if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True) and \
                                   (lo_patch_refined[2] <= up_subdomain[2]) and (up_patch_refined[2] >= lo_subdomain[2]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] + domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] + domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] + domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] + domain_shape[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2]
                                        up_subdomain_shifted[2] = up_subdomain[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the left-top edge.
                                
                                if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True) and \
                                   (lo_patch_refined[2] <= up_subdomain[2]) and (up_patch_refined[2] >= lo_subdomain[2]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] + domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] + domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] - domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] - domain_shape[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2]
                                        up_subdomain_shifted[2] = up_subdomain[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the left-back edge.
                                
                                if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (lo_subdomain[2] == lo_root_level_refined[2] - num_ghosts[2]) and (periodic_dimension[2] == True) and \
                                   (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] + domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] + domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1]
                                        up_subdomain_shifted[1] = up_subdomain[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2] + domain_shape[2]
                                        up_subdomain_shifted[2] = up_subdomain[2] + domain_shape[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the left-front edge.
                                
                                if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (up_subdomain[2] == up_root_level_refined[2] + num_ghosts[2]) and (periodic_dimension[2] == True) and \
                                   (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] + domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] + domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1]
                                        up_subdomain_shifted[1] = up_subdomain[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2] - domain_shape[2]
                                        up_subdomain_shifted[2] = up_subdomain[2] - domain_shape[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the right-bottom edge.
                                
                                if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True) and \
                                   (lo_patch_refined[2] <= up_subdomain[2]) and (up_patch_refined[2] >= lo_subdomain[2]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] - domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] - domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] + domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] + domain_shape[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2]
                                        up_subdomain_shifted[2] = up_subdomain[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the right-top edge.
                                
                                if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True) and \
                                   (lo_patch_refined[2] <= up_subdomain[2]) and (up_patch_refined[2] >= lo_subdomain[2]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] - domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] - domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] - domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] - domain_shape[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2]
                                        up_subdomain_shifted[2] = up_subdomain[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the right-back edge.
                                
                                if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (lo_subdomain[2] == lo_root_level_refined[2] - num_ghosts[2]) and (periodic_dimension[2] == True) and \
                                   (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] - domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] - domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1]
                                        up_subdomain_shifted[1] = up_subdomain[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2] + domain_shape[2]
                                        up_subdomain_shifted[2] = up_subdomain[2] + domain_shape[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the right-front edge.
                                
                                if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (up_subdomain[2] == up_root_level_refined[2] + num_ghosts[2]) and (periodic_dimension[2] == True) and \
                                   (lo_patch_refined[1] <= up_subdomain[1]) and (up_patch_refined[1] >= lo_subdomain[1]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] - domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] - domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1]
                                        up_subdomain_shifted[1] = up_subdomain[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2] - domain_shape[2]
                                        up_subdomain_shifted[2] = up_subdomain[2] - domain_shape[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the bottom-back edge.
                                
                                if (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True) and \
                                   (lo_subdomain[2] == lo_root_level_refined[2] - num_ghosts[2]) and (periodic_dimension[2] == True) and \
                                   (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0]
                                        up_subdomain_shifted[0] = up_subdomain[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] + domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] + domain_shape[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2] + domain_shape[2]
                                        up_subdomain_shifted[2] = up_subdomain[2] + domain_shape[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the bottom-front edge.
                                
                                if (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True) and \
                                   (up_subdomain[2] == up_root_level_refined[2] + num_ghosts[2]) and (periodic_dimension[2] == True) and \
                                   (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0]
                                        up_subdomain_shifted[0] = up_subdomain[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] + domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] + domain_shape[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2] - domain_shape[2]
                                        up_subdomain_shifted[2] = up_subdomain[2] - domain_shape[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the top-back edge.
                                
                                if (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True) and \
                                   (lo_subdomain[2] == lo_root_level_refined[2] - num_ghosts[2]) and (periodic_dimension[2] == True) and \
                                   (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0]
                                        up_subdomain_shifted[0] = up_subdomain[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] - domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] - domain_shape[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2] + domain_shape[2]
                                        up_subdomain_shifted[2] = up_subdomain[2] + domain_shape[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the top-front edge.
                                
                                if (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True) and \
                                   (up_subdomain[2] == up_root_level_refined[2] + num_ghosts[2]) and (periodic_dimension[2] == True) and \
                                   (lo_patch_refined[0] <= up_subdomain[0]) and (up_patch_refined[0] >= lo_subdomain[0]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0]
                                        up_subdomain_shifted[0] = up_subdomain[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] - domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] - domain_shape[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2] - domain_shape[2]
                                        up_subdomain_shifted[2] = up_subdomain[2] - domain_shape[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the left-bottom-back corner.
                                
                                if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True) and \
                                   (lo_subdomain[2] == lo_root_level_refined[2] - num_ghosts[2]) and (periodic_dimension[2] == True):
                                    if (up_patch_refined[0] + num_ghosts[0] > up_root_level_refined[0]) and \
                                       (up_patch_refined[1] + num_ghosts[1] > up_root_level_refined[1]) and \
                                       (up_patch_refined[2] + num_ghosts[2] > up_root_level_refined[2]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] + domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] + domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] + domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] + domain_shape[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2] + domain_shape[2]
                                        up_subdomain_shifted[2] = up_subdomain[2] + domain_shape[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the left-top-back corner.
                                
                                if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True) and \
                                   (lo_subdomain[2] == lo_root_level_refined[2] - num_ghosts[2]) and (periodic_dimension[2] == True):
                                    if (up_patch_refined[0] + num_ghosts[0] > up_root_level_refined[0]) and \
                                       (lo_patch_refined[1] - num_ghosts[1] < lo_root_level_refined[1]) and \
                                       (up_patch_refined[2] + num_ghosts[2] > up_root_level_refined[2]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] + domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] + domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] - domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] - domain_shape[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2] + domain_shape[2]
                                        up_subdomain_shifted[2] = up_subdomain[2] + domain_shape[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the right-bottom-back corner.
                                
                                if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True) and \
                                   (lo_subdomain[2] == lo_root_level_refined[2] - num_ghosts[2]) and (periodic_dimension[2] == True):
                                    if (lo_patch_refined[0] - num_ghosts[0] < lo_root_level_refined[0]) and \
                                       (up_patch_refined[1] + num_ghosts[1] > up_root_level_refined[1]) and \
                                       (up_patch_refined[2] + num_ghosts[2] > up_root_level_refined[2]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] - domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] - domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] + domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] + domain_shape[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2] + domain_shape[2]
                                        up_subdomain_shifted[2] = up_subdomain[2] + domain_shape[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the right-top-back corner.
                                
                                if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True) and \
                                   (lo_subdomain[2] == lo_root_level_refined[2] - num_ghosts[2]) and (periodic_dimension[2] == True):
                                    if (lo_patch_refined[0] - num_ghosts[0] < lo_root_level_refined[0]) and \
                                       (lo_patch_refined[1] - num_ghosts[1] < lo_root_level_refined[1]) and \
                                       (up_patch_refined[2] + num_ghosts[2] > up_root_level_refined[2]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] - domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] - domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] - domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] - domain_shape[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2] + domain_shape[2]
                                        up_subdomain_shifted[2] = up_subdomain[2] + domain_shape[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the left-bottom-front corner.
                                
                                if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True) and \
                                   (up_subdomain[2] == up_root_level_refined[2] + num_ghosts[2]) and (periodic_dimension[2] == True):
                                    if (up_patch_refined[0] + num_ghosts[0] > up_root_level_refined[0]) and \
                                       (up_patch_refined[1] + num_ghosts[1] > up_root_level_refined[1]) and \
                                       (lo_patch_refined[2] - num_ghosts[2] < lo_root_level_refined[2]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] + domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] + domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] + domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] + domain_shape[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2] - domain_shape[2]
                                        up_subdomain_shifted[2] = up_subdomain[2] - domain_shape[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the left-top-front corner.
                                
                                if (lo_subdomain[0] == lo_root_level_refined[0] - num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True) and \
                                   (up_subdomain[2] == up_root_level_refined[2] + num_ghosts[2]) and (periodic_dimension[2] == True):
                                    if (up_patch_refined[0] + num_ghosts[0] > up_root_level_refined[0]) and \
                                       (lo_patch_refined[1] - num_ghosts[1] < lo_root_level_refined[1]) and \
                                       (lo_patch_refined[2] - num_ghosts[2] < lo_root_level_refined[2]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] + domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] + domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] - domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] - domain_shape[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2] - domain_shape[2]
                                        up_subdomain_shifted[2] = up_subdomain[2] - domain_shape[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the right-bottom-front corner.
                                
                                if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (lo_subdomain[1] == lo_root_level_refined[1] - num_ghosts[1]) and (periodic_dimension[1] == True) and \
                                   (up_subdomain[2] == up_root_level_refined[2] + num_ghosts[2]) and (periodic_dimension[2] == True):
                                    if (lo_patch_refined[0] - num_ghosts[0] < lo_root_level_refined[0]) and \
                                       (up_patch_refined[1] + num_ghosts[1] > up_root_level_refined[1]) and \
                                       (lo_patch_refined[2] - num_ghosts[2] < lo_root_level_refined[2]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] - domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] - domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] + domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] + domain_shape[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2] - domain_shape[2]
                                        up_subdomain_shifted[2] = up_subdomain[2] - domain_shape[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                                
                                # Check the right-top-front corner.
                                
                                if (up_subdomain[0] == up_root_level_refined[0] + num_ghosts[0]) and (periodic_dimension[0] == True) and \
                                   (up_subdomain[1] == up_root_level_refined[1] + num_ghosts[1]) and (periodic_dimension[1] == True) and \
                                   (up_subdomain[2] == up_root_level_refined[2] + num_ghosts[2]) and (periodic_dimension[2] == True):
                                    if (lo_patch_refined[0] - num_ghosts[0] < lo_root_level_refined[0]) and \
                                       (lo_patch_refined[1] - num_ghosts[1] < lo_root_level_refined[1]) and \
                                       (lo_patch_refined[2] - num_ghosts[2] < lo_root_level_refined[2]):
                                        lo_subdomain_shifted[0] = lo_subdomain[0] - domain_shape[0]
                                        up_subdomain_shifted[0] = up_subdomain[0] - domain_shape[0]
                                        
                                        lo_subdomain_shifted[1] = lo_subdomain[1] - domain_shape[1]
                                        up_subdomain_shifted[1] = up_subdomain[1] - domain_shape[1]
                                        
                                        lo_subdomain_shifted[2] = lo_subdomain[2] - domain_shape[2]
                                        up_subdomain_shifted[2] = up_subdomain[2] - domain_shape[2]
                                        
                                        self.__loadDataFromPatchToSubdomain(lo_subdomain_shifted, up_subdomain_shifted, \
                                            lo_patch_refined, up_patch_refined, \
                                            level_data[var_name][level_num][component_idx, :], patch_data)
                
                f_input.close()
        
        else:
            raise RuntimeError('Problem dimension < 1 or > 3 not supported!')
        
        
        # Combine data at all levels.
       
        for var_name in var_names:
            self.__data[var_name] = level_data[var_name][0]
            
            for level_idx in range(1, num_levels):
                is_finite_idx = numpy.isfinite(level_data[var_name][level_idx])
                self.__data[var_name][is_finite_idx] = level_data[var_name][level_idx][is_finite_idx]
        
        self.__data_loaded = True
    
    
    def __loadDataFromPatchToSubdomain(self, \
            lo_subdomain, \
            up_subdomain, \
            lo_patch_refined, \
            up_patch_refined, \
            subdomain_data, \
            patch_data):
        """
        Private method to load data from patch to sub-domain.
        """
        
        dim = lo_subdomain.shape[0]
        
        if dim == 1:
            # Get the global start and end indices.
            
            x_gelf.__basic_infoobal_start_idx = max(lo_patch_refined[0], lo_subdomain[0])
            x_global_end_idx = min(up_patch_refined[0] + 1, up_subdomain[0] + 1)
            
            # Get the local start and end indices in the sub-domain.
            
            x_local_start_idx = x_global_start_idx - lo_subdomain[0]
            x_local_start_idx = max(0, x_local_start_idx)
            
            x_local_end_idx = x_global_end_idx - lo_subdomain[0]
            x_local_end_idx = min(up_subdomain[0] - lo_subdomain[0] + 1, x_local_end_idx)
            
            # Get the local start and end indices in the current patch.
            
            x_patch_start_idx = x_global_start_idx - lo_patch_refined[0]
            x_patch_start_idx = max(0, x_patch_start_idx)
            
            x_patch_end_idx = patch_data.shape[0] - ((up_patch_refined[0] + 1) - x_global_end_idx)
            x_patch_end_idx = min(patch_data.shape[0], x_patch_end_idx)
            
            if x_local_start_idx < x_local_end_idx:
                subdomain_data[x_local_start_idx:x_local_end_idx] = \
                    patch_data[x_patch_start_idx:x_patch_end_idx]
        
        elif dim == 2:
            # Get the global start and end indices.
            
            x_global_start_idx = max(lo_patch_refined[0], lo_subdomain[0])
            x_global_end_idx = min(up_patch_refined[0] + 1, up_subdomain[0] + 1)
            
            y_global_start_idx = max(lo_patch_refined[1], lo_subdomain[1])
            y_global_end_idx = min(up_patch_refined[1] + 1, up_subdomain[1] + 1)
            
            # Get the local start and end indices in the sub-domain.
            
            x_local_start_idx = x_global_start_idx - lo_subdomain[0]
            x_local_start_idx = max(0, x_local_start_idx)
            
            x_local_end_idx = x_global_end_idx - lo_subdomain[0]
            x_local_end_idx = min(up_subdomain[0] - lo_subdomain[0] + 1, x_local_end_idx)
            
            y_local_start_idx = y_global_start_idx - lo_subdomain[1]
            y_local_start_idx = max(0, y_local_start_idx)
            
            y_local_end_idx = y_global_end_idx - lo_subdomain[1]
            y_local_end_idx = min(up_subdomain[1] - lo_subdomain[1] + 1, y_local_end_idx)
            
            # Get the local start and end indices in the current patch.
            
            x_patch_start_idx = x_global_start_idx - lo_patch_refined[0]
            x_patch_start_idx = max(0, x_patch_start_idx)
            
            x_patch_end_idx = patch_data.shape[0] - ((up_patch_refined[0] + 1) - x_global_end_idx)
            x_patch_end_idx = min(patch_data.shape[0], x_patch_end_idx)
            
            y_patch_start_idx = y_global_start_idx - lo_patch_refined[1]
            y_patch_start_idx = max(0, y_patch_start_idx)
            
            y_patch_end_idx = patch_data.shape[1] - ((up_patch_refined[1] + 1) - y_global_end_idx)
            y_patch_end_idx = min(patch_data.shape[1], y_patch_end_idx)
            
            if (x_local_start_idx < x_local_end_idx) and (y_local_start_idx < y_local_end_idx):
                subdomain_data[x_local_start_idx:x_local_end_idx, y_local_start_idx:y_local_end_idx] = \
                    patch_data[x_patch_start_idx:x_patch_end_idx, y_patch_start_idx:y_patch_end_idx]
        
        elif dim == 3:
            # Get the global start and end indices.
            
            x_global_start_idx = max(lo_patch_refined[0], lo_subdomain[0])
            x_global_end_idx = min(up_patch_refined[0] + 1, up_subdomain[0] + 1)
            
            y_global_start_idx = max(lo_patch_refined[1], lo_subdomain[1])
            y_global_end_idx = min(up_patch_refined[1] + 1, up_subdomain[1] + 1)
            
            z_global_start_idx = max(lo_patch_refined[2], lo_subdomain[2])
            z_global_end_idx = min(up_patch_refined[2] + 1, up_subdomain[2] + 1)
            
            # Get the local start and end indices in the sub-domain.
            
            x_local_start_idx = x_global_start_idx - lo_subdomain[0]
            x_local_start_idx = max(0, x_local_start_idx)
            
            x_local_end_idx = x_global_end_idx - lo_subdomain[0]
            x_local_end_idx = min(up_subdomain[0] - lo_subdomain[0] + 1, x_local_end_idx)
            
            y_local_start_idx = y_global_start_idx - lo_subdomain[1]
            y_local_start_idx = max(0, y_local_start_idx)
            
            y_local_end_idx = y_global_end_idx - lo_subdomain[1]
            y_local_end_idx = min(up_subdomain[1] - lo_subdomain[1] + 1, y_local_end_idx)
            
            z_local_start_idx = z_global_start_idx - lo_subdomain[2]
            z_local_start_idx = max(0, z_local_start_idx)
            
            z_local_end_idx = z_global_end_idx - lo_subdomain[2]
            z_local_end_idx = min(up_subdomain[2] - lo_subdomain[2] + 1, z_local_end_idx)
            
            # Get the local start and end indices in the current patch.
            
            x_patch_start_idx = x_global_start_idx - lo_patch_refined[0]
            x_patch_start_idx = max(0, x_patch_start_idx)
            
            x_patch_end_idx = patch_data.shape[0] - ((up_patch_refined[0] + 1) - x_global_end_idx)
            x_patch_end_idx = min(patch_data.shape[0], x_patch_end_idx)
            
            y_patch_start_idx = y_global_start_idx - lo_patch_refined[1]
            y_patch_start_idx = max(0, y_patch_start_idx)
            
            y_patch_end_idx = patch_data.shape[1] - ((up_patch_refined[1] + 1) - y_global_end_idx)
            y_patch_end_idx = min(patch_data.shape[1], y_patch_end_idx)
            
            z_patch_start_idx = z_global_start_idx - lo_patch_refined[2]
            z_patch_start_idx = max(0, z_patch_start_idx)
            
            z_patch_end_idx = patch_data.shape[2] - ((up_patch_refined[2] + 1) - y_global_end_idx)
            z_patch_end_idx = min(patch_data.shape[2], z_patch_end_idx)
            
            if (x_local_start_idx < x_local_end_idx) and \
                (y_local_start_idx < y_local_end_idx) and \
                (z_local_start_idx < z_local_end_idx):
                subdomain_data[x_local_start_idx:x_local_end_idx, \
                                      y_local_start_idx:y_local_end_idx, \
                                      z_local_start_idx:z_local_end_idx] = \
                    patch_data[x_patch_start_idx:x_patch_end_idx, \
                               y_patch_start_idx:y_patch_end_idx, \
                               z_patch_start_idx:z_patch_end_idx]
        
        else:
            raise RuntimeError('Problem dimension < 1 or > 3 not supported!')

