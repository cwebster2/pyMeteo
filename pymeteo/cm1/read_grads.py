"""Class to read GRaDS style CM1 model output data

"""
#TODO exception throwing?  Error handling, support the u,v,w,stats files

import numpy as np
import mmap

class CM1(object):
   """Class that implements reading CM1 model data

   :param path: path to CM1 data files
   :param datasetname: the CM1 data files basename
   
   """
   nx   = 0
   """gridpoints in the x direction"""
   ny   = 0
   """gridpoints in the y direction"""
   nz   = 0
   """gridpoints in the z direction"""
   nv   = 0
   """number of variables in the dataset"""
   nt   = 0
   """number of time levels in the dataset"""
   dt   = 0
   """timestep between time levels"""
   dimX = 0
   """X grid dimension (1D)"""
   dimY = 0
   """Y grid dimension (1D)"""
   dimZ = 0
   """Z grid dimension (1D)"""
   dimT = 0
   """T grid dimension (1D)"""
   vars = 0
   n2d  = 0
   recl = 0
   dsetname = ''
   path = ''
   mem_file = 0
   dat_file = 0

#-------------------------------------------------------

   def __init__(self,__path,__datasetname):
      
      self.path = __path
      self.dsetname = __datasetname

      #We need to open the control file and read the dataset metadata 

      #create filename.  This assumes cm1 created files and the scalar file at that.
      ctl_filename = self.path + '/' + self.dsetname + '_s.ctl'

      # open the control file
      with open(ctl_filename, 'r') as ctl_file:
         print('Opened '+ ctl_filename+ ' for reading')
         #read ctl file
         self.parsectl(ctl_file)
      #close ctl file
      ctl_file.closed
      print(ctl_filename+ ' closed')

#-------------------------------------------------------
#  This function currently only works for cm1 domains with stretched x/y/z grids.

   def parsectl(self,ctl_file):

      print('  Parsing control file')
      line = ctl_file.readline()
      while line:
         # This finds the line that describes the X dimension
        if line.startswith('xdef'):
          nX = [int(s) for s in line.split() if s.isdigit()].pop()
          print('    Found X with {0} values'.format(nX))
          X = []
          for i in range(nX):
            line = ctl_file.readline()
            X.append(float(line))
          if nX != len(X):
            print('Error: Expected {0} values, read {1} values, exiting!'.format(nX,len(X)))
            quit()
          self.dimX = np.array(X)
          self.nx = nX

         # This finds the line that describes the Y dimension
        if line.startswith('ydef'):
          nY = [int(s) for s in line.split() if s.isdigit()].pop()
          print('    Found Y with {0} values'.format(nY))
          Y = []
          for i in range(nY):
            line = ctl_file.readline()
            Y.append(float(line))
          if nY != len(Y):
            print('Error: Expected {0} values, read {1} values, exiting!'.format(nY,len(Y)))
            quit()
          self.dimY = np.array(Y)
          self.ny = nY

         # This finds the line that describes the Z dimension
        if line.startswith('zdef'):
          nZ = [int(s) for s in line.split() if s.isdigit()].pop()
          print('    Found Z with {0} values'.format(nZ))
          Z = []
          for i in range(nZ):
            line = ctl_file.readline()
            Z.append(float(line))
          if nZ != len(Z):
            print('Error: Expected {0} values, read {1} values, exiting!'.format(nZ,len(Z)))
            quit()
          self.dimZ = np.array(Z)
          self.nz = nZ

         # This finds the line that describes the time dimension
        if line.startswith('tdef'):
          nT = [s for s in line.split()]
          self.nt = int(nT[1])
          self.dt = float([s for s in nT[4].split('YR')][0])
          print('    Found T with {0} steps of {1} s'.format(self.nt, self.dt))
          T = []
          for i in range(self.nt):
            T.append(float(i*self.dt))
            self.dimT = np.array(T)

         # This finds the line that begins the variable declaration
        if line.startswith('vars'):
          nV = [int(s) for s in line.split() if s.isdigit()].pop()
          print('    Found {0} variables'.format(nV))
          V = []
          for i in range(nV):
            line = ctl_file.readline()
            newV = [i+1] #varid

            vdef = ([s for s in line.split(None,3)])
            newV.append(vdef[0]) # varname
            newV.append(int(vdef[1])) #levels
            if int(vdef[1]) == 0:
              self.n2d = self.n2d + 1
            newV.append(vdef[3].rstrip()) # description
            newV = tuple(newV)

            V.append(newV)
          if nV != len(V):
            print('Error: Expected {0} values, read {1} values, exiting!'.format(nV,len(V)))
            quit()
          self.vars = V
          self.nv = nV
          #eat the endvars line
          line = ctl_file.readline()

        line = ctl_file.readline()
      
      self.recl = self.nx * self.ny * 4

#-------------------------------------------------------

   def getVarByName(self,varname):
      idx = [i for i in self.vars if i[1] == varname]
      return {'id': int(idx[0][0]), 'name': idx[0][1], 'nlevs': int(idx[0][2]), 'desc': idx[0][3]}

#-------------------------------------------------------
# This function reads count number of bytes sequentially from location loc
# from the memory region mm (which may be a memory mapped file).

   def readData(self, mm, loc, count):
      
      # Reads data starting at loc, length count
      bytes = 0
      data = b''

      # seek to location
      mm.seek(loc)

      #read until we have all of the bytes
      while (bytes < count):
         newdata = mm.read(count-bytes)
         if newdata==0:
            print('Error reading memmapped data')
            #Throw
         bytes += len(newdata)
         data += newdata

      # did we get all of the bytes?
      if count != len(data):
         print('Error reading memmapped data')
	       #throw

      return data

#-------------------------------------------------------
#  This function reads one record from the datafile.  Each
#  record is nx*ny*4 bytes and is one xy-slice of the domain

   def read3dXYSlice(self, memmap, varid, level):

      # Calculate the record number and then mult by the record
      # length to calculate the offset into the file the record begins at

      idx = (self.n2d+(varid-1-self.n2d)*(self.nz))+level # record number
      loc = idx*self.recl

      # Read the data starting at location loc, length recl

      data = self.readData(memmap, loc, self.recl)

      # did we get all of the bytes?
      if self.recl != len(data):
         print('Error reading memmapped data')
	 #throw

      # convert data to numpy array of floats and shape it into i*j
      newarray = np.fromstring(data, dtype=np.float32, count=int(self.recl/4))
      newarray = newarray.reshape((self.nx, self.ny)).T
      return newarray

#-------------------------------------------------------
# Reads a single 3d variable from the datafile

   def read3d(self, time, varname):
      """Reads a 3D variable from the dataset

      :param time: the timelevel to read
      :param varname: the variable name to read
      :reaturns: 3D array containing the variable at the time
      """
      var = self.getVarByName(varname)

      # open dat file 
      dat_filename = self.path + '/' + self.dsetname + '_{0:06d}_s.dat'.format(int(time))

      print('  Opening {0} for reading'.format(dat_filename))
      print('    Reading {0}({1}) at time {2} s'.format(varname, var['id'], time))

      # open the data file
      with open(dat_filename, 'rb') as dat_file:

        # memmap the data file
        memmap = mmap.mmap(dat_file.fileno(), 0, prot=mmap.PROT_READ)

	      # init array for the variable
        data = np.empty((self.nx, self.ny, self.nz), dtype=np.float32)

        # Read the data one record at a time until we have the full domain volume
        for z in range(self.nz):
          data[:,:,z] = self.read3dXYSlice(memmap, var['id'], z)
      
      #close memory map and the data file
      memmap.close()
      dat_file.closed
      print('  {0} closed'.format(dat_filename))

      return data

#-------------------------------------------------------
# These 3 functions seperate the opening, reading and closing
#  actions so that mutliple variables can be read from the same
#  datafile with only 1 open/close action.
      
   def read3dMultStart(self,time):

      # create datafile name
      dat_filename = self.path + '/' + self.dsetname + '_{0:06d}_s.dat'.format(int(time))
      print('  Opening {0} for reading'.format(dat_filename))
      print('    Reading time {0} s'.format(time))
   
      # open and memory map the file
      self.dat_file = open(dat_filename, 'rb') 
      self.map_file = mmap.mmap(self.dat_file.fileno(), 0, prot=mmap.PROT_READ)
      
#-------------------------------------------------------

   def read3dMultStop(self):

      # close the memory map and the data file
      self.map_file.close()
      self.dat_file.closed
      print('  File closed')

#-------------------------------------------------------

   def read3dMult(self, varname):

      # get the variable information
      var = self.getVarByName(varname)
      print('    Reading var {0}({1})'.format(varname,var['id']))

      #initialize the array to hold the data
      data = np.empty((self.nx, self.ny, self.nz), dtype=np.float32)

      # read the data one slice at a time until we have the full volume
      for z in range(self.nz):
         data[:,:,z] = self.read3dXYSlice(self.map_file, var['id'], z)

      return data

#-------------------------------------------------------

# def get var by id, get varid by name
# def var exists?
