
# cm1.py
#
# access routines for cm1 grads format output

#TODO exception throwing?  Error handling, support the u,v,w,stats files

import numpy as np
import mmap
import h5py
import glob
import re
import os

class CM1(object):
   nx   = 0
   ny   = 0
   nz   = 0
   nxp1 = 0
   nyp1 = 0
   nzp1 = 0
   nt   = 0
   dt   = 0
   dimX = 0
   dimY = 0
   dimZ = 0
   dimU = 0
   dimV = 0
   dimW = 0
   dimT = 0
   dsetname = ''
   path = ''
   datafile = 0

#-------------------------------------------------------

   def __init__(self,__path,__datasetname):

      self.path = __path
      self.dsetname = __datasetname

      #We need to open the control file and read the dataset metadata

      #read ctl file
      self.getmetadata()

#-------------------------------------------------------
#  This function currently only works for cm1 domains with stretched x/y/z grids.

   def getmetadata(self):

      # This searches files that describe the time dimension
      cm1hdf_files = glob.glob("{0}/{1}*h5".format(self.path, self.dsetname))
      cm1hdf_files.sort()
      # replaces curved90-qv14.(\d{5}).h5  with \1
      self.dimT = np.array([int(re.sub(r'{0}.(\d{{5}}).h5'.format(self.dsetname),r'\1',os.path.basename(f))) for f in cm1hdf_files])
      if len(self.dimT) == 0:
         raise ValueError('Invalid dataset specified: {0}/{1}'.format(self.path, self.dsetname))
#TODO: check dimT, if 0 then no dataset found, bail out
      self.nt   = len(self.dimT)
      #self.dt   = self.dimT[1] - self.dimT[0]
      self.dt = 15
      self.nt = 601
      # TODO!!!! FIX THIS
      print('    Found T with {0} steps of {1} s'.format(self.nt, self.dt))

      filename = "{0}/{1}.{2:05d}.h5".format(self.path, self.dsetname, self.dimT[0])
      print('  Getting metadata from {0}'.format(filename))
      cm1file = h5py.File(filename,'r')

      self.nx = int(cm1file['/grid/nx'][()])
      self.ny = int(cm1file['/grid/ny'][()])
      self.nz = int(cm1file['/grid/nz'][()])
      self.nxp1 = self.nx + 1
      self.nyp1 = self.ny + 1
      self.nzp1 = self.nz + 1

      # scale dimensions to km for existing code
      self.dimX = np.array(cm1file['/mesh/xh'])
      self.dimX[:] /= 1000.
      self.dimY = np.array(cm1file['/mesh/yh'])
      self.dimY[:] /= 1000.
      self.dimZ = np.array(cm1file['/mesh/zh'])
      self.dimZ[:] /= 1000.
      self.dimU = np.array(cm1file['/mesh/xf'])
      self.dimU[:] /= 1000.
      self.dimV = np.array(cm1file['/mesh/yf'])
      self.dimV[:] /= 1000.
      self.dimW = np.array(cm1file['/mesh/zf'])
      self.dimW[:] /= 1000.

      print('    Found X with {0} values'.format(self.nx))
      print('    Found Y with {0} values'.format(self.ny))
      print('    Found Z with {0} values'.format(self.nz))

      cm1file.close()

#-------------------------------------------------------

   def read2d_slice(self, time, ib, ie, jb, je, varname):
      filename = self.path + '/' + self.dsetname + '.{0:05d}.h5'.format(int(time))
      datafile = h5py.File(filename, 'r')
      print('    Reading {0} from ({1}:{2},{3}:{4}) at time {5} s'.format(varname, ib, ie, jb, je, time))

      nx = ie-ib
      ny = je-jb

      data = np.empty((ny,nx), dtype=np.float32)
      data = datafile[varname][jb:je,ib:ie]

      datafile.close()
      return data

#-------------------------------------------------------

   def read3d_slice(self, time, ib, ie, jb, je, kb, ke, varname):
      filename = self.path + '/' + self.dsetname + '.{0:05d}.h5'.format(int(time))
      datafile = h5py.File(filename, 'r')
      print('    Reading {0} from ({1}:{2},{3}:{4},{5}:{6}) at time {7} s'.format(varname, ib, ie, jb, je, kb, ke, time))

      nx = ie-ib
      ny = je-jb
      nz = ke-kb

      data = np.empty((nz,ny,nx), dtype=np.float32)
      data = datafile[varname][kb:ke,jb:je,ib:ie]

      datafile.close()
      return data

   
   def read3d_slice_derived(self, time, ib, ie, jb, je, kb, ke, varname3D, varname1D):
      filename = self.path + '/' + self.dsetname + '.{0:05d}.h5'.format(int(time))
      datafile = h5py.File(filename, 'r')
      print('    Reading {0} and {8} from ({1}:{2},{3}:{4},{5}:{6}) at time {7} s'.format(varname3D, ib, ie, jb, je, kb, ke, time, varname1D))

      nx = ie-ib
      ny = je-jb
      nz = ke-kb

      data = np.empty((nz, ny, nx), dtype=np.float32)
      # Have to transpose before right-broadcasting doesnt work
      data = (datafile[varname3D][kb:ke, jb:je, ib:ie].T + datafile[varname1D][kb:ke].T).T

      datafile.close()
      return data


#-------------------------------------------------------
# Reads a single 3d variable from the datafile

   def read3d(self, time, grid, varname):

      # open dat file
      filename = self.path + '/' + self.dsetname + '.{0:05d}.h5'.format(int(time))

      #print('  Opening {0} for reading'.format(filename))
      print('    Reading {0} from grid {1} at time {2} s'.format(varname, grid, time))

      # open the data file
      datafile = h5py.File(filename, 'r')

      if (grid == 'u'):
         nx = self.nxp1
      else:
         nx = self.nx
      if (grid == 'v'):
         ny = self.nyp1
      else:
         ny = self.ny
      if (grid == 'w'):
         nz = self.nzp1
      else:
         nz = self.nz

      data = np.empty((nz, ny, nx), dtype=np.float32)
      data = datafile[varname]
      print(data)

      #print('  {0} closed'.format(filename))

      datafile.close()

      return data

#-------------------------------------------------------
# These 3 functions seperate the opening, reading and closing
#  actions so that mutliple variables can be read from the same
#  datafile with only 1 open/close action.

   def read3dMultStart(self,time):

      filename = self.path + '/' + self.dsetname + '.{0:05d}.h5'.format(int(time))

      print('  Opening {0} for reading'.format(filename))
      print('    Reading {0} from grid {1} at time {2} s'.format(varname, grid, time))

      # open the data file
      self.datafile = h5py.File(filename, 'r')

#-------------------------------------------------------

   def read3dMultStop(self):

      # close the memory map and the data file
      self.datafile.close()
      print('  File closed')

#-------------------------------------------------------

   def read3dMult(self, grid, varname):

      if (grid == 'u'):
         nx = self.nxp1
      else:
         nx = self.nx
      if (grid == 'v'):
         ny = self.nyp1
      else:
         ny = self.ny
      if (grid == 'w'):
         nz = self.nzp1
      else:
         nz = self.nz

      data = np.empty((nz, ny, nx), dtype=np.float32)
      data = datafile[varname]

      return np.array(data).T

#-------------------------------------------------------

   def restrict_bounds(self,east,west,north,south,height):

       def find_nearest_idx(array,value):
           return (np.abs(array-value)).argmin()

       idx_s_e = find_nearest_idx(self.dimX, east)
       idx_u_e = find_nearest_idx(self.dimU, east)

       idx_s_w = find_nearest_idx(self.dimX, west)
       idx_u_w = find_nearest_idx(self.dimU, west)

       idx_s_n = find_nearest_idx(self.dimY, north)
       idx_v_n = find_nearest_idx(self.dimV, north)

       idx_s_s = find_nearest_idx(self.dimY, south)
       idx_v_s = find_nearest_idx(self.dimV, south)

       idx_s_z = find_nearest_idx(self.dimZ, height)
       idx_w_z = find_nearest_idx(self.dimW, height)

       return (idx_s_e, idx_s_w, idx_s_n, idx_s_s,
               idx_u_e, idx_u_w, idx_v_n, idx_v_s,
               idx_s_z, idx_w_z)

#-------------------------------------------------------
# def get var by id, get varid by name
# def var exists?
