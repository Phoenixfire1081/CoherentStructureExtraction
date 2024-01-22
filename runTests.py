import array
import numpy as np
from extractStructuresWithMC import extractStructuresMC
import unittest

#---------------------------------------------------------------------#

# Author: Abhishek Harikrishnan
# Email: abhishek.harikrishnan@fu-berlin.de
# Last updated: 06-12-2023
# Test to ensure correct working of extractStructuresSerial python module

#---------------------------------------------------------------------#

class TestBoxCount(unittest.TestCase):
	
	'''
	
	This script tests the 3D extraction/segmentation algorithm with and
	without the Marching Cubes correction.  
	
	'''
	
	# Select file to run tests on
	_filenameRead = 'testData.bin'
	
	# What indexing does it use?
	# This refers to array order
	# See here for a full description: https://docs.oracle.com/cd/E19957-01/805-4940/z400091044d0/index.html
	# Python uses C-order indexing (or row-major order). For this set _zFastest = True
	_zFastest = True
	
	# Threshold value for testing
	_threshVal = 47

	# Set additional parameters

	# Writes information relating to percolation analysis.
	_writePercolationData = False

	# Writes structure location information.
	# Useful for visualization purposes
	_writeNeighborInformation = True

	# Set data related parameters
	xlen = 200 
	ylen = 328
	zlen = 234

	# Set precision of binary data. 'f' is 32-bit floating point data.
	# For others, see here: https://docs.python.org/3/library/array.html
	precision = 'f'

	# Use array to read data. This is fast for binary files.
	data = array.array(precision)
	fr = open(_filenameRead, 'rb')
	data.fromfile(fr, (xlen*ylen*zlen))
	fr.close()

	# Convert to numpy array
	data = np.array(data, dtype = np.float32)
	
	# Print out details
	_verbose = False
	
	def test_3d_without_MC(self):
	
		_marchingCubesExt = False
		print('\nRunning 3D test case without MC...')
		
		extractStructuresObj = extractStructuresMC(self._threshVal, self.data, \
		self.xlen, self.ylen, self.zlen, self._zFastest, self._verbose, \
		self._writeNeighborInformation, self._writePercolationData, _marchingCubesExt)
		structureGrid = extractStructuresObj.extract()
		
		print('Expected number of structures at threshold 47:', 1)
		print('NOTE: structure 0 is empty space')
		
		# -1 is included because len also counts empty space as one structure
		self.assertEqual(len(np.unique(structureGrid))-1, 1)
	
	def test_3d_with_MC(self):
		
		# Use Marching Cubes neighbor correction (more computation power required)
		_marchingCubesExt = True
		print('\nRunning 3D test case with MC...')
		
		extractStructuresObj = extractStructuresMC(self._threshVal, self.data, \
		self.xlen, self.ylen, self.zlen, self._zFastest, self._verbose, \
		self._writeNeighborInformation, self._writePercolationData, _marchingCubesExt)
		structureGrid = extractStructuresObj.extract()
		
		print('Expected number of structures at threshold 47:', 31)
		print('NOTE: structure 0 is empty space')
		
		# -1 is included because len also counts empty space as one structure
		self.assertEqual(len(np.unique(structureGrid))-1, 31)

if __name__ == '__main__':
	
	# Run unit tests.
	# If everything checks out, it'll print OK.
	
    unittest.main()
