# Extraction (segmentation) of structures* from 3D scalar fields

***structure** - set of adjacent points satisfying a threshold (see Moisy & Jiménez 2004, JFM)

This code is an extension of the neighbor scanning procedure of Moisy & Jiménez, 2004 [1] to extract structures from 3D scalar fields such as vorticity magnitude, Q-criterion etc. This type of partitioninng in the computer graphics or computer vision community is termed 'segmentation'. For the classic segmentation techniques, see the review paper of Pal & Pal, 1993 [2].

## Neighbor scanning procedure with marching cubes correction

While the neighbor scanning procedure yields good results at large thresholds (where structures tend to be spaced apart), it groups together closely spaced structures at low threshold values. Upon visualization, it is easily observable that these structures are separate. When these structures are used for further post-processing, such as geometry classification with local parameters such as Shape Index and Curvedness, it will result in an inaccurate classification of the geometry. To mitigate this issue, a correction with the marching cubes algorithm (see Lorensen & Cline, 1987 [3]) is applied. For details on the implementation, please see Harikrishnan et al., 2021 [4].

## Installation

The code is available as a package from PyPI: https://pypi.org/project/extractstructuresMC/

```
pip install extractstructuresMC
```

## Example - Extracting structures from the given test data

Download testData.bin from the GitHub website first.

```
import numpy as np
import array
from extractStructuresWithMC import extractStructuresMC

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
_verbose = True

# Use Marching Cubes neighbor correction (more computation power required)
_marchingCubesExt = True

extractStructuresObj = extractStructuresMC(_threshVal, data, \
xlen, ylen, zlen, _zFastest, _verbose, \
_writeNeighborInformation, _writePercolationData, _marchingCubesExt)
structureGrid = extractStructuresObj.extract()

print('Expected number of structures at threshold 47: 31, result: ', len(np.unique(structureGrid))-1)
print('NOTE: structure 0 is empty space. len(np.unique(structureGrid)) counts it as well.')
```

## References

[1] Moisy, Frédéric, and Javier Jiménez. "Geometry and clustering of intense structures in isotropic turbulence." Journal of fluid mechanics 513 (2004): 111-133.

[2] Pal, Nikhil R., and Sankar K. Pal. "A review on image segmentation techniques." Pattern recognition 26.9 (1993): 1277-1294.

[3] Cline, Harvey E., et al. "3D reconstruction of the brain from magnetic resonance images using a connectivity algorithm." Magnetic Resonance Imaging 5.5 (1987): 345-352.

[4] Harikrishnan, Abhishek, et al. "Geometry and organization of coherent structures in stably stratified atmospheric boundary layers." arXiv preprint arXiv:2110.02253 (2021).
