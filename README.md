# Extraction (segmentation) of structures* from 3D scalar fields

***structure** - set of adjacent points satisfying a threshold (see Moisy & Jiménez 2004, JFM)

This code is an extension of the neighbor scanning procedure of Moisy & Jiménez, 2004 [1] to extract structures from 3D scalar fields such as vorticity magnitude, Q-criterion etc. This type of partitioninng in the computer graphics or computer vision community is termed 'segmentation'. For the classic segmentation techniques, see the review paper of Pal & Pal, 1993 [2].

## Neighbor scanning and marching cubes correction

While the neighbor scanning procedure yields good results at large thresholds (where structures tend to be spaced apart), it groups together closely spaced structures at low threshold values. Upon visualization, it is easily observable that these structures are separate. When these structures are used for further post-processing, such as geometry classification with local parameters such as Shape Index and Curvedness, it will result in an inaccurate classification of the geometry. To mitigate this issue, a correction with the marching cubes algorithm (see Lorensen & Cline, 1987 [3]) is applied. For details on the implementation, please see Harikrishnan et al., 2021 [4].

## References

[1] Moisy, Frédéric, and Javier Jiménez. "Geometry and clustering of intense structures in isotropic turbulence." Journal of fluid mechanics 513 (2004): 111-133.

[2] Pal, Nikhil R., and Sankar K. Pal. "A review on image segmentation techniques." Pattern recognition 26.9 (1993): 1277-1294.

[3] Cline, Harvey E., et al. "3D reconstruction of the brain from magnetic resonance images using a connectivity algorithm." Magnetic Resonance Imaging 5.5 (1987): 345-352.

[4] Harikrishnan, Abhishek, et al. "Geometry and organization of coherent structures in stably stratified atmospheric boundary layers." arXiv preprint arXiv:2110.02253 (2021).
