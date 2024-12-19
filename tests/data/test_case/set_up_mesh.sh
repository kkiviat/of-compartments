rm -r constant/polyMesh
blockMesh
mirrorMesh -dict /data/openfoam-10/run/bioreactor_cfd/configs/mirrorMeshDict.y -overwrite
mirrorMesh -dict /data/openfoam-10/run/bioreactor_cfd/configs/mirrorMeshDict.z -overwrite
transformPoints "Rz=90"
transformPoints "scale=(2710.6878788298623 10842.751515319449 2710.6878788298623)"
