#pragma once

#include "common/utils.h"

std::array<Vector3, 4> getTetVertexPositions();

std::vector<std::vector<size_t>> getTetFaceIndices();

// I wanna write a blender script that animates certain given aplitting of a cube to a bunch of tetrahedra. 
// The cube vertices are respectively
// A: (0,0,0) B: (1,0,0) C: (1,1,0) D: (0,1,0) E: (0,0,1) F: (1,0,1) G: (1,1,1) H: (0,1,1)
// ACFG, DACG, BACF, EFGA, HFCG