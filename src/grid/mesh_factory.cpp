#include "grid/mesh_factory.h"

std::array<Vector3, 4> getTetVertexPositions() {
    Vector3 p0{0, std::sqrt(6)/4, 0};
    Vector3 p1{-0.5, -std::sqrt(6)/12, -std::sqrt(3)/6};
    Vector3 p2{0.5, -std::sqrt(6)/12, -std::sqrt(3)/6};
    Vector3 p3{0, -std::sqrt(6)/12, std::sqrt(3)/3};
    return {p0, p1, p2, p3};
}


std::vector<std::vector<size_t>> getTetFaceIndices(){
    return {{0, 1, 2}, {0, 3, 1}, {0, 2, 3}, {1, 3, 2}};
}