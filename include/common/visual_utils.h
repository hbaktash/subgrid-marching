#pragma once

#include "common/utils.h"

// Visualization helpers are only available when built with the Polyscope viewer
// (SUBGRID_POLYSCOPE_VIEWER=ON). In headless builds this header is a no-op.
#ifdef HAVE_POLYSCOPE

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

// 
void initialize_polyscope();

// vector<Vector3> of boundary polygons
void visualizeBoundaryCurves(
    const std::vector<std::vector<Vector3>>& polygons,
    bool closed,
    const glm::vec3& color,
    std::string name = "boundary curves",
    double radius = 0.001
);

// Register one point cloud per tet edge showing its intersection points.
void visualizeIntersections(
    const std::array<std::vector<double>, 6> edge_isect_ts,
    const std::array<Vector3,4>& tet_positions,
    const std::string& name_suffix = ""
);

#endif // HAVE_POLYSCOPE


