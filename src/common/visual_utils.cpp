#include "common/visual_utils.h"

#include <memory>

#ifdef HAVE_POLYSCOPE

void initialize_polyscope(){
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    polyscope::options::ssaaFactor = 2; // remove for speed
    polyscope::view::upDir = polyscope::UpDir::ZUp;
    polyscope::init();
}

// vector<Vector3> of boundary polygons
void visualizeBoundaryCurves(
    const std::vector<std::vector<Vector3>>& curves,
    bool closed,
    const glm::vec3& color,
    std::string name,    
    double radius
){
    std::vector<std::array<size_t,2>> curveEdges;
    std::vector<Vector3> curve_points;
    size_t offset = 0;
    for (const auto& curve : curves){
        size_t n = curve.size();
        for (size_t i = 0; i < n; ++i){
            curve_points.push_back(curve[i]);
            if (closed || i < n - 1){
                curveEdges.push_back({offset + i, offset + ((i + 1) % n)});
            }
        }
        offset += n;
    }
    polyscope::registerCurveNetwork(name, curve_points, curveEdges)->setRadius(radius)->setColor(color)->setEnabled(true);
}

void visualizeIntersections(
    const std::array<std::vector<double>, 6> edge_isect_ts,
    const std::array<Vector3,4>& tet_positions,
    const std::string& name_suffix
){
    for (size_t e = 0; e < 6; ++e){
        std::string name = "Edge " + std::to_string(ALL_TET_PAIRS[e].first) + "-" + std::to_string(ALL_TET_PAIRS[e].second) + " Isects" + name_suffix;
        if (edge_isect_ts[e].empty()) {
            if (polyscope::hasPointCloud(name))
                polyscope::removePointCloud(name);
            continue;
        }
        std::vector<Vector3> isect_points;
        Vector3 v0 = tet_positions[ALL_TET_PAIRS[e].first];
        Vector3 v1 = tet_positions[ALL_TET_PAIRS[e].second];
        for (double t : edge_isect_ts[e]){
            isect_points.push_back(v0 * (1.0 - t) + v1 * t);
        }
        polyscope::registerPointCloud(name, isect_points)->setPointColor(glm::vec3(0.7, 0.3, 0.4))->setPointRadius(.008)->setEnabled(true);
    }
}

#endif // HAVE_POLYSCOPE
