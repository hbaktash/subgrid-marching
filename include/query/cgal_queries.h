#pragma once

#include "query/input_query_handler.h"

#include <memory>

// Exact mesh intersection queries backed by CGAL's EPECK AABB tree, selected via
// the --cgal CLI flag. Only available when built with -DSUBGRID_WITH_CGAL=ON.
//
// Predicates and intersection constructions are computed exactly and rounded to
// double only at the end, so this is the reliability reference for the fcpw-based
// MeshQueryHandler. Per-tet queries are serial; parallelization is future work.
//
// CGAL headers are kept out of this header (pImpl) so translation units that only
// need to construct the handler do not pull in CGAL/Boost/GMP.
#ifdef HAVE_CGAL

class CGALQueryHandler : public InputQueryHandler {
public:
    CGALQueryHandler(
        const std::vector<Vector3>& positions,
        const std::vector<std::vector<size_t>>& polygons
    );
    ~CGALQueryHandler() override;

    void query_intersections(
        const std::array<Vector3,4>& tet_positions,
        std::array<std::vector<double>,6>& edge_isect_ts,
        std::array<std::vector<Vector3>,6>& edge_isect_normals,
        bool useRobust = false,
        bool recordNormals = true
    ) override;

    void query_normal(const Vector3& q, Vector3& normal, bool verbose = false) override;

    bool is_sdf() const override { return false; }
    bool is_mesh() const override { return true; }
    bool has_mesh_data() const override { return true; }
    const std::vector<Vector3>& get_mesh_positions() const override { return positions; }
    const std::vector<std::vector<size_t>>& get_mesh_polygons() const override { return polygons; }

private:
    std::vector<Vector3> positions;
    std::vector<std::vector<size_t>> polygons;
    struct Impl;                    // owns the EPECK triangle soup + AABB tree
    std::unique_ptr<Impl> impl;
};

#endif // HAVE_CGAL
