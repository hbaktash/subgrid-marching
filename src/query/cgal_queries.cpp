#include "query/cgal_queries.h"

#ifdef HAVE_CGAL

#include "common/utils.h"   // ALL_TET_PAIRS, log_info

#include <algorithm>
#include <iterator>
#include <utility>
#include <vector>

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include <CGAL/AABB_tree.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/version.h>
// CGAL 6 renamed these (the pre-6 headers still work but warn as deprecated).
#if CGAL_VERSION_NR >= 1060000000
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_triangle_primitive_3.h>
#else
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#endif

namespace {

// EPECK: exact predicates *and* exact constructions (the reference kernel).
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::FT        FT;
typedef K::Point_3   Point;
typedef K::Vector_3  CVec;
typedef K::Triangle_3 Triangle;
typedef K::Segment_3  Segment;
typedef std::vector<Triangle>::iterator TriIterator;
#if CGAL_VERSION_NR >= 1060000000
typedef CGAL::AABB_triangle_primitive_3<K, TriIterator> Primitive;
typedef CGAL::AABB_traits_3<K, Primitive>               AABBTraits;
#else
typedef CGAL::AABB_triangle_primitive<K, TriIterator>   Primitive;
typedef CGAL::AABB_traits<K, Primitive>                 AABBTraits;
#endif
typedef CGAL::AABB_tree<AABBTraits>                     Tree;
typedef boost::optional<Tree::Intersection_and_primitive_id<Segment>::Type>
    IntersectionResult;

// CGAL 6 returns the intersection alternative in a std::variant; 5.x used a
// boost::variant, and the two have no common accessor.
template <typename T, typename Variant>
inline const T* intersection_get(const Variant* v) {
#if CGAL_VERSION_NR >= 1060000000
    return std::get_if<T>(v);
#else
    return boost::get<T>(v);
#endif
}

} // namespace

// Owns the triangle soup and the AABB tree built over it. Primitives store
// iterators into `triangles`, so it must outlive the tree and not reallocate
// after the build (filled fully up front).
struct CGALQueryHandler::Impl {
    std::vector<Triangle> triangles;
    std::vector<Vector3>  tri_normals;   // geometric face normal (double), parallel to triangles
    std::unique_ptr<Tree> tree;

    // Surface crossings of a single segment a->b, written into out_ts (and
    // matching face normals into out_normals when record_normals). Outputs are
    // cleared first and returned sorted ascending with exactly-coincident hits
    // deduped, all t in [0, 1]. This is the per-edge unit query_intersections
    // runs over each tet edge.
    void edge_intersections(const Vector3& a, const Vector3& b, bool record_normals,
                            std::vector<double>& out_ts,
                            std::vector<Vector3>& out_normals);
};

void CGALQueryHandler::Impl::edge_intersections(
    const Vector3& a, const Vector3& b, bool record_normals,
    std::vector<double>& out_ts, std::vector<Vector3>& out_normals
) {
    out_ts.clear();
    out_normals.clear();

    Point op(a.x, a.y, a.z), ep(b.x, b.y, b.z);
    if (op == ep) return;   // degenerate zero-length edge

    const CVec d = ep - op;
    const FT   dd = d * d;   // squared length (> 0 here)

    std::vector<IntersectionResult> hits;
    tree->all_intersections(Segment(op, ep), std::back_inserter(hits));

    // Collect (exact t, triangle index). A Segment_3 result is a coplanar overlap
    // (measure zero under grid decorrelation, topologically ambiguous as a
    // crossing) and is ignored; only Point_3 transversal hits count.
    std::vector<std::pair<FT, size_t>> found;
    for (const auto& r : hits) {
        if (!r) continue;
        if (const Point* p = intersection_get<Point>(&(r->first))) {
            const FT t = ((*p - op) * d) / dd;   // exact projection parameter
            if (t >= FT(0) && t <= FT(1)) {
                size_t idx = static_cast<size_t>(std::distance(triangles.begin(), r->second));
                found.emplace_back(t, idx);
            }
        }
    }

    std::sort(found.begin(), found.end(),
              [](const std::pair<FT,size_t>& x, const std::pair<FT,size_t>& y) {
                  return x.first < y.first;
              });

    // Exact dedup: crossings exactly on a shared triangle edge/vertex are reported
    // once per incident triangle; collapse them to one crossing.
    for (size_t k = 0; k < found.size(); ++k) {
        if (k > 0 && found[k].first == found[k - 1].first) continue;
        double td = CGAL::to_double(found[k].first);
        td = std::min(1.0, std::max(0.0, td));
        out_ts.push_back(td);
        if (record_normals)
            out_normals.push_back(tri_normals[found[k].second]);
    }
}

CGALQueryHandler::CGALQueryHandler(
    const std::vector<Vector3>& pos,
    const std::vector<std::vector<size_t>>& polys
) : positions(pos), polygons(polys), impl(std::make_unique<Impl>()) {
    for (const auto& poly : polygons) {
        if (poly.size() < 3) continue;
        // Fan-triangulate (input is already triangulated by preprocessing, but
        // this stays correct for any convex polygon).
        const Vector3& a = positions[poly[0]];
        for (size_t k = 2; k < poly.size(); ++k) {
            const Vector3& b = positions[poly[k - 1]];
            const Vector3& c = positions[poly[k]];
            Triangle tri(Point(a.x, a.y, a.z), Point(b.x, b.y, b.z), Point(c.x, c.y, c.z));
            if (tri.is_degenerate()) continue;   // AABB primitives assume non-degenerate faces
            impl->triangles.push_back(tri);
            Vector3 n = cross(b - a, c - a);
            double len = n.norm();
            impl->tri_normals.push_back(len > 0.0 ? n / len : Vector3{0.0, 0.0, 0.0});
        }
    }
    impl->tree = std::make_unique<Tree>(impl->triangles.begin(), impl->triangles.end());
    impl->tree->build();
    log_info("built CGAL (EPECK) AABB tree for mesh with " + std::to_string(positions.size()) +
             " vertices and " + std::to_string(impl->triangles.size()) + " triangles.");
}

CGALQueryHandler::~CGALQueryHandler() = default;

void CGALQueryHandler::query_intersections(
    const std::array<Vector3,4>& tet_positions,
    std::array<std::vector<double>,6>& edge_isect_ts,
    std::array<std::vector<Vector3>,6>& edge_isect_normals,
    bool /*useRobust*/,
    bool recordNormals
) {
    for (int e = 0; e < 6; ++e) {
        const Vector3& a = tet_positions[ALL_TET_PAIRS[e].first];
        const Vector3& b = tet_positions[ALL_TET_PAIRS[e].second];
        impl->edge_intersections(a, b, recordNormals, edge_isect_ts[e], edge_isect_normals[e]);
    }
}

void CGALQueryHandler::query_normal(const Vector3& q, Vector3& normal, bool /*verbose*/) {
    auto pp = impl->tree->closest_point_and_primitive(Point(q.x, q.y, q.z));
    size_t idx = static_cast<size_t>(std::distance(impl->triangles.begin(), pp.second));
    normal = impl->tri_normals[idx];
}

#endif // HAVE_CGAL
