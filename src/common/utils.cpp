#include "common/utils.h"


bool segment_triangle_intersection(
    const geometrycentral::Vector3& p0,
    const geometrycentral::Vector3& p1,
    const geometrycentral::Vector3& a,
    const geometrycentral::Vector3& b,
    const geometrycentral::Vector3& c,
    double& out_t)
{
    using namespace geometrycentral;
    const double EPS = 1e-12;

    Vector3 dir = p1 - p0;
    Vector3 edge1 = b - a;
    Vector3 edge2 = c - a;
    Vector3 pvec = cross(dir, edge2);
    double det = dot(edge1, pvec);

    if (std::abs(det) < EPS) return false;

    double invDet = 1.0 / det;
    Vector3 tvec = p0 - a;
    double u = dot(tvec, pvec) * invDet;
    if (u < -EPS || u > 1.0 + EPS) return false;

    Vector3 qvec = cross(tvec, edge1);
    double v = dot(dir, qvec) * invDet;
    if (v < -EPS || u + v > 1.0 + EPS) return false;

    double t_ray = dot(edge2, qvec) * invDet;

    double dir_len2 = dir.norm2();
    if (dir_len2 <= EPS) return false;
    double t_param = t_ray;
    if (t_param < -EPS || t_param > 1.0 + EPS) return false;

    out_t = std::max(0.0, std::min(1.0, t_param));
    return true;
}


std::vector<int> complement(std::initializer_list<int> inputs) {
    // Universe
    std::set<int> universe = {0, 1, 2, 3};
    
    // Remove provided elements
    for (int x : inputs) {
        universe.erase(x);
    }

    // Pack result
    return std::vector<int>(universe.begin(), universe.end());
}

std::vector<int> complement(const std::vector<int> &inputs){
    std::set<int> universe = {0, 1, 2, 3};
    for (int x : inputs){
        universe.erase(x);
    }
    std::vector<int> result(universe.begin(), universe.end());
    return result;
}

std::vector<std::pair<int, int>> all_pairs(std::vector<int> universe){
    std::vector<std::pair<int, int>> pairs;
    for (size_t i = 0; i < universe.size(); i++){
        for (size_t j = i+1; j < universe.size(); j++){
            pairs.push_back({universe[i], universe[j]});
        }
    }
    return pairs;
}

std::vector<std::pair<int, int>> all_pairs(std::array<int, 3> universe){
    return std::vector<std::pair<int, int>>{
        {universe[0], universe[1]},
        {universe[0], universe[2]},
        {universe[1], universe[2]}
    };
}

// normality checks
bool triangle_inequality(int a, int b, int c, 
                         int &res_a, int &res_b, int &res_c){
    int max = std::max({a, b, c});
    bool feasible = (a + b + c) >= 2 * max;
    res_a = 0; res_b = 0; res_c = 0;
    if (!feasible){
        if (a == max){
            res_a = a - b - c;
            res_b = 0; res_c = 0;
        }
        else if (b == max){
            res_b = b - a - c;
            res_a = 0; res_c = 0;
        }
        else { // c == max
            res_c = c - a - b;
            res_a = 0; res_b = 0;
        }
    }
    return feasible;
}


void center_and_normalize(SurfaceMesh& mesh, VertexPositionGeometry &geometry, double r){
    // compute center
    Vector3 center{0,0,0};
    for (Vertex v: mesh.vertices()){
        center += geometry.inputVertexPositions[v];
    }
    center /= double(mesh.nVertices());

    // translate to center
    for (Vertex v: mesh.vertices()){
        geometry.inputVertexPositions[v] -= center;
    }

    // compute scale
    double max_dist = 0.0;
    for (Vertex v: mesh.vertices()){
        double dist = geometry.inputVertexPositions[v].norm();
        if (dist > max_dist)
            max_dist = dist;
    }

    // normalize
    if (max_dist > 0.0){
        for (Vertex v: mesh.vertices()){
            geometry.inputVertexPositions[v] /= max_dist;
            geometry.inputVertexPositions[v] *= r;
        }
    }
}

// center to fit in the unit cube
void center_and_normalize(std::vector<Vector3> &points, double r){
    // compute center
    Vector3 center{0,0,0};
    for (const Vector3 &p: points){
        center += p;
    }
    center /= double(points.size());

    // translate to center
    for (Vector3 &p: points){
        p -= center;
    }

    // compute coordinate scale
    double max_coord = 0.0;
    for (const Vector3 &p: points){
        double temp_max_coord = std::max({std::abs(p.x), std::abs(p.y), std::abs(p.z)});
        if (temp_max_coord > max_coord)
            max_coord = temp_max_coord;
    }

    // normalize
    if (max_coord > 0.0){
        for (Vector3 &p: points){
            p /= max_coord;
            p *= r;
        }
    }
}

// Redraws an in-place progress bar on the current line (via '\r'); never emits a
// newline, so the caller is responsible for printing one after the final update.
void print_progress(double progress) {
    if (progress < 0.0) progress = 0.0;
    if (progress > 1.0) progress = 1.0;
    const int barWidth = 50;
    int pos = static_cast<int>(barWidth * progress);
    std::cout << "\r" << ANSI_FG_YELLOW << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << std::fixed << std::setprecision(1)
              << (progress * 100.0) << "%" << ANSI_RESET;
    std::cout.flush();
}

void triangulate_polygons(
    const std::vector<std::vector<size_t>> &input_polygons,
    std::vector<std::vector<size_t>> &output_triangles,
    bool remove_segments
){
    output_triangles.clear();
    for (const auto &poly : input_polygons){
        if (poly.size() < 3 && remove_segments) {
            continue; // skip degenerate polygons (reported in aggregate by the caller)
        }
        else if (poly.size() <= 3){
            output_triangles.push_back(poly); // already triangle
        }
        else {
            // fan triangulation
            for (size_t i = 1; i < poly.size() - 1; i++){
                output_triangles.push_back({poly[0], poly[i], poly[i+1]});
            }
        }
    }
}



std::vector<std::vector<size_t>> 
segments_to_polygons(
    const std::vector<std::pair<size_t, size_t>> &segments,
    const size_t pool_size
){
    // now build the polygons out of the pairs
    // build adjacency map
    std::map<size_t, std::vector<size_t>> cv_adjacency;
    for (const auto &pair: segments){
        size_t a = pair.first, b = pair.second;
        cv_adjacency[a].push_back(b);
        cv_adjacency[b].push_back(a);
    }
    // find polygons by walking the adjacency
    std::vector<bool> cv_visited(pool_size, false);
    std::vector<std::vector<size_t>> boundary_polygons;
    for (size_t start_idx = 0; start_idx < pool_size; start_idx++){
        if (cv_visited[start_idx]) continue;
        std::vector<size_t> polygon;
        size_t current_idx = start_idx;
        size_t previous_idx = INVALID_IND; // invalid
        while (true){
            polygon.push_back(current_idx);
            cv_visited[current_idx] = true;
            // find next
            const auto &neighbors = cv_adjacency[current_idx];
            size_t next_idx = INVALID_IND;
            for (size_t neighbor_idx : neighbors){
                if (cv_visited[neighbor_idx] != true){
                    next_idx = neighbor_idx;
                    break;
                }
            }
            if (next_idx == INVALID_IND) break; // no next found, close polygon loop
            previous_idx = current_idx;
            current_idx = next_idx;
        }
        boundary_polygons.push_back(polygon);
    }
    return boundary_polygons;
}


bool is_inside_tet(
    const Vector3 &p,
    const std::array<Vector3,4> &tet_positions,
    double eps
){
    // check if pd is inside tet
    auto signed_tet_vol = [](const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d) {
        return dot(cross((b - a),(c - a)),(d - a));
    };
    double from_pd_vol = 0.0;
    for (auto triplet: ALL_TET_TRIPLETS){
        from_pd_vol += std::abs(signed_tet_vol( tet_positions[triplet[0]], tet_positions[triplet[1]], tet_positions[triplet[2]], p));
    }
    double tet_vol = std::abs(signed_tet_vol( tet_positions[0], tet_positions[1], tet_positions[2], tet_positions[3]));
    if (from_pd_vol <= tet_vol + eps){
        return true; // inside
    }
    else return false;
}



void
pairss_to_tupless(
    const std::vector<std::vector<std::pair<size_t, size_t>>> &pairss,
    std::vector<std::vector<std::tuple<size_t, size_t>>> &tupless
){
    tupless.clear();
    for (const auto &pairs : pairss){
        std::vector<std::tuple<size_t, size_t>> tuples;
        for (const auto &p : pairs){
            tuples.push_back(std::make_tuple(p.first, p.second));
        }
        tupless.push_back(tuples);
    }
}



std::vector<size_t>
pairs_to_single_polygon(
    const std::vector<std::pair<size_t, size_t>> &segments,
    bool warn /*= false*/
){
    // Build adjacency (undirected), ignoring self-loops
    std::map<size_t, std::vector<size_t>> adjacency;
    for (const auto &seg : segments){
        size_t a = seg.first;
        size_t b = seg.second;
        if (a == b) continue; // ignore self-loops
        adjacency[a].push_back(b);
        adjacency[b].push_back(a);
    }

    // If no edges, no paths/loops exist
    if (adjacency.empty()){
        return {};
    }

    // Check degree condition: all degrees must be <= 2
    for (const auto &kv : adjacency){
        if (kv.second.size() > 2){
            throw std::runtime_error("pairs_to_single_polygon: vertex degree > 2; graph is not a collection of paths and loops.");
        }
    }

    // Find connected components, classify them as paths or loops,
    // and build an ordered list of vertices for each component.
    std::set<size_t> visited_global;
    std::vector<std::vector<size_t>> components; // only store components with size > 1

    for (const auto &kv : adjacency){
        size_t start = kv.first;
        if (visited_global.count(start)) continue;

        // Gather this connected component with DFS/BFS
        std::vector<size_t> stack;
        std::vector<size_t> comp_vertices;
        stack.push_back(start);
        visited_global.insert(start);

        while (!stack.empty()){
            size_t v = stack.back();
            stack.pop_back();
            comp_vertices.push_back(v);

            const auto &neighbors = adjacency[v];
            for (size_t nb : neighbors){
                if (!visited_global.count(nb)){
                    visited_global.insert(nb);
                    stack.push_back(nb);
                }
            }
        }

        if (comp_vertices.size() <= 1){
            // Single isolated vertex (or degenerate); ignore as requested
            continue;
        }

        // Classify component: path or loop
        int endpoint_count = 0;
        size_t endpoint_start = comp_vertices[0];
        for (size_t v : comp_vertices){
            size_t deg = adjacency[v].size();
            if (deg == 1){
                endpoint_count++;
                endpoint_start = v; // last seen endpoint; either is fine as a start
            }
            else if (deg == 0){
                // Should not occur in adjacency-based component; but guard anyway
                throw std::runtime_error("pairs_to_single_polygon: found vertex of degree 0 inside a component.");
            }
        }

        std::vector<size_t> ordered;
        ordered.reserve(comp_vertices.size());

        if (endpoint_count == 0){
            // Loop: all vertices should have degree 2
            size_t current = comp_vertices[0];
            size_t prev = std::numeric_limits<size_t>::max();
            for (size_t step = 0; step < comp_vertices.size(); ++step){
                ordered.push_back(current);
                const auto &neighbors = adjacency[current];
                if (neighbors.empty()){
                    throw std::runtime_error("pairs_to_single_polygon: loop component has vertex with degree 0.");
                }

                size_t next = std::numeric_limits<size_t>::max();
                for (size_t nb : neighbors){
                    if (nb != prev){
                        next = nb;
                        break;
                    }
                }

                if (next == std::numeric_limits<size_t>::max()){
                    // Should not happen for a proper loop
                    throw std::runtime_error("pairs_to_single_polygon: failed to advance along loop component.");
                }

                prev = current;
                current = next;
            }

            // After traversing, we should be back at the start vertex
            if (current != ordered.front()){
                throw std::runtime_error("pairs_to_single_polygon: loop component did not close properly.");
            }
        }
        else if (endpoint_count == 2){
            // Path: exactly two endpoints of degree 1
            size_t current = endpoint_start;
            size_t prev = std::numeric_limits<size_t>::max();
            for (size_t step = 0; step < comp_vertices.size(); ++step){
                ordered.push_back(current);
                const auto &neighbors = adjacency[current];

                size_t next = std::numeric_limits<size_t>::max();
                for (size_t nb : neighbors){
                    if (nb != prev){
                        next = nb;
                        break;
                    }
                }

                if (next == std::numeric_limits<size_t>::max()){
                    // Reached an endpoint; should only happen at the last step
                    if (step + 1 != comp_vertices.size()){
                        throw std::runtime_error("pairs_to_single_polygon: path component terminated prematurely.");
                    }
                    break;
                }

                prev = current;
                current = next;
            }

            if (ordered.size() != comp_vertices.size()){
                throw std::runtime_error("pairs_to_single_polygon: path component size mismatch.");
            }
        }
        else{
            // Any other number of endpoints is invalid for a simple path or loop
            throw std::runtime_error("pairs_to_single_polygon: component is neither a simple path nor a simple loop.");
        }

        components.push_back(ordered);
    }

    if (components.empty()){
        // no non-trivial paths or loops
        return {};
    }

    // If there are multiple components, warn and return the largest one
    if (components.size() > 1 && warn){
        log_warn("pairs_to_single_polygon: multiple components detected; returning the largest one.");
    }

    size_t best_idx = 0;
    size_t best_size = components[0].size();
    for (size_t i = 1; i < components.size(); ++i){
        if (components[i].size() > best_size){
            best_size = components[i].size();
            best_idx = i;
        }
    }

    return components[best_idx];
}


VoxelActivityTracker::VoxelActivityTracker(size_t grid_resolution)
    : grid_resolution_(grid_resolution) {}


size_t VoxelActivityTracker::tet_indices_to_voxel_index(const std::array<size_t, 4>& tet_indices) const {
    size_t points_per_axis = grid_resolution_ + 1;
    size_t slab = points_per_axis * points_per_axis;

    auto index_to_coords = [&](size_t idx) {
        size_t z = idx / slab;
        size_t rem = idx % slab;
        size_t y = rem / points_per_axis;
        size_t x = rem % points_per_axis;
        return std::array<size_t, 3>{x, y, z};
    };

    auto c0 = index_to_coords(tet_indices[0]);
    size_t min_x = c0[0], min_y = c0[1], min_z = c0[2];

    for (size_t i = 1; i < 4; ++i) {
        auto c = index_to_coords(tet_indices[i]);
        min_x = std::min(min_x, c[0]);
        min_y = std::min(min_y, c[1]);
        min_z = std::min(min_z, c[2]);
    }

    if (min_x >= grid_resolution_ || min_y >= grid_resolution_ || min_z >= grid_resolution_) {
        throw std::runtime_error("VoxelActivityTracker: tet indices map outside voxel grid.");
    }

    return (min_z * grid_resolution_ + min_y) * grid_resolution_ + min_x;
}


std::array<size_t, 3> VoxelActivityTracker::voxel_index_to_coords(size_t voxel_index) const {
    size_t voxels_per_axis = grid_resolution_;
    size_t plane = voxels_per_axis * voxels_per_axis;
    size_t z = voxel_index / plane;
    size_t rem = voxel_index % plane;
    size_t y = rem / voxels_per_axis;
    size_t x = rem % voxels_per_axis;
    return {x, y, z};
}


size_t VoxelActivityTracker::voxel_coords_to_index(size_t x, size_t y, size_t z) const {
    return (z * grid_resolution_ + y) * grid_resolution_ + x;
}


size_t VoxelActivityTracker::vertex_coords_to_index(size_t x, size_t y, size_t z) const {
    size_t points_per_axis = grid_resolution_ + 1;
    return (z * points_per_axis + y) * points_per_axis + x;
}


void VoxelActivityTracker::mark_voxel_vertices(size_t voxel_index, std::unordered_set<size_t>& vertex_set) {
    auto [x, y, z] = voxel_index_to_coords(voxel_index);
    for (size_t dz = 0; dz <= 1; ++dz) {
        for (size_t dy = 0; dy <= 1; ++dy) {
            for (size_t dx = 0; dx <= 1; ++dx) {
                vertex_set.insert(vertex_coords_to_index(x + dx, y + dy, z + dz));
            }
        }
    }
}


void VoxelActivityTracker::mark_active_voxel_neighbors(size_t voxel_index) {
    if (!active_voxels_.insert(voxel_index).second) {
        return;
    }

    active_or_neighbor_voxels_.insert(voxel_index);
    mark_voxel_vertices(voxel_index, active_vertices_);
    mark_voxel_vertices(voxel_index, active_or_neighbor_vertices_);

    auto [x, y, z] = voxel_index_to_coords(voxel_index);
    size_t voxels_per_axis = grid_resolution_;

    auto add_neighbor = [&](size_t nx, size_t ny, size_t nz) {
        size_t neighbor_index = voxel_coords_to_index(nx, ny, nz);
        active_or_neighbor_voxels_.insert(neighbor_index);
        mark_voxel_vertices(neighbor_index, active_or_neighbor_vertices_);
    };

    if (x > 0) add_neighbor(x - 1, y, z);
    if (x + 1 < voxels_per_axis) add_neighbor(x + 1, y, z);
    if (y > 0) add_neighbor(x, y - 1, z);
    if (y + 1 < voxels_per_axis) add_neighbor(x, y + 1, z);
    if (z > 0) add_neighbor(x, y, z - 1);
    if (z + 1 < voxels_per_axis) add_neighbor(x, y, z + 1);
}


void VoxelActivityTracker::report_tet(const std::array<size_t, 4>& tet_indices, bool has_polygon) {
    if (!has_polygon) {
        return;
    }

    size_t voxel_index = tet_indices_to_voxel_index(tet_indices);
    mark_active_voxel_neighbors(voxel_index);
}


VoxelActivityStats VoxelActivityTracker::stats() const {
    VoxelActivityStats result;
    result.active_voxels = active_voxels_.size();
    result.active_or_neighbor_voxels = active_or_neighbor_voxels_.size();
    result.active_vertices = active_vertices_.size();
    result.active_or_neighbor_vertices = active_or_neighbor_vertices_.size();
    result.active_voxel_equivalent_tets = result.active_voxels * 5;
    return result;
}


void str_to_pair(const std::string& str, size_t& a, size_t& b){
    size_t comma_pos = str.find(',');
    if (comma_pos == std::string::npos){
        throw std::runtime_error("str_to_pair: input string does not contain a comma.");
    }
    std::string first_part = str.substr(0, comma_pos);
    std::string second_part = str.substr(comma_pos + 1);
    try {
        a = std::stoul(first_part);
        b = std::stoul(second_part);
    }
    catch (const std::exception &e){
        throw std::runtime_error("str_to_pair: failed to convert string parts to size_t.");
    }
}


bool check_even_sum(
    const std::array<int, 6> e_isects
){
    // 01, 02, 03, 12, 13, 23
    // triangles: 01, 02, 12: 0,1,3
    //            01, 03, 13: 0,2,4
    //            02, 03, 23: 1,2,6
    //            12, 13, 23: 3,4,5
    if ((e_isects[0] + e_isects[1] + e_isects[3]) % 2 != 0) return false;
    if ((e_isects[0] + e_isects[2] + e_isects[4]) % 2 != 0) return false;
    if ((e_isects[1] + e_isects[2] + e_isects[5]) % 2 != 0) return false;
    if ((e_isects[3] + e_isects[4] + e_isects[5]) % 2 != 0) return false;
    return true;
}


bool check_triangle_ineq(
    const std::array<int, 6> e_isects
){
    // 01, 02, 03, 12, 13, 23
    // triangles: 01, 02, 12: 0,1,3
    //            01, 03, 13: 0,2,4
    //            02, 03, 23: 1,2,6
    //            12, 13, 23: 3,4,5
    if ((e_isects[0] + e_isects[1] + e_isects[3]) < 2*std::max({e_isects[0], e_isects[1], e_isects[3]})) return false;
    if ((e_isects[0] + e_isects[2] + e_isects[4]) < 2*std::max({e_isects[0], e_isects[2], e_isects[4]})) return false;
    if ((e_isects[1] + e_isects[2] + e_isects[5]) < 2*std::max({e_isects[1], e_isects[2], e_isects[5]})) return false;
    if ((e_isects[3] + e_isects[4] + e_isects[5]) < 2*std::max({e_isects[3], e_isects[4], e_isects[5]})) return false;
    return true;
}



vector<double> 
generate_uniform_isect_ts_single_edge(
    int edge_int, double t_max
){
    vector<double> isect_ts;
    for (int i = 0; i < edge_int; i++){
        isect_ts.push_back(static_cast<double>(i+1) * t_max / (double)(edge_int + 1));
    }
    return isect_ts;
}