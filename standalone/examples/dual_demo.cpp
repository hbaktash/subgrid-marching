// Minimal demo + global-merge snippet for the standalone dual port.
//   c++ -std=c++17 -I../cpp dual_demo.cpp -o dual_demo && ./dual_demo
#include "subgrid_mt.hpp"

#include <iostream>
#include <unordered_map>

using namespace subgrid_mt;

int main() {
    // One tet with a few intersections on its edges, plus a surface normal per
    // intersection (the dual QEF solve needs them; here we just make some up).
    array<Vec3, 4> pos = {Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}};
    array<size_t, 4> gidx = {0, 1, 2, 3};
    array<vector<double>, 6> ts;   // edges (0,1),(0,2),(0,3),(1,2),(1,3),(2,3)
    array<vector<Vec3>, 6> nrm;    // parallel to ts
    // A single loop cutting off vertex 0 (crosses the 3 edges incident to it).
    ts[0] = {0.3}; nrm[0] = {Vec3{1, 1, 1}};
    ts[1] = {0.3}; nrm[1] = {Vec3{1, 1, 1}};
    ts[2] = {0.3}; nrm[2] = {Vec3{1, 1, 1}};

    DualResult r = subgrid_dual(pos, gidx, ts, nrm, /*reg_alpha=*/1e-3, /*project_duals=*/false);
    std::cout << "boundary-loop vertices: " << r.vertices.size()
              << ", loops (faces): " << r.faces.size()
              << ", dual points: " << r.dual_positions.size()
              << ", non_normal: " << (r.non_normal ? "true" : "false") << "\n";

    // ---- Global cross-tet merge sketch (mirrors add_local_soup) ----
    // The dual boundary-loop vertices merge across tets exactly like the primal
    // ones (via signature_key). r.dual_positions[f] is the QEF point for loop f,
    // and r.faces_per_edge[f] gives the tet face {i,j,k} each loop edge crosses;
    // the downstream dual-contouring step uses those two to stitch dual points.
    std::vector<Vec3> global_v;
    std::unordered_map<std::string, size_t> key_to_gid;
    std::vector<size_t> local_to_global(r.vertices.size());
    for (size_t vi = 0; vi < r.vertices.size(); vi++) {
        std::string key = signature_key(r.signatures[vi]);
        if (key.empty()) {  // STANDALONE -> always fresh
            local_to_global[vi] = global_v.size();
            global_v.push_back(r.vertices[vi]);
        } else {
            auto it = key_to_gid.find(key);
            if (it != key_to_gid.end()) {
                local_to_global[vi] = it->second;
            } else {
                local_to_global[vi] = global_v.size();
                key_to_gid[key] = global_v.size();
                global_v.push_back(r.vertices[vi]);
            }
        }
    }
    std::cout << "merged global boundary vertices: " << global_v.size() << "\n";
    return 0;
}
