// Minimal demo + global-merge snippet for the standalone primal port.
//   c++ -std=c++17 -I../cpp primal_demo.cpp -o primal_demo && ./primal_demo
#include "subgrid_mt.hpp"

#include <iostream>
#include <unordered_map>

using namespace subgrid_mt;

int main() {
    // One tet with a few intersections on its edges.
    array<Vec3, 4> pos = {Vec3{0, 0, 0}, Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}};
    array<size_t, 4> gidx = {0, 1, 2, 3};
    array<vector<double>, 6> ts;  // edges (0,1),(0,2),(0,3),(1,2),(1,3),(2,3)
    // A single loop cutting off vertex 0 (crosses the 3 edges incident to it):
    // every face is even-sum, so this is a normal curve -> one triangle.
    ts[0] = {0.3};
    ts[1] = {0.3};
    ts[2] = {0.3};

    PrimalResult r = subgrid_primal(pos, gidx, ts);
    std::cout << "vertices: " << r.vertices.size() << ", faces: " << r.faces.size()
              << ", non_even: " << (r.non_even ? "true" : "false") << "\n";

    // ---- Global cross-tet merge sketch (mirrors add_local_soup_comb) ----
    // Accumulate PrimalResult r for every active tet, then:
    std::vector<Vec3> global_v;
    std::vector<std::vector<size_t>> global_f;
    std::unordered_map<std::string, size_t> key_to_gid;
    auto add_result = [&](const PrimalResult& res) {
        std::vector<size_t> local_to_global(res.vertices.size());
        for (size_t vi = 0; vi < res.vertices.size(); vi++) {
            std::string key = signature_key(res.signatures[vi]);
            if (key.empty()) {  // STANDALONE -> always fresh
                local_to_global[vi] = global_v.size();
                global_v.push_back(res.vertices[vi]);
            } else {
                auto it = key_to_gid.find(key);
                if (it != key_to_gid.end()) {
                    local_to_global[vi] = it->second;
                } else {
                    local_to_global[vi] = global_v.size();
                    key_to_gid[key] = global_v.size();
                    global_v.push_back(res.vertices[vi]);
                }
            }
        }
        for (const auto& f : res.faces) {
            std::vector<size_t> gf;
            for (size_t idx : f) gf.push_back(local_to_global[idx]);
            global_f.push_back(gf);
        }
    };
    add_result(r);
    std::cout << "merged global vertices: " << global_v.size() << ", faces: " << global_f.size() << "\n";
    return 0;
}
