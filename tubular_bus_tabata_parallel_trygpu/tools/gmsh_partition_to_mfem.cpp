#include "mfem.hpp"
#include <mpi.h>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <cstdint>
#include <unordered_map>
#include <unordered_set>

using namespace mfem;

// ── 64-bit 哈希 (XOR-rotate, 输入: 任意数量的 uint64_t) ──
static uint64_t mix_hash(const uint64_t *p, int n) {
    uint64_t h = 0x9e3779b97f4a7c15ULL;
    for (int i = 0; i < n; i++)
        h ^= p[i] + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    return h;
}

// ── 3 个顶点的排序坐标 → 9 double → 64-bit hash ──
static uint64_t face_hash(const double *a, const double *b, const double *c3) {
    const double *pts[3] = {a, b, c3};
    if (pts[0] > pts[1]) std::swap(pts[0], pts[1]);  // 指针排序 (按地址, 但坐标不变)
    // 按字典序 (x→y→z) 排序三个点
    std::sort(pts, pts + 3, [](const double *p, const double *q) {
        for (int d = 0; d < 3; d++) {
            if (p[d] != q[d]) return p[d] < q[d];
        }
        return false;
    });
    double buf[9];
    for (int i = 0; i < 3; i++)
        for (int d = 0; d < 3; d++)
            buf[i*3+d] = pts[i][d];
    return mix_hash((const uint64_t*)buf, 9);
}

// ── 顶点哈希: 3 个 double 直接 XOR ──
static uint64_t vert_hash(const double *c) {
    return mix_hash((const uint64_t*)c, 3);
}

// ── 边哈希: 2 个排序顶点的 6 个 double ──
static uint64_t edge_hash(const double *ca, const double *cb) {
    const double *a = ca, *b2 = cb;
    if (a > b2) { auto t = a; a = b2; b2 = t; }  // 指针排序, fallback
    double buf[6];
    // 按字典序排序两顶点
    if (a[0] != b2[0]) { if (a[0] > b2[0]) { auto t=a; a=b2; b2=t; } }
    else if (a[1] != b2[1]) { if (a[1] > b2[1]) { auto t=a; a=b2; b2=t; } }
    else { if (a[2] > b2[2]) { auto t=a; a=b2; b2=t; } }
    for (int d = 0; d < 3; d++) { buf[d] = a[d]; buf[3+d] = b2[d]; }
    return mix_hash((const uint64_t*)buf, 6);
}

// ── 每个分界面携带的信息 ──
struct SharedFace {
    uint64_t hash;
    int local_v[3]; // 该面在本进程的局部顶点索引
};

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    double t_start = MPI_Wtime();
    int my_rank, num_procs;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &num_procs);

    // 参数: <mesh1> <mesh2> ... -o <prefix>
    std::vector<std::string> files;
    std::string prefix = "mesh.part";
    bool save_vertices = true;  // 默认生成 .part.vertices 文件
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            prefix = argv[++i];
        } else {
            files.push_back(argv[i]);
        }
    }

    if ((int)files.size() != num_procs) {
        if (my_rank == 0)
            printf("Usage: mpirun -np N %s <file0> <file1> ... -o <prefix>\n", argv[0]);
        MPI_Finalize();
        return 1;
    }

    // ── 1. 读取 GMSH 分片文件 ──
    const char *my_file = files[my_rank].c_str();
    Mesh mesh(my_file, 1, 1);  // generate_edges=1, refine=0 (分片已在GMSH中加密)

    if (my_rank == 0)
        printf("[Rank %d] Loaded %d elements, %d vertices from %s\n",
               my_rank, mesh.GetNE(), mesh.GetNV(), my_file);

    // ── 解析 GMSH 2.2 $Nodes 获取全局节点编号 ──
    // MFEM Mesh 不保留 GMSH 全局 ID，需单独解析
    std::vector<int> gmsh_global_tags(mesh.GetNV());
    std::vector<double> gmsh_coords(mesh.GetNV() * 3);

    {
        std::ifstream ifs(my_file);
        std::string line;
        while (std::getline(ifs, line) && line.find("$Nodes") == std::string::npos);
        std::getline(ifs, line);
        int total_nodes = std::stoi(line);

        // 读入所有 (x,y,z) -> gtag
        struct NodeRec { double x, y, z; int gtag; };
        std::vector<NodeRec> nrecs(total_nodes);
        for (int i = 0; i < total_nodes; i++) {
            std::getline(ifs, line);
            // 手动解析: "tag x y z" → 避免 istringstream 开销
            const char *s = line.c_str();
            char *end;
            int gtag = (int)strtol(s, &end, 10);
            double x = strtod(end, &end);
            double y = strtod(end, &end);
            double z = strtod(end, &end);
            nrecs[i] = {x, y, z, gtag};
        }

        // 按坐标排序，用于二分匹配
        std::sort(nrecs.begin(), nrecs.end(), [](const NodeRec& a, const NodeRec& b) {
            if (a.x != b.x) return a.x < b.x;
            if (a.y != b.y) return a.y < b.y;
            return a.z < b.z;
        });

        // MFEM 局部顶点匹配到 GMSH 全局标签
        for (int vi = 0; vi < mesh.GetNV(); vi++) {
            const double *v = mesh.GetVertex(vi);
            auto it = std::lower_bound(nrecs.begin(), nrecs.end(), NodeRec{v[0],v[1],v[2],0},
                [](const NodeRec& a, const NodeRec& b) {
                    if (a.x != b.x) return a.x < b.x;
                    if (a.y != b.y) return a.y < b.y;
                    return a.z < b.z;
                });
            if (it != nrecs.end() && it->x == v[0] && it->y == v[1] && it->z == v[2]) {
                gmsh_global_tags[vi] = it->gtag;
                gmsh_coords[vi*3+0] = v[0];
                gmsh_coords[vi*3+1] = v[1];
                gmsh_coords[vi*3+2] = v[2];
            }
        }
    }

    // ── MPI Allgather: 交换全局节点数据 (用于生成 .part.vertices) ──
    std::vector<int> all_node_counts(num_procs);
    int my_nv = mesh.GetNV();
    MPI_Allgather(&my_nv, 1, MPI_INT, all_node_counts.data(), 1, MPI_INT, comm);

    std::vector<int> node_displs(num_procs);
    std::vector<int> coord_counts(num_procs), coord_displs(num_procs);
    int total_nv = 0, total_coord = 0;
    for (int r = 0; r < num_procs; r++) {
        node_displs[r] = total_nv; total_nv += all_node_counts[r];
        coord_displs[r] = total_coord; total_coord += all_node_counts[r] * 3;
        coord_counts[r] = all_node_counts[r] * 3;
    }

    std::vector<int> all_gtags(total_nv);
    std::vector<double> all_coords(total_coord);

    MPI_Allgatherv(gmsh_global_tags.data(), my_nv, MPI_INT,
                   all_gtags.data(), all_node_counts.data(), node_displs.data(), MPI_INT, comm);
    MPI_Allgatherv(gmsh_coords.data(), my_nv * 3, MPI_DOUBLE,
                   all_coords.data(), coord_counts.data(), coord_displs.data(), MPI_DOUBLE, comm);

    // 生成 .part.vertices: 按全局标签排序, 去重
    if (save_vertices) {
        std::vector<std::pair<int,int>> tag_idx; // (global_tag, index in all_*)
        for (int i = 0; i < total_nv; i++)
            tag_idx.push_back({all_gtags[i], i});
        std::sort(tag_idx.begin(), tag_idx.end());

        std::string vert_file = prefix + ".vertices";
        std::ofstream vf(vert_file);
        vf.precision(14);
        int sd = 3, unique = 0;
        for (size_t i = 0; i < tag_idx.size(); i++) {
            if (i > 0 && tag_idx[i].first == tag_idx[i-1].first) continue; // 去重
            unique++;
        }
        vf << unique << " " << sd << "\n";
        for (size_t i = 0; i < tag_idx.size(); i++) {
            if (i > 0 && tag_idx[i].first == tag_idx[i-1].first) continue;
            int idx = tag_idx[i].second;
            vf << all_coords[idx*3+0] << " " << all_coords[idx*3+1] << " " << all_coords[idx*3+2] << "\n";
        }
        vf.close();
        if (my_rank == 0)
            printf("[Vertices] Written %d unique vertices to %s\n", unique, vert_file.c_str());
    }

    int dim = mesh.Dimension();
    int sdim = mesh.SpaceDimension();

    // ── 2. 统计每个面的出现次数 (hash → {elem, v0, v1, v2}) ──
    struct FaceRec { int elem; int v[3]; };
    std::unordered_map<uint64_t, std::vector<FaceRec>> face_map;
    int face_idx[4][3] = {{0,1,2},{0,1,3},{0,2,3},{1,2,3}};
    for (int e = 0; e < mesh.GetNE(); e++) {
        Element *el = mesh.GetElement(e);
        if (el->GetGeometryType() != Geometry::TETRAHEDRON) continue;
        const int *v = el->GetVertices();
        for (int fi = 0; fi < 4; fi++) {
            uint64_t h = face_hash(
                mesh.GetVertex(v[face_idx[fi][0]]),
                mesh.GetVertex(v[face_idx[fi][1]]),
                mesh.GetVertex(v[face_idx[fi][2]]));
            face_map[h].push_back({e, {v[face_idx[fi][0]], v[face_idx[fi][1]], v[face_idx[fi][2]]}});
        }
    }

    // ── 3. 物理边界 hash 集合 (用属性区分) ──
    std::unordered_set<uint64_t> phys_boundary;
    for (int i = 0; i < mesh.GetNBE(); i++) {
        Element *bel = mesh.GetBdrElement(i);
        int attr = bel->GetAttribute();
        if (attr < 1001 || attr > 1005 || bel->GetNVertices() != 3) continue;
        const int *v = bel->GetVertices();
        phys_boundary.insert(face_hash(
            mesh.GetVertex(v[0]), mesh.GetVertex(v[1]), mesh.GetVertex(v[2])));
    }

    // ── 4. 分区界面 = count==1 且不在物理边界中 ──
    std::vector<SharedFace> my_part_faces;
    int n_all_bdr = 0;
    for (auto &kv : face_map) {
        if (kv.second.size() != 1) continue;
        n_all_bdr++;
        if (phys_boundary.count(kv.first)) continue;
        auto &fr = kv.second[0];
        SharedFace sf;
        sf.hash = kv.first;
        sf.local_v[0] = fr.v[0]; sf.local_v[1] = fr.v[1]; sf.local_v[2] = fr.v[2];
        my_part_faces.push_back(sf);
    }

    if (my_rank == 0)
        printf("[Rank %d] all_boundary=%d, phys_boundary=%zu, part_bdr=%zu\n",
               my_rank, n_all_bdr, phys_boundary.size(), my_part_faces.size());

    if (my_rank == 0)
        printf("[Rank %d] partition boundary faces: %zu\n", my_rank, my_part_faces.size());

    // ── 5. 从分区界面提取顶点和边(去重) ──
    struct VertInfo { uint64_t hash; int local_v; };
    struct EdgeInfo { uint64_t hash; int local_a, local_b; };
    std::vector<VertInfo> my_verts;
    std::vector<EdgeInfo> my_edges;
    {
        std::unordered_set<uint64_t> vs, es;
        for (auto &sf : my_part_faces) {
            for (int vi = 0; vi < 3; vi++) {
                int lv = sf.local_v[vi];
                uint64_t vh = vert_hash(mesh.GetVertex(lv));
                if (vs.insert(vh).second) my_verts.push_back({vh, lv});
            }
            for (int ei = 0; ei < 3; ei++) {
                int a = sf.local_v[ei], b = sf.local_v[(ei+1)%3];
                if (a > b) std::swap(a, b);
                uint64_t eh = edge_hash(mesh.GetVertex(a), mesh.GetVertex(b));
                if (es.insert(eh).second) my_edges.push_back({eh, a, b});
            }
        }
    }

    // ── 6. MPI 交换: 面hash + 顶点hash + 边hash ──
    auto exchange = [&](const auto &items, auto gethash) {
        int n = (int)items.size();
        std::vector<uint64_t> h(n);
        for (int i = 0; i < n; i++) h[i] = gethash(items[i]);
        std::vector<int> cnt(num_procs), ds(num_procs);
        MPI_Allgather(&n, 1, MPI_INT, cnt.data(), 1, MPI_INT, comm);
        int tot = 0;
        for (int r = 0; r < num_procs; r++) { ds[r] = tot; tot += cnt[r]; }
        std::vector<uint64_t> all(tot);
        MPI_Allgatherv(h.data(), n, MPI_UINT64_T, all.data(), cnt.data(), ds.data(), MPI_UINT64_T, comm);
        std::vector<std::unordered_set<uint64_t>> rs(num_procs);
        for (int r = 0; r < num_procs; r++) {
            rs[r].reserve(cnt[r]);
            for (int i = 0; i < cnt[r]; i++) rs[r].insert(all[ds[r] + i]);
        }
        return rs;
    };
    auto fhs = exchange(my_part_faces, [](const SharedFace &s){ return s.hash; });
    auto vhs = exchange(my_verts, [](const VertInfo &v){ return v.hash; });
    auto ehs = exchange(my_edges, [](const EdgeInfo &e){ return e.hash; });

    // ── 8. 匹配面 + 直接检测顶点/边属于哪些rank ──
    std::map<int, std::vector<SharedFace>> nb_faces;
    for (auto &sf : my_part_faces)
        for (int r = 0; r < num_procs; r++)
            if (r != my_rank && fhs[r].count(sf.hash))
                { nb_faces[r].push_back(sf); break; }

    std::map<int, std::set<int>> vr; // local_v -> {neighbor_ranks}
    for (auto &vi : my_verts)
        for (int r = 0; r < num_procs; r++)
            if (r != my_rank && vhs[r].count(vi.hash))
                vr[vi.local_v].insert(r);

    std::map<std::pair<int,int>, std::set<int>> er;
    for (auto &ei : my_edges)
        for (int r = 0; r < num_procs; r++)
            if (r != my_rank && ehs[r].count(ei.hash))
                er[{ei.local_a, ei.local_b}].insert(r);

    // ── 9. 构建 group: rank集合 -> SharedData ──
    struct SD { std::set<int> V; std::set<std::pair<int,int>> E; std::set<std::array<int,3>> F; };
    std::map<std::set<int>, SD> gm;

    for (auto &kv : nb_faces) {
        std::set<int> rs = {kv.first};
        for (auto &sf : kv.second)
            gm[rs].F.insert({sf.local_v[0], sf.local_v[1], sf.local_v[2]});
    }
    for (auto &kv : vr) gm[kv.second].V.insert(kv.first);
    for (auto &kv : er) gm[kv.second].E.insert(kv.first);

    // 排序输出
    std::vector<std::pair<std::set<int>, SD>> all_groups(gm.begin(), gm.end());
    std::sort(all_groups.begin(), all_groups.end(),
        [](auto &a, auto &b) { return a.first.size() < b.first.size(); });

    // ── 10. 写出 MFEM 并行格式 ──
    {
        std::ostringstream fname;
        fname << prefix << "." << std::setfill('0') << std::setw(6) << my_rank;
        std::ofstream ofs(fname.str());
        ofs.precision(14);

        ofs << "MFEM mesh v1.2\n\n#\n# MFEM Geometry Types (see mesh/geom.hpp):\n#\n"
            << "# POINT       = 0\n# SEGMENT     = 1\n# TRIANGLE    = 2\n"
            << "# SQUARE      = 3\n# TETRAHEDRON = 4\n# CUBE        = 5\n"
            << "# PRISM       = 6\n# PYRAMID     = 7\n#\n\n";

        ofs << "dimension\n" << dim << "\n\n";

        // elements
        ofs << "elements\n" << mesh.GetNE() << '\n';
        for (int e = 0; e < mesh.GetNE(); e++) {
            Element *el = mesh.GetElement(e);
            ofs << el->GetAttribute() << ' ' << el->GetGeometryType();
            const int *v = el->GetVertices();
            for (int j = 0; j < el->GetNVertices(); j++)
                ofs << ' ' << v[j];
            ofs << '\n';
        }

        // boundary
        ofs << "\nboundary\n" << mesh.GetNBE() << '\n';
        for (int i = 0; i < mesh.GetNBE(); i++) {
            Element *bel = mesh.GetBdrElement(i);
            ofs << bel->GetAttribute() << ' ' << bel->GetGeometryType();
            const int *v = bel->GetVertices();
            for (int j = 0; j < bel->GetNVertices(); j++)
                ofs << ' ' << v[j];
            ofs << '\n';
        }

        // vertices
        ofs << "\nvertices\n" << mesh.GetNV() << '\n' << sdim << '\n';
        for (int i = 0; i < mesh.GetNV(); i++) {
            const double *v = mesh.GetVertex(i);
            ofs << v[0];
            for (int d = 1; d < sdim; d++)
                ofs << ' ' << v[d];
            ofs << '\n';
        }

        // parallel info
        ofs << "\nmfem_serial_mesh_end\n";

        // communication_groups
        int ngroups = 1 + (int)all_groups.size(); // group 0 + neighbors
        ofs << "\ncommunication_groups\nnumber_of_groups " << ngroups << "\n\n";
        ofs << "# number of entities in each group, followed by ranks in group\n";
        // Group 0: 仅本进程 (ProcToLProc 要求 GetJ()[0] == MyRank)
        ofs << "1 " << my_rank << '\n';
        // Neighbor groups: 包含 my_rank + 邻居rank集合
        for (size_t gi = 0; gi < all_groups.size(); gi++) {
            // 构建完整rank集合: my_rank + 邻居集合
            std::set<int> full_set = all_groups[gi].first;
            full_set.insert(my_rank);
            ofs << full_set.size();
            for (int r : full_set) ofs << ' ' << r;
            ofs << '\n';
        }

        // totals
        int tsv = 0, tse = 0, tsf = 0;
        for (auto &g : all_groups) {
            tsv += (int)g.second.V.size();
            tse += (int)g.second.E.size();
            tsf += (int)g.second.F.size();
        }
        ofs << "\ntotal_shared_vertices " << tsv << '\n';
        ofs << "total_shared_edges " << tse << '\n';
        ofs << "total_shared_faces " << tsf << '\n';
        ofs << "\n# group 0 has no shared entities\n";

        // per-group shared entities (每个组三项必须全部列出, 即使为0)
        for (size_t gi = 0; gi < all_groups.size(); gi++) {
            auto &sd = all_groups[gi].second;
            ofs << "\n# group " << (gi + 1) << '\n';
            ofs << "shared_vertices " << sd.V.size() << '\n';
            for (int v : sd.V) ofs << v << '\n';
            ofs << "\nshared_edges " << sd.E.size() << '\n';
            for (auto &e : sd.E) ofs << e.first << ' ' << e.second << '\n';
            ofs << "\nshared_faces " << sd.F.size() << '\n';
            for (auto &f : sd.F)
                ofs << "2 " << f[0] << ' ' << f[1] << ' ' << f[2] << '\n';
        }
        ofs << "\nmfem_mesh_end\n";
        ofs.close();
    }

    if (my_rank == 0) {
        printf("\nDone. Output: %s.*\n", prefix.c_str());
        printf("Neighbors:\n");
    }
    MPI_Barrier(comm);
    for (int r = 0; r < num_procs; r++) {
        if (r == my_rank) {
            printf("  Rank %d neighbors: ", my_rank);
            for (auto &g : all_groups) {
                printf("[");
                for (int r : g.first) printf("%d ", r);
                printf("] ");
            }
            printf("\n");
        }
        MPI_Barrier(comm);
    }

    double t_elapsed = MPI_Wtime() - t_start;
    double t_max;
    MPI_Reduce(&t_elapsed, &t_max, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
    if (my_rank == 0)
        printf("\nTotal time: %.3f s\n", t_max);

    MPI_Finalize();
    return 0;
}
