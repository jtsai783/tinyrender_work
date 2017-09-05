#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "model.h"

Model::Model(const char *filename) : verts_(), faces_() {
    std::ifstream in;
    in.open (filename, std::ifstream::in);
    if (in.fail()) return;
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            Vec3f v;
            for (int i=0;i<3;i++) iss >> v.raw[i];
                // std::cout << v;

            verts_.push_back(v);
        } else if (!line.compare(0, 3, "vt ")) {
            iss >> trash >> trash;
            Vec3f v;
            for (int i=0;i<3;i++) iss >> v.raw[i];
            // std::cout << v;
            texture_verts_.push_back(v);
        } else if (!line.compare(0, 3, "vn ")) {
            iss >> trash >> trash;
            Vec3f v;
            for (int i=0;i<3;i++) iss >> v.raw[i];
            // std::cout << v;
            normal_verts_.push_back(v);
        }else if (!line.compare(0, 2, "f ")) {
            std::vector<int> f;
            // std::vector<int> ft;
            int vn, idx, vt;
            iss >> trash;
            while (iss >> idx >> trash >> vt >> trash >> vn) {
                idx--; // in wavefront obj all indices start at 1, not zero
                vt--;
                vn--;
                f.push_back(idx);
                f.push_back(vt);
                f.push_back(vn);
            }
            faces_.push_back(f);
            // texture_faces_.push_back(ft);
        }
    }
    std::cerr << "# v# " << verts_.size() << " f# "  << faces_.size() << std::endl;

    diffuse.read_tga_file("obj/african_head_diffuse.tga");
    diffuse.flip_vertically();

    normal_map.read_tga_file("obj/african_head_nm_tangent.tga");
    normal_map.flip_vertically();

    spec.read_tga_file("obj/african_head_spec.tga");
    spec.flip_vertically();
}

Model::~Model() {
}

int Model::nverts() {
    return (int)verts_.size();
}

int Model::nfaces() {
    return (int)faces_.size();
}

std::vector<int> Model::face(int idx) {
    return faces_[idx];
}

// std::vector<int> Model::texture_face(int idx) {
//     return texture_faces_[idx];
// }

Vec3f Model::vert(int i) {
    return verts_[i];
}

Vec3f Model::texture_vert(int i) {
    return texture_verts_[i];
}

Vec3f Model::normal_vert(int i) {
    return normal_verts_[i];
}
