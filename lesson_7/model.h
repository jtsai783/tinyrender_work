#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>
#include "geometry.h"
#include "tgaimage.h"

class Model {
private:
	std::vector<Vec3f> verts_;
	std::vector<Vec3f> texture_verts_;
	std::vector<Vec3f> normal_verts_;
	std::vector<std::vector<int> > faces_;
	// std::vector<std::vector<int> > texture_faces_;

public:
	Model(const char *filename);
	~Model();
	int nverts();
	int nfaces();
	Vec3f vert(int i);
	Vec3f texture_vert(int i);
	Vec3f normal_vert(int i);
	std::vector<int> face(int idx);
	TGAImage diffuse;
	TGAImage normal_map;
	TGAImage spec;
	// std::vector<int> texture_face(int idx);
};

#endif //__MODEL_H__
