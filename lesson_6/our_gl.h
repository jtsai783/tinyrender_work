#ifndef __OUR_GL_H__
#define __OUR_GL_H__

#include "geometry.h"
#include "tgaimage.h"

class Matrix;

Matrix view(Vec3f center, Vec3f camera, Vec3f up);

Matrix clip(int x, int y, int width, int height);

struct IShader {
    virtual ~IShader();
    virtual Vec3i vertex(int iface, int nthvert) = 0;
    virtual bool fragment(Vec3f bar, TGAColor &color) = 0;
};

// triangle(texture_coords, zbuffer, screen_coords, image, texture_img, light_dir, normal_coords); 
void triangle(Vec3f *texture_coords, float *zbuffer, Vec3f *pts, TGAImage &image, TGAImage &texture_img, Vec3f light_dir, Vec3f *normal_coords, int width);

#endif