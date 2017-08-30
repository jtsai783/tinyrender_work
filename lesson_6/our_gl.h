#ifndef __OUR_GL_H__
#define __OUR_GL_H__

#include "geometry.h"
#include "tgaimage.h"
#include "matrix.h"

extern Matrix Viewport;
extern Matrix ModelView;
extern Matrix Projection;

void view(Vec3f center, Vec3f camera, Vec3f up);

void clip(int x, int y, int width, int height);

void proj(float coeff);

struct IShader {
    virtual ~IShader();
    virtual Matrix vertex(Vec3f v, int nthvert) = 0;
    virtual bool fragment(Vec3f bar, TGAColor &color) = 0;
};

// triangle(texture_coords, zbuffer, screen_coords, image, texture_img, light_dir, normal_coords); 
void triangle(int *zbuffer, Matrix *pts, TGAImage &image, IShader &shader, int width);

#endif