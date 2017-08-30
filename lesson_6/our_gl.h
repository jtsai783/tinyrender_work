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
    virtual Vec3f vertex(int nface, int nthvert) = 0;
    virtual bool fragment(float *bc, TGAColor &color) = 0;
};

void triangle(int *zbuffer, Vec3f *pts, TGAImage &image, IShader &shader, int width);

#endif