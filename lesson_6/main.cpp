#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "geometry.h"
#include <iostream>
#include "model.h"
#include "matrix.h"
#include "our_gl.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
Model *model = NULL;
TGAImage *texture_img = NULL;

const int width  = 2000;
const int height = 2000;

const Vec3f camera(1, 1, 3);
const Vec3f center(0, 0, 0);
const Vec3f up(0, 1, 0);
const Vec3f light_dir(1,1,1);

Vec3f cross2(Vec3f a, Vec3f b){
    float x = a.y * b.z - a.z * b.y;
    float y = a.z * b.x - a.x * b.z;
    float z = a.x * b.y - a.y * b.x;
    return Vec3f(x,y,z);
}

Vec3f transform_v(Matrix m, Vec3f v){
    Matrix v_m = Matrix(4,1);
    v_m(0,0) = v.x;
    v_m(1,0) = v.y;
    v_m(2,0) = v.z;
    v_m(3,0) = 0.0;
    v_m = m * v_m;
    return Vec3f(v_m(0,0),v_m(1,0),v_m(2,0));
}

Vec3f transform_p(Matrix m, Vec3f v){
    Matrix v_m = Matrix(4,1);
    v_m(0,0) = v.x;
    v_m(1,0) = v.y;
    v_m(2,0) = v.z;
    v_m(3,0) = 1.0;
    v_m = m * v_m;
    v_m = v_m / v_m(3,0);
    return Vec3f(v_m(0,0),v_m(1,0),v_m(2,0));
}

struct Shader : public IShader {
    Vec3f texture_coords[3];
    Matrix cob;
    Matrix cob_IT;

    Matrix M;
    Matrix M_IT;

    Vec3f triangle_normal;
    Vec3f vertex_normal[3];
    Vec3f verts[3];
    Vec3f original_normal[3];
    Vec3f original_verts[3];

    Matrix TB;
    Matrix TBN;

    virtual Vec3f vertex(int nface, int nthvert, TGAImage &image){
        std::vector<int> face = model->face(nface);
        Vec3f v = model->vert(face[nthvert * 3]);
        Vec3f thisVert = transform_p(cob, v);
        texture_coords[nthvert] = model->texture_vert(face[nthvert * 3 + 1]);
        Vec3f normal = model->normal_vert(face[nthvert * 3 + 2]);
        vertex_normal[nthvert] = transform_v(M_IT, normal);
        verts[nthvert] = thisVert;
        original_verts[nthvert] = v;
        if(nthvert == 2){
            Vec3f tri_edge1 = verts[1] - verts[0];
            Vec3f tri_edge2 = verts[2] - verts[0];

            Vec3f orig_tri_edge1 = original_verts[1] - original_verts[0];
            Vec3f orig_tri_edge2 = original_verts[2] - original_verts[0];

            triangle_normal = cross2(orig_tri_edge1, orig_tri_edge2);
            triangle_normal = transform_v(M_IT, triangle_normal);
            triangle_normal = triangle_normal.normalize();

            Vec3f text_edge1 = texture_coords[1] - texture_coords[0];
            Vec3f text_edge2 = texture_coords[2] - texture_coords[0];

            Matrix tri_m(2,3);
            Matrix uv_m(2,2);

            tri_m(0,0) = tri_edge1.x;
            tri_m(0,1) = tri_edge1.y;
            tri_m(0,2) = tri_edge1.z;
            tri_m(1,0) = tri_edge2.x;
            tri_m(1,1) = tri_edge2.y;
            tri_m(1,2) = tri_edge2.z;

            uv_m(0,0) = text_edge1.x;
            uv_m(0,1) = text_edge1.y;
            uv_m(1,0) = text_edge2.x;
            uv_m(1,1) = text_edge2.y;

            TB = uv_m.inverse() * tri_m;

        }
        return thisVert;
    }

    virtual bool fragment(float *bc, TGAColor &color){

        //get color
        Vec3f texture = texture_coords[0] * bc[0] + texture_coords[1] * bc[1] + texture_coords[2] * bc[2];

        float texture_x = texture.x * (model->diffuse).get_width();
        float texture_y = texture.y * (model->diffuse).get_height();
        color = (model->diffuse).get(roundf(texture_x), roundf(texture_y));

        //get normal color
        TGAColor normal_color = (model->normal_map).get(roundf(texture_x), roundf(texture_y));
        Vec3f normal_map_vec = Vec3f(normal_color.r, normal_color.g, normal_color.b);
        normal_map_vec = normal_map_vec * (1.0/255.0) * 2.0 - Vec3f(1.0, 1.0, 1.0);
        normal_map_vec = normal_map_vec.normalize();

        Matrix tangent_normal_matrix(1, 3);
        tangent_normal_matrix(0,0) = normal_map_vec.x;
        tangent_normal_matrix(0,1) = normal_map_vec.y;
        tangent_normal_matrix(0,2) = normal_map_vec.z;

        //interpolate normal
        Vec3f n = vertex_normal[0] * bc[0] + vertex_normal[1] * bc[1] + vertex_normal[2] * bc[2];


        // Vec3f n = triangle_normal;
        

        Vec3f t(TB(0, 0), TB(0, 1), TB(0, 2));
        Vec3f b(TB(1, 0), TB(1, 1), TB(1, 2));

        t = t.normalize();
        b = b.normalize();
        n = n.normalize();

        TBN(0,0) = t.x;
        TBN(0,1) = t.y;
        TBN(0,2) = t.z;

        TBN(1,0) = b.x;
        TBN(1,1) = b.y;
        TBN(1,2) = b.z;

        TBN(2,0) = n.x;
        TBN(2,1) = n.y;
        TBN(2,2) = n.z;

        Matrix model_space_normal = tangent_normal_matrix * TBN;

        Vec3f model_space_normal_vec = Vec3f(model_space_normal(0,0), model_space_normal(0,1), model_space_normal(0,2));

        model_space_normal_vec = model_space_normal_vec.normalize();


        Vec3f normal_color_vec = (model_space_normal_vec + Vec3f(1.0, 1.0, 1.0)) * 0.5 * 255.0;
        TGAColor normal_vec_color(normal_color_vec.x, normal_color_vec.y, normal_color_vec.z, 255);

        // color = normal_vec_color;

        //specular
        Vec3f temp_light = transform_v(M, light_dir);
        temp_light = temp_light.normalize();

        Vec3f r = model_space_normal_vec * ( model_space_normal_vec * temp_light * 2.0 ) - temp_light;
        r = r.normalize();

        //get the spec value from spec map
        TGAColor spec_color = (model->spec).get(roundf(texture_x), roundf(texture_y));
        // float power = spec_color.val / 255.0;
        // std::cout << spec_color.val << std::endl;
        float spec = pow(std::max(r.z, 0.0f), spec_color.val);
        float diff = std::max(0.f,  model_space_normal_vec * temp_light);

        // float intensity = 1;

        color = TGAColor((float)color.r * (diff + 0.8 * spec) + 5.0, (float)color.g * (diff + 0.8 * spec) + 5.0, (float)color.b * (diff + 0.8 * spec) + 5.0, 255);

        return false;
    }
};


int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }


    TGAImage image(width, height, TGAImage::RGB);
    TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);

    // image = TGAImage(width, height, TGAImage::RGB);

    view(center, camera, up);
    proj(-1.0/(camera - center).norm());
    clip(width/8, height/8, width*3/4, height*3/4);


    Shader shader;
    shader.cob = Viewport * Projection * ModelView;
    shader.cob_IT = shader.cob.inverse().transpose();

    shader.M = Projection * ModelView;
    shader.M_IT = shader.M.inverse().transpose();

    shader.TB = Matrix(2,3);
    shader.TBN = Matrix(3,3);

    for (int i= 0; i < model->nfaces(); i++) { 
        std::vector<int> face = model->face(i);
        Vec3f screen_coords[3]; 
        for (int j=0; j<3; j++) { 
            screen_coords[j] = shader.vertex(i, j, image);
        }
        triangle(zbuffer, screen_coords, image, shader); 
    }

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    zbuffer.flip_vertically();
    image.write_tga_file("output.tga");
    zbuffer.write_tga_file("zbuffer.tga");

    delete model;
    return 0;
}

