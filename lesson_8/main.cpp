#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "geometry.h"
#include <iostream>
#include "model.h"
#include "matrix.h"
#include "our_gl.h"


TGAImage *aoimage = NULL;
const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
Model *model = NULL;
TGAImage *texture_img = NULL;

float *shadowBuffer = NULL;

const int width  = 800;
const int height = 800;

const Vec3f camera(-1, 0, 3);
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

struct DepthShader : public IShader{
    Vec3f verts[3];

    virtual Vec3f vertex(int nface, int nthvert, TGAImage &image){
        std::vector<int> face = model->face(nface);
        Vec3f v = model->vert(face[nthvert * 3]);
        Vec3f thisVert = transform_p(Viewport * Projection * ModelView, v);
        verts[nthvert] = thisVert;
        return thisVert;
    }

    virtual bool fragment(float *bc, TGAColor &color){
        Vec3f p = verts[0] * bc[0] + verts[1] * bc[1] + verts[2] * bc[2];
        // std::cout << p.z << std::endl;
        color = TGAColor(p.z, p.z, p.z, 255);
        return false;
    }
};

struct AOShader : public IShader{
    Vec3f verts[3];

    virtual Vec3f vertex(int nface, int nthvert, TGAImage &image){
        std::vector<int> face = model->face(nface);
        Vec3f v = model->vert(face[nthvert * 3]);
        Vec3f thisVert = transform_p(Viewport * Projection * ModelView, v);
        verts[nthvert] = thisVert;
        return thisVert;
    }

    virtual bool fragment(float *bc, TGAColor &color){
        color = TGAColor(0.0, 0.0, 0.0, 255.0);
        return false;
    }
};

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

    Matrix shadowTransform;

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
        Vec3f p = verts[0] * bc[0] + verts[1] * bc[1] + verts[2] * bc[2];
        Vec3f shadowP = transform_p(shadowTransform, p);
        int idx = int(shadowP.x) + int(shadowP.y) * width;
        float shadow = .3+.7*(shadowBuffer[idx] < shadowP.z + 43.34);
        // std::cout << shadow << std::endl;



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

        color = TGAColor((float)color.r * (diff + 0.8 * spec) * shadow + 5.0, (float)color.g * (diff + 0.8 * spec) * shadow + 5.0, (float)color.b * (diff + 0.8 * spec) * shadow + 5.0, 255);

        return false;
    }
};

float max_elevation_angle(float *buffer, int x, int y, float a){
    Vec2f heading(cos(a), sin(a));

    float x_f = (float)x;
    float y_f = (float)y;

    float x_f_old;
    float y_f_old;

    float max_angle = 0;

    float bufferHeight = buffer[x + y * width];
    // std::cout << bufferHeight << std::endl;
    for(float p = 0; p < 1000.0; p+=1.0){
        // std::cout << max_angle * 57.2958 << std::endl;
        x_f_old = x_f;
        y_f_old = y_f;
        x_f = x_f + heading.x;
        y_f = y_f + heading.y;
        if(x_f > width || x_f < 0 || y_f > height || y_f < 0) return max_angle;
        if((int)x_f - (int)x_f_old == 0 && (int)y_f - (int)y_f_old == 0) continue;
        //get the height
        float neighborHeight = buffer[(int)x_f + (int)y_f * width];
        // std::cout << "heighborheight" << neighborHeight << std::endl;
        if(neighborHeight > bufferHeight){

            float diff = neighborHeight - bufferHeight;

            float angle = atanf(diff/p);
            // std::cout << angle << std::endl;
            if(angle > max_angle){
                max_angle = angle;
            }
        }
    }
    return max_angle;
};

int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }


    TGAImage image(width, height, TGAImage::RGB);
    // TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);
    // TGAImage shadowBuffer(width, height, TGAImage::GRAYSCALE);
    TGAImage depth(width, height, TGAImage::RGB);
    TGAImage aoimage(width, height, TGAImage::RGB);

    TGAImage aodebugbuffer(width, height, TGAImage::GRAYSCALE);

    float *aobuffer = new float[width * height];
    float *zbuffer = new float[width*height];
    shadowBuffer   = new float[width*height];
    for (int i=width*height; --i; ) {
        aobuffer[i] = zbuffer[i] = shadowBuffer[i] = -std::numeric_limits<float>::max();
    }

    // image = TGAImage(width, height, TGAImage::RGB);

    view(center, light_dir, up);
    proj(0);
    clip(width/8, height/8, width*3/4, height*3/4);

    DepthShader depthShader;

    for (int i= 0; i < model->nfaces(); i++) { 
        std::vector<int> face = model->face(i);
        Vec3f screen_coords[3]; 
        for (int j=0; j<3; j++) { 
            screen_coords[j] = depthShader.vertex(i, j, image);
        }
        triangle(shadowBuffer, screen_coords, depth, depthShader); 
    }

    Matrix shadowView = Viewport * Projection * ModelView;

    view(center, camera, up);
    proj(-1.0/(camera - center).norm());
    clip(width/8, height/8, width*3/4, height*3/4);


    AOShader aoshader;

    for (int i= 0; i < model->nfaces(); i++) { 
        std::vector<int> face = model->face(i);
        Vec3f screen_coords[3]; 
        for (int j=0; j<3; j++) { 
            screen_coords[j] = aoshader.vertex(i, j, image);
        }
        triangle(aobuffer, screen_coords, aoimage, aoshader); 
    }

    for(int x = 0; x < width ; x++){
        for(int y = 0; y < height; y++){
            float aobufferVal = aobuffer[x + y * width];
            if(aobufferVal < 0){
                aobufferVal = 0.0;
            } else {
                aobufferVal = aobufferVal / 255.0;
            }
            aodebugbuffer.set(x, y, TGAColor(aobufferVal * 255, 255));
        }
    }

    for(int x = 0; x < width ; x++){
        for(int y = 0; y < height; y++){
    // int x = width / 2;
    // int y = height / 2;
            // if (aobuffer[x + y * width] < 0) continue;
            //emit 8 rays
            float total = 0;
            for(float a = 0; a < M_PI * 2 - 0.0001; a += M_PI/4){
            // float a = 135/57.2958;
                float max_angle = max_elevation_angle(aobuffer, x, y, a);
                // std::cout << a * 57.2958 << std::endl;
                // std::cout << max_angle << std::endl << std::endl;
                float diff = M_PI/2.0 - max_angle;
                // std::cout << "angle diff" << M_PI/2.0 << std::endl << std::endl;
                total += diff;
                // std::cout << "angle total " << total * 57.2958 << std::endl << std::endl;
            }

            total = total / 8.0;
            total = total / (M_PI/2.0);

            // std::cout << "angle total ave" << total * 57.2958 << std::endl << std::endl;
            // total = pow(total, 100.f);
            std::cout << total << std::endl;
            aoimage.set(x, y, TGAColor(total * 255, total * 255, total * 255, 255));
        }
    }



    // Shader shader;
    // shader.cob = Viewport * Projection * ModelView;
    // shader.cob_IT = shader.cob.inverse().transpose();

    // shader.M = Projection * ModelView;
    // shader.M_IT = shader.M.inverse().transpose();

    // shader.TB = Matrix(2,3);
    // shader.TBN = Matrix(3,3);

    // shader.shadowTransform = shadowView * shader.cob.inverse();

    // for (int i= 0; i < model->nfaces(); i++) { 
    //     std::vector<int> face = model->face(i);
    //     Vec3f screen_coords[3]; 
    //     for (int j=0; j<3; j++) { 
    //         screen_coords[j] = shader.vertex(i, j, image);
    //     }
    //     triangle(zbuffer, screen_coords, image, shader); 
    // }

    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    // zbuffer.flip_vertically();
    image.write_tga_file("output.tga");
    // zbuffer.write_tga_file("zbuffer.tga");

    depth.flip_vertically();
    depth.write_tga_file("shadowBuffer.tga"); 

    aoimage.flip_vertically();
    aoimage.write_tga_file("aobuffer.tga"); 

    aodebugbuffer.flip_vertically();
    aodebugbuffer.write_tga_file("aodebugbuffer.tga"); 

    delete model;
    return 0;
}

