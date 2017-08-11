#include "tgaimage.h"
#include <iostream>

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color){
		bool steep = false;

		if(std::abs(x0 - x1) < std::abs(y0 - y1)){
			std::swap(x0, y0);
			std::swap(x1, y1);
			steep = true;
		}

		int dx = x1 - x0;
		int dy = y1 - y0;

		bool reflect = false;
		
		if((dx < 0 && dy > 0) || (dy < 0 && dx > 0)){
			std::swap(y0, y1);
			reflect = true;
		}

		if(std::abs(x1) < std::abs(x0)){
			std::swap(x0, x1);
			std::swap(y0, y1);
		}


		dx = x1 - x0;
		dy = y1 - y0;

		int a = -dy;
		int b = dx;

		int delta_d_ne = a + b + a + b;
		int delta_d_e = a + a;

		//calculate initial d
		int d = a + a + b;

		int x = x0;
		int y = y0;

		//paint initial point


		for(int x = x0; x <= x1; x++){

			int actual_x;
			int actual_y;

			actual_x = x;
			actual_y = y;
			if(reflect){
				actual_y = y1 - (y - y0);
			}
			if(steep){
				std::swap(actual_y, actual_x);
			}
			image.set(actual_x, actual_y, color);

			if(d < 0){
				y++;
				d += delta_d_ne;
			} else {
				d += delta_d_e;
			}

		}
}

int main(int argc, char** argv) {
	TGAImage image(200, 200, TGAImage::RGB);
	// image.set(52, 41, red);
	line(5, 5, 150, 150, image, white);
	line(5, 5, 150, 5, image, white);
	line(5, 5, 5 ,150, image, white);
	line(7, 130, 1200 ,8, image, red);
	// line(20, 13, 40, 80, image, red); 
	// line(80, 40, 13, 20, image, red);
	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("output.tga");
	return 0;
}
