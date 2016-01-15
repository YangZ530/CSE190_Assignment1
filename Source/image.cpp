#include "image.h"
#include "bmp.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>


/**
 * Image
 **/
Image::Image (int width_, int height_)
{
    assert(width_ > 0);
    assert(height_ > 0);

    width           = width_;
    height          = height_;
    num_pixels      = width * height;
    pixels          = new Pixel[num_pixels];
    sampling_method = IMAGE_SAMPLING_POINT;

    assert(pixels != NULL);
}


Image::Image (const Image& src)
{
    width           = src.width;
    height          = src.height;
    num_pixels      = width * height;
    pixels          = new Pixel[num_pixels];
    sampling_method = IMAGE_SAMPLING_POINT;

    assert(pixels != NULL);
    memcpy(pixels, src.pixels, src.width * src.height * sizeof(Pixel));
}


Image::~Image ()
{
    delete [] pixels;
    pixels = NULL;
}

/*
void Image::AddNoise (double factor)
{

}
*/

void Image::Brighten (double factor)
{
  /* Your Work Here  (section 3.2.1 of assignment)*/
	if (factor < 0){
		printf("Error: brightness factor cannot be negative");
		return;
	}
	for (int i = 0; i < num_pixels; i++){
		pixels[i].r = ComponentScale(pixels[i].r, factor);
		pixels[i].g = ComponentScale(pixels[i].g, factor);
		pixels[i].b = ComponentScale(pixels[i].b, factor);
	}
}


void Image::ChangeContrast(double factor)
{
	/* Your Work Here (section 3.2.2) */
	Component avg_grey = 0;
	Pixel avgGrey;
	for (int i = 0; i < num_pixels; i++){
		avg_grey += pixels[i].Luminance();
	}
	avg_grey = avg_grey / num_pixels;
	avgGrey = Pixel(avg_grey, avg_grey, avg_grey);
	if (factor > 0){
		for (int i = 0; i < num_pixels; i++){
			pixels[i].r = ComponentLerp(avg_grey, pixels[i].r, factor);
			pixels[i].g = ComponentLerp(avg_grey, pixels[i].g, factor);
			pixels[i].b = ComponentLerp(avg_grey, pixels[i].b, factor);
		}
	}
	else{
		for (int i = 0; i < num_pixels; i++){
			pixels[i].r = ComponentLerp(avg_grey, 255 - pixels[i].r, -factor);
			pixels[i].g = ComponentLerp(avg_grey, 255 - pixels[i].g, -factor);
			pixels[i].b = ComponentLerp(avg_grey, 255 - pixels[i].b, -factor);
		}
	}
}


void Image::ChangeSaturation(double factor)
{
  /* Your Work Here (section 3.2.3) */
	//Pixel *copy;
	//copy = new Pixel[num_pixels];
	Pixel p;
	for (int i = 0; i < num_pixels; i++){
		//copy[i] = Pixel(pixels[i].Luminance(), pixels[i].Luminance(), pixels[i].Luminance());
		p = Pixel(pixels[i].Luminance(), pixels[i].Luminance(), pixels[i].Luminance());
		pixels[i] = PixelLerp(p, pixels[i], factor);
	}
	//delete[] copy;
}

void Image::ChangeGamma(double factor)
{
  /* Your Work Here (section 3.2.4) */
	if (factor < 0){
		printf("Error: gamma factor cannot be negative");
		return;
	}
	for (int i = 0; i < num_pixels; i++){
		pixels[i].r = floor(pow(pixels[i].r / 256.0, 1.0 / factor) * 256.0);
		pixels[i].g = floor(pow(pixels[i].g / 256.0, 1.0 / factor) * 256.0);
		pixels[i].b = floor(pow(pixels[i].b / 256.0, 1.0 / factor) * 256.0);
	}
}

Image* Image::Crop(int x, int y, int w, int h)
{
  /* Your Work Here (section 3.2.5) */
	if (x < 0 || y < 0 || w < 0 || h < 0){
		printf("Error: crop cannot take negative parameters");
		return NULL;
	}
	Image crop = Image(w, h);
	for (int i = 0; i < w; i++){
		for (int j = 0; j < h; j++){
			crop.pixels[i * j] = pixels[(i + x)*(j + y)];
		}
	}
	return &crop;
}

/*
void Image::ExtractChannel(int channel)
{
  // For extracting a channel (R,G,B) of image.  
  // Not required for the assignment
}
*/

void Image::Quantize (int nbits)
{
  /* Your Work Here (Section 3.3.1) */
	if (nbits < 1 || nbits > 8){
		printf("Error: quantization bits shoule be in range 1-8");
		return;
	}
	for (int i = 0; i < num_pixels; i++){
		pixels[i].r = floor(floor(pixels[i].r * pow(2.0, nbits) / 256.0) * 255.0 / (pow(2.0, nbits) - 1));
		pixels[i].g = floor(floor(pixels[i].g * pow(2.0, nbits) / 256.0) * 255.0 / (pow(2.0, nbits) - 1));
		pixels[i].b = floor(floor(pixels[i].b * pow(2.0, nbits) / 256.0) * 255.0 / (pow(2.0, nbits) - 1));
	}
}


void Image::RandomDither (int nbits)
{
  /* Your Work Here (Section 3.3.2) */
	if (nbits < 1 || nbits > 8){
		printf("Error: quantization bits shoule be in range 1-8");
		return;
	}
	for (int i = 0; i < num_pixels; i++){
		pixels[i].r = floor(floor(pixels[i].r * pow(2.0, nbits) / 256.0 + (rand() % 2 - 0.5)) * 255.0 / (pow(2.0, nbits) - 1));
		pixels[i].g = floor(floor(pixels[i].g * pow(2.0, nbits) / 256.0 + (rand() % 2 - 0.5)) * 255.0 / (pow(2.0, nbits) - 1));
		pixels[i].b = floor(floor(pixels[i].b * pow(2.0, nbits) / 256.0 + (rand() % 2 - 0.5)) * 255.0 / (pow(2.0, nbits) - 1));
	}
}


/* Matrix for Bayer's 4x4 pattern dither. */
/* uncomment its definition if you need it */

/*
static int Bayer4[4][4] =
{
    {15, 7, 13, 5},
    {3, 11, 1, 9},
    {12, 4, 14, 6},
    {0, 8, 2, 10}
};


void Image::OrderedDither(int nbits)
{
  // For ordered dithering
  // Not required for the assignment
}

*/

/* Error-diffusion parameters for Floyd-Steinberg*/
const double
    ALPHA = 7.0 / 16.0,
    BETA  = 3.0 / 16.0,
    GAMMA = 5.0 / 16.0,
    DELTA = 1.0 / 16.0;

void Image::FloydSteinbergDither(int nbits)
{
  /* Your Work Here (Section 3.3.3) */
	if (nbits < 1 || nbits > 8){
		printf("Error: quantization bits shoule be in range 1-8");
		return;
	}
	for (int i = 0; i < num_pixels; i++){
		Pixel prev = pixels[i];
		pixels[i].r = floor(floor(pixels[i].r * pow(2.0, nbits) / 256.0) * 255.0 / (pow(2.0, nbits) - 1));
		pixels[i].g = floor(floor(pixels[i].g * pow(2.0, nbits) / 256.0) * 255.0 / (pow(2.0, nbits) - 1));
		pixels[i].b = floor(floor(pixels[i].b * pow(2.0, nbits) / 256.0) * 255.0 / (pow(2.0, nbits) - 1));

		//caculate error
		int er = (int)prev.r - (int)pixels[i].r;
		int eg = (int)prev.g - (int)pixels[i].g;
		int eb = (int)prev.b - (int)pixels[i].b;

		//apply error
		//right
		if ((i+1) % width != 0 && i+1 <num_pixels){
			pixels[i + 1].r = floor(pixels[i + 1].r + er*ALPHA);
			pixels[i + 1].g = floor(pixels[i + 1].g + eg*ALPHA);
			pixels[i + 1].b = floor(pixels[i + 1].b + eb*ALPHA);
		}
		else if ((i + 1) % width == 0){
			pixels[i - width + 1].r += er*ALPHA;
			pixels[i - width + 1].g += eg*ALPHA;
			pixels[i - width + 1].b += eb*ALPHA;
		}
		//bottom left
		if (i + width < num_pixels && i % width != 0){
			pixels[i + width - 1].r = floor(pixels[i + width - 1].r + er*BETA);
			pixels[i + width - 1].g = floor(pixels[i + width - 1].g + eg*BETA);
			pixels[i + width - 1].b = floor(pixels[i + width - 1].b + eb*BETA);
		}
		else if (i + width >= num_pixels && i % width != 0){
			pixels[i % width - 1].r = floor(pixels[i % width - 1].r + er*BETA);
			pixels[i % width - 1].g = floor(pixels[i % width - 1].g + eg*BETA);
			pixels[i % width - 1].b = floor(pixels[i % width - 1].b + eb*BETA);
		}
		else if (i + width < num_pixels && i % width == 0){
			pixels[i + width*2 - 1].r = floor(pixels[i + width*2 - 1].r + er*BETA);
			pixels[i + width*2 - 1].g = floor(pixels[i + width*2 - 1].g + eg*BETA);
			pixels[i + width*2 - 1].b = floor(pixels[i + width*2 - 1].b + eb*BETA);
		}
		else if (i + width >= num_pixels && i % width == 0){
			pixels[i % width + width - 1].r = floor(pixels[i % width + width - 1].r + er*BETA);
			pixels[i % width + width - 1].g = floor(pixels[i % width + width - 1].g + eg*BETA);
			pixels[i % width + width - 1].b = floor(pixels[i % width + width - 1].b + eb*BETA);
		}
		//bottom
		if (i + width < num_pixels){
			pixels[i + width].r = floor(pixels[i + width].r + er*GAMMA);
			pixels[i + width].g = floor(pixels[i + width].g + eg*GAMMA);
			pixels[i + width].b = floor(pixels[i + width].b + eb*GAMMA);
		}
		else if (i + width >= num_pixels){
			pixels[i % width].r = floor(pixels[i % width].r + er*GAMMA);
			pixels[i % width].g = floor(pixels[i % width].g + eg*GAMMA);
			pixels[i % width].b = floor(pixels[i % width].b + eb*GAMMA);
		}
		//bottom right
		if ((i + 1) % width != 0 && i + width < num_pixels){
			pixels[i + width + 1].r = floor(pixels[i + width + 1].r + er*DELTA);
			pixels[i + width + 1].g = floor(pixels[i + width + 1].g + eg*DELTA);
			pixels[i + width + 1].b = floor(pixels[i + width + 1].b + eb*DELTA);
		}
		else if ((i + 1) % width == 0 && i + width < num_pixels){
			pixels[i + 1].r = floor(pixels[i + 1].r + er*DELTA);
			pixels[i + 1].g = floor(pixels[i + 1].g + eg*DELTA);
			pixels[i + 1].b = floor(pixels[i + 1].b + eb*DELTA);
		}
		else if ((i + 1) % width != 0 && i + width >= num_pixels){
			pixels[i % width + 1].r = floor(pixels[i % width + 1].r + er*DELTA);
			pixels[i % width + 1].g = floor(pixels[i % width + 1].g + eg*DELTA);
			pixels[i % width + 1].b = floor(pixels[i % width + 1].b + eb*DELTA);
		}
		else if ((i + 1) % width == 0 && i + width >= num_pixels){
			pixels[0].r = floor(pixels[0].r + er*DELTA);
			pixels[0].g = floor(pixels[0].g + eg*DELTA);
			pixels[0].b = floor(pixels[0].b + eb*DELTA);
		}
	}
	/*
	for (int i = 0; i < width; i++){
		for (int j = 0; j < height; j++){
			Pixel prev = pixels[i*j];
			pixels[i + j*width].r = floor(floor(pixels[i + j*width].r * pow(2.0, nbits) / 256.0) * 255.0 / (pow(2.0, nbits) - 1));
			pixels[i + j*width].g = floor(floor(pixels[i + j*width].g * pow(2.0, nbits) / 256.0) * 255.0 / (pow(2.0, nbits) - 1));
			pixels[i + j*width].b = floor(floor(pixels[i + j*width].b * pow(2.0, nbits) / 256.0) * 255.0 / (pow(2.0, nbits) - 1));
			Pixel curr = pixels[i*j];

			//caculate error
			Component er = prev.r - curr.r;
			Component eg = prev.g - curr.g;
			Component eb = prev.b - curr.b;
			Pixel err = Pixel(er, eg, eb);

			//apply error
			if (i + 1 <width)
				pixels[i+1 + j*width] = pixels[i + 1+j*width] + err * ALPHA;
			if (j + 1 < height){
				if (i-1 >0)
					pixels[(i - 1)+(j + 1)*width] = pixels[(i - 1)+(j + 1)*width] + err * BETA;
				pixels[i+(j + 1)*width] = pixels[i+(j + 1)*width] + err * GAMMA;
				if (i+1<width)
					pixels[(i + 1)+(j + 1)*width] = pixels[(i + 1)+(j + 1)*width] + err * DELTA;
			}
		}
		
	}
	*/
}

void ImageComposite(Image *bottom, Image *top, Image *result)
{
  // Extra Credit (Section 3.7).
  // This hook just takes the top image and bottom image, producing a result
  // You might want to define a series of compositing modes as OpenGL does
  // You will have to use the alpha channel here to create Mattes
  // One idea is to composite your face into a famous picture
}

void Image::Convolve(int *filter, int n, int normalization, int absval) {
  // This is my definition of an auxiliary function for image convolution 
  // with an integer filter of width n and certain normalization.
  // The absval param is to consider absolute values for edge detection.
  
  // It is helpful if you write an auxiliary convolve function.
  // But this form is just for guidance and is completely optional.
  // Your solution NEED NOT fill in this function at all
  // Or it can use an alternate form or definition
}

void Image::Blur(int n)
{
  /* Your Work Here (Section 3.4.1) */
}

void Image::Sharpen() 
{
  /* Your Work Here (Section 3.4.2) */
}

void Image::EdgeDetect(int threshold)
{
  /* Your Work Here (Section 3.4.3) */
}


Image* Image::Scale(int sizex, int sizey)
{
  /* Your Work Here (Section 3.5.1) */
  return NULL ;
}

void Image::Shift(double sx, double sy)
{
  /* Your Work Here (Section 3.5.2) */
}


/*
Image* Image::Rotate(double angle)
{
  // For rotation of the image
  // Not required in the assignment
  // But you can earn limited extra credit if you fill it in
  // (It isn't really that hard) 

    return NULL;
}
*/


void Image::Fun()
{
    /* Your Work Here (Section 3.6) */
}


Image* ImageMorph (Image* I0, Image* I1, int numLines, Line* L0, Line* L1, double t)
{
  /* Your Work Here (Section 3.7) */
  // This is extra credit.
  // You can modify the function definition. 
  // This definition takes two images I0 and I1, the number of lines for 
  // morphing, and a definition of corresponding line segments L0 and L1
  // t is a parameter ranging from 0 to 1.
  // For full credit, you must write a user interface to join corresponding 
  // lines.
  // As well as prepare movies 
  // An interactive slider to look at various morph positions would be good.
  // From Beier-Neely's SIGGRAPH 92 paper

    return NULL;
}


/**
 * Image Sample
 **/
void Image::SetSamplingMethod(int method)
{
  // Sets the filter to use for Scale and Shift
  // You need to implement point sampling, hat filter and mitchell

    assert((method >= 0) && (method < IMAGE_N_SAMPLING_METHODS));
    sampling_method = method;
}

Pixel Image::Sample (double u, double v, double sx, double sy)
{
  // To sample the image in scale and shift
  // This is an auxiliary function that it is not essential you fill in or 
  // you may define it differently.
  // u and v are the floating point coords of the points to be sampled.
  // sx and sy correspond to the scale values. 
  // In the assignment, it says implement MinifyX MinifyY MagnifyX MagnifyY
  // separately.  That may be a better way to do it.
  // This hook is primarily to get you thinking about that you have to have 
  // some equivalent of this function.

  if (sampling_method == IMAGE_SAMPLING_POINT) {
    // Your work here
  }

  else if (sampling_method == IMAGE_SAMPLING_HAT) {
    // Your work here
  }

  else if (sampling_method == IMAGE_SAMPLING_MITCHELL) {
    // Your work here
  }

  else {
    fprintf(stderr,"I don't understand what sampling method is used\n") ;
    exit(1) ;
  }

  return Pixel() ;
}

