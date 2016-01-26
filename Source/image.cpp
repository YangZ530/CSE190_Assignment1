#include "image.h"
#include "bmp.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <algorithm>    // std::max

/**
 * Image
 **/

# define PI 3.14159
# define WIDTH 2

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
	Image *crop = new Image(w, h);
	for (int i = 0; i < w; i++){
		for (int j = 0; j < h; j++){
			crop->pixels[i * w + j] = pixels[(i + x) * width + (j + y)];
		}
	}
	return crop;
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

		int nr = 0;
		int ng = 0;
		int nb = 0;

		//apply error
		//right
		if ((i+1) % width != 0 && i+1 <num_pixels){
			nr = floor(pixels[i + 1].r + er*ALPHA);
			ng = floor(pixels[i + 1].g + eg*ALPHA);
			nb = floor(pixels[i + 1].b + eb*ALPHA);
			pixels[i + 1].SetClamp(nr, ng, nb);
		}
		else if ((i + 1) % width == 0){
			nr = floor(pixels[i - width + 1].r + er*ALPHA);
			ng = floor(pixels[i - width + 1].g + eg*ALPHA);
			nb = floor(pixels[i - width + 1].b + eb*ALPHA);
			pixels[i - width + 1].SetClamp(nr, ng, nb);
		}
		//bottom left
		if (i + width < num_pixels && i % width != 0){
			nr = floor(pixels[i + width - 1].r + er*BETA);
			ng = floor(pixels[i + width - 1].g + eg*BETA);
			nb = floor(pixels[i + width - 1].b + eb*BETA);
			pixels[i + width - 1].SetClamp(nr, ng, nb);
		}
		else if (i + width >= num_pixels && i % width != 0){
			nr = floor(pixels[i % width - 1].r + er*BETA);
			ng = floor(pixels[i % width - 1].g + eg*BETA);
			nb = floor(pixels[i % width - 1].b + eb*BETA);
			pixels[i % width - 1].SetClamp(nr, ng, nb);
		}
		else if (i + width < num_pixels && i % width == 0){
			nr = floor(pixels[i + width*2 - 1].r + er*BETA);
			ng = floor(pixels[i + width*2 - 1].g + eg*BETA);
			nb = floor(pixels[i + width*2 - 1].b + eb*BETA);
			pixels[i + width * 2 - 1].SetClamp(nr, ng, nb);
		}
		else if (i + width >= num_pixels && i % width == 0){
			nr = floor(pixels[i % width + width - 1].r + er*BETA);
			ng = floor(pixels[i % width + width - 1].g + eg*BETA);
			nb = floor(pixels[i % width + width - 1].b + eb*BETA);
			pixels[i % width + width - 1].SetClamp(nr, ng, nb);
		}
		//bottom
		if (i + width < num_pixels){
			nr = floor(pixels[i + width].r + er*GAMMA);
			ng = floor(pixels[i + width].g + eg*GAMMA);
			nb = floor(pixels[i + width].b + eb*GAMMA);
			pixels[i + width].SetClamp(nr, ng, nb);
		}
		else if (i + width >= num_pixels){
			nr = floor(pixels[i % width].r + er*GAMMA);
			ng = floor(pixels[i % width].g + eg*GAMMA);
			nb = floor(pixels[i % width].b + eb*GAMMA);
			pixels[i % width].SetClamp(nr, ng, nb);
		}
		//bottom right
		if ((i + 1) % width != 0 && i + width < num_pixels){
			nr = floor(pixels[i + width + 1].r + er*DELTA);
			ng = floor(pixels[i + width + 1].g + eg*DELTA);
			nb = floor(pixels[i + width + 1].b + eb*DELTA);
			pixels[i + width + 1].SetClamp(nr, ng, nb);
		}
		else if ((i + 1) % width == 0 && i + width < num_pixels){
			nr = floor(pixels[i + 1].r + er*DELTA);
			ng = floor(pixels[i + 1].g + eg*DELTA);
			nb = floor(pixels[i + 1].b + eb*DELTA);
			pixels[i + 1].SetClamp(nr, ng, nb);
		}
		else if ((i + 1) % width != 0 && i + width >= num_pixels){
			nr = floor(pixels[i % width + 1].r + er*DELTA);
			ng = floor(pixels[i % width + 1].g + eg*DELTA);
			nb = floor(pixels[i % width + 1].b + eb*DELTA);
			pixels[i % width + 1].SetClamp(nr, ng, nb);
		}
		else if ((i + 1) % width == 0 && i + width >= num_pixels){
			nr = floor(pixels[0].r + er*DELTA);
			ng = floor(pixels[0].g + eg*DELTA);
			nb = floor(pixels[0].b + eb*DELTA);
			pixels[0].SetClamp(nr, ng, nb);
		}
	}
}

void ImageComposite(Image *bottom, Image *top, Image *result)
{
  // Extra Credit (Section 3.7).
  // This hook just takes the top image and bottom image, producing a result
  // You might want to define a series of compositing modes as OpenGL does
  // You will have to use the alpha channel here to create Mattes
  // One idea is to composite your face into a famous picture
}

void Image::Convolve(int *filter, int w, int normalization, int absval) {
  // This is my definition of an auxiliary function for image convolution 
  // with an integer filter of width w and certain normalization.
  // The absval param is to consider absolute values for edge detection.
  
  // It is helpful if you write an auxiliary convolve function.
  // But this form is just for guidance and is completely optional.
  // Your solution NEED NOT fill in this function at all
  // Or it can use an alternate form or definition
	
	Pixel *newPixels = new Pixel[num_pixels];
	int n = (w - 1) / 2;
	double red, green, blue;

	// iterate through pixels
	for (int b = 0; b < height; b++)
		for (int a = 0; a < width; a++ )		
		{
			red = 0;
			green = 0;
			blue = 0;

			// iterate through nearby pixels and convolve
			for (int i = b - n; i <= b + n; i++)
				for (int j = a - n; j <= a + n; j++ )				
				{
					double weight = (double)filter[(i - b + n) * w + (j - a + n)] / (double)normalization; // get normalized weight
					
					int x = j;
					int y = i;

					// reflected indexing					
					if (x < 0)
						x = -x;
					if (x >= width)
						x = width - (x - width) - 1;
					if (y < 0)
						y = -y;
					if (y >= height)
						y = height - (y - height) - 1;
						
						red += weight * (double)GetPixel(x, y).r;
						green += weight * (double)GetPixel(x, y).g;
						blue += weight * (double)GetPixel(x, y).b;
				}


			if (absval == 1)
			{
				red = abs(red);
				green = abs(green);
				blue = abs(blue);
			}

			newPixels[b*width + a].SetClamp(round(red), round(green), round(blue), GetPixel(a, b).a );
		}
	pixels = newPixels;
}

void Image::Blur(int n)
{
	if (n > 19 || n < 3 || n % 2 == 0)
	{
		printf("n must be odd and between 3 and 19");
		return;
	}

	double ratio = 0;
	double mean = (n - 1) / 2;
	double sigma = floor((float)n / 2) / 2;
	int kernel[19 * 19];
	int normalization = 0;

	ratio = (2 * sigma * sigma * PI) / exp(-(pow((float)0 - mean, 2) + pow((float)0 - mean, 2)) / (2 * pow(sigma, 2)));

	for (int y = 0; y < n; y++)
		for (int x = 0; x < n; x++)
		{
			kernel[y * n + x] = round(ratio * 1 / (2 * sigma * sigma * PI) * exp(-(pow((float)x - mean, 2) + pow((float)y - mean, 2)) / (2 * pow(sigma, 2))));
			normalization += kernel[y * n + x];
		}

	Convolve(kernel, n, normalization, 0);
}

void Image::Sharpen() 
{
	int kernel[9] = 
	{	
		-1, -2, -1,
		-2, 19, -2,
		-1, -2, -1 
	};
	Convolve( kernel, 3, 7, 0);
}

void Image::EdgeDetect(int threshold)
{
	Pixel *Gx = new Pixel[num_pixels];
	Pixel *Gy = new Pixel[num_pixels];
	Pixel *G = new Pixel[num_pixels];
	Pixel *original = new Pixel[num_pixels];

	int horizontal[9] =
	{
		-1, 0, 1,
		-2, 0, 2,
		-1, 0, 1
	};

	int vertical[9] =
	{
		1, 2, 1,
		0, 0, 0,
		-1, -2, -1
	};

	original = pixels;

	Convolve(horizontal, 3, 1, 1);
	Gx = pixels;
	pixels = original;
	Convolve(vertical, 3, 1, 1);
	Gy = pixels;
	for (int i = 0; i < num_pixels; i++)
	{
		double l;
		l = sqrt(pow(Gx[i].Luminance(), 2) + pow(Gy[i].Luminance(), 2));  // sqrt((r * 76 + g * 150 + b * 29) ^ 2 + (r * 76 + g * 150 + b * 29) ^ 2)		
		if ( l >= threshold)
			pixels[i].Set( 255, 255, 255);
		else
			pixels[i].Set(0, 0, 0);						
	}	
}


Image* Image::Scale(int sizex, int sizey)
{
	
	double sx = (double)sizex / width;
	double sy = (double)sizey / height;	
	Pixel *intermediatePixels = new Pixel[sizex * height];
	Pixel *newPixels = new Pixel[sizex * sizey];
	Image * newImg = new Image(sizex, sizey);

	for (int y = 0; y < height; y++)
		for (int x = 0; x < sizex; x++)
	{
		if (sx > 1)
			intermediatePixels[y * sizex + x] = MagnifyX( x, y, sx);
		else
			intermediatePixels[y * sizex + x] = MinifyX( x, y, sx);
	}

	pixels = intermediatePixels;
	width = sizex;

	for (int y = 0; y < sizey; y++)
		for (int x = 0; x < sizex; x++)
		{
			if (sy > 1)
				newPixels[y * sizex + x] = MagnifyY( x, y, sy);
			else
				newPixels[y * sizex + x] = MinifyY( x, y, sy);
		}
	
	newImg->pixels = newPixels;

	return newImg;
}

double hat(double x)
{
	return std::max(1 - abs(x), 0.0);
}

double mitchell(double n)
{
	double x = abs(n);
	if (x >= 0.0 && x < 1.0)
		return (7.0 * pow(x, 3.0) - 12.0 * pow(x, 2.0) + 16.0 / 3.0) / 6.0; // 1 / 6 * (7 * x ^ 3 - 12 * x ^ 2 + 16 / 3)
	else if (x >= 1.0 && x < 2.0)
		return (-7.0 / 3.0 * pow(x, 3.0) + 12.0 * pow(x, 2.0) - 20.0 * x + 32.0 / 3.0) / 6.0; // 1 / 6 * (-7 / 3 * x ^ 3 + 12 * x ^ 2 - 20 * x + 32 / 3)
	else
		return 0.0;
}

Pixel Image::MagnifyX(double a, double b, double s)
{
	double x = round(a / s);
	double y = b;
	double w = WIDTH / s;
	double normalization = 0.0;
	Pixel p;
	double red, green, blue;

	red = 0;
	green = 0;
	blue = 0;

	p.Set(0, 0, 0);

	switch (sampling_method)
	{
	case IMAGE_SAMPLING_POINT:
		return GetPixel(floor(a / s), b);
		break;
	case IMAGE_SAMPLING_HAT:
		for (double n = x - w; n <= x + w; n++)
		{
			double i = round(n);
			double weight = hat(i - a / s);

			if (i < 0)
				i = -i;
			if (i >= width)
				i = width - (i - width) - 1;

			if (weight != 0)
			{
				red += weight * GetPixel(i, y).r;
				green += weight * GetPixel(i, y).g;
				blue += weight * GetPixel(i, y).b;
				normalization += weight;
			}
		}

		p.SetClamp(round(red / normalization), round(green / normalization), round(blue / normalization));
		return p;
		break;
	case IMAGE_SAMPLING_MITCHELL:
		for (double n = x - w; n <= x + w; n++)
		{
			double i = round(n);
			double weight = mitchell(i - a / s);

			if (i < 0)
				i = -i;
			if (i >= width)
				i = width - (i - width) - 1;

			if (weight != 0)
			{
				red += weight * GetPixel(i, y).r;
				green += weight * GetPixel(i, y).g;
				blue += weight * GetPixel(i, y).b;
				normalization += weight;
			}
		}

		p.SetClamp(round(red / normalization), round(green / normalization), round(blue / normalization));
		return p;
		break;
	default:
		printf("sampling model unrecognize");
		break;
	}
}

Pixel Image::MinifyX(double a, double b, double s)
{
	double x = round(a / s);
	double y = b;
	double w = WIDTH / s;
	double normalization = 0.0;
	Pixel p;
	double red, green, blue;
	double ored, ogreen, oblue;

	ored = GetPixel(x, y).r;
	ogreen = GetPixel(x, y).g;
	oblue = GetPixel(x, y).b;

	red = 0;
	green = 0;
	blue = 0;

	p.Set(0, 0, 0);

	switch (sampling_method)
	{
	case IMAGE_SAMPLING_POINT:
		return GetPixel(round(a / s), b);
		break;
	case IMAGE_SAMPLING_HAT:		
		for (double n = x - w; n <= x + w; n++)
		{
			double i = round(n);
			double weight = hat(i * s - a);

			if (i < 0)
				i = -i;
			if (i >= width)
				i = width - (i - width) - 1;


			if (weight != 0)
			{
				red += weight * GetPixel(i, y).r;
				green += weight * GetPixel(i, y).g;
				blue += weight * GetPixel(i, y).b;
				normalization += weight;
			}
		}

		p.SetClamp(floor(red / normalization), floor(green / normalization), floor(blue / normalization));
		return p;
		break;
	case IMAGE_SAMPLING_MITCHELL:
		for (double n = x - w; n <= x + w; n++)
		{
			double i = round(n);
			double weight = mitchell(i * s - a);
			
			if (i < 0)
				i = -i;
			if (i >= width)
				i = width - (i - width) - 1;

			if (weight != 0)
			{
				red += weight * GetPixel(i, y).r;
				green += weight * GetPixel(i, y).g;
				blue += weight * GetPixel(i, y).b;
				normalization += weight;
			}
		}

		p.SetClamp(floor(red / normalization), floor(green / normalization), floor(blue / normalization));
		return p;
		break;
	default:
		printf("sampling model unrecognize");
		break;
	}
}

Pixel Image::MagnifyY(double a, double b, double s)
{
	double x = a;
	double y = round(b / s);
	double w = WIDTH / s;	
	double normalization = 0.0;
	Pixel p;
	double red, green, blue;

	red = 0;
	green = 0;
	blue = 0;

	p.Set(0, 0, 0);

	switch (sampling_method)
	{
	case IMAGE_SAMPLING_POINT:
		return GetPixel(a, floor(b / s));
		break;
	case IMAGE_SAMPLING_HAT:
		for (double n = y - w; n <= y + w; n++)
		{
			double i = round(n);
			double weight = hat(i - b / s);

			if (i < 0)
				i = -i;
			if (i >= height)
				i = height - (i - height) - 1;

			if (weight != 0)
			{
				red += weight * GetPixel(x, i).r;
				green += weight * GetPixel(x, i).g;
				blue += weight * GetPixel(x, i).b;
				normalization += weight;
			}
		}
		p.SetClamp(round(red / normalization), round(green / normalization), round(blue / normalization));
		return p;
		break;
	case IMAGE_SAMPLING_MITCHELL:
		for (double n = y - w; n <= y + w; n++)
		{
			double i = round(n);
			double weight = mitchell(i - b / s);

			if (i < 0)
				i = -i;
			if (i >= height)
				i = height - (i - height) - 1;

			if (weight != 0)
			{
				red += weight * GetPixel(x, i).r;
				green += weight * GetPixel(x, i).g;
				blue += weight * GetPixel(x, i).b;
				normalization += weight;
			}
		}
		p.SetClamp(round(red / normalization), round(green / normalization), round(blue / normalization));
		return p;
		break;
	default:
		printf("sampling model unrecognize");
		break;
	}
}

Pixel Image::MinifyY(double a, double b, double s)
{
	double x = a;
	double y = round(b / s);
	double w = WIDTH / s;
	double normalization = 0.0;
	Pixel p;
	double red, green, blue;

	red = 0;
	green = 0;
	blue = 0;

	p.Set(0, 0, 0);

	switch (sampling_method)
	{
	case IMAGE_SAMPLING_POINT:
		return GetPixel(a, round(b / s));
		break;
	case IMAGE_SAMPLING_HAT:		
		for (double n = y - w; n <= y + w; n++)
		{
			double i = round(n);
			double weight = hat(i * s - b);

			if (i < 0)
				i = -i;
			if (i >= height)
				i = height - (i - height) - 1;
			
			if (weight != 0)
			{
				red += weight * GetPixel(x, i).r;
				green += weight * GetPixel(x, i).g;
				blue += weight * GetPixel(x, i).b;
				normalization += weight;
			}
		}
		p.SetClamp(round(red / normalization), round(green / normalization), round(blue / normalization));
		return p;
		break;
	case IMAGE_SAMPLING_MITCHELL:
		
		for (double n = y - w; n <= y + w; n++)
		{
			double i = round(n);
			double weight = mitchell(i * s - b);

			if (i < 0)
				i = -i;
			if (i >= height)
				i = height - (i - height) - 1;
			
			if (weight != 0)
			{
				red += weight * GetPixel(x, i).r;
				green += weight * GetPixel(x, i).g;
				blue += weight * GetPixel(x, i).b;
				normalization += weight;
			}
		}
		p.SetClamp(round(red / normalization), round(green / normalization), round(blue / normalization));
		return p;
		
		break;
	default:
		printf("sampling model unrecognize");
		break;
	}
}

void Image::Shift(double sx, double sy)
{	
	Pixel *intermediatePixels = new Pixel[num_pixels];
	Pixel *newPixels = new Pixel[num_pixels];

	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
			intermediatePixels[y * width + x] = ShiftX(x, y, sx);

	pixels = intermediatePixels;

	for (int y = 0; y < height; y++)
		for (int x = 0; x < width; x++)
			newPixels[y * width + x] = ShiftY(x, y, sy);
	
	pixels = newPixels;
}

Pixel Image::ShiftX(double a, double b, double s)
{
	double x = round(a - s);
	double y = b;
	double w = WIDTH;
	Pixel p;
	double red, green, blue;

	red = 0;
	green = 0;
	blue = 0;

	p.Set(0, 0, 0);

	if (x < 0 || x >= width)
		return p;

	switch (sampling_method)
	{
	case IMAGE_SAMPLING_POINT:
		return GetPixel(x, y);
		break;
	case IMAGE_SAMPLING_HAT:
		for (double n = x - w; n <= x + w; n++)
		{
			double i = round(n);

			if (i < 0 || i >= width)
				continue;

			red += hat(i - a + s) * GetPixel(i, y).r;
			green += hat(i - a + s) * GetPixel(i, y).g;
			blue += hat(i - a + s) * GetPixel(i, y).b;

		}

		p.SetClamp(round(red), round(green), round(blue));
		return p;
		break;
	case IMAGE_SAMPLING_MITCHELL:
		for (double n = x - w; n <= x + w; n++)
		{
			double i = round(n);

			if (i < 0 || i >= width)
				continue;

			red += mitchell(i - a + s) * GetPixel(i, y).r;
			green += mitchell(i - a + s) * GetPixel(i, y).g;
			blue += mitchell(i - a + s) * GetPixel(i, y).b;
		}

		p.SetClamp(round(red), round(green), round(blue));
		return p;
		break;
	default:
		printf("sampling model unrecognize");
		break;
	}
}

Pixel Image::ShiftY(double a, double b, double s)
{
	double x = a;
	double y = round(b - s);
	double w = WIDTH;
	Pixel p;
	double red, green, blue;

	red = 0;
	green = 0;
	blue = 0;

	p.Set(0, 0, 0);

	if (y < 0 || y >= width)
		return p;

	switch (sampling_method)
	{
	case IMAGE_SAMPLING_POINT:
		return GetPixel(x, y);
		break;
	case IMAGE_SAMPLING_HAT:
		for (double n = y - w; n <= y + w; n++)
		{
			double i = round(n);

			if (i < 0 || i >= height)
				continue;

			red += hat(i - b + s) * GetPixel(x, i).r;
			green += hat(i - b + s) * GetPixel(x, i).g;
			blue += hat(i - b + s) * GetPixel(x, i).b;
		}
		p.SetClamp(round(red), round(green), round(blue));
		return p;
		break;
	case IMAGE_SAMPLING_MITCHELL:
		for (double n = y - w; n <= y + w; n++)
		{
			double i = round(n);

			if (i < 0 || i >= height)
				continue;

			red += mitchell(i - b + s) * GetPixel(x, i).r;
			green += mitchell(i - b + s) * GetPixel(x, i).g;
			blue += mitchell(i - b + s) * GetPixel(x, i).b;
		}
		p.SetClamp(round(red), round(green), round(blue));
		return p;
		break;
	default:
		printf("sampling model unrecognize");
		break;
	}
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

