
#include <math.h>
#include <stdio.h>
#include "inttypes.h"

#include "cube_solver.h"
#include "test_image.h"
#include "ee109-lib/lcd.h"
#include "ee109-lib/char_lcd.h"
#include "ee109-lib/camera.h"
#include "ee109-lib/pushbuttons.h"
#include "ee109-lib/colors.h"
#include "ee109-lib/vga.h"
#include <string.h>

#define LCD_RES_X  400
#define LCD_RES_Y  240
#define NUM_SCORES 4
#define THRESHOLD 0.4
#define MAX_INT 1000000
#define MIN_INT -1000000

int SOBEL_ARRAY_X[9] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
int SOBEL_ARRAY_Y[9] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};

volatile int imageCaptured = 0;
volatile void *imagePtr = (volatile void *) LCD_DEFAULT_FRONT_BUFF_BASE;
struct pixel capturedImage[LCD_RES_X*LCD_RES_Y];
short pixelData[LCD_RES_X];

// Convert short to pixel
struct pixel shortToPixel(short rawData){

	uint8_t blue = rawData & 0x1F;
	uint8_t green = (rawData & 0x7E0) >> 5;
	uint8_t red = (rawData & 0xF800) >> 11;

	struct pixel p; 
	p.blue = blue;
	p.green = green; 
	p.red = red;

	return p;

}

// Function that turns a colored image to black and white
void bwConversion(struct pixel* originalImage, double* bwImage){
	printf("bwConversion\n");
	int j;
	int i;
	for(j = 0; j < LCD_RES_Y; j++) {
		// if (j%20 == 0) printf("bwConversion row %d\0", j);
		for(i = 0; i < LCD_RES_X; i++) {
			// printf("bwConversion col %d", i);
			struct pixel currPixel = originalImage[j*LCD_RES_X + i];
			bwImage[j*LCD_RES_X + i] = sqrt(currPixel.red * currPixel.red + currPixel.blue * currPixel.blue + currPixel.green * currPixel.green) / 255;
			// printf("R: %d G: %d B: %d\n Mag: %g\n", currPixel.red, currPixel.green, currPixel.blue, sqrt(currPixel.red * currPixel.red + currPixel.blue * currPixel.blue + currPixel.green * currPixel.green));
		}
	}
}

// Function that computes the gradient matrix
void convolutionFn(double* bwImage, int* sobelArray, double* gradientArray){
	printf("convolution");
	int j;
	int i;
	for(j = 2; j < LCD_RES_Y; j++){
		if (j%20 == 0) printf("convolution row %d\0", j);
		for(i = 2; i < LCD_RES_X; i++){
			// printf("bwConversion col %d", i);
			gradientArray[j*LCD_RES_X + i] = sobelArray[0]*bwImage[j*LCD_RES_X + i] + sobelArray[1]*bwImage[j*LCD_RES_X + i-1] + sobelArray[2]*bwImage[j*LCD_RES_X + i-2] + 
			sobelArray[3]*bwImage[(j-1)*LCD_RES_X + i] + sobelArray[4]*bwImage[(j-1)*LCD_RES_X + i-1] + sobelArray[5]*bwImage[(j-1)*LCD_RES_X + i-2] + 
			sobelArray[6]*bwImage[(j-2)*LCD_RES_X + i] + sobelArray[7]*bwImage[(j-2)*LCD_RES_X + i-1] + sobelArray[8]*bwImage[(j-2)*LCD_RES_X + i-2];
		}
	}
}

// Trace of a 2x2 matrix. The trace is the sum of the elements of a 
// matrix along the main diagonal
double trace(double* mat){
	return mat[0] + mat[3];
}

double determinant(double* mat){
	return mat[0] * mat[3] - mat[1] * mat[2];
}

int computeScoresAndThreshold(double* mMatrix){
	double tr = trace(mMatrix);
	double score = fabs(determinant(mMatrix)) - tr*tr;
	if (score > THRESHOLD) {
		return 1;
	}
	return 0;
}

// Function that performs 3 tasks:
// 1. Compute M matrix
// 2. Compute pixel scores (to be done in computeScoresAndThreshold)
// 3. Compute filtered scores (to be done in computeScoresAndThreshold)
void computePixelScores(double* gradientArrayX, double* gradientArrayY, 
						int* mMatrixBinary){
		printf("Compute pixel scores\n");

	int j;
	int i;
	for (j = 2; j < LCD_RES_Y - 2; j++){
		if (j%20 == 0) printf("pixel scores row %d\0", j);
		for (i = 2; i < LCD_RES_X - 2; i++){
			double mMatrix[NUM_SCORES];
			int currIndex = j*LCD_RES_X + i;
			double mX = 0;
			double mY = 0;
			double mXY = 0;
			int k;
			int l;
			for (k = -1; k < 2; k++){
				for (l = -1; l < 2; l++){
					mX += gradientArrayX[currIndex + k*LCD_RES_X + l] * gradientArrayX[currIndex + k*LCD_RES_X + l];
					mY += gradientArrayY[currIndex + k*LCD_RES_X + l] * gradientArrayY[currIndex + k*LCD_RES_X + l];
					mXY += gradientArrayX[currIndex + k*LCD_RES_X + l] * gradientArrayY[currIndex + k*LCD_RES_X + l] * gradientArrayX[currIndex + k*LCD_RES_X + l] * gradientArrayY[currIndex + k*LCD_RES_X + l];
				}
			}
			mMatrix[0] = mX;
			mMatrix[1] = mXY;
			mMatrix[2] = mXY;
			mMatrix[3] = mY;
			int binaryScore = computeScoresAndThreshold(mMatrix);
			mMatrixBinary[currIndex] = binaryScore;
		}
	}
}

void extractCorners(int* mMatrixBinary, int* upperCorner,
					int* lowerCorner){
		printf("extract corners\n");

	int minSum = MAX_INT;
	int maxSum = MIN_INT;
	int j;
	int i;
	for (j = 5; j < LCD_RES_Y - 5; j++) {
		 printf("extract corners %d\0", j);
		for (i = 5; i < LCD_RES_X - 5; i++) {
			int currIndex = j * LCD_RES_X + i;
			if (mMatrixBinary[currIndex] && (i+j) < minSum) {
				minSum = i+j;
				upperCorner[0] = j;
				upperCorner[1] = i;
			}
			if (mMatrixBinary[currIndex] && (i+j) > maxSum) {
				maxSum = i+j;
				lowerCorner[0] = j;
				lowerCorner[1] = i;
			}
		}
	}
}

void exportToMatlab(double* bwImage, double* gradientArrayX, double* gradientArrayY, 
						int* mMatrixBinary){
	printf("export to matlab\n");


	// Export matrices for debugging
	// char *origString = "origFile.txt";
	char *bwString = "bwFile.txt";
	char *gradXString = "xGrad.txt";
	char *gradYString = "yGrad.txt";
	char *mMatrixString = "mMatrix.txt";

	// FILE *f0 = fopen(origString, "w");
	FILE *f1 = fopen(bwString, "w");
	FILE *f2 = fopen(gradXString, "w");
	FILE *f3 = fopen(gradYString, "w");
	FILE *f4 = fopen(mMatrixString, "w");

	// fprintf(f0, "%s", "orig_mat = [");
	fprintf(f1, "%s", "bw_mat = [");	
	fprintf(f2, "%s", "x_grad_mat = [");
	fprintf(f3, "%s", "y_grad_mat = [");
	fprintf(f4, "%s", "m_matrix = [");

	int j;
	int i;

	for(j = 0; j < LCD_RES_Y; j++){
		for(i = 0; i < LCD_RES_X; i++){
			if(i == LCD_RES_X - 1) {
				// fprintf(f0, "%g", originalImage[j*LCD_RES_Y + i]);
				fprintf(f1, "%g", bwImage[j*LCD_RES_X + i]);
				fprintf(f2, "%g", gradientArrayX[j*LCD_RES_X + i]);
				fprintf(f3, "%g", gradientArrayY[j*LCD_RES_X + i]);
				fprintf(f4, "%d", mMatrixBinary[j*LCD_RES_X+ i]);
			}
			else{
				// fprintf(f0, "%g,", originalImage[j*LCD_RES_Y + i]);
				fprintf(f1, "%g,", bwImage[j*LCD_RES_X + i]);
				fprintf(f2, "%g,", gradientArrayX[j*LCD_RES_X + i]);
				fprintf(f3, "%g,", gradientArrayY[j*LCD_RES_X + i]);
				fprintf(f4, "%d,", mMatrixBinary[j*LCD_RES_X+ i]);
			}
		}

		// Don't print the semicolon the last time
		if(j == LCD_RES_Y - 1) break;

		// fprintf(f0, "%s", "; ");	
		fprintf(f1, "%s", "; ");	
		fprintf(f2, "%s", "; ");
		fprintf(f3, "%s", ";");
		fprintf(f4, "%s", ";");
	}

	// fprintf(f0, "%s", "];");
	fprintf(f1, "%s", "];");	
	fprintf(f2, "%s", "];");
	fprintf(f3, "%s", "];");
	fprintf(f4, "%s", "];");

	// fclose(f0);
	fclose(f1);
	fclose(f2);
	fclose(f3);
	fclose(f4);
}

void captureImageData() {
	int rowNum = 0;
	int yOffset = 9;

	printf("C");

    for (rowNum; rowNum < LCD_RES_Y; rowNum++){
        memcpy(pixelData, 
        	  (volatile short *)imagePtr + (rowNum << yOffset), 
        	   NUM_BYTES_PER_PIXEL * LCD_RES_X);
        int index;
        for (index = 0; index < LCD_RES_X; index++){
        	short currRawData = pixelData[index];
        	struct pixel p = shortToPixel(currRawData);
        	capturedImage[rowNum * LCD_RES_X + index] = p;
        }
    }	
    char_lcd_write("Snapshot taken\0");
    imageCaptured = 1;
    printf("D");

}


void pushbuttons_isr(void* context, unsigned int id) {
	printf("A");
	uint32_t edge_addr = pushbuttons_get_edge_capture();
	printf("%x", (int) edge_addr);
	if (edge_addr == 2){
		printf("B");
		captureImageData();
	}
	pushbuttons_clear_edge_capture();
}

// Setting up and interfacing with the hardware
void setupHardware() {
	camera_enable_dma(1);
	lcd_enable_dma(1);
	//lcd_draw_rectangle(0, 0, LCD_RES_X, LCD_RES_Y, GREEN);
	pushbuttons_enable_interrupts(pushbuttons_isr);

	vga_enable_dma(1);
	vga_draw_rectangle(0, 0, LCD_RES_X, LCD_RES_Y, BLACK);
}

void displayImage(double *matrix){
	int i;
	int j;

	volatile short *vga_front_buffer = (volatile short *)VGA_DEFAULT_FRONT_BUFF_BASE;

	short max_val = 0xF;
	short pixel;
	for(j = 0; j < VGA_RES_Y; j++){
		for(i = 0; i < VGA_RES_X; i++){
			pixel = (short)(matrix[j*LCD_RES_X + i] * max_val);
			vga_front_buffer[(j << 9) + i] = pixel;
		}
	}
}

void displayColorImage(struct pixel *matrix){
	int i;
	int j;

	volatile short *vga_front_buffer = (volatile short *)VGA_DEFAULT_FRONT_BUFF_BASE;

	struct pixel p;
	short pixel;
	for(j = 0; j < VGA_RES_Y; j++){
		for(i = 0; i < VGA_RES_X; i++){
			p = matrix[j*LCD_RES_X + i];
			pixel = (((uint16_t) p.red) << 11) | (((uint16_t) p.green) << 5) | (uint16_t) p.blue;
			vga_front_buffer[(j << 9) + i] = pixel;
		}
	}
}

int main(int argc, char *argv[]){
	
	// int structPixelSize = sizeof(structPixelSize);

	setupHardware();

	printf("Waiting to acquire image\0");
	while(!imageCaptured){}
	printf("Image acquired\0");

	struct pixel* originalImage = (struct pixel*) capturedImage;

	// struct pixel originalImage[LCD_RES_X * LCD_RES_Y];

	double bwImage[LCD_RES_X * LCD_RES_Y];
	displayColorImage(originalImage);
	bwConversion(originalImage, bwImage);

	//Debug
	displayImage(bwImage);

	int mMatrixBinary[LCD_RES_X * LCD_RES_Y];
	double gradientArrayX[LCD_RES_X * LCD_RES_Y];
	double gradientArrayY[LCD_RES_X * LCD_RES_Y];

	convolutionFn(bwImage, SOBEL_ARRAY_X, gradientArrayX);
	convolutionFn(bwImage, SOBEL_ARRAY_Y, gradientArrayY);

	computePixelScores(gradientArrayX, gradientArrayY, mMatrixBinary);
	int upperCorner[2];
	int lowerCorner[2];
	extractCorners(mMatrixBinary, upperCorner, lowerCorner);
	exportToMatlab(bwImage, gradientArrayX, gradientArrayY, 
						mMatrixBinary);
	

	printf("Upper 0: %d, Upper 1: %d\n", upperCorner[0], upperCorner[1]);
	printf("Lower 0: %d, Lower 1: %d\n", lowerCorner[0], lowerCorner[1]);

	return 0;

}
