#ifndef __CUBE_SOLVER_H__
#define __CUBE_SOLVER_H__

#include "inttypes.h"

#define NUM_BYTES_PER_PIXEL 2
// Body of header goes here.

struct pixel {
	uint8_t red;
	uint8_t green;
	uint8_t blue;
};

void intFPArray(int *arr, int *fpArr, size_t numElems);

void fpShortToDoubleMatrix(short *fpMat, double *doubleMat);

// Function that turns a colored image to black and white
void bwConversion(struct pixel* originalImage, double* bwImage);
void bwConversionFixedPoint(struct pixel *originalImage, short *bwImage);

// Functioh that computes the gradient matrix
void convolutionFn(double* bwImage, int* sobelArray, double* gradientArray);
void convolutionFnFP(short* bwImage, int* sobelArray, short* gradientArray);

void productTruncate(int *product, short *trunc_product);

// Function that performs 3 tasks:
// 1. Compute M matrix
// 2. Compute pixel scores (to be done in computeScoresAndThreshold)
// 3. Compute filtered scores (to be done in computeScoresAndThreshold)
void computePixelScores(double* gradientArrayX, double* gradientArrayY, 
						int* mMatrixBinary);
void computePixelScoresFP(short* gradientArrayX, short* gradientArrayY, 
						char* mMatrixBinary);

double trace(double* mat);
short traceFP(short *mat);

double determinant(double* mat);
short determinantAbsFP(short *mat);

int computeScoresAndThreshold(double* mMatrix);
char computeScoresAndThresholdFP(short* mMatrix);

void extractCorners(int* mMatrixBinary, int* upperCorner,
					int* lowerCorner);
void extractCornersFP(char* mMatrixBinary, int* upperCorner,
					int* lowerCorner);

void exportToMatlab(double* bwImage, double* gradientArrayX, double* gradientArrayY, 
						int* mMatrixBinary);

void exportToMatlabFP(short* bwImage, short* gradientArrayX, short* gradientArrayY, 
						char* mMatrixBinary);

void displayImage(double *matrix);

void displayColorImage(struct pixel *matrix);

#endif /* __CUBE_SOLVER_H__ */
