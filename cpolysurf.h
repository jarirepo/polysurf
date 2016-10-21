#ifndef cpolysurf_h
#define cpolysurf_h

#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <string.h>
#include <mex.h>
#include <matrix.h>

#define SYNTAX "USAGE: [Gx,Gy,Gz] = cpolysurf(P [,surfu=16, surfv=16, tighten=0])"
#define DEBUG 0

//define inputs and outputs
#define P_IN 		prhs[0]
#define SURFU_IN	prhs[1]
#define SURFV_IN	prhs[2]
#define TIGHTEN_IN	prhs[3]

#define GX_OUT		plhs[0]
#define GY_OUT		plhs[1]
#define GZ_OUT		plhs[2]

#define MAX_SEGM	99			// max. value for SURFU and SURFV
#define MAX_PTS		999			// max. number of polyline pts.
#define ZERO 		1E-12
#define PARAM_SCALE 65535

typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned long uint32;

typedef struct {
	double **pp;		//d-by-n array with polyline pts.
	int d;				//dimensionality of pts.
	int n;				//number of pts.
	//double *t;			//parameter 0<=t[...]<=1
	int *idx;			//mapped parameters t[...] to line segment indices			
	//double *h;			//local line parameter	
	uint16 *t;
	uint16 *h;
} Polyline;

//globals
static char str[256] = {0};		//used for general output messages
static char tmp[32] = {0};
//double *pu, *pv;
//uint16 *pu, *pv;

//function prototypes
void parameterize(Polyline* S);
void injectParams(const Polyline* S, double*, int &);
void mapParams(Polyline* S, const double*, int);

//functions supporting 16-bit fixed point parameterization
void parameterize16(Polyline* S);
void injectParams(const Polyline* S, uint16* pu, int &nu);
void mapParams(Polyline* S, const uint16* pu, int nu);
void radixSort(uint16* x, int n, int order=1);
void generateSurface(double* px, double* py, double* pz, 
	const Polyline* S, const uint16* u, int nu, const uint16* v, int nv);
#endif
