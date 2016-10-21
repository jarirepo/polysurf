/* cpolysurf.cpp

	MEX implementation of the polysurf function to generate a bilinear 
	surface from 4 connected 3-D polylines.
	
	USAGE: [Gx,Gy,Gz] = cpolysurf(P,surfu,surfv [,tighten=0])
	
	INPUTS
	P	:	a 4-element cell array containing the definition points of the 
			four polylines. Each polyline is a sequence of points that 
			represent a boundary edge/curve G(u,v) of the final bilinear 
			surface G(u,v) = [x(u,v),y(u,v),z(u,v)]
			
			The organization of P is as follows.
			P{1}: curve G(u,0)
			P{2}: curve G(u,1)
			P{3}: curve G(0,v)
			P{4}: curve G(1,v)
				
	surfu	:	number of surface segments in the U-direction (16)
	surfv	:	number of surface segments in the V-direction (16)
	
	tighten	: 	set to 1 to include the polyline control points in the
				generated surface to generate tight surface boundaries. 
				The number of segments in the u- and v-directions may 
				then become larger than specified by surfu and surfv. 
				On the other hand, the resolution will also be slightly 
				improved near the definition points. Default is tight=0

	OUTPUTS
	Gx,Gy,Gz	:	matrices containing the X-,Y- and Z-values of the generated 
					surface. These can be used directly as input to MATLAB-
					functions surf() and contour() to visualize the surface 
					properties.
					
	Date	:	2016-10-19
	Author	:	Jari Repo, University West, jari.repo@hv.se
	
	Change log
	
	2016-10-20, JRE
	-The parameterization of the surface boundaries is now based on 16-bit 
	fixed point numbers for faster execution which invloves parameter 
	injection and sorting of the parameter arrays u and v.
	Some decimal precision is lost due to this but hardly noticeable.
	
	Refs.
	* https://en.wikipedia.org/wiki/Coons_patch
	* https://se.mathworks.com/help/matlab/cc-mx-matrix-library.html
	* https://se.mathworks.com/help/matlab/mex-library.html	
*/

#include "cpolysurf.h"

//Gateway function
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
	if (nlhs<2 || nlhs>3 || nrhs<1 || nrhs>4) {
		mexErrMsgTxt(SYNTAX);
	}
	
	/* check on input P
	P must be a 4-element cell array containing the control points of the four 
	polylines that define the surface boundaries in the u- and v-irections:	
		P{1}: G(u,0)
		P{2}: G(u,1)
		P{3}: G(0,v)
		P{4}: G(1,v)
	*/
	if (!mxIsCell(P_IN)) {
		mexErrMsgTxt("Input P must be a cell array");
	}
	
	/* check individual elements in P */
	if (mxGetM(P_IN)!=1 || mxGetN(P_IN)!=4) {
		mexErrMsgTxt("Input P must be a 1-by-4 cell array");
	}
	
	/** Check on the individual elements in P
		p[4]  	: pointers to the individual cells in P
		Np[4] 	: number of polyline control pts.
		d 		: polyline dimensionality (must be either 2 or 3)
	**/
	mxArray *p[4] = {0};
	int i, j, d, 
		Np[4] = {0},
		surfu =  16, 
		surfv = 16, 
		tighten = 0;	
	
	for (i=0; i<4; i++) {
		p[i] = mxGetCell(P_IN, (mwIndex) i);	//get ptr. to array i
		//check data type
		if (mxIsNumeric(p[i])) {
			//polyline dimensionality check (must be consistent; either 2 or 3)
			if (i==0) {
				d = (int) mxGetM(p[i]);
				if (d<2 || d>3) {
					mexErrMsgTxt("Polyline dimensions must be either 2 or 3");
				}
			}else {
				if (mxGetM(p[i])!=d) {
					mexErrMsgTxt("Polyline dimensions must be identical");
				}
			}
			//check on number of polyline pts. Must be at least 2 for a single linear segment
			Np[i] = (int) mxGetN(p[i]);
			if (Np[i]<2 || Np[i]>MAX_PTS) {
				sprintf(str,"Number of polyline points must be between 2 and %d\n",MAX_PTS);
				mexErrMsgTxt(str);
			}				
		}			
	}
	
	//check to see that the dimensionality d matches the number of outputs
	if (d!=nlhs) {
		mexErrMsgTxt("Polyline dimensionality must match the number of outputs");
	}
			
	if (nrhs>1) {
		//check on input SURFU
		if (!mxIsNumeric(SURFU_IN)) {
			mexErrMsgTxt("Input SURFU must be a real positive integer");
		}
		surfu = (int) mxGetScalar(SURFU_IN);
		if (surfu<1 || surfu>MAX_SEGM) {
			sprintf(str,"Input SURFU must be within 1 and %d\n",MAX_SEGM);
			mexErrMsgTxt(str);
		}
	}
	
	if (nrhs>2) {
		// check on input SURFV
		if (!mxIsNumeric(SURFV_IN)) {
			mexErrMsgTxt("Input SURFV must be a real positive integer");
		}
		surfv = (int) mxGetScalar(SURFV_IN);
		if (surfv<1 || surfv>MAX_SEGM) {
			sprintf(str,"Input SURFV must be within 1 and %d\n",MAX_SEGM);
			mexErrMsgTxt(str);
		}
	}

	if (nrhs>3) {
		//check on input TIGHTEN
		if (!mxIsNumeric(TIGHTEN_IN)) {
			mexErrMsgTxt("Input TIGHTEN_IN must a real integer");
		}
		tighten = (int) mxGetScalar(TIGHTEN_IN);		
		if (!(tighten==0 || tighten==1)) {
			mexErrMsgTxt("Input TIGHTEN must be either 0 or 1");
		}
	}
	
	#if DEBUG
	mexPrintf("%s\n\tsurfu: %d\n\tsurfv: %d\n\ttighten: %d\n",
		SYNTAX,surfu,surfv,tighten);	
	#endif
		
	/* 	Allocate buffers for the surface parameter grids u and v
		Make these larger than surfu+1 and surfv+1 to allow injection
		of line parameters just in case tighten=1, to avoid dynamic
		memory re-allocation. */
	int nu = surfu+1,
		nv = surfv+1;
		
	//double* pu = (double *) mxCalloc((mwSize) (nu+Np[0]+Np[1]), sizeof(double));	
	//double* pv = (double *) mxCalloc((mwSize) (nv+Np[2]+Np[3]), sizeof(double));	
	uint16* pu = (uint16 *) mxCalloc((mwSize) (nu+Np[0]+Np[1]), sizeof(uint16));	
	uint16* pv = (uint16 *) mxCalloc((mwSize) (nv+Np[2]+Np[3]), sizeof(uint16));
	
	// genearte u- and v-grids
	for (i=0; i<nu; i++) {
		//pu[i] = (double) i/surfu;
		//sprintf(str,"%-16.6e\n",pu[i]); mexPrintf(str);
		pu[i] = (uint16) (65535*i/surfu+0.5f);
	}
	for (j=0; j<nv; j++) {
		//pv[j] = (double) j/surfv;
		pv[j] = (uint16) (65535*j/surfv+0.5f);
	}
	
	// allocate 4 Polyline structs
	Polyline S[4];
	
	// parameterize the polylines (choord length parametrization)
	for (i=0; i<4; i++) {
		//map the polyline point buffer p[i] to an matrix with d rows and Np[i] columns
		S[i].n = Np[i];
		S[i].d = d;
		//S[i].t = new double[Np[i]];
		S[i].t = new uint16[Np[i]];
		
		/*
		S[i].pp = new double*[d];
		S[i].pp[0] = (double *) mxGetPr(p[i]);
		for (j=1; j<d; j++) {
			S[i].pp[j] = S[i].pp[j-1] + Np[i];
		}			
		*/
		
		// MATLAB stores M-by-N matrices in transposed order in memory
		S[i].pp = new double*[ S[i].n ];
		S[i].pp[0] = (double *) mxGetPr( p[i] );
		for (j=1; j<S[i].n; j++) {
			S[i].pp[j] = S[i].pp[j-1] + d;
		}
		
		// do the chordal parameterization
		//parameterize( &S[i] );
		parameterize16( &S[i] );
		
		#if DEBUG		
		mexPrintf("Polyline parameterization\nt: ");
		for (j=0; j<S[i].n; j++) {
			//sprintf(tmp,"%-16.6e\n",S[i].t[j]);
			sprintf(tmp,"%d",S[i].t[j]);
			if (j<S[i].n-1) strcat(tmp,",");
			mexPrintf(tmp);
		}
		mexPrintf("\n");
		#endif
	}
	
	#if DEBUG
	//output polyline 0
	for (i=0; i<S[0].n; i++) {
		*str = 0;
		for (j=0; j<S[0].d; j++) {
			sprintf(tmp,"%.2f",S[0].pp[i][j]);
			if (j<S[0].d-1)
				strcat(tmp,", ");
			strcat(str,tmp);
		}
		strcat(str,"\n");
		mexPrintf(str);
	}
	#endif
		
	if (tighten) {
		/** Inject polyline parameters into the u- and v-grid to ensure a tight surface boundary.
		This will update pu and nu. Note that the parameter injection is done prior to the sorting 
		to reduce the computational overhead. **/
		injectParams(&S[0],pu,nu);
		injectParams(&S[1],pu,nu);
		injectParams(&S[2],pv,nv);
		injectParams(&S[3],pv,nv);
		
		#if DEBUG		
		mexPrintf("u = ");
		for (i=0; i<nu; i++) {
			//sprintf(tmp,"%.4f",pu[i]);
			sprintf(tmp,"%d",pu[i]);
			if (i<nu-1) strcat(tmp,", ");
			mexPrintf(tmp);
		}
		mexPrintf("\nv = ");
		for (j=0; j<nv; j++) {
			//sprintf(tmp,"%.4f",pv[j]);
			sprintf(tmp,"%d",pv[j]);
			if (j<nv-1) strcat(tmp,", ");
			mexPrintf(tmp);
		}
		//sprintf(str,"\n\n\tnu = %d, nv = %d\n", nu, nv);
		//mexPrintf(str);
		mexPrintf("\n");
		#endif
		
		// sort the u- and v-arrays	in ascending order using a fast base-2 radix sort algorithm (16-bit)
		radixSort(pu,nu);
		radixSort(pv,nv);

		#if DEBUG	
		mexPrintf("Sorted parameter vectors:\n");
		mexPrintf("u = ");
		for (i=0; i<nu; i++) {
			//sprintf(tmp,"%.4f",pu[i]);
			sprintf(tmp,"%d",pu[i]);
			if (i<nu-1) strcat(tmp,", ");
			mexPrintf(tmp);
		}
		mexPrintf("\nv = ");
		for (j=0; j<nv; j++) {
			//sprintf(tmp,"%.4f",pv[j]);
			sprintf(tmp,"%d",pv[j]);
			if (j<nv-1) strcat(tmp,", ");
			mexPrintf(tmp);
		}
		sprintf(str,"\n\n\tnu = %d, nv = %d\n", nu, nv);
		mexPrintf(str);
		#endif
	}
	
	/** Map the parameter values (u,v) to boundary curve segment numbers to speed up 
	later curve interpolation. This will update S[...].idx[...] which containts 
	the line segment numbers 0,1,2,..., and h[...] which contains the local line 
	parameters to improve the performance in the final bilinear interpolation. **/
	S[0].idx = new int[nu];
	S[1].idx = new int[nu];
	S[2].idx = new int[nv];
	S[3].idx = new int[nv];
	
	/** 
	S[0].h = new double[nu];
	S[1].h = new double[nu];
	S[2].h = new double[nv];
	S[3].h = new double[nv];
	**/
	
	S[0].h = new uint16[nu];
	S[1].h = new uint16[nu];
	S[2].h = new uint16[nv];
	S[3].h = new uint16[nv];
	
	mapParams(&S[0],pu,nu);
	mapParams(&S[1],pu,nu);
	mapParams(&S[2],pv,nv);
	mapParams(&S[3],pv,nv);
	
	/** Generate the bilinear surface
	The surface points G(u,v)=[x(u,v),y(u,v),z(u,v)] are written directly to the 
	output buffers **/	
	GX_OUT = mxCreateNumericMatrix((mwSize)nv,(mwSize)nu,mxDOUBLE_CLASS,mxREAL);
	GY_OUT = mxCreateNumericMatrix((mwSize)nv,(mwSize)nu,mxDOUBLE_CLASS,mxREAL);
	if (d==3 && nlhs>2) {
		GZ_OUT = mxCreateNumericMatrix((mwSize)nv,(mwSize)nu,mxDOUBLE_CLASS,mxREAL);
	}
	
	generateSurface((double *)mxGetPr(GX_OUT), (double *)mxGetPr(GY_OUT), (double *)mxGetPr(GZ_OUT),
		S, pu, nu, pv, nv);
		
	/** free memory and terminate **/
	mxFree(pu);
	mxFree(pv);	
	for (i=0; i<4; i++) {
		delete [] S[i].pp;
		delete [] S[i].t;
		delete [] S[i].idx;
		delete [] S[i].h;
	}
}

/** parameterize(...)

**/
void parameterize(Polyline *S) {	
	int i,j,n=S->n,d=S->d;
	double a, L = 0;
	double **pp = S->pp;	
	
	S->t[0] = 0;	
	for (i=1; i<n; i++) {
		a = 0;
		for (j=0; j<d; j++) {
			a += pow(pp[i][j]-pp[i-1][j], (double)2);
		}
		L += sqrt(a);	// cumulative chord. length
		S->t[i] = L;
	}	
	for(i=1;i<n-1; i++) {
		S->t[i] /= L;
	}	
	S->t[n-1] = 1;
}

/** injectParams(...)

**/
void injectParams(const Polyline* S, double *u, int &nu) {
	//injection of line parameters S->t into the parameter grid u,
	//avoiding dublicate parameter values,
	//updates u and nu
	static const double tol = 1e-11;
	int i,j=1,k,n=S->n,nu0=nu;	
	bool found;
	double ti;
	
	for (i=1; i<n-1; i++) {
		ti = S->t[i];
		//set starting index for the search
		k = j;
		//we don't have to scan u for longer than its original number
		//of elements nu0-1
		found = false;
		//while (k<(nu0-1) && ti<u[k] && (!found)) {
		while (k<(nu0-1) && (!found)) {
			found = (abs(u[k++]-ti)<tol)?true:false;
		}
		if (!found) {
			//insert t[i] into u and update nu
			u[nu++] = ti;
		}else {
			//ti was found in u[k-1],
			//since t is monotonically increasing, the next search 
			//can be started at index k
			j = k;
		}
	}	
}

/** mapParams(...)

**/
void mapParams(Polyline* S, const double *u, int nu) {
	int i,k,n=S->n;
	S->idx[0] = 0;
	S->h[0] = 0;
	i = 1;	
	for (k=1; k<nu; k++) {
		while (u[k] >= S->t[i] && i<n-1) {
			i++;
		}
		i--;
		S->idx[k] = i;
		//local line parameter
		S->h[k] = (u[k]-S->t[i])/(S->t[i+1]-S->t[i]);
	}
}

/** uint16 version **/
void parameterize16(Polyline *S) {	
	int i,j,n=S->n,d=S->d;
	double a, L;
	double **pp = S->pp;	
	double *s = new double[n];
	s[0] = 0;
	for (i=1; i<n; i++) {
		a = 0;
		for (j=0; j<d; j++) {
			a += pow(pp[i][j]-pp[i-1][j], (double)2);
		}
		s[i] = s[i-1]+sqrt(a);	//cum. chord length
	}
	L = s[n-1];
	S->t[0] = 0;	
	for(i=1;i<n-1; i++) {
		S->t[i] = (uint16)(65535*s[i]/L+0.5f);	
	}	
	S->t[n-1] = 65535;
	delete [] s;
}

void injectParams(const Polyline* S, uint16 *u, int &nu) {
	//injection of line parameters S->t into the parameter grid u,
	//avoiding dublicate parameter values,
	//updates u and nu
	int i,j=1,k,n=S->n,nu0=nu;	
	bool found;
	uint16 ti;
	
	for (i=1; i<n-1; i++) {
		ti = S->t[i];
		//set starting index for the search
		k = j;
		//we don't have to scan u for longer than its original number
		//of elements nu0-1
		found = false;
		//while (k<(nu0-1) && ti<u[k] && (!found)) {
		while (k<(nu0-1) && (!found)) {
			found = (u[k++]==ti)?true:false;	//integer comparison
		}
		if (!found) {
			//insert t[i] into u and update nu
			u[nu++] = ti;
		}else {
			//ti was found in u[k-1],
			//since t is monotonically increasing, the next search 
			//can be started at index k
			j = k;
		}
	}	
}

void mapParams(Polyline* S, const uint16 *u, int nu) {
	int i,k,n=S->n;
	S->idx[0] = 0;
	S->h[0] = 0;
	i = 1;	
	for (k=1; k<nu; k++) {
		while (u[k] >= S->t[i] && i<n-1) {
			i++;
		}
		i--;
		S->idx[k] = i;
		//local line parameter
		S->h[k] = (uint16)(65535*(u[k]-S->t[i])/(S->t[i+1]-S->t[i])+0.5f);
	}
	#if DEBUG
	mexPrintf("Local line parameters:\n");
	for (k=0; k<nu; k++) {
		sprintf(tmp,"%d",S->h[k]);
		if (k<nu-1) strcat(tmp,", ");
		mexPrintf(tmp);
	}
	mexPrintf("\n");
	#endif
}

/** radixSort(...)

	Base 2 radix sort algorithm for 16-bit numbers
**/	
void radixSort(uint16 *x, int n, int order) {
	uint16 *x0 = new uint16[n],
		   *x1 = new uint16[n];		   
    int n0,n1,i,r;	
	uint32 m = 1;	
    for (r=0; r<16; r++, m<<=1) {
        n0 = n1 = 0;
        for (i=0; i<n; i++) {
            if (x[i] & m)
                x1[n1++] = x[i];
            else
                x0[n0++] = x[i];
        }
		if (order==1) {	//ascending order
			for (i=0; i<n0; i++)
				x[i] = x0[i];
			for (i=0; i<n1; i++)
				x[n0+i] = x1[i];
		}else {	//descending order		
			for (i=0; i<n1; i++)
				x[i] = x1[i];
			for (i=0; i<n0; i++)
				x[n1+i] = x0[i];
		}
    }	
    delete [] x0;
	delete [] x1;
}

/** generateSurface(...)

	Generate surface based on Coon's bilinear interpolation
	Remarks: Gx,Gy,Gz are all nv-by-nu arrays	
**/	
void generateSurface(double *px, double *py, double *pz, 
	const Polyline* S, const uint16* u, int nu, const uint16* v, int nv) {
					
	//get dimensionality of the polyline control points
	int d = S[0].d;
	
	//get surface "corners" i.e. the polyline endpoints	
	double 	*G00 = S[0].pp[0],
			*G10 = S[3].pp[0],
			*G01 = S[1].pp[0],
			*G11 = S[1].pp[ S[1].n-1 ];
	
	#if DEBUG
	mexPrintf("Surface corners:\n");
	if (d==3) {
		sprintf(str,"G00=(%.4f,%.4f,%.4f)\n",G00[0],G00[1],G00[2]);	mexPrintf(str);	
		sprintf(str,"G10=(%.4f,%.4f,%.4f)\n",G10[0],G10[1],G10[2]);	mexPrintf(str);	
		sprintf(str,"G01=(%.4f,%.4f,%.4f)\n",G01[0],G01[1],G01[2]);	mexPrintf(str);	
		sprintf(str,"G11=(%.4f,%.4f,%.4f)\n",G11[0],G11[1],G11[2]);	mexPrintf(str);
	}else {
		sprintf(str,"G00=(%.4f,%.4f)\n",G00[0],G00[1]);	mexPrintf(str);	
		sprintf(str,"G10=(%.4f,%.4f)\n",G10[0],G10[1]);	mexPrintf(str);	
		sprintf(str,"G01=(%.4f,%.4f)\n",G01[0],G01[1]);	mexPrintf(str);	
		sprintf(str,"G11=(%.4f,%.4f)\n",G11[0],G11[1]);	mexPrintf(str);	
	}
	#endif
		
	double P0v[3],P1v[3],Pu0[3],Pu1[3],Puv[3];
	int i,j,k,a,b;
	long offs = 0L;
	
	for (j=0; j<nv; j++) {	//v-direction
		a = S[2].idx[j];
		b = S[3].idx[j];

		for (k=0; k<d; k++) {	//x,y,z
			P0v[k] = S[2].pp[a][k] + (S[2].pp[a+1][k] - S[2].pp[a][k]) * (double)S[2].h[j]/65535;							
			P1v[k] = S[3].pp[b][k] + (S[3].pp[b+1][k] - S[3].pp[b][k]) * (double)S[3].h[j]/65535;							
		}
		
		for (i=0; i<nu; i++) {	//u-direction
			a = S[0].idx[i];
			b = S[1].idx[i];
			
			for (k=0; k<d; k++) {	//x,y,z
				
				Pu0[k] = S[0].pp[a][k] + (S[0].pp[a+1][k] - S[0].pp[a][k]) * (double)S[0].h[i]/65535;							
				Pu1[k] = S[1].pp[b][k] + (S[1].pp[b+1][k] - S[1].pp[b][k]) * (double)S[1].h[i]/65535;
				
				//apply Coon's bilinear interpolation formula
				//We need to divide by 65535 since we use 16-bit integers for all the parameter values
				
				/**
				Puv[k] = (Pu0[k]*(65535-v[j]) + Pu1[k]*v[j] + P0v[k]*(65535-u[i]) + P1v[k]*u[i])/65535 - 
					(G00[k]*(65535-u[i])*(65535-v[j]) +
					G01[k]*v[j]*(65535-u[i]) +
					G10[k]*u[i]*(65535-v[j]) + 
					G11[k]*u[i]*v[j])/4294836225L;					
				**/
				
				Puv[k] = (Pu0[k]*(65535-v[j]) + Pu1[k]*v[j] + P0v[k]*(65535-u[i]) + P1v[k]*u[i] - 
						(G00[k]*(65535-u[i])*(65535-v[j]) + G01[k]*v[j]*(65535-u[i]) +
						G10[k]*u[i]*(65535-v[j]) + G11[k]*u[i]*v[j])/65535)/65535;				
			}
			
			/** store surface point.
			Again, MATLAB wants them transposed **/			
			*(px + i*nv + j) = Puv[0];
			*(py + i*nv + j) = Puv[1];
			if (d==3) {
				*(pz + i*nv + j) = Puv[2];
			}
		}
	}
}
