/* legendreFunctionsmx.c
 * Description: compute legendrefunctions (4pi-normalized) up to given degree
 *              and order at given specific co-latitude.
 *              Normalization given for example in 
 *              (cf (3.92), Torge, 2012, Geodesy):
 *              \bar{P}_{lm}(t) =
 *                sqrt(k(2l+1)((l-m)!/(l+m)!))P_{lm}(t),
 *              (k = 1 for m == 0, k = 2 for m != 0).
 * Syntax: P = legendreFunctionsmx(theta, max_deg)
 * @input theta: scalar double co-latitude
 * @input max_degree: scalar int maximum degree
 * @output P: matrix with dimensions (max_degree x max_degree).
 *            rows - degree, columns - order.
 *            legendre functions in lower triangular matrix.
 *            upper triangular matrix = 0.
 */
#include "mex.h"  
#include <math.h>
#include <gsl/gsl_sf.h>

void mexFunction(int nlhs, mxArray *plhs[],        // output variables
                 int nrhs, const mxArray *prhs[])  // input variables
{
	#define P_OUT plhs[0]
	#define THETA_IN prhs[0]
	#define MAXDEG_IN prhs[1]

	double *P;
	double theta, cos_theta;
	size_t max_deg;
	int D;
	int m_sum = 0;
	int m, n;

	theta = mxGetScalar(THETA_IN);
	max_deg = (size_t) mxGetScalar(MAXDEG_IN);
	cos_theta = cos(theta);

	// get array of associated legendre polynomials for lower 
	// triangular matrix of P
	double legendre_array[gsl_sf_legendre_array_n(max_deg)];
	if (gsl_sf_legendre_array_e(GSL_SF_LEGENDRE_SCHMIDT,
	                            max_deg, 
	                            cos_theta, 
	                            1.0,
	                            legendre_array))
		mexErrMsgTxt("could not compute legendre polynomials.");

	// create uninitialized P matrix (max_deg+1 x mag_deg+1)
	D = max_deg + 1;
	P_OUT = mxCreateDoubleMatrix(0, 0, mxREAL);
	mxSetM(P_OUT, D);
	mxSetN(P_OUT, D);
	mxSetData(P_OUT, mxMalloc(sizeof(double)*D*D));
	P = mxGetPr(P_OUT);

	// write output matrix P
	for (m = 0; m < D; m++) {
		m_sum += m;
		// lower triangular matrix
		for (n = 0; n < m+1; n++) 
			// get 4-pi normalization from multiplying with
			// sqrt(2*m + 1)
			P[m + D*n] = legendre_array[m_sum + n] * sqrt(2*m + 1);
		// upper triangular matrix
		for (n = m+1; n < D; n++)
			P[m + D*n] = 0;
	}
	return;
}
