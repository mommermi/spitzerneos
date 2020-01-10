/* Michael Mommert, DLR, 26. Apr. 2010
 * michael.mommert@dlr.de
*/

//-------------------------------------- random number generator from Numerical Recipes

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

static long idum = (-1);

float ran2(long *idum)

{
        int j;
        long k;
        static long idum2=123456789;
        static long iy=0;
        static long iv[NTAB];
        float temp;

        if (*idum <= 0) {
                if (-(*idum) < 1) *idum=1;
                else *idum = -(*idum);
                idum2=(*idum);
                for (j=NTAB+7;j>=0;j--) {
                        k=(*idum)/IQ1;
                        *idum=IA1*(*idum-k*IQ1)-k*IR1;
                        if (*idum < 0) *idum += IM1;
                        if (j < NTAB) iv[j] = *idum;
                }
                iy=iv[0];
        }
        k=(*idum)/IQ1;
        *idum=IA1*(*idum-k*IQ1)-k*IR1;
        if (*idum < 0) *idum += IM1;
        k=idum2/IQ2;
        idum2=IA2*(idum2-k*IQ2)-k*IR2;
        if (idum2 < 0) idum2 += IM2;
        j=iy/NDIV;
        iy=iv[j]-idum2;
        iv[j] = *idum;
        if (iy < 1) iy += IMM1;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;
}

// --------------------------------------------------------------------------
// generate uniform random number between a and b
double randnum (double a, double b)
{
	return (double)ran2(&idum)*(b-a)+a;
}


// generate gaussian random number around mean with standard deviation sigma
double gaussrand (double mean, double sigma)
{
	//algorithm derived from Numerical Recipes for C, chapter 7-2, routine 'gasdev'
	double uni1, uni2, gauss1, gauss2, radius = 0.;
	
	while ( (radius>1.) || (radius==0.) )
	{
		//create two uniform random numbers in square [-1:1][-1:1]
		uni1 = randnum(-1.,1.);
		uni2 = randnum(-1.,1.);
		
		//calculate radius from origin (0,0)
		radius = uni1*uni1 + uni2*uni2;
	}
	
	gauss1 = uni1 * sqrt(-2.0*log(radius)/radius) * sigma + mean;
	gauss2 = uni2 * sqrt(-2.0*log(radius)/radius) * sigma + mean;
	
	return gauss1;
}


// restricts gaussrand to a given range
double gaussrand_restricted (double mean, double sigma, double a, double b)
{
	double rand = a-1.;
	
	while ( (rand < a) || (rand > b) )
		rand = gaussrand (mean, sigma);
	
	return rand;
}
