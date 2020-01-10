 /*
 *      NEOSURVEY_MODEL - modelling software umbrella routine for neosurvey
 *
 *      Feb 24, 2015 michael.mommert@nau.edu
 *      
 *      command line arguments: 
 *      model input model [-e value] [-s]
 *        input           data input file
 *        model           models to be used for calculation
 *                        1: ref. STM, 2: FRM, 3: NEATM (floating eta), 
 *                        4: NEATM (fixed eta) [needs d/w/f]
 *                        in case of '4': d: Delbo 2004, w: Wolters 2007, f: fix eta,
 *                        (if none of the above: 4f is used)
 *        [-tpm fix/rand] perform TPM calculations with fix/random parameters
 *        [-eta value]    fix eta value/initial guess
 *        [-b value]      fix beta value
 *        [-fc]           flux correction (for reflected sunlight)
 *        [-fr min max]   thermal flux ratio filter
 *        [-rr value]     solar reflectance ratio for all wavelengths
 *        [-ncc]           no color correction (input fluxes are thermal fluxes)
 *        [-er min max]   eta range
 *        [-montecarlo]   perform error estimation
 *        [-chi2montecarlo] perform error estimation with flux error adjustment
 *        [-s]            calculate each measurement separately
 *        [-tno]          uses Brucker-relation and TNO D-pV-H    
 *        [-Herror]       read out additional column in input file
 *        [-p]            plot data
 *      
 *      Flux conversion:
 *        F[mJy] = 3.335640952e14 * wavelength^2 * F[W m^-2 um^-1]
 *      
 *      
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include "gnuplot_i.h"		//to use gnuplot interface
#include "statusbar.h"

#define INPUT_LINE_LENGTH 2048    //max length for inputfile lines
#define INPUT_COLUMNS 35          //number of columns in inputfile
#define EXP_NUMBER_OF_OBJECTS 10000 //expected number of objects
#define JMAX 25                   //total number of iterations in integration

//constants
// use PI from math.h library
#define PI   M_PI
#define PI2  M_PI_2
#define EMISSIVITY 0.9                 //emissivity
#define AU          149597870.691        //Astronomical Units in kilometers
#define HC2         5.95521967e7         //h*c^2 in W m^2 (=/1e-24)
#define HCK         14387.6866           //h*c/k using um in c
#define SOLARCONST  1367                 //solar constant in W m^-2 for 1 AU distance
#define STEFANBOLTZMANN   5.67e-8        //Stefan Boltzmann constant
#define CLIGHT      299792458            //speed of light


//string containing model to be applied
static char *model;
/* flag for model currently used
 * 0: no model
 * 1: STM
 * 2: FRM
 * 3: NEATM (floating eta)
 * 4: NEATM (fixed eta)
*/

//static double eta_max = 2.7; // for TNOs
static double eta_max = 3.14;
static double eta_min = 0.6;

static int TNO = 0; //=1: use Brucker et al. 2009 p-q relation and 1342.6 instead of 1329 in D-pV-H relation

static int colcorr = 1; //apply color correction

static int colcorr_error = 0; //turns to one if color correction not possible for at least one object

//define ASTEROID structure containing parameters
typedef struct
{	char   *name;
	double wavelength, eta, chi2, modelparameter;
	double geodist, heliodist, phaseangle, reflectanceratio;
	double absmag, slopepar, absmag_error, slopepar_error, reflectanceratio_error;
	double flux, fluxerr, colorcorrection, flux_raw, diameter, pv;
	double diameter_sigma, pv_sigma, eta_sigma, modelparameter_sigma;
	double diameter_mean, pv_mean, eta_mean, modelparameter_mean;
	double diameter_median, pv_median, eta_median, modelparameter_median;
	double diameter_lower, pv_lower, eta_lower, modelparameter_lower;
	double diameter_upper, pv_upper, eta_upper, modelparameter_upper;
        double diameter_lower3sig, pv_lower3sig, eta_lower3sig;
        double diameter_upper3sig, pv_upper3sig, eta_upper3sig;

	double airmass; //only for flux estimations
	
	//'next_measurement' points to index of the next ASTEROID object with the same name
	//if there are no more measurements -> 0
	int next_measurement;
	
	//'fitparflag' is set to 1 in case of an unrealistic high/low eta/beta value, otherwise it is 0
	int fitparflag;
	//'badflux' is set to 1 in case of flux ratio breaking check condition
	int badflux;
} ASTEROID;


//initialize some functions
void fluxcorrection ( ASTEROID *obj, int start_asteroid, int number_of_objects);
void output_on_screen ( ASTEROID *obj, int number_of_objects, int fitting );
void runmodel(ASTEROID *obj, int start_asteroid, int number_of_objects, int fitting, double epsilon, double eta, char* resultsfile, int display, char *notes);
int long_output_in_file ( ASTEROID obj, char *name, char *notes );
int short_output_in_file ( ASTEROID obj, char *name, char *notes );
void plotmodel ( ASTEROID *obj, int i, double epsilon, double (*planck_model)(ASTEROID, double) );



double pv_from_D (double D, double H)
{
	if ( TNO == 0 )
		return (1329*1329)*pow(10,(-.4)*H)/(D*D);
	else
//		return (1342.6*1342.6)*pow(10,(-.4)*H)/(D*D);
		return (1330.3*1330.3)*pow(10,(-.4)*H)/(D*D);  //using m_solar=-26.76 (provided by Esa)
}



double D_from_pv (double pv, double H)
{
	if ( TNO == 0 )
		return 1329/sqrt(pv)*pow(10,-H/5.);
	else
//		return 1342.6/sqrt(pv)*pow(10,-H/5.);
		return 1330.3/sqrt(pv)*pow(10,-H/5.);  //using m_solar=-26.76 (provided by Esa)
}

double bondalbedo (double G, double pv)
{
	double A;

	if ( TNO == 0 )
		A = (0.29 + 0.684*G)*pv;
	else
		A = (0.479 + 0.336*pv)*pv; //Brucker et al. 2009
	if (A >= 1.0)
	{
		printf("WARNING: A > 1 (%f)\n", A);
		return 0.99999;
	}
	else
		return A;
}

// Romberg Integration for model integration in one dimension
double romberg ( double (*integration_model)(double, double, ASTEROID, double),
			double a, double b, double eps, ASTEROID obj, double T0 )
{
	double I[JMAX+1][JMAX+1];
	double h[JMAX+1];
	double sum;
	int n, i, j;
		
	//calculate I_{1,1} = trapezoidal rule
	h[1] = (b-a);
	I[1][1] = 0.5*h[1]*(integration_model(a, eps, obj, T0) 
	                    + integration_model(b, eps, obj, T0));
	//printf("      mit I[1][1] = %.3e\n", I[1][1]);
	for ( n = 2; n <= JMAX; n++ )
	{	
		h[n] = h[n-1]*0.5;
		
		//calculate I_{n,1} using less calculation steps
		for ( i = 1, sum = 0; i <= pow(2,(n-2)); i++ )
			sum += integration_model(a + (2*i-1) * h[n], eps, obj, T0);
		I[n][1] = 0.5*I[n-1][1] + h[n]*sum;
		
      		//calculate I_{n,1} using existing results
		int m = n, k = 1;
		for ( j = 1; j < n; j++ )
		{
			k++;
			m--;
			
			//calculate I_{m,k}
			I[m][k] = I[m+1][k-1] + (I[m+1][k-1] - I[m][k-1])
					/(pow(2,(2*k-2)) - 1);
		}
		
		//check for nan
		if ( isnan(I[1][n]) )
		{
			printf("ERROR in Romberg integration: NaN (rerun with slight offset)\n");
			return I[1][n];
		}
				
		//do at least 5 iterations
		if (n > 5)
			//check for convergence inbetween epsilon and special case
			//of I[1][j] == 0 for (0<j<6)
			if ( (fabs(I[1][n]-I[1][n-1])/I[1][n] < eps) ||
			      (I[1][n] == 0 && I[1][n-1] == 0 ) )
				return I[1][n];
	}
	
	printf("WARNING: Romberg Integration did not converge! Return %.3e\n", I[1][n-1]);
	return I[1][n-1];
}



/*read in masterlist input file and put data into ASTEROID objects
  returns number of objects read from file = number of lines*/
int readinput ( char *filename, ASTEROID *obj, double eta, double diameter, double pv, double modelparameter, 
                int corrflux, int fitting, double reflectanceratio, int inputfilestyle )
{
	FILE *inputdata = fopen (filename, "rt");
	
	//auxiliary variables
	char *line = malloc(sizeof(char)*INPUT_LINE_LENGTH); 
	char rubbish_string[1000];
	int j, k;

	if ( inputdata == NULL )
	{
		printf("ERROR: can't open file %s\n", filename);
		return 1;
	}
	
	//-------------------------------------------------------------------------- NEOSurvey input file
	//input format: name geodist heliodist phaseangle absmag absmag_error slopepar wavelength flux fluxerr
	else
	{
		j = 0;  //line number and asteroid index
		while ( !feof(inputdata) && (j < EXP_NUMBER_OF_OBJECTS) ) 
		{
			sprintf(line, " "); //clear line to prevent multiple reading of the last line
			fgets(line, INPUT_LINE_LENGTH, inputdata);
			
			//if the line contains a '#' or it is too short, ignore it
			if ( (strpbrk(line, "#") == NULL) && (strlen(line) > 5) &&
			     (strstr(line, "ignoreme") == NULL) )
			{
				if ( j > EXP_NUMBER_OF_OBJECTS)
				{
					printf("ERROR: expected number of objects exceeded!\n");
					abort();
				}
				
				//write data into ASTEROID struct
				if ( inputfilestyle == 0)
					sscanf(line, "%s %lg %lg %lg %lg %lg %lg %lg %lg", 
					       obj[j].name, 
					       &obj[j].geodist,
					       &obj[j].heliodist,
					       &obj[j].phaseangle,
					       &obj[j].absmag,
					       &obj[j].slopepar,
					       &obj[j].wavelength,
					       &obj[j].flux,
					       &obj[j].fluxerr);
				else
					sscanf(line, "%s %lg %lg %lg %lg %lg %lg %lg %lg %lg", 
					       obj[j].name, 
					       &obj[j].geodist,
					       &obj[j].heliodist,
					       &obj[j].phaseangle,
					       &obj[j].absmag,
					       &obj[j].absmag_error,
					       &obj[j].slopepar,
					       &obj[j].wavelength,
					       &obj[j].flux,
					       &obj[j].fluxerr);

				// add calibration uncertainty
				obj[j].fluxerr = sqrt(obj[j].fluxerr*obj[j].fluxerr + 0.05*obj[j].flux*0.05*obj[j].flux);
				
				obj[j].phaseangle *= PI2/90.;
				
				obj[j].colorcorrection = 1.;
				obj[j].reflectanceratio = reflectanceratio;
				obj[j].modelparameter = modelparameter;
				obj[j].diameter = diameter;
				obj[j].pv       = pv;
				
				//convert mJy in W/m²/um
				obj[j].flux  = obj[j].flux_raw = obj[j].flux/(3.335640952e14*obj[j].wavelength*obj[j].wavelength);
				obj[j].fluxerr /= (3.335640952e14*obj[j].wavelength*obj[j].wavelength);
				
				//errors
				obj[j].slopepar_error   = 0.0;	
				obj[j].reflectanceratio_error = 0.3;
				
				//assign eta values
				if ( strpbrk(model, "1") != NULL )  //STM
					obj[j].eta = eta;
				if ( strpbrk(model, "2") != NULL )  //FRM
					obj[j].eta = PI;
				if ( strpbrk(model, "3") != NULL )  //NEATM (floating eta)
					obj[j].eta = eta;  //as an initial value
				if ( strpbrk(model, "4") != NULL )  //NEATM
				{
					if ( strpbrk(model, "d") != NULL ) //use Delbo 2004 method
						obj[j].eta = (0.013*obj[j].phaseangle/PI2*90.)+0.89+eta;
					else if ( strpbrk(model, "w") != NULL ) //use Wolters 2007 method
						obj[j].eta = (0.013*obj[j].phaseangle/PI2*90.)+0.91+eta;
					else if ( strpbrk(model, "m") != NULL ) //use Mainzer 2011 method
						obj[j].eta = (0.00963*obj[j].phaseangle/PI2*90.)+0.761+eta;
					else if ( strpbrk(model, "n") != NULL ) //use neosurvey method
						obj[j].eta = (0.01*obj[j].phaseangle/PI2*90.)+0.87+eta;
					else
						obj[j].eta = eta;  //fixed eta
				}
				
				obj[j].next_measurement = 0; //in case of first appearence
				
				obj[j].fitparflag = 0;
				
				j++;
			}
		}
		obj[j-1].next_measurement = 0; //make sure, last object points to no other object
	}



	
	//free(line); free(rubbish_string);
	
	//flux correction (correction for reflected light)
	if ( corrflux == 1 )
	{
		fluxcorrection ( obj, 0, j );
	}
	
	fclose(inputdata);
	
	return j;
}


//include model specific functions
#include "stm.h"
#include "frm.h"
#include "neatm.h"
#include "random.h"
#include "errors.h"
#include "colorcorr.h"

//---------------------
//calculation functions


//perform flux corrections for reflected light (Bowell et al., Asteroids II)
void fluxcorrection ( ASTEROID *obj, int start_asteroid, int number_of_objects )
{
	int i;
	
	for ( i=start_asteroid; i<number_of_objects; i++)
	{
		//flux emitted by sun; adapted from Gueymard 2004 for Spitzer wavelengths
		//and Rieke 2008 (http://www.iop.org/EJ/article/1538-3881/135/6/2245/aj_135_6_2245.tables.html)
		//fluxes in W m^-2 um^-1
		double sunflux;	
		

		//printf("--------------%s\n", obj[i].name);
		
		//calculate apparent magnitude of the asteroid as seen by the observer
		double v_alpha = obj[i].absmag - 2.5 
					* log10( (1-obj[i].slopepar)* exp( -3.33*pow(tan(obj[i].phaseangle/2),0.63))
				                  + obj[i].slopepar * exp( -1.87*pow(tan(obj[i].phaseangle/2),1.22)) );
		double m_v = v_alpha + 5*log10(obj[i].heliodist*obj[i].geodist);
		
		//printf("H=%f V=%f\n", obj[i].absmag, m_v);

		// my approach
		
		//use fix solar fluxes for each wavelength
		//wavelengths for Spitzer observations
		if ( obj[i].wavelength == 3.545 || obj[i].wavelength == 3.55 || obj[i].wavelength == 3.6 )
		{
			sunflux = 1.384e1; //3.551 um from Rieke 2008
			//sunflux = 1.356e1; //for 3.55 um from Gueymard 2004
			//obj[i].colorcorrection = 1.16388; //default NEO color correction factors (Trilling 2008?)
		}
		else if ( obj[i].wavelength == 4.493 || obj[i].wavelength == 4.5 || obj[i].wavelength == 4.49)
		{
			sunflux = 5.4586e0; //linear approximation for 4.493 um from Rieke 2008
			//sunflux = 5.059e0;  //linear approximation for 4.493 um from Gueymard 2004
			//obj[i].colorcorrection = 1.08555;//default NEO color correction factors (Trilling 2008?)
		}
		else if ( obj[i].wavelength == 5.731 || obj[i].wavelength == 5.8 )
		{
			sunflux = 2.167e0;  //5.731 um from Rieke 2008
			//sunflux = 1.96e0;  //linear approximation for 5.731 um from Gueymard 2004
			//obj[i].colorcorrection = 1.04;//default NEO color correction factors (Trilling 2008?)
		}
		else if ( obj[i].wavelength == 7.872 || obj[i].wavelength == 8.0 )
		{
			sunflux = 6.3608e-1;  //linear approximation for 7.872 um from Rieke 2008
			//sunflux = 6.03e-1;  //linear approximation for 7.872 um from Gueymard 2004
			//obj[i].colorcorrection = 1.01;//default NEO color correction factors (Trilling 2008?)
		}
		else if ( obj[i].wavelength == 3.3526 )
		{
			sunflux = 1.718e1; //3.3510 1.718E-02 um from Rieke 2008
		}
		else if ( obj[i].wavelength == 4.6028)
		{
		        sunflux = 4.882e0; //4.6010 4.882E-03 from Rieke 2008
		}
		else if ( obj[i].wavelength == 11.5608 || obj[i].wavelength == 11.0984 )
		{
			sunflux = 1.400e-1;  //11.5610 1.400E-04 from Rieke 2008
		}
		else if ( obj[i].wavelength == 16 )
		{
			sunflux = 3.858e-2;  //16.0010 3.858E-05 from Rieke 2008
		}
		else if ( obj[i].wavelength == 22.0883 || obj[i].wavelength == 22.6405 )
		{
			sunflux = 1.075e-2;  //22.0860 1.075E-05 from Rieke 2008
		}

		else
		{
			obj[i].colorcorrection = 1.0;
		        printf("WARNING: no flux correction possible! Wavelength (%f um) unknown\n", obj[i].wavelength);
			//abort();
		}
		
		//calculate solar flux at specific wavelength using planck function and
		//normalization at 3.6 um (5.5464*10^16 mJy = 13.2020 W m^-2 um^-1) in latter dimension
		sunflux = 7483.465532/(obj[i].wavelength*obj[i].wavelength*obj[i].wavelength*obj[i].wavelength*obj[i].wavelength
					*(exp(2.4701724/obj[i].wavelength)-1));

		
		//fraction of sunlight flux reflected by asteroid is ...
		double refl_fraction = pow(10, -0.4*(m_v-(-26.73)));
		

		//spectral reflectivity at lambda >= 3.6um is reflectanceratio*reflectivity(V)
		// if reflected solar light is not all the measured flux
		if (sunflux*refl_fraction*obj[i].reflectanceratio < 0.99*obj[i].flux)
		  refl_fraction *= obj[i].reflectanceratio;
		else
		  refl_fraction = 0.99*obj[i].flux/sunflux;
		
		if (obj[i].reflectanceratio < 0)
		  refl_fraction = 0;

		// this line replaces the reflectance ratio with the fraction of the measured flux that is reflected solar light
		obj[i].reflectanceratio = sunflux * refl_fraction/obj[i].flux;

		//thermal flux
		obj[i].flux = (obj[i].flux_raw - (sunflux * refl_fraction));

      	}
}


//determine chi^2 sum of one object (using fitting or not)
double chi2sum_object ( ASTEROID *obj, int i, double diameter, double epsilon, int fitting,
                        double (*planck_model)(ASTEROID, double) )
{
	double chi2_sum = 0;
	
	// color correction
	double pv = pv_from_D(diameter, obj[i].absmag);
	double T0 = pow((1-bondalbedo(obj[i].slopepar, pv))*SOLARCONST/(obj[i].heliodist*obj[i].heliodist*obj[i].eta*EMISSIVITY*STEFANBOLTZMANN), 0.25);

	// printf("T_mean: %f\n", T0/pow(2,0.25));

	if ( fitting == 1)   //use chi^2 fitting
	{
		int all_done = 0;
		while ( all_done == 0) 
		{
			obj[i].diameter = diameter;

			// -----------------------------------------------------------
			
			double calcflux = planck_model(obj[i], epsilon);
			
			obj[i].chi2 = (obj[i].flux-calcflux)*(obj[i].flux-calcflux)
					/(obj[i].fluxerr*obj[i].fluxerr);

			
			chi2_sum += obj[i].chi2;
			
				
			//breakout if all measurements concerning this
			//object are processed
			if ( obj[i].next_measurement != 0)
				i = obj[i].next_measurement;
			else 
				all_done = 1;
		}
	} 
	if ( fitting == 0)        //don't use chi^2 fitting
	{
		obj[i].diameter = diameter;
		
		double calcflux = planck_model(obj[i], epsilon);
		
		//chi2 determination
		/* obj[i].chi2 = (obj[i].flux/obj[i].colorcorrection-calcflux) */
		/* 		*(obj[i].flux/obj[i].colorcorrection-calcflux) */
		/* 		/(obj[i].fluxerr/obj[i].colorcorrection */
		/* 		  *obj[i].fluxerr/obj[i].colorcorrection); */

		obj[i].chi2 = (obj[i].flux-calcflux)*(obj[i].flux-calcflux)
			  /(obj[i].fluxerr*obj[i].fluxerr);

		
		/*printf("diameter: %.3f, flux_meas: %.6e, flux_calc: %.6e, flux_err: %.6e, chi2: %.6e\n", obj[i].diameter,
				obj[i].flux/obj[i].colorcorrection, calcflux, obj[i].fluxerr, obj[i].chi2);*/
		
		chi2_sum = obj[i].chi2;
	}

	return chi2_sum;
}

/* find best chi^2 fit for model and determine diameter and pv
 * results are written into ASTEROID struct
 * returns chi^2-value of best fit
*/
double find_diameter ( ASTEROID *obj, int j, double epsilon, int fitting,
                        double (*planck_model)(ASTEROID, double))
{
	//calculate diameter limit estimations
	//... for pv = 0.01 as upper limit for diameter 
	//0.3 is an empirical value
	double diameter_max = 0.3*D_from_pv(0.01, obj[j].absmag);
	//... for A = 0.9 as lower limit for diameter
	double diameter_min = D_from_pv(0.9/(0.29+0.684*obj[j].slopepar), obj[j].absmag);
	
	double chi2_sum,  chi2_lastsum = 1e6,
	       diameter, diameter_stepsize = (diameter_max - diameter_min)/10;
	
	/*--------------
	 * coarse diameter estimation over the physically possible range
	 * by calculating chi^2 for 10 equidistant diameters
	*/
	double chi2_min = -1, diameter_bestguess;
	for ( diameter = diameter_min; 
			diameter <= diameter_max;
			diameter += diameter_stepsize )
	{
		chi2_sum = chi2sum_object ( obj, j, diameter, epsilon, fitting, planck_model );
		
		if ( (chi2_sum < chi2_min) || (chi2_min == -1) )
		{
			chi2_min = chi2_sum;
			diameter_bestguess = diameter;
		}
		//printf("coarse  d=%.5f, chi2=%.5e, chi2_last=%.5e eta=%.3f (b/z)eta=%.3f\n", obj[j].diameter, chi2_sum, chi2_lastsum, obj[j].eta, obj[j].modelparameter);
	}

	/*--------------
	 * enhance first diameter estimation by minimization of chi^2
	 * by increasing the diameter value using a predefined stepsize
	 * starting at the best guess diameter; if chi^2 decreases keep on
	 * moving, if chi^2 increases negate and half stepsize
	 * repeat until chi2 converges 
	*/	
	//assign initial guesses to global variables and calculate chi2_sums
	diameter_stepsize = 0.1*diameter_bestguess;
	diameter = diameter_bestguess - diameter_stepsize;
	chi2_sum = chi2sum_object ( obj, j, diameter, epsilon, fitting, planck_model );
	//while ( fabs(chi2_sum - chi2_lastsum) > epsilon*epsilon )  very old criterion
	while ( fabs(chi2_sum - chi2_lastsum)/chi2_sum > epsilon ) 
	{
		//do step
		diameter += diameter_stepsize;
		
		chi2_lastsum = chi2_sum;
		chi2_sum = chi2sum_object ( obj, j, diameter, epsilon, fitting, planck_model );
		
		//printf("fine  d=%.5f, chi2=%.5e, chi2_last=%.5e eta=%.3f (b/z)eta=%.3f\n", obj[j].diameter, chi2_sum, chi2_lastsum, obj[j].eta, obj[j].modelparameter);
		
		//if chi2 increases with diameter, turn around with smaller stepsize
		if ( chi2_sum > chi2_lastsum )
			diameter_stepsize /= (-2);
	}
	
	//determine number of measurements of this object
	int number_this_object = 0;
	if ( fitting == 1 )
	{
		int all_done = 0;
		int i = j;              //start with j, since j has to be the first object of this kind
		while ( all_done == 0) 
		{
			number_this_object++;
			//breakout if all measurements concerning this
			//object are processed
			if ( obj[i].next_measurement != 0)
				i = obj[i].next_measurement;
			else 
				all_done = 1;
		}
	} 
	else number_this_object = 1;

	//disable fitting if there is only measurement of this object
	if ( number_this_object == 1 )
		fitting = 0;

	//assign calculated diameter and pv to each object included in the calculations
	//determine chi2_red and assign it, too
	if ( fitting == 1 )
	{
		// SWITCH HERE FROM CHI2 TO REDUCED CHI2
		if ( strpbrk(model, "3") != NULL )  //NEATM (floating eta methods)
			if ( number_this_object > 2 )
				chi2_sum = chi2_sum/(number_this_object-2);	
			else chi2_sum = -chi2_sum;
		else //NEATM (fixed eta methods)
			if ( number_this_object > 1 )
				chi2_sum = chi2_sum/(number_this_object-1);	
			else chi2_sum = -chi2_sum;

		int all_done = 0;
		int i = j;              //start with j, since j has to be the first object of this kind
		while ( all_done == 0) 
		{
			obj[i].diameter = diameter;
			obj[i].pv = pv_from_D(diameter, obj[j].absmag);

			obj[i].chi2 = chi2_sum;

			//breakout if all measurements concerning this
			//object are processed
			if ( obj[i].next_measurement != 0)
				i = obj[i].next_measurement;
			else 
				all_done = 1;
		}
	}
	else 
	{
		obj[j].diameter = diameter;
		obj[j].pv = pv_from_D(diameter, obj[j].absmag);

		// SWITCH HERE FROM CHI2 TO REDUCED CHI2
		// since there is only measurement for this object, Chi2 has no meaning
		obj[j].chi2 = 0;
	}
	
	// return sqrt(chi2_sum);   why sqrt????
	return chi2_sum;

}


/* find best chi^2 fit for a set of measurements and determine eta
 * results are written into ASTEROID struct
 * returns chi^2-value of best fit
 * fitting == 1 is mandatory
*/
double fit_eta ( ASTEROID *obj, int j, double epsilon, double eta_guess, 
                        double (*planck_model)(ASTEROID, double) )
{
	double chi2_sum = 1e6, chi2_lastsum = 0, lasteta = 1e6;
	double eta;
	//initialize stepsize
	int k = 0;  //number of iteration steps performed
	double eta_stepsize = 0.1;

	/*--------------
	 * coarse eta estimation over the physically possible range
	 * by calculating chi^2 for 10 equidistant diameters
	*/
	double chi2_min = -1;
	for ( eta = 0.7; eta <= 2.5; eta += eta_stepsize )
	{
		//update eta values of ASTEROIDS
		int all_done = 0;
		int i = j;              //start with j, since j has to be the first object of this kind
		while ( all_done == 0) 
		{
			obj[i].eta = eta;
			
			//breakout if all measurements concerning this
			//object are processed
			if ( obj[i].next_measurement != 0)
				i = obj[i].next_measurement;
			else 
				all_done = 1;
		}

		chi2_sum = fabs(find_diameter(obj, j, epsilon, 1, planck_model));
		
		if ( (chi2_sum < chi2_min) || (chi2_min == -1) )
		{
			chi2_min = chi2_sum;
			eta_guess = eta;
		}
		//printf("coarse  d=%.5f, chi2=%.5e, chi2_last=%.5e eta=%.3f (b/z)eta=%.3f\n", obj[j].diameter, chi2_sum, chi2_lastsum, obj[j].eta, obj[j].modelparameter);
	}

	eta = eta_guess;

	/*improve eta estimate
        ---------------------------------------------*/

	//while ( ((fabs(chi2_sum - chi2_lastsum)/chi2_sum > epsilon) || (k <= (int)((eta_max-eta_guess)/eta_stepsize))) && (k < 100) ) //do at least as many steps to reach the upper eta limit
	while ( (fabs(eta-lasteta)/eta > epsilon) || (k <= (int)((eta_max-eta_guess)/eta_stepsize)) )
	{
		chi2_lastsum = chi2_sum;
		lasteta = eta;

		//do step
		eta += eta_stepsize;
		
		//update eta values of ASTEROIDS
		int all_done = 0;
		int i = j;              //start with j, since j has to be the first object of this kind
		while ( all_done == 0) 
		{
			obj[i].eta = eta;
			
			//breakout if all measurements concerning this
			//object are processed
			if ( obj[i].next_measurement != 0)
				i = obj[i].next_measurement;
			else 
				all_done = 1;
		}

		chi2_sum = fabs(find_diameter(obj, j, epsilon, 1, planck_model));

		//update chi2 values of ASTEROIDS
		all_done = 0;
		i = j;              //start with j, since j has to be the first object of this kind
		while ( all_done == 0) 
		{
			obj[i].chi2 = chi2_sum;
			
			//breakout if all measurements concerning this
			//object are processed
			if ( obj[i].next_measurement != 0)
				i = obj[i].next_measurement;
			else 
				all_done = 1;
		}

		//check for realistic eta values
		if ( (eta > eta_max) || (eta < eta_min) )
		{
			printf("\b  ERROR: unrealistic eta value (%.2f)! \n", eta);
			
			all_done = 0;
			i = j;              //start with j, since j has to be the first object of this kind
			while ( all_done == 0) 
			{
				chi2_sum           = 1e9;
				obj[i].chi2       = chi2_sum;
				obj[i].fitparflag = 1;
				
				//breakout if all measurements concerning this
				//object are processed
				if ( obj[i].next_measurement != 0)
					i = obj[i].next_measurement;
				else 
					all_done = 1;
			}
			
			break;
		}
		
		//printf("   eta:  d=%.5f, chi2=%.5e, eta=%.3f modelpar=%.3f\n", obj[j].diameter, chi2_sum, eta, obj[j].modelparameter /* /PI*180. */);
		
		//if chi2 increases with eta, turn around with smaller stepsize
		if ( chi2_sum > chi2_lastsum )
			eta_stepsize /= (-2);
				
		k++;
		
		printf("\b%c", rotate_statusbar());
		fflush(stdout);
	}
	printf("\b");
	
	return chi2_sum;
}

/* find best chi^2 fit for a set of measurements and determine modelparameter
 * results are written into ASTEROID struct
 * returns chi^2-value of best fit
 * fitting == 1 is mandatory
*/
double fit_modelparameter ( ASTEROID *obj, int j, double epsilon, double (*planck_model)(ASTEROID, double) )
{
	double modelparameter, lastmodelparameter = 1e6;
	double chi2_sum[21];
	
	double max_modelparameter;
	if ( strpbrk(model, "56") != NULL ) max_modelparameter = PI2;
	if ( strpbrk(model, "9") != NULL ) max_modelparameter = 1.;
	
	//--------------
	// coarse modelparameter estimation
	// by calculating chi^2 for 20 equidistant modelparameters

	int m, n, number_of_minima = 0;
	double modelparameter_bestguess[5], eta_bestguess[5];
	for (m=0; m<5; m++) { modelparameter_bestguess[m]=0; eta_bestguess[m]=0; }
	for ( modelparameter = 0, m=0; modelparameter < max_modelparameter+max_modelparameter/20/2; modelparameter += max_modelparameter/20, m++ )
	{
		//update modelparameter values of ASTEROIDS
		int all_done = 0;
		int i = j;              //start with j, since j has to be the first object of this kind
		while ( all_done == 0) 
		{
			obj[i].modelparameter = modelparameter;
			obj[i].chi2 = 0;
			
			obj[i].eta = 1;   //only in case no eta fitting is applied, otherwise eta will be determined by fitting
			
			//breakout if all measurements concerning this
			//object are processed
			if ( obj[i].next_measurement != 0)
				i = obj[i].next_measurement;
			else 
				all_done = 1;
		}
		
		//chi2_sum[m] = find_diameter ( obj, j, epsilon, 1, planck_model );  //do not use eta fitting
		chi2_sum[m] = fit_eta ( obj, j, epsilon, 1.5, planck_model );   //do use eta fitting
		if ( (chi2_sum[m]+1) == chi2_sum[m] )
			chi2_sum[m] = 1e9;
		
		printf("coarse sslat: (%d) d=%.5f, pv=%.5f, chi2=%.5e, eta=%.3f sslat=%.3f\n", m, obj[j].diameter, obj[j].pv, chi2_sum[m], obj[j].eta, obj[j].modelparameter );
		
		printf("\b%c", rotate_statusbar());
		fflush(stdout);
	}

	printf("\n\n\n");

	//find minima in chi2 distribution
	for (m=0, n=0; m<=20; m++)
	{
		printf ("%f\n", chi2_sum[m]);
		if ( ( (m == 0) && (chi2_sum[1] > chi2_sum[0]) ) ||          //minimum at lowest modelparameter value 
		     ( (m == 19) && (chi2_sum[18] > chi2_sum[19]) ) ||       //minimum at highest modelparameter value
		     ( (chi2_sum[m-1] > chi2_sum[m]) && (chi2_sum[m+1] > chi2_sum[m]) ) ) //minimum inbetween modelparameter range
		{
			printf("min @ %f\n", max_modelparameter/20*m);
			modelparameter_bestguess[n] = max_modelparameter/20*m;
			n++; number_of_minima++;
		}
	}

	printf ("%d minima found: modelparameter = ", number_of_minima);
	for (m=0; m<number_of_minima; m++)
		printf ("%f ", modelparameter_bestguess[m]);
	printf("\n");

	printf("\b# ");

	for (m=0; m<number_of_minima; m++)
	{
		printf("examine minimum no. %d\n", m+1);

		//initialize stepsize
		int k = 0;  //number of iteration steps performed
		double modelparameter_stepsize = max_modelparameter/(2*10);
		modelparameter = modelparameter_bestguess[m];
		double chi2_sum = 1e6, chi2_lastsum = 0;
		//----------------
		// fine tuning of modelparameter
		// 
		//while ( (fabs(chi2_sum - chi2_lastsum)/chi2_sum > epsilon) || (k<=5) ) //do at least 5 steps
		while ( ((fabs(modelparameter-lastmodelparameter)/modelparameter > epsilon) || (k<=5)) && (k<100) ) //do at least 5 steps
		{
			chi2_lastsum = chi2_sum;
			lastmodelparameter = modelparameter;			

			//check for modelparameter being outside realistic range	
			//if limit is reached, turn around with smaller stepsize
			if ( modelparameter < 1e-6 )
			{
				modelparameter = 0;
				modelparameter_stepsize = fabs(modelparameter_stepsize)/2;
			}
			if ( modelparameter > (max_modelparameter-1e-6) )
			{
				modelparameter = max_modelparameter;
				modelparameter_stepsize = -abs(modelparameter_stepsize)/2;
			}
				
			//update modelparameter values of ASTEROIDS
			int all_done = 0;
			int i = j;              //start with j, since j has to be the first object of this kind
			while ( all_done == 0) 
			{
				obj[i].modelparameter = modelparameter;
				obj[i].chi2 = chi2_sum;
				
				//breakout if all measurements concerning this
				//object are processed
				if ( obj[i].next_measurement != 0)
					i = obj[i].next_measurement;
				else 
					all_done = 1;
			}
			
			//chi2_sum = find_diameter ( obj, j, epsilon, 1, planck_model );  //do not use eta fitting
			chi2_sum = fit_eta ( obj, j, epsilon, 1, planck_model );   //do use eta fitting		

			printf("fine sslat: d=%.5f, pv=%.5f, chi2=%.5e, chi2_last=%.5e eta=%.3e sslat=%.3e\n", obj[j].diameter, obj[j].pv, chi2_sum, chi2_lastsum, obj[j].eta, obj[j].modelparameter );

			//if chi2 increases with modelparameter, turn around with smaller stepsize
			if ( chi2_sum > chi2_lastsum )
				modelparameter_stepsize /= (-2);
			
			//do step
			k++;
			modelparameter += modelparameter_stepsize;
			
			
			printf("\b%c", rotate_statusbar());
			fflush(stdout);
		}
		printf("\b\b");

		printf("solution for minimum %d: diameter %f, albedo %f, chi2 %e\n", m, obj[j].diameter, obj[j].pv, obj[j].chi2);
	}


	return 0;
}

//checks for occurring flags in data and distributes them over all data of this object
void checkforflags (ASTEROID *obj, int j)
{
	int badflux = 0, fitparflag = 0;
	
	//search for flags in object data
	int all_done = 0;
	int i = j;              //start with j, since j has to be the first object of this kind
	while ( all_done == 0) 
	{
		if ( obj[i].fitparflag == 1 )
			fitparflag = 1;
		if ( obj[i].badflux == 1 )
			badflux = 1;
		
		//breakout if all measurements concerning this
		//object are processed
		if ( obj[i].next_measurement != 0)
			i = obj[i].next_measurement;
		else 
			all_done = 1;
	}
	
	//mark all data of this object
	all_done = 0;
	i = j;              //start with j, since j has to be the first object of this kind
	while ( all_done == 0) 
	{
		if ( fitparflag == 1 )
			obj[i].fitparflag = 1;
		if ( badflux == 1 )
			obj[i].badflux = 1;
		
		//breakout if all measurements concerning this
		//object are processed
		if ( obj[i].next_measurement != 0)
			i = obj[i].next_measurement;
		else 
			all_done = 1;
	}
}

//------------------------------------------------------
//output functions

void output_on_screen ( ASTEROID *obj, int number_of_objects, int fitting )
{
	int i;
	
	for ( i = 0; i < number_of_objects ; i++ )
		printf( "%25.25s %.5f %.4f %.3f %.1f | %.6f %.6f %.6f | %.6f %.6f\n",
			obj[i].name, obj[i].diameter, obj[i].pv, obj[i].eta, obj[i].wavelength, 
			obj[i].flux_raw*3.335640952e14*obj[i].wavelength*obj[i].wavelength,
			obj[i].flux*3.335640952e14*obj[i].wavelength*obj[i].wavelength,
			obj[i].flux/obj[i].colorcorrection*3.335640952e14*obj[i].wavelength*obj[i].wavelength, 
			obj[i].fluxerr*3.335640952e14*obj[i].wavelength*obj[i].wavelength,
			obj[i].fluxerr/obj[i].colorcorrection*3.335640952e14*obj[i].wavelength*obj[i].wavelength);
}

int long_output_in_file ( ASTEROID obj, char *filename, char *notes )
{
	//check whether file already exists
	int exist = 0;
	FILE *test = fopen (filename, "r");
	if (test) { exist = 1; fclose(test); }

	FILE *outputresults = fopen (filename, "at");
	
	if ( outputresults == NULL )
	{	printf("ERROR: can't open file results\n");
		return 1;
	}
	
	if (exist == 0)
	  fprintf( outputresults, "#                                                       name    diameter     d_sig      d_low      d_up    d_low3sig   d_up3sig       pv       pv_sig     pv_low     pv_up   pv_low3sig  pv_up3sig     eta      eta_sig    eta_low     eta_up  eta_low3sig eta_up3sig   sslat    sslat_sig  sslat_low   sslat_up   geodist   heliodist     alpha     absmag   absmag_err  slopepar slopepar_err reflratio refl_error  lambda    flux_raw    fluxerr     flux      colcorr chi2_red flags_notes\n");

	char flags[100];
	strcpy (flags, "");
	
	if ( obj.flux < 0 )
		sprintf( flags, "#Flux<0_");
	if ( obj.fitparflag == 1 )
		sprintf( flags, "#FitUnrealistic_");
	if ( obj.badflux == 1 )
		sprintf( flags, "#BadFluxRatio_");
	
	
	fprintf( outputresults, "%60.60s  %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f"
	                        " %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f"
	                        " %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5e %s%s\n",
	obj.name, obj.diameter, obj.diameter_sigma, obj.diameter_lower, obj.diameter_upper, obj.diameter_lower3sig, obj.diameter_upper3sig, 
	obj.pv, obj.pv_sigma, obj.pv_lower, obj.pv_upper, obj.pv_lower3sig, obj.pv_upper3sig, 
	obj.eta, obj.eta_sigma, obj.eta_lower, obj.eta_upper, obj.eta_lower3sig, obj.eta_upper3sig,
	obj.modelparameter/PI*180, obj.modelparameter_sigma/PI*180, obj.modelparameter_lower/PI*180, obj.modelparameter_upper/PI*180,
	obj.geodist, obj.heliodist, obj.phaseangle/PI2*90, obj.absmag, obj.absmag_error, obj.slopepar, obj.slopepar_error,
	obj.reflectanceratio, obj.reflectanceratio_error, obj.wavelength,
	obj.flux_raw*3.335640952e14*obj.wavelength*obj.wavelength,
	obj.fluxerr*3.335640952e14*obj.wavelength*obj.wavelength, 
	obj.flux/obj.colorcorrection*3.335640952e14*obj.wavelength*obj.wavelength, 
	obj.colorcorrection, obj.chi2, flags, notes );
	
	fclose(outputresults);
	
	return 0;
}

int short_output_in_file ( ASTEROID obj, char *filename, char *notes )
{
	//check whether file already exists
	int exist = 0;
	FILE *test = fopen (filename, "r");
	if (test) { exist = 1; fclose(test); }

	FILE *outputresults = fopen (filename, "at");
	
	if ( outputresults == NULL )
	{	printf("ERROR: can't open file results\n");
		return 1;
	}

	if (exist == 0)
		fprintf( outputresults, "#                                                       name   diameter    pv     eta    sslat   "
                         "d_sig   pv_sig eta_sig sslat_sig alpha   chi2_red  flags_notes\n");
	
	char flags[100];
	strcpy (flags, "");
	
	if ( obj.flux < 0 )
		sprintf( flags, "#Flux<0_");
	if ( obj.fitparflag == 1 )
		sprintf( flags, "#FitUnrealistic_");
	if ( obj.badflux == 1 )
		sprintf( flags, "#BadFluxRatio_");
	
	fprintf( outputresults, "%60.60s  %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5e %s%s\n",
	obj.name, obj.diameter, obj.pv, obj.eta, obj.modelparameter/PI*180, 
	obj.diameter_sigma, obj.pv_sigma, obj.eta_sigma, obj.modelparameter_sigma,
	obj.phaseangle/PI2*90, obj.chi2, flags, notes );
	
	fclose(outputresults);
	
	return 0;
}




/* create files containing flux measurements and model spectral energy distributions of *ONE* asteroid
*/
void plotmodel ( ASTEROID *obj, int i, double epsilon,
		 double (*planck_model)(ASTEROID, double) ) // 'i' denotes the first object to be plotted from array 'obj' 
{
	//lower and upper limit for wavelength
	double lambda_min = 1e6, lambda_max = 0;
	
    	//create normalization model asteroid using the parameters of the 'i'-th object in array 'obj'
        ASTEROID normal = obj[i];

	//-----------------------------------
	//create file containing measurements
	char measurementsfilename[250];
	sprintf(measurementsfilename, "%s.meas", obj[i].name);
	FILE *measurements = fopen(measurementsfilename, "wt");
	int all_done = 0;
	int j = i;
	while ( all_done == 0 ) 
	{
		// determine flux correction factor:
	        // fluxcorrection = flux of object 'j' / flux of normal object ('i')
	        normal.wavelength = obj[j].wavelength;
		double fluxcorrection = planck_model(obj[j], epsilon)/planck_model(normal, epsilon);

		/* measurement file output:
		 *   wavelength [microns]
		 *   flux       [W m^-2 microns^-1]
		 *   fluxerror  [W m^-2 microns^-1]
		*/			
		fprintf(measurements, "%.2e   %.5e   %.5e\n",
		        obj[j].wavelength, obj[j].flux/obj[j].colorcorrection*3.335640952e14*obj[j].wavelength*obj[j].wavelength/fluxcorrection, 
			obj[j].fluxerr/obj[j].colorcorrection*3.335640952e14*obj[j].wavelength*obj[j].wavelength/fluxcorrection);

  		if ( obj[j].wavelength < lambda_min )
			lambda_min = obj[j].wavelength;

		if ( obj[j].wavelength > lambda_max )
			lambda_max = obj[j].wavelength;
		
		//replace actual heliocentric and geocentric distance by normalised ones
		obj[j].heliodist = normal.heliodist;
		obj[j].geodist   = normal.geodist;

		//breakout if all measurements concerning this
		//object are processed
		if ( obj[j].next_measurement != 0)
			j = obj[j].next_measurement;
		else 
			all_done = 1;
	}
	fclose(measurements);

	//adjust wavelength limits
	//lambda_min = lambda_min - (lambda_max-lambda_min)*.1;
	lambda_min *= .5;
	if ( lambda_min < 1 ) 
		lambda_min = 1;
	lambda_max = lambda_max + (lambda_max-lambda_min)*3.;

	//for one-wavelength data
	if ( lambda_min == lambda_max )
	{
		lambda_min *=0.5;
		lambda_max *=2;
	}


	printf("plot SED from %.2f to %.2f\n", lambda_min, lambda_max);
	
	//set wavelength range manually
	lambda_min = 1;
	lambda_max = 25;

	char fluxfilename[250];
	char fluxfilenameSI[250];
	sprintf(fluxfilename, "%s.flux", obj[j].name);
	sprintf(fluxfilenameSI, "%s.fluxSI", obj[j].name);

	//-----------------------------------
	//create files containing flux data
	if ( obj != NULL ) 
	{
		FILE *fluxdata = fopen(fluxfilename, "wt");
		FILE *fluxdataSI = fopen(fluxfilenameSI, "wt");
		
		double lambda;
		for ( lambda = lambda_min; lambda <= lambda_max; lambda+=(lambda_max-lambda_min)/1000 )
		{
			/* fluxdata file output:
			 *   wavelength [microns]
			 *   flux       [mJy]
			*/
			obj[i].wavelength = lambda;
			fprintf(fluxdata, "%.5e   %.5e\n", lambda, planck_model(obj[i], epsilon)*3.335640952e14*obj[i].wavelength*obj[i].wavelength);
			fprintf(fluxdataSI, "%.5e   %.5e\n", lambda, planck_model(obj[i], epsilon));
		}
		fclose(fluxdata);
		fclose(fluxdataSI);
	}

}


/* create datapoints for surface plot of chi2(diameter, eta)
 * can be displayed in gnuplot using:
 *       splot 'eta_surfaceplot' using 1:2:3 with lines
 * adapt zrange to resolve minimum
*/ 
void eta_surfaceplot ( ASTEROID *obj, int j, double eta_min, double eta_max, 
		double eta_step, double diameter_min, double diameter_max,
		double diameter_step, double epsilon, int fitting,
		double (*planck_model)(ASTEROID, double) )
{
	FILE *eta_surfaceplot = fopen("eta_surfaceplot", "wt");
	
	double eta, diameter;
	for ( eta = eta_min; eta <= eta_max; eta += eta_step )
	{
		for ( diameter = diameter_min; diameter <= diameter_max; 
			diameter += diameter_step )
		{
			if ( fitting == 1 )
			{
				int all_done = 0, i = j;
				while ( all_done == 0 ) 
				{
					obj[i].eta = eta;
					
					//breakout if all measurements concerning this
					//object are processed
					if ( obj[i].next_measurement != 0)
						i = obj[i].next_measurement;
					else 
						all_done = 1;
				}
			}
			fprintf(eta_surfaceplot, "%.5e   %.5e   %.5e\n", diameter, eta, 
				chi2sum_object(obj, j, diameter, epsilon, fitting, planck_model));
		}
		fprintf(eta_surfaceplot, "\n");
	}
	
	fclose(eta_surfaceplot);
}

//-------------------------------------------------------------------------------

void runmodel(ASTEROID *obj, int start_asteroid, int number_of_objects, int fitting, double epsilon, double eta, char *resultsfile, int display, char *notes)
{
	int j;
	
	//string containing names of all objects processed so far
	char *objects_processed_so_far = malloc(sizeof(char)*200*EXP_NUMBER_OF_OBJECTS);
	
	//create filename for long results-files
	char *longresults = malloc(sizeof(char)*200);
	strcpy(longresults, "");  //clear string
	strcat(longresults, resultsfile);
	strcat(longresults, "_long");
	
	//string, which is appended to results-file line
	char catstring[INPUT_LINE_LENGTH*2];
	strcpy(catstring, "");  //clear string
	
	//-------
	//STM calculations
	if ( strpbrk(model, "1") != NULL )  //use refined STM
	{
		for ( j = start_asteroid; j < start_asteroid+number_of_objects; j++ )
		{
			if ( strstr (objects_processed_so_far, obj[j].name) == NULL)
			{
				strcpy(catstring, "STM");
				strcat(catstring, notes);
				
				find_diameter ( obj, j, epsilon, fitting, planck_stm );
				if ( fitting == 1 )
				{
					strcat (objects_processed_so_far, obj[j].name);
					checkforflags(obj, j);
				}
				if ( display == 1 )
					printf("%s (STM, fixed eta): d=%.3f, pv=%.3f, eta=%.3f, chi2=%.3e\n", obj[j].name, 
					obj[j].diameter, obj[j].pv, obj[j].eta, obj[j].chi2);
				
				if ( strstr (resultsfile, "no_output") == NULL )
				  {
				    if ( short_output_in_file(obj[j], resultsfile, catstring) == 1 )
				      printf("ERROR: no output to file possible!\n");
				  }
			}
			
			//tabular results in file
			if ( strstr (resultsfile, "no_output") == NULL )
			  {
			    if ( long_output_in_file(obj[j], longresults, catstring) == 1 )
			      printf("ERROR: no output to file possible!\n");
			  }
		}
	}
	
	strcpy(catstring, "");                 //clear string
	strcpy(objects_processed_so_far, "");  //clear string
	
	//-------
	//FRM calculations
	if ( strpbrk(model, "2") != NULL )  //use FRM
	{
		for ( j = start_asteroid; j < start_asteroid+number_of_objects; j++ )
		{
			if ( strstr (objects_processed_so_far, obj[j].name) == NULL)
			{
				strcpy(catstring, "FRM");
				strcat(catstring, notes);
				
				find_diameter ( obj, j, epsilon, fitting, planck_frm );
				if ( fitting == 1 )
				{
					strcat (objects_processed_so_far, obj[j].name);
					checkforflags(obj, j);
				}
				if ( display == 1 )
					printf("%s (FRM, fixed eta): d=%.3f, pv=%.3f, eta=%.3f, chi2=%.3e\n", obj[j].name, 
					obj[j].diameter, obj[j].pv, obj[j].eta, obj[j].chi2);
				
				if ( strstr (resultsfile, "no_output") == NULL )
				  {
				    if ( short_output_in_file(obj[j], resultsfile, catstring) == 1 )
				      printf("ERROR: no output to file possible!\n");
				  }
			}
			
			//tabular results in file
			if ( strstr (resultsfile, "no_output") == NULL )
			  {
			    if ( long_output_in_file(obj[j], longresults, catstring) == 1 )
			      printf("ERROR: no output to file possible!\n");
			  }
		}
	}
	
	strcpy(catstring, "");                 //clear string
	strcpy(objects_processed_so_far, "");  //clear string
	
	//-------
	//NEATM calculations (floating eta)
	if ( (strpbrk(model, "3") != NULL) )  //use NEATM (floating eta)
	{
		for ( j = start_asteroid; j < start_asteroid+number_of_objects; j++ )
		{
			if ( strstr (objects_processed_so_far, obj[j].name) == NULL )
			{
				strcpy(catstring, "NEATM(floating_eta)");
				strcat(catstring, notes);
				
				fit_eta(obj, j, epsilon, 1.0, planck_neatm);
				if ( fitting == 1 )
				{
					strcat (objects_processed_so_far, obj[j].name);
					checkforflags(obj, j);
				}
				if ( display == 1 )
					printf("%s (NEATM, floating eta): d=%.3f, pv=%.3f, eta=%.3f, chi2=%.3e\n", obj[j].name, 
					obj[j].diameter, obj[j].pv, obj[j].eta, obj[j].chi2);
					
				if ( strstr (resultsfile, "no_output") == NULL )
				  {
				    if ( short_output_in_file(obj[j], resultsfile, catstring) == 1 )
				      printf("ERROR: no output to file possible!\n");
				  }
			}
			
			//tabular results in file
			if ( strstr (resultsfile, "no_output") == NULL )
			  {
			    if ( long_output_in_file(obj[j], longresults, catstring) == 1 )
			      printf("ERROR: no output to file possible!\n");
			  }
		}
	}
	
	strcpy(catstring, "");                 //clear string
	strcpy(objects_processed_so_far, "");  //clear string
	
	//-------
	//NEATM calculations (fixed eta)
	if ( (strpbrk(model, "4") != NULL) )  //use NEATM (fixed eta)
	{
		for ( j = start_asteroid; j < start_asteroid+number_of_objects; j++ )
		{
		  char needle[100];
		  sprintf(needle, "%s|", obj[j].name);
			if ( strstr (objects_processed_so_far, needle) == NULL)
			{
				strcpy(catstring, "NEATM(fixed_eta)");
				strcat(catstring, notes);
				
				find_diameter(obj, j, epsilon, fitting, planck_neatm);
				if ( fitting == 1 )
				{
					strcat (objects_processed_so_far, needle);
					checkforflags(obj, j);
				}
				if ( display == 1 )
					printf("%s (NEATM, fixed eta): d=%.3f, pv=%.3f, eta=%.3f, chi2=%.3e\n", obj[j].name, 
					obj[j].diameter, obj[j].pv, obj[j].eta, obj[j].chi2);
				
				if ( strstr (resultsfile, "no_output") == NULL )
				  {
				    if ( short_output_in_file(obj[j], resultsfile, catstring) == 1 )
				      printf("ERROR: no output to file possible!\n");
				  }
			}

			//tabular results in file
			if ( strstr (resultsfile, "no_output") == NULL )
			  {
			    if ( long_output_in_file(obj[j], longresults, catstring) == 1 )
			      printf("ERROR: no output to file possible!\n");
			  }
		}
	}
	
	strcpy(catstring, "");                 //clear string
	strcpy(objects_processed_so_far, "");  //clear string
	
	
	strcpy(catstring, "");                 //clear string
	strcpy(objects_processed_so_far, "");  //clear string
	
	//-------
	//combo calculations (STM, FRM, NEATM)
	if ( (strpbrk(model, "c") != NULL) ) 
	{
		printf("\ncombined calculations:\n----------------------\n");
		
		//------------------------- STM
		model = "1";
		for ( j = start_asteroid; j < start_asteroid+number_of_objects; j++ )
			obj[j].eta = eta;
		
		runmodel(obj, start_asteroid, number_of_objects, fitting, epsilon, 0, "results_STM", 1, notes);
		
		
		//------------------------- STM
		model = "2";
		for ( j = start_asteroid; j < start_asteroid+number_of_objects; j++ )
			obj[j].eta = PI;
		
		runmodel(obj, start_asteroid, number_of_objects, fitting, epsilon, 0, "results_FRM", 1, notes);
		
		
		//------------------------- NEATM
		model = "3";
		
		runmodel(obj, start_asteroid, number_of_objects, fitting, epsilon, 0, "results_NEATM", 1, notes);
		
		model = "c";
	}
	
	free (objects_processed_so_far);
}

/*estimates expected flux for observations
*/
void estimate_flux (ASTEROID *obj, double number_of_objects, double epsilon)
{
	int i;
	// remove existing latex and model input files
	remove("latex");
	remove("model.input");

	FILE *latexfile = fopen ("latex", "at");
	FILE *modelfile = fopen ("model.input", "at");
	
	fprintf(latexfile, "\\centering\n\\begin{tabular}{r c c c c c c c c}\n  \\hline\\hline \\\\ \n"
	       "  Object &  r (AU)& $\\Delta$ (AU)& $\\alpha$ (deg) & V mag & $\\lambda$ ($\\mu$m) & Flux (mJy) & SNratio \\\\ \n"
	       "  \\hline\\\\ \n");
	
	for (i=0; i<number_of_objects; i++)
	{
		//calculate thermal flux
		if ( (strpbrk(model, "1") != NULL) ) //STM
			obj[i].flux_raw = planck_stm(obj[i], epsilon);
		if ( (strpbrk(model, "2") != NULL) ) //FRM
			obj[i].flux_raw = planck_frm(obj[i], epsilon);
		if ( (strpbrk(model, "34") != NULL) ) //NEATM
			obj[i].flux_raw = planck_neatm(obj[i], epsilon);

		// color correction
		double pv = pv_from_D(obj[i].diameter, obj[i].absmag);
		double T0 = pow((1-bondalbedo(obj[i].slopepar, pv))*SOLARCONST/(obj[i].heliodist*obj[i].heliodist*obj[i].eta*EMISSIVITY*STEFANBOLTZMANN), 0.25);

		obj[i].colorcorrection = 1.;

		// PACS COLOR CORRECTION
		if ( obj[i].wavelength == 70 )
		  {
		    T0 = T0*pow(0.5,0.25); // TNOs are Cool uses mean surface temperature
		    obj[i].colorcorrection = 4.71879e-8*T0*T0*T0*T0-1.19493e-5*T0*T0*T0+0.00113054*T0*T0-0.0473228*T0+1.71981;
		  }
		else if ( obj[i].wavelength == 100 )
		  {
		    T0 = T0*pow(0.5,0.25); // TNOs are Cool uses mean surface temperature
		    obj[i].colorcorrection = 1.66534e-8*T0*T0*T0*T0-4.07308e-6*T0*T0*T0+0.000363925*T0*T0-0.0135295*T0+1.1566;
		  }
		else if ( obj[i].wavelength == 160 )
		  {
		    T0 = T0*pow(0.5,0.25); // TNOs are Cool uses mean surface temperature
		    obj[i].colorcorrection = 1.03868e-9*T0*T0*T0*T0-3.04287e-8*T0*T0*T0-3.10479e-5*T0*T0+0.00404142*T0+0.882475;	
		  }
		
		// IRAC COLOR CORRECTION
		// ExploreNEOs uses subsolar temperature
		else if ( (obj[i].wavelength == 3.55) || (obj[i].wavelength == 3.545) || (obj[i].wavelength == 3.6))
		  obj[i].colorcorrection = -7.6839e-08*T0*T0*T0 + 8.14039e-05*T0*T0 - 0.0297487*T0 + 4.85495;
		//obj[i].colorcorrection = IRACcolorcorr(obj, i, 1, epsilon);
		else if ( (obj[i].wavelength == 4.493) || (obj[i].wavelength == 4.5) )
		  obj[i].colorcorrection = -4.93676e-08*T0*T0*T0 + 5.15348e-05*T0*T0 - 0.0185212*T0 + 3.33937;
		//obj[i].colorcorrection = IRACcolorcorr(obj, i, 2, epsilon);
		
		// MIPS COLOR CORRECTION
		// uses subsolar temperature
		else if ( (obj[i].wavelength == 24) || (obj[i].wavelength == 23.68) || (obj[i].wavelength == 23.675) )
		  obj[i].colorcorrection = 6.48111e+07/(T0*T0*T0*T0*T0) - 5.45981e+06/(T0*T0*T0*T0) + 206351/(T0*T0*T0) - 2567.6/(T0*T0) + 11.3144/T0 + 0.93459;
		//obj[i].colorcorrection = MIPScolorcorr(obj, i, 24, epsilon);
		else if ( (obj[i].wavelength == 71.42) || (obj[i].wavelength == 71.44) )
		  obj[i].colorcorrection = 1824.15/(T0*T0*T0) + 95.5894/(T0*T0) - 8.51701/T0 + 4.01837e-10*T0*T0 - 8.5828e-06*T0 + 1.01135;
		//obj[i].colorcorrection = MIPScolorcorr(obj, i, 70, epsilon);
		
		// IRS PEAK--UP COLOR CORRECTION (only blue channel, approximate values derived from manual)
		// uses subsolar temperature
		else if ( obj[i].wavelength == 16 )
		  obj[i].colorcorrection = 1.30223e+10/(T0*T0*T0*T0*T0) - 4.50356e+08/(T0*T0*T0*T0) + 5.86753e+06/(T0*T0*T0) - 28405.3/(T0*T0) + 31.223/T0 + 1.02542;
		
		// WISE COLOR CORRECTION
		// uses subsolar temperature
		else if ( (obj[i].wavelength == 3.3526) )
		  obj[i].colorcorrection =  7.94861e+09/(T0*T0*T0*T0) - 9.41457e+07/(T0*T0*T0) + 6012.43/T0 - 3.73172e-11*T0*T0*T0*T0 + 1.27526e-07*T0*T0*T0 - 0.000170268*T0*T0 + 0.113257*T0 - 38.0054; 
		else if ( (obj[i].wavelength == 4.6028) )
		  obj[i].colorcorrection  = 6.25545e+06/(T0*T0*T0) - 71263.6/(T0*T0) + 535.888/T0 + 8.03351e-10*T0*T0*T0 - 2.56492e-06*T0*T0 + 0.0031402*T0 - 0.867562;
		else if ( (obj[i].wavelength == 11.5608) )
		  obj[i].colorcorrection  = 3.63093e+06/(T0*T0*T0) - 36438.8/(T0*T0) + 279.812/T0 + 6.43141e-10*T0*T0*T0 - 2.02289e-06*T0*T0 + 0.00240654*T0 - 0.347453;
		else if ( (obj[i].wavelength == 22.0883) )
		  obj[i].colorcorrection  = 1075.28/(T0*T0) - 12.6009/T0 + 1.09716e-08*T0*T0 - 2.68914e-05*T0 + 1.02422;  
		

		// reverse color correction (despite this line, fluxes in the fluxest.results file are still monochromatic!)
		obj[i].flux = obj[i].colorcorrection * obj[i].flux_raw;

				
		//flux error estimation = 10%
		obj[i].fluxerr = obj[i].flux_raw*.1;
		
		//calculate vmag
		double v_alpha = obj[i].absmag - 2.5 
					* log10( (1-obj[i].slopepar)* exp( -3.33*pow(tan(obj[i].phaseangle/2),0.63))
				                  + obj[i].slopepar * exp( -1.87*pow(tan(obj[i].phaseangle/2),1.22)) );
		double vmag = v_alpha + 5*log10(obj[i].heliodist*obj[i].geodist);
		
		//determine SN-ratio
		double snratio = 0;
		if ( obj[i].wavelength == 8.7 ) //MIRSI 8.7um
			snratio = obj[i].flux*3.335640952e14*obj[i].wavelength*obj[i].wavelength/14.7;
		else if ( obj[i].wavelength == 11.6 ) //MIRSI 11.6um
			snratio = obj[i].flux*3.335640952e14*obj[i].wavelength*obj[i].wavelength/22.0;
		else if ( obj[i].wavelength == 18.4 ) //MIRSI 18.4um
			snratio = obj[i].flux*3.335640952e14*obj[i].wavelength*obj[i].wavelength/57.3;
		else
			snratio = 0; //sqrt(obj[i].flux*3.335640952e14*obj[i].wavelength*obj[i].wavelength);
		
		
		//write results into file
		if ( long_output_in_file(obj[i], "fluxest.results", "") == 1 )
				printf("ERROR: no output to file possible!\n");
		
		fprintf(latexfile, "  %s & %.3f & %.3f & %.1f & %.2f & %.1f & %.1f & %.1f \\\\\n", 
			obj[i].name, obj[i].heliodist, obj[i].geodist, obj[i].phaseangle/PI2*90, vmag,
			obj[i].wavelength, obj[i].flux*3.335640952e14*obj[i].wavelength*obj[i].wavelength, snratio );

		//model input file output
		fprintf(modelfile, "%15.15s   %.3f   %.3f   %.1f   %.2f   %.2f   %.2f   %.3f  0.1 \n", obj[i].name, 
			obj[i].geodist, obj[i].heliodist, obj[i].phaseangle/PI2*90, obj[i].absmag,
			obj[i].slopepar, obj[i].wavelength, obj[i].flux*3.335640952e14*obj[i].wavelength*obj[i].wavelength);
	}
	
	//latex screen output
	fprintf(latexfile, "  \\hline\\hline \\\\ \n\\end{tabular}\n");

	fclose(latexfile);
	fclose(modelfile);
}




//-------------------------------------------------------------------------------

int main(int argc, char** argv)
{
	//definition of internal variables
	char *inputfile;
	double diameter_global, pv_global, eta_global, modelparameter_global;
	ASTEROID neo[EXP_NUMBER_OF_OBJECTS];
	int number_of_objects;
	int dothetpm = 0;           //=1: perform TPM calculations with fixed parameters
	                            //=2: perform TPM calculations with random parameters
	int fitting = 1;            //=1: use all measurements for chi^2-fitting
				    //=0: process each measurement separately
	int plotting = 0;           //plot measurements and best fit results
	int corrflux = 1;           //apply flux correction
	int dotheerror = 0;          //perform error estimation
	double reflectanceratio = 1.6; // global solar reflectance ratio (Mainzer et al. 2011)
	int chi2montecarlo = 0;	    //adjust flux errors in monte carlo method
	int inputfilestyle = 0;     //0: nothing special about input file
	int i, j;
	
	double epsilon = 1e-3;  //desired relative accuracy of flux calculations
	
	//delete existing output files
	remove ("results");
		
	
	//allocate memory for ASTEROID names
	for ( i = 0; i < EXP_NUMBER_OF_OBJECTS; i++ )
	{
		neo[i].name = malloc(sizeof(char)*INPUT_LINE_LENGTH);
		neo[i].diameter_sigma = neo[i].pv_sigma = neo[i].eta_sigma = neo[i].modelparameter_sigma = 0;
		neo[i].diameter_median = neo[i].pv_median = neo[i].eta_median = neo[i].modelparameter_median = 0;
		neo[i].diameter_lower = neo[i].pv_lower = neo[i].eta_lower = neo[i].modelparameter_lower = 0;
		neo[i].diameter_upper = neo[i].pv_upper = neo[i].eta_upper = neo[i].modelparameter_upper = 0;
		neo[i].diameter_lower3sig = neo[i].pv_lower3sig = neo[i].eta_lower3sig = 0;
		neo[i].diameter_upper3sig = neo[i].pv_upper3sig = neo[i].eta_upper3sig = 0;
	}
	
	//check for command line arguments
	if ( argc < 3 ) 
	{	printf( "ERROR: wrong number of arguments (%d)!\ntry this syntax:\n"
			"model input models [-e value] [-s]\n"
			"  input          data input file\n"
			"  models          models to be used for calculation\n"
			"                  1: ref. STM, 2: FRM, 3: NEATM (floating eta),\n"
			"                  4: NEATM (fixed eta) [needs d/w/f]\n"
			"                  in case of '4': f: fixed eta, d: Delbo 2004, w: Wolters 2007, m: Mainzer 2011, n: neosurvey\n"
			"                  (if none of the above: 4f is used)\n"
			"  [-tpm fix/rand] perform TPM calculations with fix/random parameters\n"
			"  [-eta value]    fix eta value\n"
			"  [-sslat value]  fix sslat value\n"
			"  [-fc]           flux correction (for reflected sunlight)\n"
			"  [-rr value]     solar reflectance ratio for all wavelengths\n"
			"  [-ncc]          no color correction (input fluxes are thermal fluxes)\n"
			"  [-montecarlo]   perform error estimation\n"
			"  [-chi2montecarlo] perform error estimation and adjust flux errors\n"
			"  [-er min max]   eta range\n"
			"  [-s]            calculate each measurement separately\n"
			"  [-p]            plot data\n"
                        "  [-tno]          uses Brucker-relation and TNO D-pV-H\n",
			(argc-1));
		return 0;
	}
	
	//initializing variables from command line arguments
	//-----------------------------------
	inputfile = argv[1];   //filename of input file
	model = argv[2];       //model to be applied
	pv_global = 0.01;      //arbitrary albedo value
	diameter_global = 0.1; //arbitrary diameter
	eta_global = 0.;       //arbitrary delta_eta
	modelparameter_global = 0;//preset modelparameter=0 (means: sun and observer in equator plane)
		
	//check for optional arguments
	for ( i = 2 ; i < argc ; i++)
	{
		if ( strcmp(argv[i], "-tpm") == 0 )
		{
			if ( strpbrk(argv[i+1], "f") != 0 )
				dothetpm = 1;
			if ( strpbrk(argv[i+1], "r") != 0 )
				dothetpm = 2;
		}
		
		if ( strcmp(argv[i], "-eta") == 0 )
			eta_global = atof(argv[i+1]);
		
		if ( strcmp(argv[i], "-sslat") == 0 )
			modelparameter_global = atof(argv[i+1]);
		
		if ( strcmp(argv[i], "-p") == 0 ) 
			plotting = 1;
		
		if ( strcmp(argv[i], "-s") == 0 ) 
			fitting = 0;
		
		/* if ( strcmp(argv[i], "-fc") == 0 ) */
		/* 	corrflux = 1; */
		
		if ( strcmp(argv[i], "-er") == 0 )
		{
			eta_min = atof(argv[i+1]);
			eta_max = atof(argv[i+2]);
		}
		
		if ( strcmp(argv[i], "-ncc") == 0 )
			colcorr = 0;
			
		if ( strcmp(argv[i], "-montecarlo") == 0 )
			dotheerror = 1;

		if ( strcmp(argv[i], "-chi2montecarlo") == 0 )
		{
			dotheerror = 1;
			chi2montecarlo = 1;
		}		

		if ( strcmp(argv[i], "-tno") == 0)
			TNO = 1;

		inputfilestyle = 1;

		if ( strcmp(argv[i], "-rr") == 0 )
			reflectanceratio = atof(argv[i+1]);
	}
	
	//-------------------------------------------------------------------------------------------------------------------------
	
	//read in input file 
	//-----------------------------------
	//for usual calculations
	if ( dothetpm == 0 )
	{
	  //printf("read in data from file %s\n", inputfile);
		number_of_objects = readinput ( inputfile, neo, eta_global, diameter_global, pv_global, modelparameter_global, 
		                                corrflux, fitting, reflectanceratio, inputfilestyle );
	}
		
	//check for correct read-in process
	if ( number_of_objects == 0 )
		return EXIT_FAILURE;
	
	
	//------------------------------------------------------------------
	//run calculations

	if ( strstr(inputfile, "fluxest") )
		estimate_flux(neo, number_of_objects, epsilon);
	else if ( dotheerror == 1 )
		error_estimation(neo, number_of_objects, 10000, 1., fitting, epsilon, corrflux, chi2montecarlo, plotting);
	else
		if ( dothetpm == 0 )
			runmodel(neo, 0, number_of_objects, fitting, epsilon, eta_global, "results", 1, "");

	//-----------------------------------
	//present results
	
	//string containing names of all objects processed so far
	char *objects_processed_so_far = malloc(sizeof(char)*40*EXP_NUMBER_OF_OBJECTS);
	
	//plot fits
	if ( plotting == 1 )
	{
		printf("\nplotting data:\n");
		for ( j = 0; j < number_of_objects; j++ )
			{
				if ( fitting == 1 )
				{
					if ( strstr (objects_processed_so_far, neo[j].name) == NULL)
					{
						if ( (strpbrk(model, "1") != NULL) )
							plotmodel(neo, j, epsilon, planck_stm);
						if ( (strpbrk(model, "2") != NULL) )
							plotmodel(neo, j, epsilon, planck_frm);
						if ( (strpbrk(model, "34") != NULL) )
							plotmodel(neo, j, epsilon, planck_neatm);
						
						strcat (objects_processed_so_far, neo[j].name);
					}
				} else
				{
					if ( (strpbrk(model, "1") != NULL) )
						plotmodel(neo, j, epsilon, planck_stm);
					if ( (strpbrk(model, "2") != NULL) )
						plotmodel(neo, j, epsilon, planck_frm);
					if ( (strpbrk(model, "34") != NULL) )
						plotmodel(neo, j, epsilon, planck_neatm);
				}
			}
	}
	
	//clear neo names and next_measurements
	for ( i = 0; i < EXP_NUMBER_OF_OBJECTS; i++ )
	{
		strcpy(neo[i].name, "");
		neo[i].next_measurement = 0;
	}
	
	
	return 0;
}
