/*
 *      FRM model routines
 *      Version 1.2 - Jan 28, 2010
 *      
 *      Michael Mommert, DLR
 *      michael.mommert@dlr.de
 * 
 * 
 *      expected relative flux errors (=imprecision) of flux integration 
 *        and fitting = flux error of measurement
 */

//Planck function integrand
//if function only one-dimensional, use only latitude
//no use for 'parameter'
double integrand_frm (double latitude, double eps, ASTEROID obj, double T0)
{
	double expfct = exp(HCK/(obj.wavelength*T0*pow(cos(latitude),0.25)));
	
	if ( isnan(expfct) )    //NaN results from division by zero in expfct
		expfct = 1e50;  //to avoid division by zero
	
	return (cos(latitude)*cos(latitude)/(expfct-1));
}


//complete flux calculation in W m^-2 um^-1 using romberg integration
//no use for 'parameter'
double planck_frm( ASTEROID obj, double epsilon )
{
	double pv = pv_from_D(obj.diameter, obj.absmag);
	double T0 = pow((1-bondalbedo (obj.slopepar, pv))*SOLARCONST/(obj.heliodist*obj.heliodist
			*PI*EMISSIVITY*STEFANBOLTZMANN), 0.25);
	
	return (2*EMISSIVITY*obj.diameter/AU*obj.diameter/AU*HC2
			/(obj.geodist*obj.geodist*obj.wavelength*obj.wavelength*obj.wavelength
			  *obj.wavelength*obj.wavelength)  //um conversion to be found in HC2 
			*romberg(integrand_frm, 0, PI2, epsilon, obj, T0 )); 
}


//create a temperature distribution plot for an ASTEROID object
void tempplot_frm (ASTEROID obj)
{
	FILE *tempfile = fopen("temp", "wt");
	double phi, theta, temp;
	
	double pv = pv_from_D(obj.diameter, obj.absmag);
	double T0 = pow((1-bondalbedo (obj.slopepar, pv))*SOLARCONST/(obj.heliodist*obj.heliodist
			*PI*EMISSIVITY*STEFANBOLTZMANN), 0.25);
	
	//coordinates in body centered system
	for (theta=-PI2; theta<=PI2; theta+=PI/180)
	{
		temp = pow(fabs(cos(theta)),0.25);
		
		for (phi=0; phi<=2*PI; phi+=PI/180)
			fprintf(tempfile, "%.3f   %.3f   %.3f\n", phi, theta, T0*temp); 
		
		fprintf(tempfile, "\n");
		
		//printf("theta %.3e, cos(theta) %.3e, pow %.3e\n", theta, cos(theta), temp);
	}
	
	fclose(tempfile);
}
