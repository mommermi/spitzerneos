/*
 *      STM model routines
 *      Version 1.2 - Jan 28, 2010
 *      
 *      Michael Mommert, DLR
 *      michael.mommert@dlr.de
 * 
 * 
 *      expected relative flux errors (=imprecision) of flux integration 
 *        and fitting = flux error of measurement
 */

//empirical phase correction factor for STM
double phase_correction_factor( double phase )
{
	return (pow(10, -0.01*fabs(phase/PI2*90)*0.4));

} 

//Planck function integrand
//if function only one-dimensional, use only latitude
//no use for 'parameter'
double integrand_stm (double latitude, double eps, ASTEROID obj, double T0)
{
	double expfct = exp(HCK/(obj.wavelength*T0*pow(cos(latitude),0.25)));
	
	if ( isnan(expfct) )    //NaN results from division by zero in expfct
		expfct = 1e30;  //to avoid division by zero
	
	return (cos(latitude)*sin(latitude)/(expfct-1));
}


//complete flux calculation in W m^-2 um^-1 using romberg integration
//no use for 'parameter'
double planck_stm(ASTEROID obj, double epsilon)
{
	//calculate maximum temperature
	// NEO: double pv = (1329*1329)*pow(10,(-.4)*obj.absmag)
	//			/(obj.diameter*obj.diameter);
	// TNO:
	double pv = (1342.6*1342.6)*pow(10,(-.4)*obj.absmag)
				/(obj.diameter*obj.diameter);
	double bond_albedo = (0.29 + 0.684*obj.slopepar)*pv;
	double T0 = pow((1-bond_albedo)*SOLARCONST/(obj.heliodist*obj.heliodist
			*obj.eta*EMISSIVITY*STEFANBOLTZMANN), 0.25);
		
	//printf("starte romberg(integrand_stm,0,PI2...)\n");
	return (EMISSIVITY*obj.diameter/AU*obj.diameter/AU*HC2*PI    //um conversion to be found in HC2 
			/(obj.geodist*obj.geodist*obj.wavelength*obj.wavelength*obj.wavelength*obj.wavelength*obj.wavelength)
			*phase_correction_factor( obj.phaseangle )  //phase correction for STM
			*romberg(integrand_stm, 0, PI2, epsilon, obj, T0)); 
}

//create a temperature distribution plot for an ASTEROID object
void tempplot_stm (ASTEROID obj)
{
	FILE *tempfile = fopen("temp", "wt");
	double phi, theta, temp;
	
	//calculate maximum temperature
	double pv = (1329*1329)*pow(10,(-.4)*obj.absmag)/(obj.diameter*obj.diameter);
	double bond_albedo = (0.29 + 0.684*obj.slopepar)*pv;
	double T0 = pow((1-bond_albedo)*SOLARCONST/(obj.heliodist*obj.heliodist
		*obj.eta*EMISSIVITY*STEFANBOLTZMANN), 0.25);
	
	printf("%.3f mit %.3f km\n", T0, obj.diameter);
	
	for (theta=-PI2; theta<=PI2; theta+=PI/180)
	{
		for (phi=0; phi<=2*PI; phi+=PI/180)
		{
			if ( sin(theta) > 0 )
				temp =  T0*pow(sin(theta),0.25);
			else
				temp = 0;
			fprintf(tempfile, "%.3f   %.3f   %.3f\n", phi, theta, temp);
		
		}
		fprintf(tempfile, "\n");
	}
	
	fclose(tempfile);
}

