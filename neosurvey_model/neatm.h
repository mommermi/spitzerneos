/*
 *      NEATM model routines
 *      Version 1.2 - Jan 28, 2010
 *      
 *      Michael Mommert, DLR
 *      michael.mommert@dlr.de
 * 
 * 
 *      expected relative flux errors (=imprecision) of flux integration 
 *        and fitting = flux error of measurement
 */

static double latitude_static;

static double longitude_upperlimit, longitude_lowerlimit;


/* if using double-romberg
   use then in planck_neatm: return (EMISSIVIT...romberg(aux_integrand_neatm, 0, PI2, eps..)
*/
double integrand_neatm (double longitude, double eps, ASTEROID obj, double T0)
{
        double local_T = T0 * pow(cos(latitude_static)*cos(longitude),0.25);

	double expfct = exp(HCK/(obj.wavelength*T0*pow(cos(latitude_static)*cos(longitude),0.25)));

	double ccc = 1;

	if (colcorr == 1) 
	  {
	    // PACS COLOR CORRECTION
	    
	    if ( obj.wavelength == 70 )
		ccc = 4.71879e-8*local_T*local_T*local_T*local_T-1.19493e-5*local_T*local_T*local_T+0.00113054*local_T*local_T-0.0473228*local_T+1.71981;
	    else if ( obj.wavelength == 100 )
		ccc = 1.66534e-8*local_T*local_T*local_T*local_T-4.07308e-6*local_T*local_T*local_T+0.000363925*local_T*local_T-0.0135295*local_T+1.1566;
	    else if ( obj.wavelength == 160 )
		ccc = 1.03868e-9*local_T*local_T*local_T*local_T-3.04287e-8*local_T*local_T*local_T-3.10479e-5*local_T*local_T+0.00404142*local_T+0.882475;	

	    // IRAC COLOR CORRECTION
	    // ExploreNEOs uses subsolar temperature
	    else if ( (obj.wavelength == 3.55) || (obj.wavelength == 3.545) || (obj.wavelength == 3.6))
	      ccc = -7.6839e-08*local_T*local_T*local_T + 8.14039e-05*local_T*local_T - 0.0297487*local_T + 4.85495;
	    //ccc = IRACcolorcorr(obj, i, 1, epsilon);
	    else if ( (obj.wavelength == 4.493) || (obj.wavelength == 4.5) )
	      ccc = -4.93676e-08*local_T*local_T*local_T + 5.15348e-05*local_T*local_T - 0.0185212*local_T + 3.33937;
	    //ccc = IRACcolorcorr(obj, i, 2, epsilon);
	    
	    
	    // MIPS COLOR CORRECTION
	    // uses subsolar temperature
	    else if ( (obj.wavelength == 24) || (obj.wavelength == 23.68) || (obj.wavelength == 23.675) )
	      ccc = 6.48111e+07/(local_T*local_T*local_T*local_T*local_T) - 5.45981e+06/(local_T*local_T*local_T*local_T) + 206351/(local_T*local_T*local_T) - 2567.6/(local_T*local_T) + 11.3144/local_T + 0.93459;
	    //ccc = MIPScolorcorr(obj, i, 24, epsilon);
	    else if ( (obj.wavelength == 71.42) || (obj.wavelength == 71.44) )
	      ccc = 1824.15/(local_T*local_T*local_T) + 95.5894/(local_T*local_T) - 8.51701/local_T + 4.01837e-10*local_T*local_T - 8.5828e-06*local_T + 1.01135;
	    //ccc = MIPScolorcorr(obj, i, 70, epsilon);
	    
	    // IRS PEAK--UP COLOR CORRECTION (only blue channel, approximate values derived from manual)
	    // uses subsolar temperature
	    else if ( obj.wavelength == 16 )
	      ccc = 1.30223e+10/(local_T*local_T*local_T*local_T*local_T) - 4.50356e+08/(local_T*local_T*local_T*local_T) + 5.86753e+06/(local_T*local_T*local_T) - 28405.3/(local_T*local_T) + 31.223/local_T + 1.02542;
	    
	    
	    // WISE COLOR CORRECTION
	    // uses subsolar temperature
	    else if ( (obj.wavelength == 3.3526) )
	      ccc =  7.94861e+09/(local_T*local_T*local_T*local_T) - 9.41457e+07/(local_T*local_T*local_T) + 6012.43/local_T - 3.73172e-11*local_T*local_T*local_T*local_T + 1.27526e-07*local_T*local_T*local_T - 0.000170268*local_T*local_T + 0.113257*local_T - 38.0054; 
	    else if ( (obj.wavelength == 4.6028) )
	      ccc  = 6.25545e+06/(local_T*local_T*local_T) - 71263.6/(local_T*local_T) + 535.888/local_T + 8.03351e-10*local_T*local_T*local_T - 2.56492e-06*local_T*local_T + 0.0031402*local_T - 0.867562;
	    else if ( (obj.wavelength == 11.5608) || (obj.wavelength == 11.0984) )
	      ccc  = 3.63093e+06/(local_T*local_T*local_T) - 36438.8/(local_T*local_T) + 279.812/local_T + 6.43141e-10*local_T*local_T*local_T - 2.02289e-06*local_T*local_T + 0.00240654*local_T - 0.347453;
	    else if ( (obj.wavelength == 22.0883) || (obj.wavelength == 22.6405) )
	      ccc  = 1075.28/(local_T*local_T) - 12.6009/local_T + 1.09716e-08*local_T*local_T - 2.68914e-05*local_T + 1.02422;  
	    
	    // APEX COLOR CORRECTION (assume uniform transmission)
	    else if ( (obj.wavelength == 6.001) )
	      ccc  = 2.84433e+08/(local_T*local_T*local_T) - 4.00502e+06/(local_T*local_T) + 24256.1/local_T + 3.05563e-08*local_T*local_T*local_T - 9.6426e-05*local_T*local_T + 0.121391*local_T - 71.5251;
	    
	    
	    else 
	      {
		if (colcorr_error == 0) 
		  {
		    //printf ("WARNING: no color correction possible for %s at %f micron\n", obj[i].name, obj.wavelength);
		    colcorr_error = 1;
		  }
		//abort();
	      }
	  }



	
	if ( isnan(expfct) )    //NaN results from division by zero in expfct
		expfct = 1e30;  //to avoid division by zero
	
	return ccc*(cos(latitude_static)*cos(latitude_static)*cos(longitude-obj.phaseangle)/(expfct-1));
}

double aux_integrand_neatm ( double latitude, double eps, ASTEROID obj, double T0 )
{
	latitude_static = latitude;
	
	//printf(" aux_integrand_neatm mit lat=%.3f -> romberg(%.3f, %.3f)\n",
	//       latitude_static, longitude1(phase), longitude2(phase));
	double rombi = romberg(integrand_neatm, longitude_lowerlimit, longitude_upperlimit, eps, obj, T0);
	if (isnan(rombi))
	{
          printf("%f %f %f %f\n", obj.diameter, pv_from_D(obj.diameter, obj.absmag), T0, bondalbedo(obj.slopepar, pv_from_D(obj.diameter, obj.absmag)));
	}

}


//complete flux calculation in W m^-2 um^-1 using 2d-Romberg integration
//no use for 'parameter'
double planck_neatm(ASTEROID obj, double epsilon)
{
	//define upper and lower integration limit for longitude
	if ( obj.phaseangle >= 0 )
	{
		longitude_lowerlimit = obj.phaseangle - PI2;
		longitude_upperlimit = PI2;
	}
	else
	{
		longitude_lowerlimit = -PI2;
		longitude_upperlimit = obj.phaseangle + PI2;
	}
	
	double pv = pv_from_D(obj.diameter, obj.absmag);
	double T0 = pow((1-bondalbedo(obj.slopepar, pv))*SOLARCONST/(obj.heliodist*obj.heliodist
			*obj.eta*EMISSIVITY*STEFANBOLTZMANN), 0.25);

	//display mean surface temperature
	//printf("mean T = %f\n", T0*pow(2,-0.25));


	double integration = romberg(aux_integrand_neatm, 0, PI2, epsilon, obj, T0);
	if (isnan(integration))
	{
	   double a = 0;
	   double b = PI2;
	   while (isnan(integration))
           {	   
	      a += 1e-5;
	      b -= 1e-5;
	      integration = romberg(aux_integrand_neatm, 0, PI2, epsilon, obj, T0);
	   }	 
	}

	return (EMISSIVITY*obj.diameter*obj.diameter*HC2
			 /(obj.geodist*AU*obj.geodist*AU*obj.wavelength*obj.wavelength*
			  obj.wavelength*obj.wavelength*obj.wavelength) //um conversion to be found in HC2
			 *integration); 

}


//complete flux calculation in W m^-2 um^-1 using 2d-Romberg integration
//no use for 'parameter'
double planck_neatm_T0(ASTEROID obj, double epsilon, double T0)
{
	//define upper and lower integration limit for longitude
	if ( obj.phaseangle >= 0 )
	{
		longitude_lowerlimit = obj.phaseangle - PI2;
		longitude_upperlimit = PI2;
	}
	else
	{
		longitude_lowerlimit = -PI2;
		longitude_upperlimit = obj.phaseangle + PI2;
	}
	
	//display mean surface temperature
	//printf("mean T = %f\n", T0*pow(2,-0.25));

	return (EMISSIVITY*obj.diameter*obj.diameter*HC2
			 /(obj.geodist*AU*obj.geodist*AU*obj.wavelength*obj.wavelength*
			  obj.wavelength*obj.wavelength*obj.wavelength) //um conversion to be found in HC2
			 *romberg(aux_integrand_neatm, 0, PI2, epsilon, obj, T0)); 

}



/*
//Planck function integrand
//if function only one-dimensional, use only latitude
//no use for 'parameter'
double integrand_neatm (double longitude, double epsilon, ASTEROID obj, double T0)
{
	double expfct, longsum;
	
	//integrate over longitude using trapezoidal rule
	double longstep = (longitude_upperlimit-longitude_lowerlimit)/100.;
	for (longitude=longitude_lowerlimit, longsum=0; longitude<=longitude_upperlimit; longitude+=longstep)
	{
		expfct = exp(HCK/(obj.wavelength*T0*pow(cos(latitude_static)*cos(longitude),0.25)));
		
		if ( isnan(expfct) )    //NaN results from division by zero in longsum
			expfct = 1e30;  //to avoid division by zero
		
		if ( (longitude==longitude_lowerlimit) || (fabs(longitude-longitude_upperlimit)<longstep) )
			longsum += 0.5*cos(latitude_static)*cos(latitude_static)*cos(longitude-obj.phaseangle)/(expfct-1);
		else
			longsum += cos(latitude_static)*cos(latitude_static)*cos(longitude-obj.phaseangle)/(expfct-1);
	}
	
	return longsum*longstep;
}
*/
