/*
 *      ERRORS - error determination for MODEL
 *      Version 1.1 - Oct 12, 2010
 *      
 *      Michael Mommert, DLR
 *      michael.mommert@dlr.de
 *
 *   V 1.0 - determine mean and sigma from Monte Carlo simulation for one or a number of objects
 *   V 1.1 - determine median and 68.2% limits from Monte Carlo simulation, number of runs adjusted so that really n simulation are available
 *   V 1.2 - implemented to determine 68.2% around the best fit model - not in use!
 */


void error_estimation (ASTEROID *obj, int number_of_objects, int number_of_labrats, double error_sigma, int fitting, double epsilon, int corrflux, int blowfluxerrors, int plotting)
{
	ASTEROID labrat[EXP_NUMBER_OF_OBJECTS]; //array size is maximum of single flux measurements (the same as for neo)
	//counter_labrat: current number of labrats for this object (this index_obj)
	int index_obj, index_labrat, counter_labrat, number_of_measurements;
	double reflectanceratio_mean, absmag_mean, slopepar_mean, flux_raw_mean;
	double reflectanceratio_sigma, absmag_sigma, slopepar_sigma, flux_raw_sigma;
	double diameter_mean, pv_mean, eta_mean, modelparameter_mean;
	double diameter_sigma, pv_sigma, eta_sigma, modelparameter_sigma;
	double k = 1;
	
	//arrays for sorting 
	double diameter[EXP_NUMBER_OF_OBJECTS];
	double pv[EXP_NUMBER_OF_OBJECTS];
	double eta[EXP_NUMBER_OF_OBJECTS];
	double modelparameter[EXP_NUMBER_OF_OBJECTS];
	
	//create result-file name
	char filename[100]; 
	//sprintf(filename, "error_single.results");
	sprintf(filename, "no_output");

	// disable eta-clipping during Monte Carlo method
	eta_max = 100;
	eta_min = 0;

	printf("perform Monte-Carlo error estimation:\n");
	
	//model all objects once in order to derive the best fit model and the multiplication factor k = S_min = sqrt(Chi2_red) later
	runmodel(obj, 0, number_of_objects, fitting, epsilon, 0, filename, 0, "_MC_preparation");


	for (index_obj=0; index_obj<number_of_objects; index_obj++)
	{
		printf(" object: %s\n", obj[index_obj].name);

		//rescale flux errors in order to achieve a good error estimation
		if (blowfluxerrors == 1)
		{
			//determine multiplication factor k = S_min = sqrt(Chi2_red)
			k = sqrt(obj[index_obj].chi2);
			if (k < 1) k = 1;
			printf("   increase flux errors by a factor of %.2f, i.e. \n", k);
				
			//display new flux uncertainties for this object
			int all_done = 0;
			int i = index_obj;
			while ( all_done == 0) 
			{
				printf("    %.3f micron flux: (%.3f +- %.3f) mJy  (instead +- %.3f)\n", obj[i].wavelength, 
						obj[i].flux/obj[i].colorcorrection*3.335640952e14*obj[i].wavelength*obj[i].wavelength,
						k*obj[i].fluxerr*3.335640952e14*obj[i].wavelength*obj[i].wavelength,obj[i].fluxerr*3.335640952e14*obj[i].wavelength*obj[i].wavelength);			

				if ( obj[i].next_measurement != 0)
					i = obj[i].next_measurement;
				else 
					all_done = 1;
			}
		}

		int rubbish_labrats = 0; //number of measurements that don't pass the quality check
		
		//reset each object
		reflectanceratio_mean = 0; absmag_mean = 0; slopepar_mean = 0; flux_raw_mean = 0;
		reflectanceratio_sigma = 0; absmag_sigma = 0; slopepar_sigma = 0; flux_raw_sigma = 0;
		diameter_mean = 0; pv_mean = 0; eta_mean = 0; modelparameter_mean = 0;
		diameter_sigma = 0; pv_sigma = 0; eta_sigma = 0; modelparameter_sigma = 0;
		
		/* source asteroid data
		printf( "source %17.17s  %.3f %.4f %.2f %.3f<%.3f<%.3f %.3f %.3f\n",
		obj[index_obj].name, obj[index_obj].heliodist, obj[index_obj].geodist, obj[index_obj].phaseangle/PI2*90,
		(obj[index_obj].flux_raw-obj[index_obj].fluxerr)*3.335640952e14*obj[index_obj].wavelength*obj[index_obj].wavelength, 
		obj[index_obj].flux_raw*3.335640952e14*obj[index_obj].wavelength*obj[index_obj].wavelength, 
		(obj[index_obj].flux_raw+obj[index_obj].fluxerr)*3.335640952e14*obj[index_obj].wavelength*obj[index_obj].wavelength, 
		obj[index_obj].flux/obj[index_obj].colorcorrection*3.335640952e14*obj[index_obj].wavelength*obj[index_obj].wavelength, obj[index_obj].reflectanceratio );
		*/
		
		//create and model synthetic asteroids one at a time, number_of_labrats+rubbish_labrats times
		for (counter_labrat=0; counter_labrat<number_of_labrats+rubbish_labrats; counter_labrat++)
		{
			//reset index_labrat to momentary obj index
			index_labrat = index_obj;
			
			if ( fitting == 1 )
			{
				int all_done = 0;
				number_of_measurements = 0;
				
				//randomize each measurement for this object
				while ( all_done == 0 ) 
				{
					number_of_measurements++;
					
					//adjust name
					char tmpstring[100];
					sprintf(tmpstring, "%s_%d", obj[index_obj].name, counter_labrat);
					labrat[index_labrat].name              = tmpstring;
					//no changes here
					labrat[index_labrat].wavelength        = obj[index_labrat].wavelength;
					labrat[index_labrat].eta               = obj[index_labrat].eta;
					labrat[index_labrat].modelparameter    = obj[index_labrat].modelparameter;
					labrat[index_labrat].geodist           = obj[index_labrat].geodist;
					labrat[index_labrat].heliodist         = obj[index_labrat].heliodist;
					labrat[index_labrat].reflectanceratio  = obj[index_labrat].reflectanceratio;
					labrat[index_labrat].reflectanceratio_error = obj[index_labrat].reflectanceratio_error;
					labrat[index_labrat].absmag_error      = obj[index_labrat].absmag_error;
					labrat[index_labrat].slopepar_error    = obj[index_labrat].slopepar_error;
					labrat[index_labrat].eta_sigma         = obj[index_labrat].eta_sigma;	
					labrat[index_labrat].phaseangle        = obj[index_labrat].phaseangle;
					labrat[index_labrat].colorcorrection   = obj[index_labrat].colorcorrection;
					labrat[index_labrat].fluxerr           = k*obj[index_labrat].fluxerr;
					labrat[index_labrat].next_measurement  = obj[index_labrat].next_measurement;
					
					//------------------------------------------------------------------------------------------
					// APPLY NO ERRORS AT ALL
					labrat[index_labrat].reflectanceratio  = obj[index_labrat].reflectanceratio;
					labrat[index_labrat].absmag            = obj[index_labrat].absmag;
					labrat[index_labrat].slopepar          = obj[index_labrat].slopepar;
					labrat[index_labrat].flux_raw          = obj[index_labrat].flux_raw;
					labrat[index_labrat].flux              = labrat[index_labrat].flux_raw;
					labrat[index_labrat].eta = obj[index_labrat].eta;
					//------------------------------------------------------------------------------------------
					
					//------------------------------------------------------------------------------------------
					// APPLY ERRORS (uncomment to apply)
					
					//randomize the following parameters
					
					
					//----------------------- REFLECTANCE RATIO IR/V
					labrat[index_labrat].reflectanceratio  = 
						gaussrand(obj[index_labrat].reflectanceratio,
						          error_sigma*obj[index_labrat].reflectanceratio_error);
					
					
					
					
					//----------------------- ABSOLUTE MAGNITUDE H
					labrat[index_labrat].absmag            = 
					        gaussrand(obj[index_labrat].absmag, error_sigma*obj[index_labrat].absmag_error);
					
					
					
					
					//----------------------- SLOPE PARAMETER G
					//use gaussian distribution in a given range (here 0.05 < G < 0.5
					/* labrat[index_labrat].slopepar          =  */
					/* 		gaussrand_restricted(obj[index_labrat].slopepar, error_sigma*obj[index_labrat].slopepar_error, 0.05, 0.5); */
					/* //use complete gaussian distribution around G */
					/* labrat[index_labrat].slopepar          =  */
					/* 	gaussrand(obj[index_labrat].slopepar, error_sigma*obj[index_labrat].slopepar_error); */
					
					
					
					//------------------------ FLUX
					//randomize flux_raw and then later recalculate final flux (or not)
					labrat[index_labrat].flux_raw          = 
					  gaussrand(obj[index_labrat].flux_raw, k*error_sigma*obj[index_labrat].fluxerr);
					while (labrat[index_labrat].flux_raw <= 0.0)
					  labrat[index_labrat].flux_raw          = 
					    gaussrand(obj[index_labrat].flux_raw, k*error_sigma*obj[index_labrat].fluxerr);
					labrat[index_labrat].flux              = labrat[index_labrat].flux_raw;

					//------------------------- ETA
					//eta errors for fixed eta models
					labrat[index_labrat].eta = obj[index_labrat].eta + exp(log(0.655)+0.628*gaussrand(0, 1))-0.675;
					while ((labrat[index_labrat].eta > 3.14) || (labrat[index_labrat].eta < 0.6))
					  labrat[index_labrat].eta = obj[index_labrat].eta + exp(log(0.655)+0.628*gaussrand(0, 1))-0.675;  

					
					//------------------------------------------------------------------------------------------
					
					labrat[index_labrat].fitparflag           = 0;
					

					printf("flux: %f", labrat[index_labrat].flux);
					
					//printf ("obj[%d]->labrat[%d] = labrat no %d\n", index_obj, index_labrat, counter_labrat);
					
					//breakout if all data sets concerning this
					//object are processed
					if ( labrat[index_labrat].next_measurement != 0)
						index_labrat = labrat[index_labrat].next_measurement;
					else 
						all_done = 1;
				}
			}
			else
			{
				number_of_measurements = 1;
				
				//adjust name
				char tmpstring[100];
				sprintf(tmpstring, "%s_%d", obj[index_obj].name, counter_labrat);
				labrat[index_labrat].name              = tmpstring;
				//no changes here
				labrat[index_labrat].wavelength        = obj[index_obj].wavelength;
				labrat[index_labrat].eta               = obj[index_obj].eta;
				labrat[index_labrat].modelparameter    = obj[index_obj].modelparameter;
				labrat[index_labrat].geodist           = obj[index_obj].geodist;
				labrat[index_labrat].heliodist         = obj[index_obj].heliodist;
				labrat[index_labrat].reflectanceratio  = obj[index_obj].reflectanceratio;
				labrat[index_labrat].absmag_error      = obj[index_obj].absmag_error;
				labrat[index_labrat].slopepar_error    = obj[index_obj].slopepar_error;
				labrat[index_labrat].reflectanceratio_error = obj[index_obj].reflectanceratio_error;
				labrat[index_labrat].phaseangle        = obj[index_obj].phaseangle;
				labrat[index_labrat].colorcorrection   = obj[index_obj].colorcorrection;
				labrat[index_labrat].fluxerr           = obj[index_obj].fluxerr;
				
				
				//----------------------- REFLECTANCE RATIO IR/V
				labrat[index_labrat].reflectanceratio  = 
				  gaussrand(obj[index_labrat].reflectanceratio,
					    error_sigma*obj[index_labrat].reflectanceratio_error);
					
					
					
				
				//----------------------- ABSOLUTE MAGNITUDE H
				labrat[index_labrat].absmag            = 
				  gaussrand(obj[index_labrat].absmag, error_sigma*obj[index_labrat].absmag_error);
				
					
					
				
				//----------------------- SLOPE PARAMETER G
				//use gaussian distribution in a given range (here 0.05 < G < 0.5
				/* labrat[index_labrat].slopepar          =  */
				/* 		gaussrand_restricted(obj[index_labrat].slopepar, error_sigma*obj[index_labrat].slopepar_error, 0.05, 0.5); */
				/* //use complete gaussian distribution around G */
				/* labrat[index_labrat].slopepar          =  */
				/* 	gaussrand(obj[index_labrat].slopepar, error_sigma*obj[index_labrat].slopepar_error); */
				
					
					
	                        //------------------------ FLUX
                                //randomize flux_raw and then later recalculate final flux (or not)
                                labrat[index_labrat].flux_raw          =
                                    gaussrand(obj[index_labrat].flux_raw, k*error_sigma*obj[index_labrat].fluxerr);
                                while (labrat[index_labrat].flux_raw <= 0.0)
                                     labrat[index_labrat].flux_raw          =
                                        gaussrand(obj[index_labrat].flux_raw, k*error_sigma*obj[index_labrat].fluxerr);
                                labrat[index_labrat].flux              = labrat[index_labrat].flux_raw;

                                //------------------------- ETA
                                //eta errors for fixed eta models
                                labrat[index_labrat].eta = obj[index_labrat].eta + exp(log(0.655)+0.628*gaussrand(0, 1))-0.675;
                                while ((labrat[index_labrat].eta > 3.14) || (labrat[index_labrat].eta < 0.6))
                                     labrat[index_labrat].eta = obj[index_labrat].eta + exp(log(0.655)+0.628*gaussrand(0, 1))-0.675;	



				/* //randomize the following parameters */
				/* labrat[index_labrat].reflectanceratio  =  */
				/* 	gaussrand(obj[index_obj].reflectanceratio, */
				/* 	          error_sigma*obj[index_obj].reflectanceratio_error); */
				/* labrat[index_labrat].absmag            =  */
				/*         gaussrand(obj[index_obj].absmag, error_sigma*obj[index_obj].absmag_error); */
				/* labrat[index_labrat].slopepar          =  */
				/* 	gaussrand(obj[index_obj].slopepar, error_sigma*obj[index_obj].slopepar_error); */
				
				/* //randomize flux_raw */
				/* labrat[index_labrat].flux_raw          =  */
				/* 	gaussrand(obj[index_obj].flux_raw, error_sigma*obj[index_obj].fluxerr); */
				/* labrat[index_labrat].flux              = labrat[index_labrat].flux_raw; */
				
				/* labrat[index_labrat].fitparflag           = 0; */
				
				/* //eta errors for fixed eta models */
				/* if ( strstr(model, "4d") != NULL ) //use Delbo 2004 method */
				/* 	labrat[index_labrat].eta = gaussrand(0.013, 0.004)*labrat[index_labrat].phaseangle/PI2*90 */
				/* 				   + gaussrand(0.89, 0.17);  */
				/* else if ( strstr(model, "4w") != NULL )  //use Wolters 2007 method */
				/*   labrat[index_labrat].eta = gaussrand(0.013, 0.004)*labrat[index_labrat].phaseangle/PI2*90 */
				/* 				   + gaussrand(0.91, 0.17);  */
				/* else if ( (strstr(model, "4f") != NULL) ) */
				/*   labrat[index_labrat].eta            =  */
				/*     gaussrand(obj[index_labrat].eta, 0.001);	 */
				
				//printf ("obj[%d]->labrat[%d] = labrat no %d\n", index_obj, index_labrat, counter_labrat);
			}
	
			//printf("    model (%d-%d)!\n", index_obj, index_obj+number_of_measurements);
			
			//apply flux correction
			if ( corrflux == 1 )
				fluxcorrection ( labrat, index_obj, index_obj+number_of_measurements );

			// check if corrected fluxes are not too small (code breaks for fluxes < 0.001mJy)
			// also check for bad eta values
			int i, negflux = 0; // =1 if one flux value of this object is negative
			int fitparflag = 0; // =1 if one eta value of this object is out of range
			for (i=index_labrat; i<index_labrat+number_of_measurements; i++)
			{
				if ( labrat[i].flux < 0.001/(3.335640952e14 * labrat[index_labrat].wavelength*labrat[index_labrat].wavelength) )
				{
				  printf("  WARNING: low flux (%fmJy)! Replace with lower limit\n", labrat[i].flux*(3.335640952e14 * labrat[index_labrat].wavelength*labrat[index_labrat].wavelength));
				  labrat[i].flux = 0.0001/(3.335640952e14 * labrat[index_labrat].wavelength*labrat[index_labrat].wavelength);
				  negflux = 1;
				}
				if ( labrat[i].fitparflag == 1 )
					fitparflag = 1;
			}

			//if (negflux != 0 || fitparflag != 0) counter_labrat--;
			if (negflux != 0 || fitparflag != 0) rubbish_labrats++;

			//display progress status
			//printf(" (%d/%d)  ", counter_labrat-rubbish_labrats+1, number_of_labrats);
			fflush(stdout);

			//--------------------------- model momentary labrat asteroid(s)
			runmodel(labrat, index_obj, number_of_measurements, fitting, epsilon, 0, filename, 0, "_error_estimation");
			//start modeling here with object index_obj, since index_labrat points at the last created measurement of this object and index_obj points at the first measurement

			if ( negflux == 0 && fitparflag == 0 )
			{
				//------------calculate mean values and sigmas predecessors of randomized parameters and modelling results------------------------
				reflectanceratio_mean += labrat[index_labrat].reflectanceratio;
				reflectanceratio_sigma += labrat[index_labrat].reflectanceratio*labrat[index_labrat].reflectanceratio;
				absmag_mean += labrat[index_labrat].absmag;
				absmag_sigma += labrat[index_labrat].absmag*labrat[index_labrat].absmag;
				slopepar_mean += labrat[index_labrat].slopepar;
				slopepar_sigma += labrat[index_labrat].slopepar*labrat[index_labrat].slopepar;
				flux_raw_mean += labrat[index_labrat].flux_raw;
				flux_raw_sigma += labrat[index_labrat].flux_raw*labrat[index_labrat].flux_raw*1e40;  //1e40 to avoid truncation errors
				
				diameter_mean += labrat[index_labrat].diameter;
				diameter_sigma += labrat[index_labrat].diameter*labrat[index_labrat].diameter;
				pv_mean += labrat[index_labrat].pv;
				pv_sigma += labrat[index_labrat].pv*labrat[index_labrat].pv;
				eta_mean += labrat[index_labrat].eta;
				eta_sigma += labrat[index_labrat].eta*labrat[index_labrat].eta;
				modelparameter_mean += labrat[index_labrat].modelparameter;
				modelparameter_sigma += labrat[index_labrat].modelparameter*labrat[index_labrat].modelparameter;
				
				diameter[counter_labrat-rubbish_labrats] = labrat[index_labrat].diameter;
				pv[counter_labrat-rubbish_labrats] = labrat[index_labrat].pv;
				eta[counter_labrat-rubbish_labrats] = labrat[index_labrat].eta;
				modelparameter[counter_labrat-rubbish_labrats] = labrat[index_labrat].modelparameter;
				
				printf("  done!\n");
			}
			
			/*
			printf( "%d %24.24s  %.3f %.4f %.3f %.3f %.3f %.3f %.3f \n", index_labrat,
			labrat[index_labrat].name, labrat[index_labrat].heliodist, labrat[index_labrat].geodist, labrat[index_labrat].phaseangle/PI2*90,
			labrat[index_labrat].flux_raw*3.335640952e14*labrat[index_labrat].wavelength*labrat[index_labrat].wavelength,
			labrat[index_labrat].flux/labrat[index_labrat].colorcorrection*3.335640952e14*labrat[index_labrat].wavelength*labrat[index_labrat].wavelength, 
			labrat[index_labrat].reflectanceratio, labrat[index_labrat].flux/labrat[index_labrat].fluxerr );
			*/ 
			
			//plot labrat
			if (plotting == 1)
				plotmodel(labrat, index_obj, epsilon, planck_neatm);
			
		}
		
		
		//------------------ calculate statistics results -----------------------
		
		//printf("number of clean labrats: %d\n", counter_labrat-rubbish_labrats);


		//randomized input parameter
		reflectanceratio_mean /= number_of_labrats;
		absmag_mean /= number_of_labrats;
		slopepar_mean /= number_of_labrats;
		flux_raw_mean /= number_of_labrats;
		
		diameter_mean /= number_of_labrats;
		pv_mean /= number_of_labrats;
		eta_mean /= number_of_labrats;
		modelparameter_mean /= number_of_labrats;
		
		reflectanceratio_sigma = sqrt( (reflectanceratio_sigma/number_of_labrats - reflectanceratio_mean*reflectanceratio_mean) );
		absmag_sigma = sqrt( (absmag_sigma/number_of_labrats - absmag_mean*absmag_mean) );
		slopepar_sigma = sqrt( (slopepar_sigma/number_of_labrats - slopepar_mean*slopepar_mean) );
		flux_raw_sigma = sqrt( (flux_raw_sigma/number_of_labrats*1e-40 - flux_raw_mean*flux_raw_mean) );
				
		diameter_sigma = sqrt( (diameter_sigma/number_of_labrats - diameter_mean*diameter_mean) );
		pv_sigma = sqrt( (pv_sigma/number_of_labrats - pv_mean*pv_mean) );
		eta_sigma = sqrt( (eta_sigma/number_of_labrats - eta_mean*eta_mean) );
		modelparameter_sigma = sqrt( (modelparameter_sigma/number_of_labrats - modelparameter_mean*modelparameter_mean) );
	
		//sort arrays
		for (counter_labrat=0; counter_labrat<number_of_labrats; counter_labrat++)
			for (index_labrat=0; index_labrat<number_of_labrats; index_labrat++)
			{
				if (diameter[index_labrat] > diameter[counter_labrat])
				{
					double swap = diameter[index_labrat];
					diameter[index_labrat] = diameter[counter_labrat];
					diameter[counter_labrat] = swap;
				}
				
				if (pv[index_labrat] > pv[counter_labrat])
				{
					double swap = pv[index_labrat];
					pv[index_labrat] = pv[counter_labrat];
					pv[counter_labrat] = swap;
				}
				
				if (eta[index_labrat] > eta[counter_labrat])
				{
					double swap = eta[index_labrat];
					eta[index_labrat] = eta[counter_labrat];
					eta[counter_labrat] = swap;
				}
				
				if (modelparameter[index_labrat] > modelparameter[counter_labrat])
				{
					double swap = modelparameter[index_labrat];
					modelparameter[index_labrat] = modelparameter[counter_labrat];
					modelparameter[counter_labrat] = swap;
				}
			}
		
		
		//determine median, upper and lower 68.2% and 99.7% limits
		double diameter_median, pv_median, eta_median, modelparameter_median;
		double diameter_lower, pv_lower, eta_lower, modelparameter_lower;
		double diameter_upper, pv_upper, eta_upper, modelparameter_upper;
		double diameter_lower3sig, pv_lower3sig, eta_lower3sig, modelparameter_lower3sig;
		double diameter_upper3sig, pv_upper3sig, eta_upper3sig, modelparameter_upper3sig;
		int lower, upper, lower3sig, upper3sig; 
		if ( ((int)(number_of_labrats/2) - (double)number_of_labrats/2) == 0 )
		{
			//even number of labrats [labrats 0...(number_of_labrats-1); average labrats around (number_of_labrats-1)/2]
			diameter_median = 0.5*(diameter[(number_of_labrats-2)/2]+diameter[(number_of_labrats)/2]);
			pv_median = 0.5*(pv[(number_of_labrats-2)/2]+pv[(number_of_labrats)/2]);
			eta_median = 0.5*(eta[(number_of_labrats-2)/2]+eta[(number_of_labrats)/2]);
			modelparameter_median = 0.5*(modelparameter[(number_of_labrats-2)/2]+modelparameter[(number_of_labrats)/2]);
			lower = (int)((number_of_labrats-2)/2 - number_of_labrats*0.341);
			upper = (int)((number_of_labrats-2)/2 + number_of_labrats*0.341);
			lower3sig = (int)((number_of_labrats-2)/2 - number_of_labrats*0.4985);
			upper3sig = (int)((number_of_labrats-2)/2 + number_of_labrats*0.4985);
		}
		else
		{
			//uneven number of labrats [labrats 0...(number_of_labrats-1);use labrat (number_of_labrats-1)/2]
			diameter_median = diameter[(number_of_labrats-1)/2];
			pv_median = pv[(number_of_labrats-1)/2];
			eta_median = eta[(number_of_labrats-1)/2];
			modelparameter_median = modelparameter[(number_of_labrats-1)/2];
			lower = (int)((number_of_labrats-1)/2 - number_of_labrats*0.341);
			upper = (int)((number_of_labrats-1)/2 + number_of_labrats*0.341);
			lower3sig = (int)((number_of_labrats-1)/2 - number_of_labrats*0.4985);
			upper3sig = (int)((number_of_labrats-1)/2 + number_of_labrats*0.4985);
		}

		//printf("%d %d\n", lower3sig, upper3sig);

		diameter_lower       = diameter[lower] - diameter_median;
		diameter_lower3sig   = diameter[lower3sig] - diameter_median;		
		pv_lower             = pv[lower] - pv_median;
		pv_lower3sig         = pv[lower3sig] - pv_median;
		eta_lower            = eta[lower] - eta_median;
		eta_lower3sig        = eta[lower3sig] - eta_median;
		modelparameter_lower = modelparameter[lower] - modelparameter_median;
		modelparameter_lower3sig = modelparameter[lower3sig] - modelparameter_median;
		diameter_upper       = diameter[upper] - diameter_median;
		diameter_upper3sig   = diameter[upper3sig] - diameter_median;
		pv_upper             = pv[upper] - pv_median;
		pv_upper3sig         = pv[upper3sig] - pv_median;
		eta_upper            = eta[upper] - eta_median;
		eta_upper3sig        = eta[upper3sig] - eta_median;
		modelparameter_upper = modelparameter[upper] - modelparameter_median;


		//printf ("diam: median %.5f, lower %.5f, upper %.5f\n", diameter_median, diameter_lower, diameter_upper);
		
		//for (counter_labrat=0; counter_labrat<number_of_labrats; counter_labrat++)
		//	printf("%d %.5f %.5f %.5f\n", counter_labrat, diameter[counter_labrat], pv[counter_labrat], eta[counter_labrat]);
		
		/*
		printf("source    mean   sigma\n %.3f  %.3f+-%.3f  reflectanceratio\n %.3f  %.3f+-%.3f  absmag\n %.3f  %.3f+-%.3f  slopepar\n %.3f  %.3f+-%.3f  flux_raw\n", 
		       obj[index_obj].reflectanceratio, reflectanceratio_mean, reflectanceratio_sigma,
		       obj[index_obj].absmag, absmag_mean, absmag_sigma,
		       obj[index_obj].slopepar, slopepar_mean, slopepar_sigma,
		       obj[index_obj].flux_raw*3.335640952e14*obj[index_obj].wavelength*obj[index_obj].wavelength, 
		       flux_raw_mean*3.335640952e14*obj[index_obj].wavelength*obj[index_obj].wavelength,
		       flux_raw_sigma*3.335640952e14*obj[index_obj].wavelength*obj[index_obj].wavelength);
		*/

		/*
		//determine 68% around the best model fit parameters
		//find index of the closest datapoints to the best fit model values of d, pv and eta
		double diameter_residual, pv_residual, eta_residual;
		diameter_residual = pv_residual = eta_residual = 10000 ;
		int diameter_bestfit_index, pv_bestfit_index, eta_bestfit_index;
		int i;
		printf("best fit values: %f  %f   %f\n", obj[index_obj].diameter, obj[index_obj].pv, obj[index_obj].eta);
		for (i=0; i<number_of_labrats; i++)
		  {
		    if (fabs(diameter[i]-obj[index_obj].diameter) < diameter_residual)
		      {
			diameter_residual = fabs(diameter[i]-obj[index_obj].diameter);
			diameter_bestfit_index = i;
		      }
		    if (fabs(pv[i]-obj[index_obj].pv) < pv_residual)
		      {
			pv_residual = fabs(pv[i]-obj[index_obj].pv);
			pv_bestfit_index = i;
		      }
		    if (fabs(eta[i]-obj[index_obj].eta) < eta_residual)
		      {
			eta_residual = fabs(eta[i]-obj[index_obj].eta);
			eta_bestfit_index = i;
		      }
		  }
		
	     	printf("best fit matches: %f  %f   %f\n", diameter[diameter_bestfit_index], pv[pv_bestfit_index], eta[eta_bestfit_index]);

		//check whether limits are inside range
		int diameter_low_index = diameter_bestfit_index-(int)(number_of_labrats*0.341);
		int diameter_up_index = diameter_bestfit_index+(int)(number_of_labrats*0.341);
		int pv_low_index = pv_bestfit_index-(int)(number_of_labrats*0.341);
		int pv_up_index = pv_bestfit_index+(int)(number_of_labrats*0.341);
		int eta_low_index = eta_bestfit_index-(int)(number_of_labrats*0.341);
		int eta_up_index = eta_bestfit_index+(int)(number_of_labrats*0.341);
		if (diameter_low_index < 0)
		  {
		    diameter_low_index = 0;
		    printf ("WARNING: lower diameter limit not sampled; apply lowest value as uncertainty\n");
		  }
		if (diameter_up_index >= number_of_labrats)
		  {
		    diameter_up_index = number_of_labrats-1;
		    printf ("WARNING: upper diameter limit not sampled; apply highest value as uncertainty\n");
		  }

		if (pv_low_index < 0)
		  {
		    pv_low_index = 0;
		    printf ("WARNING: lower pv limit not sampled; apply lowest value as uncertainty\n");
		  }
		if (pv_up_index >= number_of_labrats)
		  {
		    pv_up_index = number_of_labrats-1;
		    printf ("WARNING: upper pv limit not sampled; apply highest value as uncertainty\n");
		  }

		if (eta_low_index < 0)
		  {
		    eta_low_index = 0;
		    printf ("WARNING: lower eta limit not sampled; apply lowest value as uncertainty\n");
		  }
		if (eta_up_index >= number_of_labrats)
		  {
		    eta_up_index = number_of_labrats-1;
		    printf ("WARNING: upper eta limit not sampled; apply highest value as uncertainty\n");
		  }
		*/
		

		if ( fitting == 1 )
		{
			int all_done = 0;
			int i = index_obj;
			while ( all_done == 0) 
			{
			  obj[i].diameter               = diameter_median;             
			  obj[i].diameter_mean          = diameter_mean;
			  obj[i].diameter_sigma         = diameter_sigma;
			  obj[i].diameter_median        = diameter_median;
			  obj[i].diameter_lower         = diameter_lower;
			  obj[i].diameter_upper         = diameter_upper;
			  obj[i].diameter_lower3sig     = diameter_lower3sig;
			  obj[i].diameter_upper3sig     = diameter_upper3sig;
			  obj[i].pv                     = pv_median;                               
			  obj[i].pv_mean                = pv_mean;
			  obj[i].pv_sigma               = pv_sigma;
			  obj[i].pv_median              = pv_median;
			  obj[i].pv_lower               = pv_lower;
			  obj[i].pv_upper               = pv_upper;
			  obj[i].pv_lower3sig           = pv_lower3sig;
			  obj[i].pv_upper3sig           = pv_upper3sig;
			  obj[i].eta                    = eta_median; 
			  obj[i].eta_mean               = eta_mean;
			  obj[i].eta_sigma              = eta_sigma;
			  obj[i].eta_median             = eta_median;
			  obj[i].eta_lower              = eta_lower;
			  obj[i].eta_upper              = eta_upper;
			  obj[i].eta_lower3sig          = eta_lower3sig;
			  obj[i].eta_upper3sig          = eta_upper3sig;
			  //obj[i].modelparameter         = modelparameter_median;
			  //obj[i].modelparameter_mean    = modelparameter_mean;
			  //obj[i].modelparameter_sigma   = modelparameter_sigma;
			  //obj[i].modelparameter_median  = modelparameter_median;
			  //obj[i].modelparameter_lower   = modelparameter_lower;
			  //obj[i].modelparameter_upper   = modelparameter_upper;
			  obj[i].reflectanceratio       = reflectanceratio_mean;
			  obj[i].reflectanceratio_error = reflectanceratio_sigma;
			  obj[i].absmag                 = absmag_mean;
			  obj[i].absmag_error           = absmag_sigma;
			  obj[i].slopepar               = slopepar_mean;
			  obj[i].slopepar_error         = slopepar_sigma;
			  obj[i].flux_raw               = flux_raw_mean;
			  obj[i].fluxerr                = flux_raw_sigma;
			  obj[i].chi2                   = 0;
				
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
		  obj[index_obj].diameter               = diameter_median;             
		  obj[index_obj].diameter_mean          = diameter_mean;
		  obj[index_obj].diameter_sigma         = diameter_sigma;
		  obj[index_obj].diameter_median        = diameter_median;
		  obj[index_obj].diameter_lower         = diameter_lower;
		  obj[index_obj].diameter_upper         = diameter_upper;
		  obj[index_obj].diameter_lower3sig     = diameter_lower3sig;
		  obj[index_obj].diameter_upper3sig     = diameter_upper3sig;
		  obj[index_obj].pv                     = pv_median;                               
		  obj[index_obj].pv_mean                = pv_mean;
		  obj[index_obj].pv_sigma               = pv_sigma;
		  obj[index_obj].pv_median              = pv_median;
		  obj[index_obj].pv_lower               = pv_lower;
		  obj[index_obj].pv_upper               = pv_upper;
		  obj[index_obj].pv_lower3sig           = pv_lower3sig;
		  obj[index_obj].pv_upper3sig           = pv_upper3sig;
		  obj[index_obj].eta                    = eta_median; 
		  obj[index_obj].eta_mean               = eta_mean;
		  obj[index_obj].eta_sigma              = eta_sigma;
		  obj[index_obj].eta_median             = eta_median;
		  obj[index_obj].eta_lower              = eta_lower;
		  obj[index_obj].eta_upper              = eta_upper;
		  obj[index_obj].eta_lower3sig          = eta_lower3sig;
		  obj[index_obj].eta_upper3sig          = eta_upper3sig;
		  //obj[index_obj].modelparameter         = modelparameter_median;
		  //obj[index_obj].modelparameter_mean    = modelparameter_mean;
		  //obj[index_obj].modelparameter_sigma   = modelparameter_sigma;
		  //obj[index_obj].modelparameter_median  = modelparameter_median;
		  //obj[index_obj].modelparameter_lower   = modelparameter_lower;
		  //obj[index_obj].modelparameter_upper   = modelparameter_upper;
		  obj[index_obj].reflectanceratio       = reflectanceratio_mean;
		  obj[index_obj].reflectanceratio_error = reflectanceratio_sigma;
		  obj[index_obj].absmag                 = absmag_mean;
		  obj[index_obj].absmag_error           = absmag_sigma;
		  obj[index_obj].slopepar               = slopepar_mean;
		  obj[index_obj].slopepar_error         = slopepar_sigma;
		  obj[index_obj].flux_raw               = flux_raw_mean;
		  obj[index_obj].fluxerr                = flux_raw_sigma;
		  obj[index_obj].chi2                   = 0;
		  
		  
			/* obj[index_obj].diameter       = diameter_mean; */
			/* obj[index_obj].diameter_sigma = diameter_sigma; */
			/* obj[index_obj].pv             = pv_mean; */
			/* obj[index_obj].pv_sigma       = pv_sigma; */
			/* obj[index_obj].eta            = eta_mean; */
			/* obj[index_obj].eta_sigma      = eta_sigma; */
			/* obj[index_obj].modelparameter  = modelparameter_mean; */
			/* obj[index_obj].modelparameter_sigma = modelparameter_sigma; */
			/* obj[index_obj].reflectanceratio = reflectanceratio_mean; */
			/* obj[index_obj].reflectanceratio_error = reflectanceratio_sigma; */
			/* obj[index_obj].absmag         = absmag_mean; */
			/* obj[index_obj].absmag_error   = absmag_sigma; */
			/* obj[index_obj].slopepar       = slopepar_mean; */
			/* obj[index_obj].slopepar_error = slopepar_sigma; */
			/* obj[index_obj].flux_raw       = flux_raw_mean; */
			/* obj[index_obj].fluxerr        = flux_raw_sigma; */
			/* obj[index_obj].chi2 = 0; */
		}
		
		//apply flux correction
		if ( corrflux == 1 )
			fluxcorrection ( labrat, index_obj, index_obj+number_of_measurements );
		
		printf(". %s (error estimation): d=%.3f (%+.3f/%+.3f, 3sig:%+.3f/%+.3f), pv=%.3f (%+.3f/%+.3f, 3sig: %+.3f/%+.3f), eta=%.3f (%+.3f/%+.3f), sslat=%.3f (%+.3f/%+.3f)\n",
					obj[index_obj].name,
					obj[index_obj].diameter, obj[index_obj].diameter_lower, obj[index_obj].diameter_upper,
		                        obj[index_obj].diameter_lower3sig, obj[index_obj].diameter_upper3sig,
					obj[index_obj].pv, obj[index_obj].pv_lower, obj[index_obj].pv_upper,
		                        obj[index_obj].pv_lower3sig, obj[index_obj].pv_upper3sig,
					obj[index_obj].eta, obj[index_obj].eta_lower, obj[index_obj].eta_upper,
					obj[index_obj].modelparameter, obj[index_obj].modelparameter_lower, obj[index_obj].modelparameter_upper);
		
		//write result to file
		if ( long_output_in_file(obj[index_obj], "error_results.long", "") == 1 )
				printf("ERROR: no output to file possible!\n");
		if ( short_output_in_file(obj[index_obj], "error_results", "") == 1 )
				printf("ERROR: no output to file possible!\n");
		
		// crap? ...
		//set index_obj to index_labrat+1 to indicate that this object already has been processed
		//if ( fitting == 1 ) index_obj = index_labrat+1;

		//find index_obj of next object
		int all_done = 0;
		while ( all_done == 0) 
		{
			if ( obj[index_obj].next_measurement != 0)
				index_obj = obj[index_obj].next_measurement;
			else 
				all_done = 1;
		}

	}
}


#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
