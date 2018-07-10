/* 
	Copyright (c) 2018, Paolo Parotto and Debora Mroczek, 
	Department of Physics, University of Houston, Houston, TX 77204, US.
*/

/*
	This file produces an EoS matching Lattice QCD at muB=0, and containing a critical point in the 3D Ising 
		universality class, in a parametrized form.
 	It allows for different choices of constraints on the shape of the critical line, which reduce the number of parameters.
*/ 

#define NRANSI

/* Import standard libraries */
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
//#include <direct.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* Import libraries from Numerical Recipes */
#include "nrD.h"
#include "nrutilD.h"

/* Import additional library for many variables used in the program */
#include "Variables.h"
#include "Functions.h"

/* Define PI and the number of variables (2) in the root finding procedure */
#define PI 3.141592653589
#define N 2


/* The main body of the program. */
int main(void)
{
	/* Define strings for filenames; used for exporting files.*/
	char nameCoords[100],nameFolder[100],nameMainfolder[100],
			nameChisIsing[100],nameChisNoIsing[100],namedChisNoIsingdT[100],named2ChisNoIsingdT2[100],
			namePressIsing3D[100],namedPdTIsing3D[100],namedPdmuBIsing3D[100],
			named2PdT2Ising3D[100],named2PdmuB2Ising3D[100],named2PdTdmuBIsing3D[100],
			namePressNoIsing3D[100],namedPdTNoIsing3D[100],namedPdmuBNoIsing3D[100],
			named2PdT2NoIsing3D[100],named2PdmuB2NoIsing3D[100],named2PdTdmuBNoIsing3D[100],
			namePressTotHRG3D[100],namedPdTTotHRG3D[100],namedPdmuBTotHRG3D[100],
			named2PdT2TotHRG3D[100],named2PdmuB2TotHRG3D[100],named2PdTdmuBTotHRG3D[100],			
			namePressIsingPlusNoIsing3D[100],namedPdTIsingPlusNoIsing3D[100],namedPdmuBIsingPlusNoIsing3D[100],
			named2PdT2IsingPlusNoIsing3D[100],named2PdmuB2IsingPlusNoIsing3D[100],named2PdTdmuBIsingPlusNoIsing3D[100],
			namePressFinal3D[100],nameEnerDensFinal3D[100],nameBarDensFinal3D[100],nameEntrFinal3D[100],
			nameSpsoundFinal3D[100],nameChi2Final3D[100];

	/* Define the name of the main folder where the program lives and the files we wish to import are located.*/
	strcpy(nameMainfolder,"../");
	
	/* Define time variables to measure time elapsed creating the coordinates tables.*/
	time_t start, stop;

	/* Define char variables to be used in selecting the parameter input mode.*/	
	char MODE, MODESTR[10];

	/* All vectors, matrices and tensors are initialized. */
	/* Vectors for the root finding. */
	x=vector(1,N);	f=vector(1,N);
	/* Matrix for Jacobian. */
	JJ=matrix(1,2,1,2);
	/* Rank 3 tensor for coordinates. */
	Coords=f3tensor(5,821,0,1202,1,2);
	/* Vectors for Chi's. */
	Chi0LatVec=vector(5,821);			Chi2LatVec=vector(5,821);			Chi4LatVec=vector(5,821);			
	Chi0IsingVec=vector(5,821);			Chi2IsingVec=vector(5,821);			Chi4IsingVec=vector(5,821);			
	Chi0NoIsingVec=vector(5,821);		Chi2NoIsingVec=vector(5,821);		Chi4NoIsingVec=vector(5,821);		
	dChi0NoIsingdTVec=vector(5,821);	dChi2NoIsingdTVec=vector(5,821);	dChi4NoIsingdTVec=vector(5,821);	
	d2Chi0NoIsingdT2Vec=vector(5,821);	d2Chi2NoIsingdT2Vec=vector(5,821);	d2Chi4NoIsingdT2Vec=vector(5,821);	
	/* Matrices for thermodynamic functions over the phase diagram. */
	PressNoIsingMat=matrix(5,821,0,601);			dPressNoIsingdTMat=matrix(5,821,0,601);				dPressNoIsingdmuBMat=matrix(5,821,0,601);	
	d2PressNoIsingdT2Mat=matrix(5,821,0,601);		d2PressNoIsingdmuB2Mat=matrix(5,821,0,601);			d2PressNoIsingdTdmuBMat=matrix(5,821,0,601);

	PressNoIsingFilterMat=matrix(5,821,0,601);		dPressNoIsingFilterdTMat=matrix(5,821,0,601);		dPressNoIsingFilterdmuBMat=matrix(5,821,0,601);	
	d2PressNoIsingFilterdT2Mat=matrix(5,821,0,601);	d2PressNoIsingFilterdmuB2Mat=matrix(5,821,0,601);	d2PressNoIsingFilterdTdmuBMat=matrix(5,821,0,601);

	PressIsingMat=matrix(5,821,0,601);				dPressIsingdTMat=matrix(5,821,0,601);				dPressIsingdmuBMat=matrix(5,821,0,601);	
	d2PressIsingdT2Mat=matrix(5,821,0,601);			d2PressIsingdmuB2Mat=matrix(5,821,0,601);			d2PressIsingdTdmuBMat=matrix(5,821,0,601);

	PressTotMat=matrix(5,821,0,601);				dPressTotdTMat=matrix(5,821,0,601);					dPressTotdmuBMat=matrix(5,821,0,601);	
	d2PressTotdT2Mat=matrix(5,821,0,601);			d2PressTotdmuB2Mat=matrix(5,821,0,601);				d2PressTotdTdmuBMat=matrix(5,821,0,601);

	PressHRGMat=matrix(5,821,0,601);				dPressHRGdTMat=matrix(5,821,0,601);					dPressHRGdmuBMat=matrix(5,821,0,601);	
	d2PressHRGdT2Mat=matrix(5,821,0,601);			d2PressHRGdmuB2Mat=matrix(5,821,0,601);				d2PressHRGdTdmuBMat=matrix(5,821,0,601);

	PressTotHRGMat=matrix(5,821,0,601);				dPressTotHRGdTMat=matrix(5,821,0,601);				dPressTotHRGdmuBMat=matrix(5,821,0,601);	
	d2PressTotHRGdT2Mat=matrix(5,821,0,601);		d2PressTotHRGdmuB2Mat=matrix(5,821,0,601);			d2PressTotHRGdTdmuBMat=matrix(5,821,0,601);

	PressFinalMat=matrix(5,821,0,601);				EntropyFinalMat=matrix(5,821,0,601);				BarDensityFinalMat=matrix(5,821,0,601);
	EnerDensityFinalMat=matrix(5,821,0,601);		SpSoundFinalMat=matrix(5,821,0,601);				Chi2FinalMat=matrix(5,821,0,601);			



	/* Chi's from Lattice are imported and stored in a vector. */
	printf("Storing Lattice Data... \n \n");
	FILE *LatIn = fopen("EOS_tables/Lattice_Data_5_821_dT1.dat","r");
	if (LatIn == 0){
    	fprintf(stderr, "failed to open Lattice Data\n");
    	exit(1);
  	}
  	FILE *Latout = fopen("EOS_tables/Chis_Lat_5_821_dT1.dat","w");
	for(i=5;fscanf(LatIn,"%lf %lf %lf %lf",&xIn1,&xIn2,&xIn3,&xIn4) !=EOF;i++){
    	Chi0LatVec[i] = xIn2*pow(i,4);
    	Chi2LatVec[i] = xIn3*pow(i,4);
    	Chi4LatVec[i] = xIn4*pow(i,4);
 
    	fprintf(Latout,"%3.1f %12.16f	%12.16f	%12.16f\n",(double) i,Chi0LatVec[i],Chi2LatVec[i],Chi4LatVec[i]);
  	}
  	fclose(LatIn);
  	fclose(Latout);
  	printf("Successfully stored Lattice Data \n \n");


	/* Parameters are read from the file parameters.dat*/
  	FILE *ParametersIn = fopen("parameters.dat", "r");
  	if (ParametersIn == 0){
  		fprintf(stderr,"failed to open paremeters file\n");
  		exit(1);
  	}
  	for(i=5;fscanf(ParametersIn,"%s %lf %lf %lf %lf %lf %lf", v0, &v1, &v2, &v3, &v4, &v5, &v6) !=EOF; i++){

  		const char * MODE = v0;
  		printf("MODE: %s\n", MODE);

  		if(strncmp(MODE,"FRE",3) == 0){
  			strcpy(MODESTR, "FREE");
  			TC = v1;
  			muBC = v2;
  			angle1 = v3;
  			angle2 = v4;
  			ww = v5;
  			rho = v6;

  			printf("The parameters you entered are:\nTC = %f\nmuBC = %f\nangle1 = %f\nangle2 = %f\nw = %f\nrho = %f\n",TC,muBC,angle1,angle2,ww,rho);
  		}

  		if(strncmp(MODE,"PAR",3) == 0){
  			strcpy(MODESTR, "PAR");
  			T0 = v1;
  			kappa = v2;
  			muBC = v3;
  			anglediff = v4;
  			ww = v5;
  			rho = v6;

  			TC = T0 + kappa/T0 * muBC * muBC; angle1 = 180/PI*fabs(atan(-2.0*kappa/T0*muBC)); angle2 = angle1 + anglediff;
  			printf("PARABOLA\n");
  			printf("The parameters you entered are:\nTC = %f\nmuBC = %f\nangle1 = %f\nangle2 = %f\nw = %f\nrho = %f\n",TC,muBC,angle1,angle2,ww,rho);
  		}

  		if(strncmp(MODE,"STR",3) == 0){
  			strcpy(MODESTR, "STR");
  			T0 = v1;
  			muBC = v2;
  			angle1 = v3;
  			angle2 = v4;
  			ww = v5;
  			rho = v6;
  			TC = T0 - atan(angle1*PI/180)*muBC;

  			printf("The parameters you entered are:\nTC = %f\nmuBC = %f\nangle1 = %f\nangle2 = %f\nw = %f\nrho = %f\n",TC,muBC,angle1,angle2,ww,rho);
  		}

  	break;
  	}
  	fclose(ParametersIn);

	/* Determine dTC and dmuBC from TC, muBC, ww and rho. */
	dTC = ww*TC; 
	dmuBC = dTC*rho;

	/* Determine the filenames according to the parameter choice made. */
	sprintf(nameFolder, "Files_%s_%d_%d_%d_%d_%d_%d",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);

	sprintf(nameCoords, "Coords_%s_%d_%d_%d_%d_%d_%d.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameChisIsing, "Chis_Ising_%s_%d_%d_%d_%d_%d_%d_muB0.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameChisNoIsing, "Chis_No_Ising_%s_%d_%d_%d_%d_%d_%d_muB0.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedChisNoIsingdT, "dChis_No_Ising_dT_%s_%d_%d_%d_%d_%d_%d_muB0.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2ChisNoIsingdT2, "d2Chis_No_Ising_dT2%s_%d_%d_%d_%d_%d_%d_muB0.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namePressIsing3D, "Press_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdTIsing3D, "dPdT_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdmuBIsing3D, "dPdmuB_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdT2Ising3D, "d2PdT2_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdmuB2Ising3D, "d2PdmuB2_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdTdmuBIsing3D, "d2PdTdmuB_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namePressNoIsing3D, "Press_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdTNoIsing3D, "dPdT_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdmuBNoIsing3D, "dPdmuB_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdT2NoIsing3D, "d2PdT2_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdmuB2NoIsing3D, "d2PdmuB2_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdTdmuBNoIsing3D, "d2PdTdmuB_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namePressIsingPlusNoIsing3D, "Press_Ising_Plus_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdTIsingPlusNoIsing3D, "dPdT_Ising_Plus_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdmuBIsingPlusNoIsing3D, "dPdmuB_Ising_Plus_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdT2IsingPlusNoIsing3D, "d2PdT2_Ising_Plus_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdmuB2IsingPlusNoIsing3D, "d2PdmuB2_Ising_Plus_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdTdmuBIsingPlusNoIsing3D, "d2PdTdmuB_Ising_Plus_No_Ising_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namePressTotHRG3D, "Press_Tot_HRG_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdTTotHRG3D, "dPdT_Tot_HRG_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namedPdmuBTotHRG3D, "dPdmuB_Tot_HRG_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdT2TotHRG3D, "d2PdT2_Tot_HRG_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdmuB2TotHRG3D, "d2PdmuB2_Tot_HRG_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(named2PdTdmuBTotHRG3D, "d2PdTdmuB_Tot_HRG_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(namePressFinal3D, "Press_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameEntrFinal3D, "Entr_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameBarDensFinal3D, "BarDens_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameEnerDensFinal3D, "EnerDens_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameSpsoundFinal3D, "SpSound_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);
	sprintf(nameChi2Final3D, "Chi2_Final_%s_%d_%d_%d_%d_%d_%d_3D.dat",MODESTR,(int)TC,(int)muBC,(int)angle1,(int)angle2,(int) dTC,(int) dmuBC);


	/* Create folder with filename, and move to the folder.*/
	mkdir(nameFolder, S_IRWXU | S_IRWXG | S_IRWXO);
	chdir(nameFolder);


	/* For calculations hereafter, angles are expressed in radians. */
	angle1 = angle1*PI/180;
	angle2 = angle2*PI/180;


	/* Once all the parameters are set, the Jacobian matrix of the (r,h) -> (T,muB) map is generated and
			a value is obtained for drdmuB and dhdmuB.  */
	Jacobian(JJ,dTC,dmuBC,angle1,angle2);
	drdmuB = JJ[1][2];
	dhdmuB = JJ[2][2];
	drdT = JJ[1][1];
	dhdT = JJ[2][1];
	printf("\nThe Jacobian calculation gave:\ndrdmuB = %12.16f\ndhdT = %12.16f\ndrdT = %12.16f\ndhmuB = %12.16f",drdmuB,dhdmuB,drdT,dhdT);


	/* Now everything is set, the program checks if the coordinates file corresponding to the current choice of parameters 
		is present in the working folder. If yes, it skips; if not, it generates it, assigning it the name stored in 
		the string "namecoords".  */
	FILE *FileCoords=fopen(nameCoords,"w");
	printf("\n\nCreating the coordinates file...\n");
	/* Time is taken at start and end to measure how long it took to generate the coordinates file.*/
	time(&start);
	for(Ti=5.0;Ti<=821.0;Ti+=1.0){
		for(muBi=-601.0;muBi<=601.0;muBi+=1.0){
			/* Condition for choosing init condition in root finding. */
			if(muBi<=0.0){x[1]=20.0;x[2]=0.0;}
			else{
				if(Ti>=condition(muBi)){x[1]=25.0;x[2]=0.5;}
				if(Ti<condition(muBi)){x[1]=25.0;x[2]=-0.5;}			
			}
			/* The actual root finding routine. */
			newt(x,N,&check,funcv);
			funcv(N,x,f);
			if (check) printf("Convergence problems.\n");
			/* Print to file the coordinates {T,muB,R,Theta} (in this order). */
			fprintf(FileCoords,"%3.1f %3.1f %2.16f %2.16f\n",Ti,muBi,x[1],x[2]);
			//	printf("%12.16f %12.16f %12.16f %12.16f\n",Ti,muBi,x[1],x[2]);
		}
	}
	time(&stop);	
	/* A control on whether the loops reached the end, and did not get stuck, otherwise file is eliminated. */	
	if(Ti == 822 && muBi == 602)	printf("\nThe file %s was created successfully in %d seconds\n",
														 nameCoords, (int) difftime(stop,start));
	else{
		printf("\nThe file %s was NOT created successfully, removing it..\n",nameCoords);
		remove(nameCoords);
	} 								
	fclose(FileCoords);


	/* Coordinates are imported and stored in a Rank 3 tensor Coords. */
	printf("\nImporting coordinates from file %s \n",nameCoords);
	FileCoords=fopen(nameCoords,"r");
	for(j=0, x1int=4;fscanf(FileCoords,"%lf %lf %lf %lf\n", &xIn1, &xIn2, &xIn3, &xIn4) != EOF; j++ ){
		x2int = j % 1203;
		if (x2int == 0) x1int++;
  		Coords[x1int][x2int][1] = xIn3;
		Coords[x1int][x2int][2] = xIn4; 
	}
	fclose(FileCoords);
    printf("\nImporting coordinates successful\n");


    /* The coordinates corresponding to muB=0 are used to generate the critical Chi's at different temperatures., and exported. */
	printf("\nGenerating critical chis at muB=0\n");
	FILE *FileChiIsing=fopen(nameChisIsing,"w");
	for(i=5;i<=821;i++){
		Tval = (double) i;
		Chi0IsingVec[i] = -G(Coords[i][601][1],Coords[i][601][2])*pow(TC,4);
		Chi2IsingVec[i] = -Tval*Tval*d2GdmuB2ConT(Coords[i][601][1],Coords[i][601][2])*pow(TC,4);
		Chi4IsingVec[i] = -Tval*Tval*Tval*Tval*d4GdmuB4ConT(Coords[i][601][1],Coords[i][601][2])*pow(TC,4);

		fprintf(FileChiIsing,"%3.1f %12.16f %12.16f %12.16f\n",(double) i,Chi0IsingVec[i],Chi2IsingVec[i],Chi4IsingVec[i]);
	}
	fclose(FileChiIsing);
	printf("\nGenerating chis successful\n");


	/* The pressure is generated over the whole range of coordinates.
		NOTE: the pressure is symmetrized around muB=0 to ensure all odd order derivatives vanish at muB=0. */
	/*NOTE: file output is currently disabled. (Un)comment file opening and closing, as well as the export line in the loop,
		to turn on/off the file export.*/
	printf("\nGenerating critical pressure in 3D\n");
	//FILE *FilePressIsing=fopen(namePressIsing3D,"w");			FILE *FiledPdTIsing=fopen(namedPdTIsing3D,"w");
	//FILE *FiledPdmuBIsing=fopen(namedPdmuBIsing3D,"w");		FILE *Filed2PdT2Ising=fopen(named2PdT2Ising3D,"w");
	//FILE *Filed2PdmuB2Ising=fopen(named2PdmuB2Ising3D,"w");	FILE *Filed2PdTdmuBIsing=fopen(named2PdTdmuBIsing3D,"w");
	//Calculate pressure in 3D
	for(i=5;i<=820;i++){
		for(j=0;j<=601;j++){
			k = j + 601;
			if(j==0){
				PressIsingMat[i][j] = - G(Coords[i][k][1],Coords[i][k][2])*pow(TC,4);
				dPressIsingdTMat[i][j] = - dGdTConmuB(Coords[i][k][1],Coords[i][k][2])*pow(TC,4);
				dPressIsingdmuBMat[i][j] = - dGdmuBConT(Coords[i][k][1],Coords[i][k][2])*pow(TC,4);
				d2PressIsingdT2Mat[i][j] = - d2GdT2ConmuB(Coords[i][k][1],Coords[i][k][2])*pow(TC,4);
				d2PressIsingdmuB2Mat[i][j] = - d2GdmuB2ConT(Coords[i][k][1],Coords[i][k][2])*pow(TC,4);
				d2PressIsingdTdmuBMat[i][j] = - d2GdTdmuB(Coords[i][k][1],Coords[i][k][2])*pow(TC,4);
			}
	 		else{
	 			PressIsingMat[i][j] = - (G(Coords[i][k][1],Coords[i][k][2]) + G(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);
	 			dPressIsingdTMat[i][j] = - (dGdTConmuB(Coords[i][k][1],Coords[i][k][2]) + dGdTConmuB(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);
	 			dPressIsingdmuBMat[i][j] = - (dGdmuBConT(Coords[i][k][1],Coords[i][k][2]) - dGdmuBConT(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);
	 			d2PressIsingdT2Mat[i][j] = - (d2GdT2ConmuB(Coords[i][k][1],Coords[i][k][2]) + d2GdT2ConmuB(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);
	 			d2PressIsingdmuB2Mat[i][j] = - (d2GdmuB2ConT(Coords[i][k][1],Coords[i][k][2]) + d2GdmuB2ConT(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);
	 			d2PressIsingdTdmuBMat[i][j] = - (d2GdTdmuB(Coords[i][k][1],Coords[i][k][2]) - d2GdTdmuB(Coords[i][k-2*j][1],Coords[i][k-2*j][2]))/2.0*pow(TC,4);
	 		}
			
	 		//fprintf(FilePressIsing,"%3.1f	%3.1f	%12.16f\n",(double) j,(double) i, PressIsingMat[i][j]);
	 		//fprintf(FiledPdTIsing,"%3.1f	%3.1f	%12.16f\n",(double) j,(double) i, dPressIsingdTMat[i][j]);
	 		//fprintf(FiledPdmuBIsing,"%3.1f	%3.1f	%12.16f\n",(double) j,(double) i, dPressIsingdmuBMat[i][j]);
	 		//fprintf(Filed2PdT2Ising,"%3.1f	%3.1f	%12.16f\n",(double) j,(double) i, d2PressIsingdT2Mat[i][j]);
	 		//fprintf(Filed2PdmuB2Ising,"%3.1f	%3.1f	%12.16f\n",(double) j,(double) i, d2PressIsingdmuB2Mat[i][j]);
	 		//fprintf(Filed2PdTdmuBIsing,"%3.1f	%3.1f	%12.16f\n",(double) j,(double) i, d2PressIsingdTdmuBMat[i][j]);
	 	}
	}
	//fclose(FilePressIsing); 	fclose(FiledPdTIsing); 			fclose(FiledPdmuBIsing); 
	//fclose(Filed2PdT2Ising);	fclose(Filed2PdmuB2Ising);		fclose(Filed2PdTdmuBIsing);


	/* Non-Ising Chi's are calculated, exported and stored. */
	/*NOTE: file output is currently disabled. (Un)comment file opening and closing, as well as the export line in the loop,
		to turn on/off the file export.*/
	printf("\nGenerating non-Ising Chi's at muB = 0\n");
	FILE *FileChiNoIsing = fopen(nameChisNoIsing,"w");
	for (i=5; i<=821; i++){
	    Tval = (double) i;
	    Chi0NoIsingVec[i] = Chi0LatVec[i] - Chi0IsingVec[i];
	    Chi2NoIsingVec[i] = Chi2LatVec[i] - Chi2IsingVec[i];
	    Chi4NoIsingVec[i] = Chi4LatVec[i] - Chi4IsingVec[i];
        
	   	fprintf(FileChiNoIsing,"%3.1f %12.16f %12.16f %12.16f\n", Tval, Chi0NoIsingVec[i],Chi2NoIsingVec[i],Chi4NoIsingVec[i]);
	}
	fclose(FileChiNoIsing);


	/* The non-Ising pressure in 3D is calculated from Taylor expansion, stored and exported. */ 
	/*NOTE: file output is currently disabled. (Un)comment file opening and closing, as well as the export line in the loop,
		to turn on/off the file export.*/
   	printf("\nCalculating non-Ising Pressure in 3D.\n");
	//FILE *FilePressNoIsing=fopen(namePressNoIsing3D,"w");			
	for (i=5; i<=821; i++){
		for(j=0; j<=601; j++){
		    muBval = (double) j;
		    Tval = (double) i;
       
    		PressNoIsingMat[i][j] = Chi0NoIsingVec[i]+(1.0/2.0)*Chi2NoIsingVec[i]*pow(muBval/Tval,2)
    												+(1.0/24.0)*Chi4NoIsingVec[i]*pow(muBval/Tval,4);

            //fprintf(FilePressNoIsing,"%3.1f	%3.1f	%12.16f\n",muBval,Tval,PressNoIsingMat[i][j]);
   		}
	}
	//fclose(FilePressNoIsing); 


	/* Filtering of regular pressure here. Widths of the filter are also set here. Result is stored and exported. */
	/*NOTE: file output is currently disabled. (Un)comment file opening and closing, as well as the export line in the loop,
		to turn on/off the file export.*/
	printf("\nFiltering non-Ising Pressure in 3D\n");
	sigmax = 18.0;	sigmay = 1.0;
	//FILE *FilePressNoIsingFilter = fopen("NoIsingPressFilter.dat","w");
	for (i=5; i<=821; i++){
		for(j=0; j<=601; j++){
			muBval = (double) j;
			Tval = (double) i;

 			kmin = i-sigmax;	kmax = i+sigmax;
 			lmin = j-sigmay;	lmax = j+sigmay;	
 
			sum = 0.0; norm = 0.0;
			for(k=kmin;k<=kmax;k++){
				for(l=lmin;l<=lmax;l++){
					if(k<=5)		kint = 5;
					else if(k>=821)	kint = 821;	
					else		 	kint = k;
					
					if(l<=0)  		lint = -l;
					else if(l>=601)	lint = 601;		
					else			lint = l;
									
					sum+=GausFunc(kint-i,lint-j,sigmax,sigmay)*PressNoIsingMat[kint][lint];
					norm+=GausFunc(kint-i,lint-j,sigmax,sigmay);
				}
			}
			PressNoIsingFilterMat[i][j] = sum/norm;
			//fprintf(FilePressNoIsingFilter, "%3.1f %3.1f %12.16f\n", muBval,Tval,PressNoIsingFilterMat[i][j]);
		}
	}
	//fclose(FilePressNoIsingFilter);


	/* Derivatives (numerical) of the filtered Non-Ising pressure and all thermodynamic quantities of interest are calculated,
		and stored, over the phase diagram. */
	printf("\nCalculating (numerical) derivatives of the filtered Non-Ising Pressure\n");

	for (i=5; i<=820; i++)  for (j=0;j<=601; j++)	dPressNoIsingFilterdTMat[i][j] = PressNoIsingFilterMat[i+1][j] - PressNoIsingFilterMat[i][j]; 
	
	for (i=5; i<=821; i++)	for (j=0;j<=600; j++)	dPressNoIsingFilterdmuBMat[i][j] = PressNoIsingFilterMat[i][j+1] - PressNoIsingFilterMat[i][j];
	
	for (i=5; i<=819; i++)  for (j=0;j<=601; j++)	d2PressNoIsingFilterdT2Mat[i][j] = dPressNoIsingFilterdTMat[i+1][j] - dPressNoIsingFilterdTMat[i][j]; 
	
	for (i=5; i<=821; i++)  for (j=0;j<=599; j++)	d2PressNoIsingFilterdmuB2Mat[i][j] = dPressNoIsingFilterdmuBMat[i][j+1] - dPressNoIsingFilterdmuBMat[i][j];

	for (i=5; i<=820; i++)  for (j=0;j<=600; j++)	d2PressNoIsingFilterdTdmuBMat[i][j] = dPressNoIsingFilterdTMat[i][j+1] - dPressNoIsingFilterdTMat[i][j];
	


  	/*	Calculate pressure and derivatives summing regular pressure and critical pressure, then normalize. Store and export. */
	/*NOTE: file output is currently disabled. (Un)comment file opening and closing, as well as the export line in the loop,
		to turn on/off the file export.*/
    printf("\nCalculating other derivatives in 3D \n");
  	//FILE *FilePressTot = fopen(namePressIsingPlusNoIsing3D, "w");					FILE *FiledPdTIsingPlusNoIsing=fopen(namedPdTIsingPlusNoIsing3D,"w");
	//FILE *FiledPdmuBIsingPlusNoIsing=fopen(namedPdmuBIsingPlusNoIsing3D,"w");		FILE *Filed2PdT2IsingPlusNoIsing=fopen(named2PdT2IsingPlusNoIsing3D,"w");
	//FILE *Filed2PdmuB2IsingPlusNoIsing=fopen(named2PdmuB2IsingPlusNoIsing3D,"w");	FILE *Filed2PdTdmuBIsingPlusNoIsing=fopen(named2PdTdmuBIsingPlusNoIsing3D,"w");
	for (i=5; i<=819; i++){
	    for(j=0;j<=599; j++){
		    Tval = (double) i;
		    muBval = (double) j;
		    if(muBval == 0) cure=0.0;
		    else cure=1.0;

			PressTotMat[i][j] = PressIsingMat[i][j] + PressNoIsingFilterMat[i][j];
    		dPressTotdTMat[i][j] = dPressIsingdTMat[i][j] + dPressNoIsingFilterdTMat[i][j];
    		dPressTotdmuBMat[i][j] = (dPressIsingdmuBMat[i][j] + dPressNoIsingFilterdmuBMat[i][j])*cure;
    		d2PressTotdT2Mat[i][j] = d2PressIsingdT2Mat[i][j] + d2PressNoIsingFilterdT2Mat[i][j];
    		d2PressTotdmuB2Mat[i][j] = d2PressIsingdmuB2Mat[i][j] + d2PressNoIsingFilterdmuB2Mat[i][j];  		
    		d2PressTotdTdmuBMat[i][j] = d2PressIsingdTdmuBMat[i][j] + d2PressNoIsingFilterdTdmuBMat[i][j];

    		//fprintf(FilePressTot,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, PressTotMat[i][j]);
	 		//fprintf(FiledPdTIsingPlusNoIsing,"%3.1f	%3.1f	%12.16f\n",muBval,Tval,dPressTotdTMat[i][j]);
	 		//fprintf(FiledPdmuBIsingPlusNoIsing,"%3.1f	%3.1f	%12.16f\n",muBval,Tval,dPressTotdmuBMat[i][j]);
	 		//fprintf(Filed2PdT2IsingPlusNoIsing,"%3.1f	%3.1f	%12.16f\n",muBval,Tval,d2PressTotdT2Mat[i][j]);
	 		//fprintf(Filed2PdmuB2IsingPlusNoIsing,"%3.1f	%3.1f	%12.16f\n",muBval,Tval,d2PressTotdmuB2Mat[i][j]);
	 		//fprintf(Filed2PdTdmuBIsingPlusNoIsing,"%3.1f	%3.1f	%12.16f\n",muBval,Tval,d2PressTotdTdmuBMat[i][j]);
    	}
	}
	//fclose(FilePressTot);					fclose(FiledPdTIsingPlusNoIsing); 		fclose(FiledPdmuBIsingPlusNoIsing); 
	//fclose(Filed2PdT2IsingPlusNoIsing);	fclose(Filed2PdmuB2IsingPlusNoIsing);	fclose(Filed2PdTdmuBIsingPlusNoIsing);



//---------------MERGING WITH HRG BEYOND THIS POINT---------------------------//

	chdir(nameMainfolder);

  	/* Now, the smooth joining with the HRG pressure starts. Import HRG pressure and store. */
	printf("\nImporting HRG Pressure \n");
	FILE *FilePressHRG = fopen("EOS_tables/Press_HRG_MUB000601_T005300_dT1.dat", "r");
	if (FilePressHRG == 0){
	    fprintf(stderr, "failed to open HRG Pressure \n");
	    exit(1);
	}
	for(i=0, x2int=0;fscanf(FilePressHRG,"%lf %lf %lf\n",&xIn1,&xIn2,&xIn3) !=EOF;i++){
	    x1int = (i % 817) + 5;
	    PressHRGMat[x1int][x2int] = xIn3*pow(x1int,4);

   	    if(x1int == 821) x2int++;
	}
	fclose(FilePressHRG); 

	chdir(nameFolder);

  	/*	Calculate derivatives of the HRG pressure. Store and export. */
	/*NOTE: file output is currently disabled. (Un)comment file opening and closing, as well as the export line in the loop,
		to turn on/off the file export.*/
    printf("\nCalculating derivatives of HRG pressure in 3D \n");
    //FILE *FiledPdTHRG=fopen("dPdT_HRG_3D.dat","w");
	//FILE *FiledPdmuBHRG=fopen("dPdmuB_HRG_3D.dat","w");		FILE *Filed2PdT2HRG=fopen("d2PdT2_HRG_3D.dat","w");
	//FILE *Filed2PdmuB2HRG=fopen("d2PdmuB2_HRG_3D.dat","w");	FILE *Filed2PdTdmuBHRG=fopen("d2PdTdmuB_HRG_3D.dat","w");
	for (i=5; i<=820; i++){
	    for(j=0;j<=600; j++){
		    Tval = (double) i;
		    muBval = (double) j;

    		dPressHRGdTMat[i][j] = PressHRGMat[i+1][j] - PressHRGMat[i][j];
    		dPressHRGdmuBMat[i][j] = PressHRGMat[i][j+1] - PressHRGMat[i][j];

	 		//fprintf(FiledPdTHRG,"%3.1f	%3.1f	%12.16f\n",muBval,Tval,dPressHRGdTMat[i][j]);
	 		//fprintf(FiledPdmuBHRG,"%3.1f	%3.1f	%12.16f\n",muBval,Tval,dPressHRGdmuBMat[i][j]);
    	}
	}
	for (i=5; i<=819; i++){
	    for(j=0;j<=599; j++){
		    Tval = (double) i;
		    muBval = (double) j;

	    	d2PressHRGdT2Mat[i][j] = dPressHRGdTMat[i+1][j] - dPressHRGdTMat[i][j];
    		d2PressHRGdmuB2Mat[i][j] = dPressHRGdmuBMat[i][j+1] - dPressHRGdmuBMat[i][j];  		
    		d2PressHRGdTdmuBMat[i][j] = dPressHRGdmuBMat[i+1][j] - dPressHRGdmuBMat[i][j];
    		
    		//fprintf(Filed2PdT2HRG,"%3.1f	%3.1f	%12.16f\n",muBval,Tval,d2PressHRGdT2Mat[i][j]);
	 		//fprintf(Filed2PdmuB2HRG,"%3.1f	%3.1f	%12.16f\n",muBval,Tval,d2PressHRGdmuB2Mat[i][j]);
	 		//fprintf(Filed2PdTdmuBHRG,"%3.1f	%3.1f	%12.16f\n",muBval,Tval,d2PressHRGdTdmuBMat[i][j]);
		}
	}
	//fclose(FiledPdTHRG); 		fclose(FiledPdmuBHRG); 
	//fclose(Filed2PdT2HRG); 	fclose(Filed2PdmuB2HRG);	fclose(Filed2PdTdmuBHRG);


	/* Final pressure, joined with the HRG one is calculated, stored and exported. The parameter deltaT for merging with HRG is also set here. */
	/*NOTE: file output is currently disabled. (Un)comment file opening and closing, as well as the export line in the loop,
		to turn on/off the file export.*/
	printf("\nMerging with HRG data.\n");
  	//FILE *FilePressTotHRG = fopen(namePressTotHRG3D, "w");		FILE *FiledPdTTotHRG=fopen(namedPdTTotHRG3D,"w");
	//FILE *FiledPdmuBTotHRG=fopen(namedPdmuBTotHRG3D,"w");		FILE *Filed2PdT2TotHRG=fopen(named2PdT2TotHRG3D,"w");
	//FILE *Filed2PdmuB2TotHRG=fopen(named2PdmuB2TotHRG3D,"w");	FILE *Filed2PdTdmuBTotHRG=fopen(named2PdTdmuBTotHRG3D,"w");
	deltaT = 17.0;
	for (i=5; i<=819; i++){
	    for(j=0;j<=599; j++){
		    Tval = (double) i;
		    muBval = (double) j;
    		Targum = (Tval - (T0 + kappa/T0*muBval*muBval - 23.0))/deltaT;

    		if(Targum >= 4.5)    	PressTotHRGMat[i][j] = PressTotMat[i][j];
    		if(Targum <= -4.5)      PressTotHRGMat[i][j] = PressHRGMat[i][j];
    		else                    PressTotHRGMat[i][j] = PressTotMat[i][j]*0.5*(1.0 + tanh(Targum)) 
                                            + PressHRGMat[i][j]*0.5*(1.0 - tanh(Targum));
	
            if(Targum >= 4.5)       dPressTotHRGdTMat[i][j] = dPressTotdTMat[i][j];
    		if(Targum <= -4.5)      dPressTotHRGdTMat[i][j] = dPressHRGdTMat[i][j];
    		else                    dPressTotHRGdTMat[i][j] = dPressTotdTMat[i][j]*0.5*(1.0 + tanh(Targum))
                                            + dPressHRGdTMat[i][j]*0.5*(1.0 - tanh(Targum))
                                        	- PressTotMat[i][j]*0.5/deltaT*1.0/pow(cosh(Targum),2)
                                        	+ PressHRGMat[i][j]*0.5/deltaT*1.0/pow(cosh(Targum),2);       

            if(Targum >= 4.5)       dPressTotHRGdmuBMat[i][j] = dPressTotdmuBMat[i][j];
   			if(Targum <= -4.5)      dPressTotHRGdmuBMat[i][j] = dPressHRGdmuBMat[i][j];
    		else                    dPressTotHRGdmuBMat[i][j] = dPressTotdmuBMat[i][j]*0.5*(1.0 + tanh(Targum))
                                            + dPressHRGdmuBMat[i][j]*0.5*(1.0 - tanh(Targum))
                                        	- PressTotMat[i][j]*0.5/pow(cosh(Targum),2)*(2*kappa/T0*muBval)/deltaT
                                        	+ PressHRGMat[i][j]*0.5/pow(cosh(Targum),2)*(2*kappa/T0*muBval)/deltaT;

            if(Targum >= 4.5)       d2PressTotHRGdT2Mat[i][j] = d2PressTotdT2Mat[i][j];
    		if(Targum <= -4.5)      d2PressTotHRGdT2Mat[i][j] = d2PressHRGdT2Mat[i][j];
    		else                    d2PressTotHRGdT2Mat[i][j] = d2PressTotdT2Mat[i][j]*0.5*(1.0 + tanh(Targum)) 
                                            + d2PressHRGdT2Mat[i][j]*0.5*(1.0 - tanh(Targum))
                                	        + dPressTotdTMat[i][j]/deltaT*1.0/pow(cosh(Targum),2)
                                        	- dPressHRGdTMat[i][j]/deltaT*1.0/pow(cosh(Targum),2)
                                        	- PressTotMat[i][j]/pow(cosh(Targum),2)*tanh(Targum)/pow(deltaT,2) 
                                        	+ PressHRGMat[i][j]/pow(cosh(Targum),2)*tanh(Targum)/pow(deltaT,2);   

            if(Targum >= 4.5)       d2PressTotHRGdmuB2Mat[i][j] = d2PressTotdmuB2Mat[i][j];
    		if(Targum <= -4.5)      d2PressTotHRGdmuB2Mat[i][j] = d2PressHRGdmuB2Mat[i][j];
    		else                    d2PressTotHRGdmuB2Mat[i][j] = d2PressTotdmuB2Mat[i][j]*0.5*(1.0 + tanh(Targum))
                                            + d2PressHRGdmuB2Mat[i][j]*0.5*(1.0 - tanh(Targum))
                                	        - dPressTotdmuBMat[i][j]/pow(cosh(Targum),2)*(2*kappa/T0*muBval)/deltaT
                                        	+ dPressHRGdmuBMat[i][j]/pow(cosh(Targum),2)*(2*kappa/T0*muBval)/deltaT
                                        	- PressTotMat[i][j]/pow(cosh(Targum),2)*(2*tanh(Targum)
                                        		*(2*kappa/T0*muBval)/deltaT
                                        		*(2*kappa/T0*muBval)/deltaT
                                        		+ 2*kappa/T0/deltaT)
                                        	+ PressHRGMat[i][j]/pow(cosh(Targum),2)*(2*tanh(Targum)
                                        		*(2*kappa/T0*muBval)/deltaT
                                        		*(2*kappa/T0*muBval)/deltaT
                                        		+ 2*kappa/T0/deltaT);

            if(Targum >= 4.5)       d2PressTotHRGdTdmuBMat[i][j] = d2PressTotdTdmuBMat[i][j];
    		if(Targum <= -4.5)      d2PressTotHRGdTdmuBMat[i][j] = d2PressHRGdTdmuBMat[i][j];
    		else                    d2PressTotHRGdTdmuBMat[i][j] = d2PressTotdTdmuBMat[i][j]*0.5*(1.0 + tanh(Targum))
                                            + d2PressHRGdTdmuBMat[i][j]*0.5*(1.0 - tanh(Targum))
                                        	+ dPressTotdmuBMat[i][j]*0.5/deltaT*1.0/pow(cosh(Targum),2)
                                        	- dPressHRGdmuBMat[i][j]*0.5/deltaT*1.0/pow(cosh(Targum),2) 
                             	            + dPressTotdTMat[i][j]*0.5/deltaT*(2*kappa/T0*muBval)/pow(cosh(Targum),2)
                                        	- dPressHRGdTMat[i][j]*0.5/deltaT*(2*kappa/T0*muBval)/pow(cosh(Targum),2)
                                        	+ PressTotMat[i][j]/pow(cosh(Targum),2)/pow(deltaT,2)
                                        		*tanh(Targum)*(2*kappa/T0*muBval)
                                        	- PressHRGMat[i][j]/pow(cosh(Targum),2)/pow(deltaT,2)
                                        		*tanh(Targum)*(2*kappa/T0*muBval);    

       		//fprintf(FilePressTotHRG,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, PressTotHRGMat[i][j]);
    		//fprintf(FiledPdTTotHRG,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, dPressTotHRGdTMat[i][j]);
    		//fprintf(FiledPdmuBTotHRG,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, dPressTotHRGdmuBMat[i][j]);
    		//fprintf(Filed2PdT2TotHRG,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, d2PressTotHRGdT2Mat[i][j]);
 		  	//fprintf(Filed2PdmuB2TotHRG,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, d2PressTotHRGdmuB2Mat[i][j]);
	  		//fprintf(Filed2PdTdmuBTotHRG,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, d2PressTotHRGdTdmuBMat[i][j]);                        	
   		}
	}
	//fclose(FilePressTotHRG);	fclose(FiledPdTTotHRG); 		fclose(FiledPdmuBTotHRG); 
	//fclose(Filed2PdT2TotHRG);	fclose(Filed2PdmuB2TotHRG);		fclose(Filed2PdTdmuBTotHRG);


	/* Now all the derivatives have been merged with the HRG, we can finish with the thermodynamics and export the quantities we want,
			in particular we will normalize with suitable powers of the temperature to obtain dimensionless quantities.*/
	printf("\nMerging with HRG complete. \n \nCalculating thermodynamics quantities. \n");
	FILE *FilePressFinal = fopen(namePressFinal3D, "w");			FILE *FileEntrFinal = fopen(nameEntrFinal3D, "w");
	FILE *FileBarDensFinal = fopen(nameBarDensFinal3D, "w");		FILE *FileEnerDensFinal = fopen(nameEnerDensFinal3D, "w");
	FILE *FileSpSoundFinal = fopen(nameSpsoundFinal3D, "w");		FILE *FileChi2Final = fopen(nameChi2Final3D, "w");
	for (i=30; i<=800; i++){
	    for(j=0;j<=450; j++){
		    Tval = (double) i;
		    muBval = (double) j;
		   	if(muBval == 0) cure=0.0;
		    else cure=1.0;

    		PressFinalMat[i][j] = PressTotHRGMat[i][j]/pow(Tval,4);
    		EntropyFinalMat[i][j] = dPressTotHRGdTMat[i][j]/pow(Tval,3);
    		BarDensityFinalMat[i][j] = (dPressTotHRGdmuBMat[i][j]/pow(Tval,3))*cure;
    		EnerDensityFinalMat[i][j] = (dPressTotHRGdTMat[i][j]*Tval - PressTotHRGMat[i][j] + muBval*dPressTotHRGdmuBMat[i][j])/pow(Tval,4);
    		SpSoundFinalMat[i][j] = (dPressTotHRGdmuBMat[i][j]*dPressTotHRGdmuBMat[i][j]*d2PressTotHRGdT2Mat[i][j] - 2.0*dPressTotHRGdTMat[i][j]
    								*dPressTotHRGdmuBMat[i][j]*d2PressTotHRGdTdmuBMat[i][j] + dPressTotHRGdTMat[i][j]
    								*dPressTotHRGdTMat[i][j]*d2PressTotHRGdmuB2Mat[i][j])
    								*1.0/(dPressTotHRGdTMat[i][j]*Tval + muBval*dPressTotHRGdmuBMat[i][j])
    								*1.0/(d2PressTotHRGdT2Mat[i][j]*d2PressTotHRGdmuB2Mat[i][j]-d2PressTotHRGdTdmuBMat[i][j]
    								*d2PressTotHRGdTdmuBMat[i][j]);     								
	      	Chi2FinalMat[i][j] = d2PressTotHRGdmuB2Mat[i][j]/pow(Tval,2);

    		fprintf(FilePressFinal,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, PressFinalMat[i][j]);
    		fprintf(FileEntrFinal,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, EntropyFinalMat[i][j]);
    		fprintf(FileBarDensFinal,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, BarDensityFinalMat[i][j]);
    		fprintf(FileEnerDensFinal,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, EnerDensityFinalMat[i][j]);
 		  	fprintf(FileSpSoundFinal,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, SpSoundFinalMat[i][j]);
	  		fprintf(FileChi2Final,"%3.16f	%3.1f	%12.16f \n", muBval, Tval, Chi2FinalMat[i][j]);
    	}
	}
	fclose(FilePressFinal);			fclose(FileEntrFinal);			fclose(FileBarDensFinal);
	fclose(FileEnerDensFinal);		fclose(FileSpSoundFinal);		fclose(FileChi2Final);	

	printf("\nProcedure completed.\n");

	/* Free memory allocations for vectors and matrices. */

	free_vector(x,1,N);
	free_vector(f,1,N);

	free_matrix(JJ,1,2,1,2);

	free_f3tensor(Coords,5,821,0,1202,1,2);

	free_vector(Chi0LatVec,5,821);		free_vector(Chi2LatVec,5,821);		free_vector(Chi4LatVec,5,821);		
	free_vector(Chi0IsingVec,5,821);	free_vector(Chi2IsingVec,5,821);	free_vector(Chi4IsingVec,5,821);	
	free_vector(Chi0NoIsingVec,5,821);	free_vector(Chi2NoIsingVec,5,821);	free_vector(Chi4NoIsingVec,5,821);	
	
	free_matrix(PressNoIsingMat,5,821,0,601);		free_matrix(dPressNoIsingdTMat,5,821,0,601);		free_matrix(dPressNoIsingdmuBMat,5,821,0,601);
	free_matrix(d2PressNoIsingdT2Mat,5,821,0,601);	free_matrix(d2PressNoIsingdmuB2Mat,5,821,0,601);	free_matrix(d2PressNoIsingdTdmuBMat,5,821,0,601);

	free_matrix(PressNoIsingFilterMat,5,821,0,601);			free_matrix(dPressNoIsingFilterdTMat,5,821,0,601);		free_matrix(dPressNoIsingFilterdmuBMat,5,821,0,601);
	free_matrix(d2PressNoIsingFilterdT2Mat,5,821,0,601);	free_matrix(d2PressNoIsingFilterdmuB2Mat,5,821,0,601);	free_matrix(d2PressNoIsingFilterdTdmuBMat,5,821,0,601);

	free_matrix(PressIsingMat,5,821,0,601);			free_matrix(dPressIsingdTMat,5,821,0,601);			free_matrix(dPressIsingdmuBMat,5,821,0,601);
	free_matrix(d2PressIsingdT2Mat,5,821,0,601);	free_matrix(d2PressIsingdmuB2Mat,5,821,0,601);		free_matrix(d2PressIsingdTdmuBMat,5,821,0,601);
	
	free_matrix(PressTotMat,5,821,0,601);			free_matrix(dPressTotdTMat,5,821,0,601);			free_matrix(dPressTotdmuBMat,5,821,0,601);
	free_matrix(d2PressTotdT2Mat,5,821,0,601);		free_matrix(d2PressTotdmuB2Mat,5,821,0,601);		free_matrix(d2PressTotdTdmuBMat,5,821,0,601);
	
	free_matrix(PressHRGMat,5,821,0,601);			free_matrix(dPressHRGdTMat,5,821,0,601);			free_matrix(dPressHRGdmuBMat,5,821,0,601);
	free_matrix(d2PressHRGdT2Mat,5,821,0,601);		free_matrix(d2PressHRGdmuB2Mat,5,821,0,601);		free_matrix(d2PressHRGdTdmuBMat,5,821,0,601);

	free_matrix(PressTotHRGMat,5,821,0,601);		free_matrix(dPressTotHRGdTMat,5,821,0,601);			free_matrix(dPressTotHRGdmuBMat,5,821,0,601);
	free_matrix(d2PressTotHRGdT2Mat,5,821,0,601);	free_matrix(d2PressTotHRGdmuB2Mat,5,821,0,601);		free_matrix(d2PressTotHRGdTdmuBMat,5,821,0,601);
	
	free_matrix(PressFinalMat,5,821,0,601);			free_matrix(EntropyFinalMat,5,821,0,601);			free_matrix(BarDensityFinalMat,5,821,0,601);
	free_matrix(EnerDensityFinalMat,5,821,0,601);	free_matrix(SpSoundFinalMat,5,821,0,601);			free_matrix(Chi2FinalMat,5,821,0,601);	


	return 0;
}


#undef NRANSI
