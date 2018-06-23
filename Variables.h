//  Variables.h

#ifndef Variables_h
#define Variables_h

/* Define constants appearing in the parametrization (R,Theta) -> (r,h), critical exponents, 
		and combinations thereof (the c's coefficients) */
#define A -0.76201
#define B 0.00804
#define beta 0.326
#define delta 4.8
#define alpha 0.11
#define M0 0.605
#define h0 0.364
#define c0 0.0424369
#define c1 0.321329
#define c2 -1.20375
#define c3 -0.00126032
#define bede 1.5648

/* Dummy Variables for parameter input*/
double v1,v2,v3,v4,v5,v6;
char v0[10];

/* Parameters of the Ising -> QCD map */
double dTC,dmuBC,TC,muBC,angle1,angle2,ww,rho;

/* Parameters used in case the choice of parametrization is such that the CEP lies on a straight line or a parabola. */
double T0, kappa, anglediff;

/* Parameters that take on the values obtained from the inverse Jacobian of the Ising -> QCD map. */
double drdmuB,dhdmuB,drdT,dhdT;

/* Dummy variable for data i/o */
double xIn1,xIn2,xIn3,xIn4,xIn5;

/* Integers that will be used as indices. */
int i,j,k,l,x1int,x2int,kint,lint;

/* Vectors and integer used in the root finding routine. */
double *x,*f;
int check;

/* Variables used for the Gauss filter. */
int kmin,kmax,lmin,lmax;
double sigmax,sigmay,sum,norm;

/* Define the dummy variables for the root finding "for" loop. */ 
double Ti, muBi;

/* Variables used to export/use T and muB as doubles. */
double Tval,muBval,cure;

/* Variables used in the merging with the HRG pressure. */
double Tprime,Targum,deltaT;

/* These variables are defined to be used inside the expression for the G and derivatives thereof, to speed up the process. */
double gg0Num, gg1Num, gg2Num, htil0Num, htil1Num, htil2Num;
double d2GdR2Num, d2GdTheta2Num, d2GdRdThetaNum, dGdRNum, dGdThetaNum;
double dGdrConhNum, dGdhConrNum;
double dThetadRConhNum, dRdThetaConhNum, dThetadRConrNum, dRdThetaConrNum;
double drdRConhNum, drdThetaConhNum, dhdRConrNum, dhdThetaConrNum, d2rdR2ConhNum, d2rdTheta2ConhNum, d2hdR2ConrNum, d2hdTheta2ConrNum;
double dRdrConhNum, dThetadrConhNum, dRdhConrNum, dThetadhConrNum, d2Rdr2ConhNum, d2Thetadr2ConhNum, d2Rdh2ConrNum, d2Thetadh2ConrNum;
double d2RdrdhNum, d2ThetadrdhNum;
double d2Gdr2ConhNum, d2Gdh2ConrNum, d2GdrdhNum;
double g0,htilnum,htilprimenum,htilsecondnum,m1Th2,ggnum,ggprimenum,ggsecondnum,
		tau,ggamma,omega,omegap,ggammap,zeta,zetap,Th2,R1malpha,R1momegap; 

/* These vectors are used to store the expansion coefficients of the pressure at different T's, 
		as well as derivatives of the non-Ising one wrt T. */
double *Chi0LatVec, *Chi2LatVec, *Chi4LatVec, *Chi6LatVec;
double *Chi0IsingVec, *Chi2IsingVec, *Chi4IsingVec, *Chi6IsingVec;
double *Chi0NoIsingVec, *Chi2NoIsingVec, *Chi4NoIsingVec, *Chi6NoIsingVec;
double *dChi0NoIsingdTVec, *dChi2NoIsingdTVec, *dChi4NoIsingdTVec, *dChi6NoIsingdTVec;
double *d2Chi0NoIsingdT2Vec, *d2Chi2NoIsingdT2Vec, *d2Chi4NoIsingdT2Vec, *d2Chi6NoIsingdT2Vec;

/* The variables are used to export the Chi's from the Ising model, and the 3D critical pressure. */
double Chi0Ising, Chi2Ising, Chi4Ising, Pij, dPdTij, dPdmuBij, d2PdT2ij, d2PdmuB2ij, d2PdTdmuBij;

/* Matrix used for the inverse Jacobian of the Ising -> QCD map. */
double **JJ;

/* Matrices to store/export pressure(s) and other thermodynamic functions over the phase diagram. */
double **PressNoIsingMat,**dPressNoIsingdTMat,**d2PressNoIsingdT2Mat,**dPressNoIsingdmuBMat,**d2PressNoIsingdmuB2Mat,
			**d2PressNoIsingdTdmuBMat;
double **PressNoIsingFilterMat,**dPressNoIsingFilterdTMat,**d2PressNoIsingFilterdT2Mat,**dPressNoIsingFilterdmuBMat,
			**d2PressNoIsingFilterdmuB2Mat,**d2PressNoIsingFilterdTdmuBMat;
double **PressIsingMat,**dPressIsingdTMat,**d2PressIsingdT2Mat,**dPressIsingdmuBMat,**d2PressIsingdmuB2Mat,
			**d2PressIsingdTdmuBMat;
double **PressTotMat,**dPressTotdTMat,**d2PressTotdT2Mat,**dPressTotdmuBMat,**d2PressTotdmuB2Mat,
			**d2PressTotdTdmuBMat;
double **PressHRGMat,**dPressHRGdTMat,**d2PressHRGdT2Mat,**dPressHRGdmuBMat,**d2PressHRGdmuB2Mat,
			**d2PressHRGdTdmuBMat;
double **PressTotHRGMat,**dPressTotHRGdTMat,**d2PressTotHRGdT2Mat,**dPressTotHRGdmuBMat,**d2PressTotHRGdmuB2Mat,
			**d2PressTotHRGdTdmuBMat;
double **PressFinalMat,**dPdTMat,**dPdmuBMat,**d2PdT2Mat,**d2PdmuB2Mat,
			**d2PdTdmuBMat,**EntropyFinalMat, **BarDensityFinalMat, **EnerDensityFinalMat, **SpSoundFinalMat, **Chi2FinalMat;

/* A Rank=3 tensor to store the coordinates. */
double ***Coords;


#endif /* Variables_h */
