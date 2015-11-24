/*=================================================================
 *
 * gait2de.c
 *
 * Explicit differential equation for 2D musculoskeletal model : dx/dt = f(x,u)
 * This is the source code for the MEX function gait2de.mexw32
 * The musculoskeletal model is documented in the file gait2de_reference.odt
 * The function documentation is in gait2de.m
 *
 * Copyright 2009-2011 Orchard Kinetics LLC
 *
 *=================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mex.h"
#include "gait2de.h"

// Compilation settings
#define VERBOSE 0					// set this to 1 for debugging

// size of the model (some other size constants are in gait2d.h)
#define NMUS 16						// number of muscles 
#define NSTATES 2*NDOF+2*NMUS		// number of system state variables, 2*NDOF + 2*NMUS

// M_PI is known in gcc but not in Visual C++ 2008
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

// define struct that holds muscle properties
typedef struct {
	char *name;				// name
	double Lceopt; 			// Optimal length of CE (m)
	double Width;			// Width of CE force-length relationship (m)
	double Fmax;			// Maximal isometric force of CE (N)
	double SEEslack;		// Slack length of the SEE, in meters
	double PEEslack;		// Slack length of the PEE, relative to Lceopt
	double L0;				// Muscle+tendon length when all DOFs are zero
	double MA[NMOM];		// Moment arms (positive when muscle causes ccw rotation of distal joint, model facing +X)
	// dependent parameters (calculated from others during preprocessing)
	double kSEE;			// Stiffness parameter of SEE, in Fmax/m^2		
} muscleprop;

// global variables
static int initialized = 0;
static muscleprop muscles[NMUS];		// contains all muscle properties
static param_struct param;				// contains other model parameters

// Other parameters	
static double par_JointK1, par_JointK2, par_JointB;					// passive joint stiffness and damping parameters
static double par_MinAngle[NMOM], par_MaxAngle[NMOM];	// passive joint range of motion
	
// General muscle parameters
static double par_HillA;					// Normalized Hill parameter a/Fmax for F-v relationship (usually 0.25)
static double par_Gmax;					// Maximum eccentric muscle force, relative to Fmax
static double par_Vmax;					// Max. contraction velocity in Lceopt/s
static double par_Tact, par_Tdeact;			// activation and deactivation time constants
static double par_gmax;					// maximum eccentric force
static double par_umax;					// strain of SEE at Fmax load
static double par_kPEE;					// Stiffness parameter of PEE, in Fmax/Lceopt^2
static double par_muscledamping;		// Muscle damping

// ===================================================================================
// SetParameters: set all model parameters
// ===================================================================================
void SetParameters() {
	int i,j,k;
	
	// multibody model parameters, from Winter book for 75 kg body mass and 1.8 m body height
	param.TrunkMass 	= 50.8500;
	param.ThighMass 	= 7.5000;
	param.ShankMass 	= 3.4875;
	param.FootMass 		= 1.0875;

	param.TrunkInertia 	= 3.1777;
	param.ThighInertia 	= 0.1522;
	param.ShankInertia 	= 0.0624;
	param.FootInertia 	= 0.0184;

	param.TrunkCMy 		= 0.3155;
	param.ThighCMy 		= -0.1910;
	param.ShankCMy 		= -0.1917;
	param.FootCMx		= 0.0768;
	param.FootCMy 		= -0.0351;
	param.ThighLen 		= 0.4410;
	param.ShankLen 		= 0.4428;
	
	// contact model parameters
	param.ContactHeelX	= -0.06;			// X coordinate of heel contact point
	param.ContactToeX	= 0.15;				// X coordinate of toe contact point
	param.ContactY		= -0.07;			// Y coordinate of both contact points
	param.ContactStiff	= 5e7;				// ground contact stiffness, N/m^3
	param.ContactDamp	= 0.85;				// ground contact damping, s/m
	param.ContactFric	= 1.0;         		// friction coefficient
	param.ContactV0		= 0.01;				// velocity constant (m/s), for |ve| = vc -> |fx| = 0.4621*c*fy
	
	// passive joint moment parameters
	par_JointK1			= 1;				// overall joint stiffness (Nm/rad)
	par_JointK2			= 10000;			// stiffness at joint limits (Nm/rad^2)
	par_JointB			=  1;				// joint damping (Nms/rad), exists always
	par_MinAngle[0] 	=  -30;				// Rhip
	par_MaxAngle[0] 	=  160;		
	par_MinAngle[1] 	= -160;				// Rknee
	par_MaxAngle[1] 	=   -5;
	par_MinAngle[2] 	=  -60;				// Rankle
	par_MaxAngle[2] 	=   60;
	// copy the right side range of motion into the left side
	for (i=0; i<NMOM/2; i++) {
		j = i + NMOM/2;
		par_MinAngle[j] = par_MinAngle[i];
		par_MaxAngle[j] = par_MaxAngle[i];
	}

	// muscle general parameters
	par_HillA			= 0.25;				// normalized Hill constant a/Fmax
	par_Gmax			= 1.8;				// max. eccentric force relative to Fmax
	par_kPEE			= 1.0;
	par_umax			= 0.04;
	par_Vmax			= 10.0;				// max shortening velocity (m/s)
	par_muscledamping	= 0.001;			// muscle damping (Fmax per Lceopt per sec)
	par_Tact			= 0.01;				// activation time constant (s), from Opensim gait model
	par_Tdeact			= 0.04;				// deactivation time constant (s), from Opensim gait model

	// set all moment arms to zero (non-zero values will be defined later)
	for (i=0; i<NMUS; i++) {
		for (j=0; j<NMOM; j++) {
			muscles[i].MA[j] = 0.0;
		}
	}
	
	// muscle-specific parameters
	i=0;
	muscles[i].name 	= "R.Iliopsoas";
	muscles[i].Fmax 	= 1500;
	muscles[i].Lceopt 	= 0.102;
	muscles[i].Width 	= 1.298;
	muscles[i].SEEslack	= 0.142;
	muscles[i].PEEslack = 1.2;
	muscles[i].L0		= 0.248;
	muscles[i].MA[0] 	= 0.05;	

	i++;
	muscles[i].name 	= "R.Glutei";
	muscles[i].Fmax 	= 3000;
	muscles[i].Lceopt 	= 0.200;
	muscles[i].Width 	= 0.625;
	muscles[i].SEEslack	= 0.157;
	muscles[i].PEEslack = 1.2;
	muscles[i].L0		= 0.271;
	muscles[i].MA[0] 	= -0.062;	
	
	i++;
	muscles[i].name 	= "R.Hamstrings";
	muscles[i].Fmax 	= 3000;
	muscles[i].Lceopt 	= 0.104;
	muscles[i].Width 	= 1.197;
	muscles[i].SEEslack	= 0.334;
	muscles[i].PEEslack = 1.2;
	muscles[i].L0		= 0.383;
	muscles[i].MA[0] 	= -0.072;		// hip moment arm
	muscles[i].MA[1]	= -0.034;		// knee moment arm
	
	i++;
	muscles[i].name 	= "R.Rectus";
	muscles[i].Fmax 	= 1200;
	muscles[i].Lceopt 	= 0.081;
	muscles[i].Width 	= 1.443;
	muscles[i].SEEslack	= 0.398;
	muscles[i].PEEslack = 1.4;
	muscles[i].L0		= 0.474;
	muscles[i].MA[0] 	= 0.034;		// hip moment arm
	muscles[i].MA[1]	= 0.050;		// knee moment arm
	
	i++;
	muscles[i].name 	= "R.Vasti";
	muscles[i].Fmax 	= 7000;
	muscles[i].Lceopt 	= 0.093;
	muscles[i].Width 	= 0.627;
	muscles[i].SEEslack	= 0.223;
	muscles[i].PEEslack = 1.4;
	muscles[i].L0		= 0.271;
	muscles[i].MA[1]	= 0.042;		// knee moment arm
	
	i++;
	muscles[i].name 	= "R.Gastrocnemius";
	muscles[i].Fmax 	= 3000;
	muscles[i].Lceopt 	= 0.055;
	muscles[i].Width 	= 0.888;
	muscles[i].SEEslack	= 0.420;
	muscles[i].PEEslack = 1.2;
	muscles[i].L0		= 0.487;
	muscles[i].MA[1]	= -0.020;		// knee moment arm
	muscles[i].MA[2]	= -0.053;		// ankle moment arm
	
	i++;
	muscles[i].name 	= "R.Soleus";
	muscles[i].Fmax 	= 4000;
	muscles[i].Lceopt 	= 0.055;
	muscles[i].Width 	= 1.039;
	muscles[i].SEEslack	= 0.245;
	muscles[i].PEEslack = 1.2;
	muscles[i].L0		= 0.284;
	muscles[i].MA[2]	= -0.053;		// ankle moment arm
	
	i++;
	muscles[i].name 	= "R.TibialisAnterior";
	muscles[i].Fmax 	= 2500;
	muscles[i].Lceopt 	= 0.082;
	muscles[i].Width 	= 0.442;
	muscles[i].SEEslack	= 0.317;
	muscles[i].PEEslack = 1.2;
	muscles[i].L0		= 0.381;
	muscles[i].MA[2]	= 0.037;		// ankle moment arm
	
	// make the left side muscles by copying the right side muscles
	for (i=0; i<NMUS/2; i++) {
		j = i + NMUS/2;
		muscles[j].name = (char *) calloc(strlen(muscles[i].name)+1, sizeof(char));
		strcpy(muscles[j].name, muscles[i].name);
		muscles[j].name[0] = 'L';
		muscles[j].Fmax 	= muscles[i].Fmax;
		muscles[j].Lceopt 	= muscles[i].Lceopt;
		muscles[j].Width 	= muscles[i].Width;
		muscles[j].SEEslack	= muscles[i].SEEslack;
		muscles[j].PEEslack = muscles[i].PEEslack;
		muscles[j].L0		= muscles[i].L0;
		for (k=0; k<NMOM/2; k++) {
			muscles[j].MA[k+3] = muscles[i].MA[k];	
		}
	}
	
	// printf("MOMENT ARM MATRIX:\n");
	// for (i=0; i<NMUS; i++) {
		// for (j=0; j<NMOM; j++) printf("%8.4f ",muscles[i].MA[j]);
		// printf("\n");
	// }
	
// Preprocessing and error checking
	for (i=0; i<NMOM; i++) {
		par_MinAngle[i] = par_MinAngle[i]*M_PI/180.0;
		par_MaxAngle[i] = par_MaxAngle[i]*M_PI/180.0;
		if (par_MinAngle[i] > par_MaxAngle[i]) {
			printf("Error in joint %d\n", i);
			mexErrMsgTxt("Max angle must be greater than min angle.");
		}
	}
	for (i=0; i<NMUS; i++) {
		// compute kSEE in units of Fmax m^-2
		muscles[i].kSEE = 1.0/pow(par_umax * muscles[i].SEEslack, 2);
		// In this model we use the theoretical Width of 0.56 (Walker & Schrodt paper)	
		muscles[i].Width = 0.56;
	}
}

// ===================================================================================
// MuscleDynamics: the explicit muscle dynamics, returns Lcedot = f(a,Lce,Lm) and also computes muscle force
// ===================================================================================
double MuscleDynamics(muscleprop *m, double Act, double Lce, double Lm, double *Force) {

	// Lce is given relative to Lceopt and forces are all computed relative to Fmax

	double Fsee, Fpee, Fce, x, Fiso, Lcedot, Vmax, cc;
	double lambda,a,b,c;
		
	// Determine force in SEE from current SEE length
	x = Lm - Lce*m->Lceopt - m->SEEslack;		// elongation of SEE, in meters
	Fsee = x / m->Fmax;							// start with 1 N/m bidirectional spring
	if (x > 0) {
		Fsee = Fsee + m->kSEE * x * x;			// quadratic term when elongated past slack
	}
	*Force = m->Fmax * Fsee;					// this is the SEE force in Newtons

	// subtract PEE force to get CE force
	x = Lce - m->PEEslack;				// elongation of PEE, in Lceopt units
	Fpee = x * m->Lceopt / m->Fmax; 	// start with 1 N/m bidirectional spring
    if (x > 0) {
		Fpee = Fpee + par_kPEE * x * x;
	}
	Fce = Fsee - Fpee;
	// Note: it is OK if Fce becomes negative (compressive), the force-velocity model allows it and responds correctly.
	
	// Fiso is the normalized isometric force-length relationship at max activation, normalized to Fmax
	x = (Lce - 1.0)/m->Width;
	Fiso = exp(-x*x);

	// activation dependent scaling of Vmax according to Chow & Darling (1999)
    lambda = 1 - exp(-3.82 * Act) + Act * exp(-3.82);
	Vmax = lambda * par_Vmax;
		
	// calculate Lcedot from Fce using Hill-Katz force-velocity model with additional damping
	if (Fce < Act*Fiso)	{					// concentric
		a = -par_muscledamping/par_HillA;
		b = (Act*Fiso + par_muscledamping*Vmax + Fce/par_HillA);
		c = Vmax*(Act*Fiso-Fce);
	}
	else {									// eccentric
		// c parameter for continuity between the two hyperbolas
		cc = Vmax*par_HillA*(par_Gmax-1)/(par_HillA+1);
		a = par_muscledamping;
		b = (par_Gmax*Act*Fiso + par_muscledamping*cc - Fce);
		c = cc*(Act*Fiso-Fce);
	}
		
	return (-b + sqrt(b*b-4*a*c))/(2*a);

}

// =========================================================================
// mexFunction: this is the actual MEX function interface
// =========================================================================
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	// working variables
	int i, j, k, nrows, ncols, iframe, nframes;
	double d, ang, angvel;
	
	// muscle variables
	double Lm[NMUS];				// muscle-tendon lengths
	double adot[NMUS];
	double Lcedot[NMUS];
	double force[NMUS];

	// multibody dynamics variables
	double *q, *qd;
	double qdd[NDOF];
	double mom[NDOF];
	double stick[2*NSTICK];
	double grf[6];
	
	// MEX function pointers to inputs and outputs
	double *x, *xdot, *u, *M;
	double *mforces_out, *mom_out;
	double *grf_out, *stick_out;
	
	// Initialize the model, if needed
	if (initialized != 1959) {
		printf("********************************************************************\n");
		printf("* GAIT2DE -- A Musculoskeletal Dynamics Model for Gait and Posture *\n");  
		printf("*         Version 1.0 -- Compiled: "__DATE__" "__TIME__"            *\n"); 
		printf("*              Licensed for non-commercial use only                *\n");
		printf("*                  (c) 2011 Orchard Kinetics LLC                   *\n");
		printf("*                                                                  *\n"); 
		printf("*  This software, and all associated files, may be distributed     *\n"); 
		printf("*  freely without restriction, provided that this text is not      *\n"); 
		printf("*  altered.  The latest version may be downloaded from             *\n");
        printf("*  http://www.orchardkinetics.com.                                 *\n"); 
		printf("*                                                                  *\n"); 
		printf("* THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED    *\n"); 
		printf("* BY APPLICABLE LAW. EXCEPT WHEN OTHERWISE STATED IN WRITING THE   *\n"); 
		printf("* COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM       *\n");
		printf("* \"AS IS\" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR        *\n");
		printf("* IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES   *\n");
		printf("* OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE     *\n");
		printf("* ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS  *\n");
		printf("* WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE     *\n");
		printf("* COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.           *\n");	
		printf("********************************************************************\n");
		printf("Initializing model...\n");
		SetParameters();
		initialized = 1959;
	}
	
	if (nrhs<2 || nrhs>3) {
		mexErrMsgTxt("gait2d: must have 2 or 3 inputs.");
	}

	// State of model is the first input (required)
	nrows = mxGetM(prhs[0]);
	ncols = mxGetN(prhs[0]);
	if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) )
		mexErrMsgTxt("gait2de: Incorrect type for state vector x, must be double.");
	if (nrows != NSTATES)
		mexErrMsgTxt("gait2de: Incorrect size for state vector x, must have 50 rows.");
	x = mxGetPr(prhs[0]);
	nframes = ncols;
		
	// Controls are second input (required)
	nrows = mxGetM(prhs[1]);
	ncols = mxGetN(prhs[1]);
	if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ) mexErrMsgTxt("gait2d: Incorrect type for u, must be double.");
	if (nrows != NMUS) mexErrMsgTxt("gait2de: Incorrect size for u, must have 16 rows.");
	if (ncols != nframes) mexErrMsgTxt("gait2de: Controls u must have same number of columns as state x.");
	u = mxGetPr(prhs[1]);
	
	// see if additional actuation was given
	if (nrhs == 3) {
		nrows = mxGetM(prhs[2]);
		ncols = mxGetN(prhs[2]);
		if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ) mexErrMsgTxt("gait2d: Incorrect type for M, must be double.");
		if (nrows != NDOF) mexErrMsgTxt("gait2de: Incorrect size for actuation M, must have 9 rows.");
		if (ncols != nframes) mexErrMsgTxt("gait2de: Actuation M must have same number of columns as state x.");
		M = mxGetPr(prhs[2]);
	}
	
	// Create matrix for the xdot output of the MEX function
	plhs[0] = mxCreateDoubleMatrix(NSTATES, nframes, mxREAL);
	xdot = mxGetPr(plhs[0]);

	// Create matrix for the optional GRF output of the MEX function
	if (nlhs > 1) {
		plhs[1] = mxCreateDoubleMatrix(6, nframes, mxREAL);
		grf_out = mxGetPr(plhs[1]);	
	}
	
	// Create matrix for the optional stick output of the MEX function
	if (nlhs > 2) {
		plhs[2] = mxCreateDoubleMatrix(NSTICK*2, nframes, mxREAL);
		stick_out = mxGetPr(plhs[2]);	
	}
	
	// Create matrix for the optional muscle forces output of the MEX function
	if (nlhs > 3) {
		plhs[3] = mxCreateDoubleMatrix(NMUS, nframes, mxREAL);
		mforces_out = mxGetPr(plhs[3]);	
	}

	// Create matrix for the optional joint moment output of the MEX function
	if (nlhs > 4) {
		plhs[4] = mxCreateDoubleMatrix(NMOM, nframes, mxREAL);
		mom_out = mxGetPr(plhs[4]);
	}
	
	for (iframe=0; iframe < nframes; iframe++) {

		// Compute the muscle Lcedot and muscle forces
		for(i=0; i<NMUS; i++) {	
			// Calculate da/dt from muscle activation model
			adot[i] = (u[i] - x[2*NDOF+NMUS+i]) * (u[i]/par_Tact + (1-u[i])/par_Tdeact);
		
			// Calculate muscle length Lm and derivatives dLm/dq from generalized coordinates in x
			Lm[i] = muscles[i].L0;			// muscle-tendon length when all angles are zero
			// Add contributions from joint angles, which are x[3+j]
			for (j=0; j<NMOM; j++) Lm[i] = Lm[i] - muscles[i].MA[j]*x[3+j];
			
			// Calculate Lcedot for the muscle
			Lcedot[i] = MuscleDynamics(&muscles[i],
				x[2*NDOF+NMUS+i],		// active state of muscle i
				x[2*NDOF+i],			// Lce of muscle i
				Lm[i],					// muscle length
				&force[i]);				// muscle force (output)
				
			// copy muscle forces to the MEX function output variable, if needed
			if (nlhs > 3) *(mforces_out++) = force[i];
		}

		// Compute the generalized forces for input to Autolev code
		for (i=0; i<NDOF; i++) {
			// if requested, use the extra actuation given as inputs
			if (nrhs == 3) {
				mom[i] = M[i];
			}
			else {
				mom[i] = 0.0;
			}
			
			// if this DOF is a joint angle, apply passive joint moments and muscle moments
			if (i>2) {
				ang = x[i];											// joint angle is one of the state variables
				angvel = x[NDOF+i];									// the corresponding angular velocity
				
				// is angle above upper limit of ROM?
				d = ang - par_MaxAngle[i-3];
				if (d > 0.0) {	
					mom[i] = -par_JointK2 * d*d;
				}
				
				// is angle below lower limit of ROM?
				d = ang - par_MinAngle[i-3];
				if (d < 0.0) {
					mom[i] = par_JointK2 * d*d;
				}
				
				// add a small amount of damping and overall stiffness
				mom[i] = mom[i] - par_JointB * angvel - par_JointK1 * ang;
					
				// add the muscle moments
				for (j=0; j<NMUS; j++) {
					mom[i] = mom[i] + muscles[j].MA[i-3] * force[j];
				}
				
				// copy joint moments to the MEX function output variable, if needed
				if (nlhs > 4) *(mom_out++) = mom[i];
			}			
		}
		
		// Call the C function that was generated by Autolev
		q = &x[0];
		qd = &x[NDOF];	
		// for (i=0; i<NDOF; i++) printf("mom[%d] = %f\n", i, mom[i]);
		gait2d_al(&param, q, qd, qdd, mom, grf, stick);
		
		// copy ground reaction forces to the MEX function output variable, if needed
		if (nlhs > 1) for (i=0; i<6; i++) *(grf_out++) = grf[i];

		// copy stick figure data to the MEX function output variable, if needed
		if (nlhs > 2) for (i=0; i<2*NSTICK; i++) *(stick_out++) = stick[i];

		// copy the elements of xdot into the MEX function output, columnwise
		for (i=0; i<NDOF; i++) *(xdot++) = x[NDOF+i];		// the first NDOF rows are: qdot = dq/dt
		for (i=0; i<NDOF; i++) *(xdot++) = qdd[i];			// the next NDOF rows are the equations of motion from Autolev
		for (i=0; i<NMUS; i++) *(xdot++) = Lcedot[i];		// the next NMUS rows are the muscle Lcedot values
		for (i=0; i<NMUS; i++) *(xdot++) = adot[i];			// the final NMUS rows are the muscle activation dynamics: da/dt = ...
	
		// advance the input pointers to the next frame
		x = x + NSTATES;
		u = u + NMUS;
		M = M + NDOF;
		
		
	
	}
	
	return;
	
}
