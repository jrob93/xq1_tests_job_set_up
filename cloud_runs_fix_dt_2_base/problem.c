/**
 * Gravitational Collapse of particle cloud in the Shearing sheet (Hill's approximation)
 *
 *This is intended to be the base code for integrating collapse of shearing cloud
 *allows the reloading and restarting of simulations
 *Note: rerunning will not neccessarily produce the same results (due to collisions)
 *
 *For the intial run (iteration 0): choose seed position text file, or select new (to generate new positions)
 *	Provide the number of particles and inflation, rotation and random velocity factors
 *Or load a binary file: iteration counter, iter= 1,2,3...
 *	Select the file to load from (0 will effetcively be a repeat, or load file i where file i is a file from the previous iteration in order to restart)
 *
 *_num_out determines how frequently a .txt and .bin is saved
 */

//general simualtion parameters
double t_0= -1.0; //placeholder for reading initial time from input file
static int pos_i=0; //file naming variable
double _num_out=18046*0.25*0.25*(3.0/2.0);
static double t_max=1e2*365*24*60*60;

//RUN PARAMETERS
char ins[64]="../seed_pos/seed0_norm_100000.txt";

//run parameters, required when generating new run (otherwise loaded from binary)
double Ntot=1e5; //number of particles, needed for new conditions, and setting f!!!
//double fac=0.0; //random velocity factor
double X=0.5;
double _f=1;

//When reloading from binary: iter>0 and give file number to restart from
int iter=0; //interation number of simulation, used for reloading/restarting
int restart=0; //file number to restart from. Will load dat<restart>_<iter>.bin

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <omp.h>
#include <time.h>
#include <sys/time.h>
#include "rebound.h"
#include "../reb_funcs/reb_func.h"
#include <unistd.h>
#include <string.h>
#include <libgen.h>

void heartbeat(struct reb_simulation* const r);

int main(int argc, char* argv[]) {

	//timing:
	struct timeval tim;
	gettimeofday(&tim, NULL);
	double timing1 = tim.tv_sec+(tim.tv_usec/1000000.0);

	// Get the number of processors
	int np = omp_get_num_procs();
	printf("Number of cores: %d \n",np);
	//int np_max_jakita=16,np_max_kelvin=20; //maximum number of processors to use, 16 on jakita, 20 on kelvin
	// int np_max_jakita=16,np_max_kelvin=np; //maximum number of processors to use, 16 on jakita, 20 on kelvin
	// int np_max=find_proc_num(np_max_kelvin,np_max_jakita,	"kelvin"); //test if kelvin is in base name and set np_max accordingly
	// if (np>10){ //use the maximum number of cores available, unless we are in jakita or kelvin
	// 	np=np_max;
	// }
omp_set_num_threads(np);

	// Setup constants and routines
	struct reb_simulation* r;
	r = NULL; //set r to null to avoid compilation errors

  struct run_params * rp;
  rp = malloc(sizeof(*rp)); //assign memory for run_params, OR INITIALISE VARIABLES?

	if (iter==0) { //This is the first iteration, set up a new simulation

		r=reb_create_simulation();
		//set rebound simulation parameters
		r->opening_angle2	= .5;					// This determines the precision of the tree code gravity calculation.
		r->integrator			= REB_INTEGRATOR_SEI;
		r->boundary			= REB_BOUNDARY_OPEN;
		r->gravity			= REB_GRAVITY_TREE;
		r->collision			= REB_COLLISION_TREE;
		r->collision_resolve    = reb_collision_resolve_merge_out;          // Choose merger collision routine.
		r->G 				= G;		// N / kg^2 m^2
		r->usleep = -1.0; //set -ve to disable visualization
		r->heartbeat			= heartbeat;	// function pointer for heartbeat

		double R_c_fac=0.6; //cloud radius factor
		double v_col=30; // expected collision velocity, m/s

		//set_Ntot(rp,ins); //set Ntot as input above, or from number of particles in seed file

		//set run parameter values, as we have found Ntot we should have correct values of f and dt
		rp->Ntot=Ntot;
		rp->a=30*AU;
rp->R_eq=100000.0;
		rp->rho=1e3;
		rp->M_tot=(4.0/3.0)*M_PI*rp->rho*pow(rp->R_eq,3.0);
		rp->OMEGA=sqrt((r->G)*M_sun/pow(rp->a,3));
		rp->R_hill=pow(r->G*rp->M_tot/(3*rp->OMEGA*rp->OMEGA),(1.0/3.0));
		rp->R_c=R_c_fac*rp->R_hill;
		rp->X=X;
		rp->OM_circ=rp->X*sqrt(r->G*rp->M_tot/pow(rp->R_c,3));
		//double _N=1e5; //Nesvorny number of particles
		//rp->f=pow(rp->Ntot/_N,-1.0/6.0)*_f; //f IS SCALED FOR N?!
		rp->f=_f;
rp->dt=rp->f*rp->R_eq/(v_col*pow(rp->Ntot,1.0/3.0))*(2.0/3.0);
		r->dt=rp->dt; //set time step according to particle size and estimated collisional velocity
		r->ri_sei.OMEGA 		= rp->OMEGA; //Set SEI rotation
		r->ri_sei.OMEGAZ=r->ri_sei.OMEGA; //setting OMEGAZ explicitly from start

    //set up simulation box
		double boxsize 			= 2*10.0*rp->R_c;			// m. set the box size
		reb_configure_box(r, boxsize, 1, 1, 1);

		//set initial conditions, either new or load a text file
		char new[12]="new";
		if(!strcmp(ins,new)){ //When new conditions are required, new positions and velocities are calculated
			printf("\n--------------\nnew initial conditions\n--------------\n");
			new_initial_conditions(r,rp);
			}
		else{
			printf("\n--------------\nload %s as initial conditions\n--------------\n",ins);
			file_input_normfix(r,rp,ins); //Positions are loaded from the seed file, and then velocities are calculated
			}

	}

	else if (iter>0) { //Load the previous iteration of the simulation

		//load an entire simulation from binary files
		pos_i=restart; //file number to reload from. 0 to redo, or some other file number to restart from
		char fname[32];
		sprintf(fname, "dat%07d_%d.bin", pos_i,iter-1); // reload the file of index pos_i from the previous run (iter-1)
		printf("iteration %d, load from rebound binary: %s \n",iter,fname);
		r = reb_create_simulation_from_binary(fname);

		//LOADED BINARY, reset function pointers and t_0
		r->heartbeat			= heartbeat;	// function pointer for heartbeat
		r->collision_resolve    = reb_collision_resolve_merge_out; // Choose merger collision routine.
		t_0=r->t;

		//load_run_params(rp,iter); //loads the original run_params, initial Ntot etc FIGURE OUT HOW TO FUNCTIONALISE THIS BIT?
		char ignore_char;
	  double ignore_float;
	  int ignore_int;
	  char param_in[64];
	  sprintf(param_in,"run_params_%d.txt",iter-1);
	  FILE *fp;
	  printf("load %s\n",param_in );
	  fp=fopen(param_in,"r");
	  fscanf(fp,"%s",&ignore_char);
	  fscanf(fp,"%lf",&rp->Ntot); //reset Ntot from run_params
	  fscanf(fp,"%lf",&rp->a);
	  fscanf(fp,"%lf",&rp->R_eq);
	  fscanf(fp,"%lf",&rp->rho);
	  fscanf(fp,"%lf",&rp->M_tot);
	  fscanf(fp,"%lf",&rp->R_c);
	  fscanf(fp,"%lf",&rp->OMEGA);
	  fscanf(fp,"%lf",&rp->X);
	  fscanf(fp,"%lf",&rp->OM_circ);
	  fscanf(fp,"%lf",&rp->f);
	  fscanf(fp,"%lf",&rp->dt);
	  fscanf(fp,"%lf",&ignore_float);
	  fscanf(fp,"%d",&ignore_int);
	  fscanf(fp,"%s",&ignore_char);
	  fclose(fp);

		if (r->dt!=rp->dt) {
			printf("\ntimestep loading error\n");
		}
		printf("Integrate from time %e until %e \n",t_0,t_max);
		printf("Timesteps: %e \n",r->dt);

	}

	else{
		printf("\nLoading error\n");
	}

 	print_constants(r,rp);

	char coll_out[64];
	sprintf(coll_out,"collisions_%d.txt",iter); // collision file is created with restart numbers
	fopen(coll_out,"w"); //wipes the collision file from previous runs

  //Generate output, before integration with run time = 0, in case the simulation is stopped early
  gen_run_params(rp,iter,0.0,np,ins);  //Generate output, before integration with run_time = 0, in case the simulation is stopped early

	//set output interval
	printf("total timesteps = %e\n",t_max/r->dt);
	printf("Save files every %e steps\n", _num_out);
	printf("Save a total of %e timesteps\n",t_max/r->dt/_num_out);
	printf("Sample timescale = %e s, or %e days\n\n",_num_out*r->dt,_num_out*r->dt/(24.0*60*60));
	//----------------------------------------------------------------------------
	srand(0);
	reb_integrate(r, t_max); //integrate up to t_max
	//----------------------------------------------------------------------------

	//finish timing
	gettimeofday(&tim, NULL);
	double timing2 = tim.tv_sec+(tim.tv_usec/1000000.0);
	double run_time=(timing2-timing1);
	printf("\ntiming:%e s\n",run_time);

	//Generate output
  gen_run_params(rp,iter,run_time,np,ins);

}
//------------------------------------------------------------------------------
void heartbeat(struct reb_simulation* const r){
	char buf[32];
	if (reb_output_check(r, _num_out*r->dt)){ //print out every nth timestep
		reb_output_timing(r, t_max); //did this slow everything down???
		//Output the txt file
		sprintf(buf, "dat%07d_%d.txt", pos_i,iter); // puts string into buffer
		reb_output_ascii_data(r,buf); //txt file output: time header then x,y,z,vx,vy,vz,m,r for each particle
		//Output the binary file
		sprintf(buf, "dat%07d_%d.bin", pos_i,iter); // puts string into buffer
		reb_output_binary(r,buf); //rebound function to save whole simulation to binary file?
		pos_i++;
	}
}
