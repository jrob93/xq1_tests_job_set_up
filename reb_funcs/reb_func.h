#include "vectmath.h"
//------------------------------------------------------------------------------
//Constants
//------------------------------------------------------------------------------

//Note: should these units be more accurate, can we change them now that simulations have been run?
double G=6.67428e-11;
double M_sun=1.98855e30; //kg
double AU=1.496e11; //m

//------------------------------------------------------------------------------
//initialise variables
double m_min; //define m_min as a global variable
double m_damp_lim;

//------------------------------------------------------------------------------
//vector operations

//------------------------------------------------------------------------------
//parameter structure
struct run_params {
  double Ntot;
  double a;
  double R_eq;
  double rho;
  double M_tot;
  double R_c;
  double OMEGA;
  double X;
  double OM_circ;
  double f;
  double dt;
  double R_hill;
  //cwd,Ntot,a,R_eq,rho,M_tot,R_c,OMEGA,X,OM_circ,f,r->dt,run_time,np,ins);
};

//------------------------------------------------------------------------------
//FUNCTIONS
//------------------------------------------------------------------------------
//Function that points to 3 floats x,y,z and creates the corresponding unit 3d vector
void rand_vec(double* x, double* y, double *z){
  *x = reb_random_normal(1.0);
  *y = reb_random_normal(1.0);
  *z = reb_random_normal(1.0);
  double mag=pow(*x**x+*y**y+*z**z,0.5);
  *x=*x/mag;
  *y=*y/mag;
  *z=*z/mag;
  //printf("%.18e\t%.18e\t%.18e\n",*x,*y,*z);
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

/*Function to output a text file with position, velocity, mass and radius of each particle
Use this to output when doing binary file loading, we do not want to add t_0 to time*/
void reb_output_ascii_data(struct reb_simulation* r, char* filename){
    const int N = r->N;
    FILE* of = fopen(filename,"w"); //writes to file, not append!
    if (of==NULL){
        reb_exit("Can not open file.");
    }
		fprintf(of,"%.18e\n",r->t); //put time at start of file
    for (int i=0;i<N;i++){
        struct reb_particle p = r->particles[i];
        fprintf(of,"%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\t%.18e\n",p.x,p.y,p.z,p.vx,p.vy,p.vz,p.m,p.r);
    }
    fclose(of);
}

//------------------------------------------------------------------------------
//define a new collision function that will output txt file
int reb_collision_resolve_merge_out(struct reb_simulation* const r, struct reb_collision c){
	if (r->particles[c.p1].lastcollision==r->t || r->particles[c.p2].lastcollision==r->t) return 0;
    // Every collision will cause two callbacks (with p1/p2 interchanged).
    // Always remove particle with larger index and merge into lower index particle.
    // This will keep N_active meaningful even after mergers.
    int swap = 0;
    int i = c.p1;
    int j = c.p2;   //want j to be removed particle
    if (j<i){
        swap = 1;
        i = c.p2;
        j = c.p1;
    }
    struct reb_particle* pi = &(r->particles[i]);
    struct reb_particle* pj = &(r->particles[j]);
    double invmass = 1.0/(pi->m + pj->m);
		const double t = r->t;
    char coll_out[64];
    sprintf(coll_out,"collisions_%d.txt",iter); // collision file is created with restart numbers
		FILE* of = fopen(coll_out,"a+");                // open file for collision output

    fprintf(of, "%.18e\t", t);                                 // time
    fprintf(of,"%d\t%d\t",i,j); //particle index
		fprintf(of, "%.18e\t%.18e\t",pi->m,pj->m);  // m1 and m2. need to write m1 and m2 BEFORE the masses are reassigned
		fprintf(of, "%.18e\t%.18e\t%.18e\t",pi->x,pi->y,pi->z);
		fprintf(of, "%.18e\t%.18e\t%.18e\t",pj->x,pj->y,pj->z);
		fprintf(of, "%.18e\t%.18e\t%.18e\t",pi->vx,pi->vy,pi->vz);
		fprintf(of, "%.18e\t%.18e\t%.18e\t",pj->vx,pj->vy,pj->vz);
    /*
    fprintf(of, "%.18e\t", t);                                 // time
    fprintf(of,"%d\t",i); //particle index
    fprintf(of, "%.18e\t",pj->m);  // m1. need to write m1 and m2 BEFORE the masses are reassigned
    fprintf(of, "%.18e\t%.18e\t%.18e\t",pi->x,pi->y,pi->z);
    fprintf(of, "%.18e\t%.18e\t%.18e\t",pi->vx,pi->vy,pi->vz);

    fprintf(of,"%d\t",j); //particle index
    fprintf(of, "%.18e\t",pj->m);  // m2. need to write m1 and m2 BEFORE the masses are reassigned
    fprintf(of, "%.18e\t%.18e\t%.18e\t",pj->x,pj->y,pj->z);
    fprintf(of, "%.18e\t%.18e\t%.18e\t",pj->vx,pj->vy,pj->vz);
    */
		vector v_rel;
		v_rel[0]=pi->vx-pj->vx;
		v_rel[1]=pi->vy-pj->vy;
		v_rel[2]=pi->vz-pj->vz;
		double abs_v_rel;
		ABSV(abs_v_rel,v_rel);
		fprintf(of, "%.18e\t", abs_v_rel);                                 // relative velocity
	//this position is just the average position between the two, should print each particle pos instead!!!
		//fprintf(of, "%.18e\t", (pi->x+pj->x)/2.);  // x position
		//fprintf(of, "%.18e\t", (pi->y+pj->y)/2.);  // y position
		//fprintf(of, "%.18e\t", (pi->z+pj->z)/2.);  // z position
		fprintf(of, "\n");
		fclose(of);

    printf("\nCollision detected. i: %d, j: %d\n",i,j);
		//printf("%.18e\t%.18e\t%.18e\t\n",pi->x,pi->y,pi->z);
		//printf("%.18e\t%.18e\t%.18e\t\n",pj->x,pj->y,pj->z);
		//printf("%.18e\t%.18e\t%.18e\t\n",pi->vx,pi->vy,pi->vz);
		//printf("%.18e\t%.18e\t%.18e\t\n",pj->vx,pj->vy,pj->vz);
		//printf("%.18e\t%.18e\t%.18e\t\n",v_rel[0],v_rel[1],v_rel[2]);

    // Merge by conserving mass, volume and momentum
    pi->vx = (pi->vx*pi->m + pj->vx*pj->m)*invmass;
    pi->vy = (pi->vy*pi->m + pj->vy*pj->m)*invmass;
    pi->vz = (pi->vz*pi->m + pj->vz*pj->m)*invmass;
    pi->x  = (pi->x*pi->m + pj->x*pj->m)*invmass;
    pi->y  = (pi->y*pi->m + pj->y*pj->m)*invmass;
    pi->z  = (pi->z*pi->m + pj->z*pj->m)*invmass;
    pi->m  = pi->m + pj->m;
    pi->r  = pow(pow(pi->r,3.)+pow(pj->r,3.),1./3.);
    pi->lastcollision = r->t;

    // If hermes calculate energy offset in global - hasn't been removed from global yet
    if (r->ri_hermes.global){
        if(r->ri_hermes.global->ri_hermes.mini_active){
            r->ri_hermes.global->ri_hermes.collision_this_global_dt = 1;
        }
    }

    return swap?1:2; // Remove particle p2 from simulation
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//define a new collision function that will output txt file with momentum dampening
int reb_collision_resolve_merge_out_damp(struct reb_simulation* const r, struct reb_collision c){
	if (r->particles[c.p1].lastcollision==r->t || r->particles[c.p2].lastcollision==r->t) return 0;
    // Every collision will cause two callbacks (with p1/p2 interchanged).
    // Always remove particle with larger index and merge into lower index particle.
    // This will keep N_active meaningful even after mergers.
    int swap = 0;
    int i = c.p1;
    int j = c.p2;   //want j to be removed particle
    if (j<i){
        swap = 1;
        i = c.p2;
        j = c.p1;
    }
    struct reb_particle* pi = &(r->particles[i]);
    struct reb_particle* pj = &(r->particles[j]);
    double invmass = 1.0/(pi->m + pj->m);

    double N_boulder = pj->m/m_min;
    printf("m1 = %e\t m2 = %e\t m_min = %e\n",pi->m,pi->m,m_min);
    printf("Number of boulders = %e\n",N_boulder);
    double inv_sqrt_N = 1.0/sqrt(N_boulder);

		const double t = r->t;
    char coll_out[64];
    sprintf(coll_out,"collisions_%d.txt",iter); // collision file is created with restart numbers
		FILE* of = fopen(coll_out,"a+");                // open file for collision output

    fprintf(of, "%.18e\t", t);                                 // time
    fprintf(of,"%d\t%d\t",i,j); //particle index
		fprintf(of, "%.18e\t%.18e\t",pi->m,pj->m);  // m1 and m2. need to write m1 and m2 BEFORE the masses are reassigned
		fprintf(of, "%.18e\t%.18e\t%.18e\t",pi->x,pi->y,pi->z);
		fprintf(of, "%.18e\t%.18e\t%.18e\t",pj->x,pj->y,pj->z);
		fprintf(of, "%.18e\t%.18e\t%.18e\t",pi->vx,pi->vy,pi->vz);
		fprintf(of, "%.18e\t%.18e\t%.18e\t",pj->vx,pj->vy,pj->vz);
    /*
    fprintf(of, "%.18e\t", t);                                 // time
    fprintf(of,"%d\t",i); //particle index
    fprintf(of, "%.18e\t",pj->m);  // m1. need to write m1 and m2 BEFORE the masses are reassigned
    fprintf(of, "%.18e\t%.18e\t%.18e\t",pi->x,pi->y,pi->z);
    fprintf(of, "%.18e\t%.18e\t%.18e\t",pi->vx,pi->vy,pi->vz);

    fprintf(of,"%d\t",j); //particle index
    fprintf(of, "%.18e\t",pj->m);  // m2. need to write m1 and m2 BEFORE the masses are reassigned
    fprintf(of, "%.18e\t%.18e\t%.18e\t",pj->x,pj->y,pj->z);
    fprintf(of, "%.18e\t%.18e\t%.18e\t",pj->vx,pj->vy,pj->vz);
    */
		vector v_rel;
		v_rel[0]=pi->vx-pj->vx;
		v_rel[1]=pi->vy-pj->vy;
		v_rel[2]=pi->vz-pj->vz;
		double abs_v_rel;
		ABSV(abs_v_rel,v_rel);
		fprintf(of, "%.18e\t", abs_v_rel);                                 // relative velocity
	//this position is just the average position between the two, should print each particle pos instead!!!
		//fprintf(of, "%.18e\t", (pi->x+pj->x)/2.);  // x position
		//fprintf(of, "%.18e\t", (pi->y+pj->y)/2.);  // y position
		//fprintf(of, "%.18e\t", (pi->z+pj->z)/2.);  // z position
		fprintf(of, "\n");
		fclose(of);

    printf("\nCollision detected. i: %d, j: %d\n",i,j);
		//printf("%.18e\t%.18e\t%.18e\t\n",pi->x,pi->y,pi->z);
		//printf("%.18e\t%.18e\t%.18e\t\n",pj->x,pj->y,pj->z);
		//printf("%.18e\t%.18e\t%.18e\t\n",pi->vx,pi->vy,pi->vz);
		//printf("%.18e\t%.18e\t%.18e\t\n",pj->vx,pj->vy,pj->vz);
		//printf("%.18e\t%.18e\t%.18e\t\n",v_rel[0],v_rel[1],v_rel[2]);

    printf("Damped collision\n");

    // Merge by conserving mass and volume, momentum is damped
    pi->vx = (pi->vx*pi->m + pj->vx*pj->m*inv_sqrt_N)*invmass;
    pi->vy = (pi->vy*pi->m + pj->vy*pj->m*inv_sqrt_N)*invmass;
    pi->vz = (pi->vz*pi->m + pj->vz*pj->m*inv_sqrt_N)*invmass;
    pi->x  = (pi->x*pi->m + pj->x*pj->m)*invmass;
    pi->y  = (pi->y*pi->m + pj->y*pj->m)*invmass;
    pi->z  = (pi->z*pi->m + pj->z*pj->m)*invmass;
    pi->m  = pi->m + pj->m;
    pi->r  = pow(pow(pi->r,3.)+pow(pj->r,3.),1./3.);
    pi->lastcollision = r->t;

    // If hermes calculate energy offset in global - hasn't been removed from global yet
    if (r->ri_hermes.global){
        if(r->ri_hermes.global->ri_hermes.mini_active){
            r->ri_hermes.global->ri_hermes.collision_this_global_dt = 1;
        }
    }

    return swap?1:2; // Remove particle p2 from simulation
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//define a new collision function that will output txt file with momentum dampening, for only particle masses below a certain limit
int reb_collision_resolve_merge_out_damp_cut(struct reb_simulation* const r, struct reb_collision c){
	if (r->particles[c.p1].lastcollision==r->t || r->particles[c.p2].lastcollision==r->t) return 0;
    // Every collision will cause two callbacks (with p1/p2 interchanged).
    // Always remove particle with larger index and merge into lower index particle.
    // This will keep N_active meaningful even after mergers.
    int swap = 0;
    int i = c.p1;
    int j = c.p2;   //want j to be removed particle
  //  printf("\np1 = %i, p2 = %i\n",i,j);

    if (j<i){
      //printf("swap\n");
        swap = 1;
        i = c.p2;
        j = c.p1;
        //printf("p1 = %i, p2 = %i\n",j,i);

    }
    struct reb_particle* pi = &(r->particles[i]);
    struct reb_particle* pj = &(r->particles[j]);

    //printf("pi = %i, pj = %i\n",i,j);
    //printf("mi = %e\t mj = %e\n",pi->m,pj->m);

    int mass_swap=0;
    if (pj->m > pi->m) {
      //printf("\n*m2 is greater than m1!*\nswap values\n");
      pi = &(r->particles[j]);
      pj = &(r->particles[i]);
      //printf("pi = %i, pj = %i\n",j,i);
      //printf("mi = %e\t mj = %e\n",pi->m,pj->m);
      mass_swap=1;
    }
    //printf("mass_swap = %i\n", mass_swap);


    double invmass = 1.0/(pi->m + pj->m);
    double N_boulder = 1.0;

		const double t = r->t;
    char coll_out[64];
    sprintf(coll_out,"collisions_%d.txt",iter); // collision file is created with restart numbers
		FILE* of = fopen(coll_out,"a+");                // open file for collision output

    fprintf(of, "%.18e\t", t);                                 // time
    fprintf(of,"%d\t%d\t",i,j); //particle index
		fprintf(of, "%.18e\t%.18e\t",pi->m,pj->m);  // m1 and m2. need to write m1 and m2 BEFORE the masses are reassigned
    fprintf(of, "%.18e\t%.18e\t",pi->r,pj->r);
		fprintf(of, "%.18e\t%.18e\t%.18e\t",pi->x,pi->y,pi->z);
		fprintf(of, "%.18e\t%.18e\t%.18e\t",pj->x,pj->y,pj->z);
		fprintf(of, "%.18e\t%.18e\t%.18e\t",pi->vx,pi->vy,pi->vz);
		fprintf(of, "%.18e\t%.18e\t%.18e\t",pj->vx,pj->vy,pj->vz);

		/*vector v_rel;
		v_rel[0]=pi->vx-pj->vx;
		v_rel[1]=pi->vy-pj->vy;
		v_rel[2]=pi->vz-pj->vz;
		double abs_v_rel;
		ABSV(abs_v_rel,v_rel);
		fprintf(of, "%.18e\t", abs_v_rel); */                               // relative velocity

    //printf("Collision detected. i: %d, j: %d\n",i,j);
    //printf("mi = %e\t mj = %e m_damp_lim = %e\n",pi->m,pj->m,m_damp_lim);

    if (pj->m <= m_damp_lim) {
      //do a momentum damped collision
      //printf("Damped collision\n");
      printf("Collision detected. i: %d, j: %d, damped\n",i,j);

      N_boulder = pj->m/m_min;
      printf("m_min = %e \tNumber of boulders = %e\n",m_min,N_boulder);
      double inv_sqrt_N = 1.0/sqrt(N_boulder);

      // Merge by conserving mass and volume, momentum is damped
      pi->vx = (pi->vx*pi->m + pj->vx*pj->m*inv_sqrt_N)*invmass;
      pi->vy = (pi->vy*pi->m + pj->vy*pj->m*inv_sqrt_N)*invmass;
      pi->vz = (pi->vz*pi->m + pj->vz*pj->m*inv_sqrt_N)*invmass;
      pi->x  = (pi->x*pi->m + pj->x*pj->m)*invmass;
      pi->y  = (pi->y*pi->m + pj->y*pj->m)*invmass;
      pi->z  = (pi->z*pi->m + pj->z*pj->m)*invmass;
      pi->m  = pi->m + pj->m;
      pi->r  = pow(pow(pi->r,3.)+pow(pj->r,3.),1./3.);
      pi->lastcollision = r->t;

    }
    else {
      //printf("inelastic merger collision\n");
      printf("Collision detected. i: %d, j: %d, inelastic\n",i,j);

      //normal collsion routine
      // Merge by conserving mass, volume and momentum
      pi->vx = (pi->vx*pi->m + pj->vx*pj->m)*invmass;
      pi->vy = (pi->vy*pi->m + pj->vy*pj->m)*invmass;
      pi->vz = (pi->vz*pi->m + pj->vz*pj->m)*invmass;
      pi->x  = (pi->x*pi->m + pj->x*pj->m)*invmass;
      pi->y  = (pi->y*pi->m + pj->y*pj->m)*invmass;
      pi->z  = (pi->z*pi->m + pj->z*pj->m)*invmass;
      pi->m  = pi->m + pj->m;
      pi->r  = pow(pow(pi->r,3.)+pow(pj->r,3.),1./3.);
      pi->lastcollision = r->t;

    }

    // If hermes calculate energy offset in global - hasn't been removed from global yet
    if (r->ri_hermes.global){
        if(r->ri_hermes.global->ri_hermes.mini_active){
            r->ri_hermes.global->ri_hermes.collision_this_global_dt = 1;
        }
    }

    fprintf(of, "%.18e\t",pi->m);
    fprintf(of, "%.18e\t",pi->r);
		fprintf(of, "%.18e\t%.18e\t%.18e\t",pi->x,pi->y,pi->z);
		fprintf(of, "%.18e\t%.18e\t%.18e\t",pi->vx,pi->vy,pi->vz);

    fprintf(of, "%.18e\t", N_boulder);                                 // number of boulders
    fprintf(of, "\n");
		fclose(of);

    if (mass_swap==1) {
      //the particles have been swapped. now store pi (which will be removed) in pj
      //printf("store pi in pj\n");
      pj->vx = pi->vx;
      pj->vy = pi->vy;
      pj->vz = pi->vz;
      pj->x  = pi->x;
      pj->y  = pi->y;
      pj->z  = pi->z;
      pj->m  = pi->m;
      pj->r  = pi->r;
      pj->lastcollision = pi->lastcollision;

    }
    //printf("mi = %e\t mj = %e\n",pi->m,pj->m);
    //printf("remove particle p%i\n\n",swap?1:2);

    return swap?1:2; // Remove particle p2 from simulation
}
//------------------------------------------------------------------------------
void file_input(struct reb_simulation* const r, struct run_params * rp, char *fname){
  //Function to input positions only from file. Recalculates velocities
  //Now finds total mass
  //printf("\nload initial conditons with random velocities of mag %e*v_circ\n",fac);
  printf("\nload initial conditons with no random velocity \n");
	double mass=0.0;
	//Input from file
	FILE *fp;
	printf("read file: %s\n",fname);
	fp=fopen(fname,"r");
	//Read number of lines
	int ch=0;
	int lines=0;
		while(!feof(fp))
	{
		ch = fgetc(fp);
		if(ch == '\n')
		{
			lines++;
		}
	}
	rewind(fp); //reset the file pointer
	printf("%s %d\n","number of lines: ",lines);
  if (lines-1!=rp->Ntot){
    printf("\nerror in particle number?!!\n");
  }
  //rp->Ntot=lines-1; //reset number of lines
	//read in data
	double x=0.0,y=0.0,z=0.0,vx=0.0,vy=0.0,vz=0.0,_m=0.0,_rad=0.0; //initialise variables, including dummy variables for mass and radius
	fscanf(fp,"%lf",&t_0); //read in the first line
	printf("%s %e\n","t_0: ",t_0 );

  //set particle mass and radius
  double m = rp->M_tot/rp->Ntot;
	double rad 	= pow((3*m/(4*M_PI*rp->rho)),(1.0/3.0));
	rad=rp->f*rad; //inflate the radius by factor f
	r->softening 			= rad;			// m set the softening to be equivalent to particle size

	for(int j=0;j<(lines-1);j++){ //read all subsequent lines
		struct reb_particle pt;
		fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf",&x,&y,&z,&vx,&vy,&vz,&_m,&_rad); //read in the text file line by line
		pt.z = z;
		pt.x = x;
		pt.y = y;
		pt.m 		= m; 	// kg
		mass+=m;
		pt.r 		= rad;						// m. Note that the radius may aleady be inflated by factor f
    vector pos={pt.x,pt.y,pt.z};
		vector vel;
		//shear velocity
		vector v_s;
		v_s[0] 		= 0;
		v_s[1] 		= -1.5*pt.x*rp->OMEGA;
		v_s[2] 		= 0;
		//uniform swarm rotation (Nesvorny)
		double om_circ[3]={0,0,rp->OM_circ};
		vector v_r;
		CROSSVP(v_r,om_circ,pos); //circular velocity at particle radius
		//initial velocity is sum of shear and circular
		ADDV(vel,v_s,v_r);
    /*
    //Find the magnitude of the random velocity: 0<vrand<vcirc
    double vcirc=pow(v_r[0]*v_r[0]+v_r[1]*v_r[1]+v_r[2]*v_r[2],0.5);
    double vrand=reb_random_uniform(0,fac*vcirc);
    //Calculate a random velocity vector
    double v_ranx,v_rany,v_ranz;
    rand_vec(&v_ranx,&v_rany,&v_ranz);
    vector v_ran;
    v_ran[0]=(v_ranx)*vrand;
    v_ran[1]=(v_rany)*vrand;
    v_ran[2]=(v_ranz)*vrand;
    //printf("vx=%e vy=%e vz=%e mag=%e \n",v_ran[0],v_ran[1],v_ran[2],pow(v_ran[0]*v_ran[0]+v_ran[1]*v_ran[1]+v_ran[2]*v_ran[2],0.5));
    //printf("velx=%e vely=%e velz=%e \n",vel[0],vel[1],vel[2]);
    ADDV(vel,vel,v_ran);
    */
    //printf("velx=%e vely=%e velz=%e \n",vel[0],vel[1],vel[2]);
		pt.vx = vel[0];
		pt.vy = vel[1];
		pt.vz = vel[2];
		reb_add(r, pt);
	}
	printf("total mass: %.18e\n",mass);
	fclose(fp);
	//r->softening 			= rad;			// m set the softening to be equivalent to particle size (UNINFLATED?)
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void file_input_normfix(struct reb_simulation* const r, struct run_params * rp, char *fname){
  //Function to input positions only from file. Recalculates velocities
  //Now finds total mass
  //RESCALES SEED POSITIONS!

  printf("\nload initial conditons with no random velocity \n");
	double mass=0.0;
	//Input from file
	FILE *fp;
	printf("read file: %s\n",fname);
	fp=fopen(fname,"r");
	//Read number of lines
	int ch=0;
	int lines=0;
		while(!feof(fp))
	{
		ch = fgetc(fp);
		if(ch == '\n')
		{
			lines++;
		}
	}
	rewind(fp); //reset the file pointer
	printf("%s %d\n","number of lines: ",lines);
  if (lines-1!=rp->Ntot){
    printf("\nerror in particle number?!! reset Ntot\n");
    rp->Ntot=lines-1;
  }
  //rp->Ntot=lines-1; //reset number of lines
	//read in data
	double x=0.0,y=0.0,z=0.0,vx=0.0,vy=0.0,vz=0.0,_m=0.0,_rad=0.0; //initialise variables, including dummy variables for mass and radius
	fscanf(fp,"%lf",&t_0); //read in the first line
	printf("%s %e\n","t_0: ",t_0 );

  //set particle mass and radius
  double m = rp->M_tot/rp->Ntot;
	double rad 	= pow((3*m/(4*M_PI*rp->rho)),(1.0/3.0));
	rad=rp->f*rad; //inflate the radius by factor f
	r->softening 			= rad;			// m set the softening to be equivalent to particle size

	for(int j=0;j<(lines-1);j++){ //read all subsequent lines
		struct reb_particle pt;
		fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf",&x,&y,&z,&vx,&vy,&vz,&_m,&_rad); //read in the text file line by line
		pt.z = z*rp->R_c;
		pt.x = x*rp->R_c;
		pt.y = y*rp->R_c;
		pt.m 		= m; 	// kg
		mass+=m;
		pt.r 		= rad;						// m. Note that the radius may aleady be inflated by factor f
    vector pos={pt.x,pt.y,pt.z};
		vector vel;
		//shear velocity
		vector v_s;
		v_s[0] 		= 0;
		v_s[1] 		= -1.5*pt.x*rp->OMEGA;
		v_s[2] 		= 0;
		//uniform swarm rotation (Nesvorny)
		double om_circ[3]={0,0,rp->OM_circ};
		vector v_r;
		CROSSVP(v_r,om_circ,pos); //circular velocity at particle radius
		//initial velocity is sum of shear and circular
		ADDV(vel,v_s,v_r);
    //printf("velx=%e vely=%e velz=%e \n",vel[0],vel[1],vel[2]);
		pt.vx = vel[0];
		pt.vy = vel[1];
		pt.vz = vel[2];
		reb_add(r, pt);
	}
	printf("total mass: %.18e\n",mass);
	fclose(fp);
	//r->softening 			= rad;			// m set the softening to be equivalent to particle size (UNINFLATED?)
}

//------------------------------------------------------------------------------
void file_input_normfix2(struct reb_simulation* const r, struct run_params * rp, char *fname){
  //Function to input positions only from file. Recalculates velocities
  //Now finds total mass
  //RESCALES SEED POSITIONS!

  printf("\nload initial conditons with no random velocity \n");
	double mass=0.0;
	//Input from file
	FILE *fp;
	printf("read file: %s\n",fname);
	fp=fopen(fname,"r");
	//Read number of lines
	int ch=0;
	int lines=0;
		while(!feof(fp))
	{
		ch = fgetc(fp);
		if(ch == '\n')
		{
			lines++;
		}
	}
	rewind(fp); //reset the file pointer
	printf("%s %d\n","number of lines: ",lines);
  if (lines-1!=rp->Ntot){
    printf("\nerror in particle number?!! reset Ntot\n");
    rp->Ntot=lines-1;
  }
  //rp->Ntot=lines-1; //reset number of lines
	//read in data
	double x=0.0,y=0.0,z=0.0,vx=0.0,vy=0.0,vz=0.0,_m=0.0,_rad=0.0; //initialise variables, including dummy variables for mass and radius
	fscanf(fp,"%lf",&t_0); //read in the first line
	printf("%s %e\n","t_0: ",t_0 );

  //set particle mass and radius
  double m = rp->M_tot/rp->Ntot;
	double rad 	= pow((3*m/(4*M_PI*rp->rho)),(1.0/3.0));
	rad=rp->f*rad; //inflate the radius by factor f
	r->softening 			= rad;			// m set the softening to be equivalent to particle size

	for(int j=0;j<(lines-1);j++){ //read all subsequent lines
		struct reb_particle pt;
		fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf",&x,&y,&z,&vx,&vy,&vz,&_m,&_rad); //read in the text file line by line
		pt.z = z*rp->R_c;
		pt.x = x*rp->R_c;
		pt.y = y*rp->R_c;
		pt.m 		= m; 	// kg
		mass+=m;
		pt.r 		= rad;						// m. Note that the radius may aleady be inflated by factor f
    vector pos={pt.x,pt.y,pt.z};
		vector vel;
		//shear velocity
		vector v_s;
		v_s[0] 		= 0;
		v_s[1] 		= 1.5*pt.x*rp->OMEGA;
		v_s[2] 		= 0;
    printf("add shear\n");
		//uniform swarm rotation (Nesvorny)
		double om_circ[3]={0,0,rp->OM_circ};
		vector v_r;
		CROSSVP(v_r,om_circ,pos); //circular velocity at particle radius
		//initial velocity is sum of shear and circular
		ADDV(vel,v_s,v_r);
    //printf("velx=%e vely=%e velz=%e \n",vel[0],vel[1],vel[2]);
		pt.vx = vel[0];
		pt.vy = vel[1];
		pt.vz = vel[2];
		reb_add(r, pt);
	}
	printf("total mass: %.18e\n",mass);
	fclose(fp);
	//r->softening 			= rad;			// m set the softening to be equivalent to particle size (UNINFLATED?)
}

//------------------------------------------------------------------------------
void file_input_normfix3(struct reb_simulation* const r, struct run_params * rp, char *fname){
  //Function to input positions only from file. Recalculates velocities
  //Now finds total mass
  //RESCALES SEED POSITIONS!

  printf("\nload initial conditons with no random velocity \n");
	double mass=0.0;
	//Input from file
	FILE *fp;
	printf("read file: %s\n",fname);
	fp=fopen(fname,"r");
	//Read number of lines
	int ch=0;
	int lines=0;
		while(!feof(fp))
	{
		ch = fgetc(fp);
		if(ch == '\n')
		{
			lines++;
		}
	}
	rewind(fp); //reset the file pointer
	printf("%s %d\n","number of lines: ",lines);
  if (lines-1!=rp->Ntot){
    printf("\nerror in particle number?!! reset Ntot\n");
    rp->Ntot=lines-1;
  }
  //rp->Ntot=lines-1; //reset number of lines
	//read in data
	double x=0.0,y=0.0,z=0.0,vx=0.0,vy=0.0,vz=0.0,_m=0.0,_rad=0.0; //initialise variables, including dummy variables for mass and radius
	fscanf(fp,"%lf",&t_0); //read in the first line
	printf("%s %e\n","t_0: ",t_0 );

  //set particle mass and radius
  double m = rp->M_tot/rp->Ntot;
	double rad 	= pow((3*m/(4*M_PI*rp->rho)),(1.0/3.0));
	rad=rp->f*rad; //inflate the radius by factor f
	r->softening 			= rad;			// m set the softening to be equivalent to particle size

	for(int j=0;j<(lines-1);j++){ //read all subsequent lines
		struct reb_particle pt;
		fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf",&x,&y,&z,&vx,&vy,&vz,&_m,&_rad); //read in the text file line by line
		pt.z = z*rp->R_c;
		pt.x = x*rp->R_c;
		pt.y = y*rp->R_c;
		pt.m 		= m; 	// kg
		mass+=m;
		pt.r 		= rad;						// m. Note that the radius may aleady be inflated by factor f
    vector pos={pt.x,pt.y,pt.z};
		vector vel;
		//shear velocity
		vector v_s;
		v_s[0] 		= 0;
		v_s[1] 		= 0;
		v_s[2] 		= 0;
    printf("add shear\n");
		//uniform swarm rotation (Nesvorny)
		double om_circ[3]={0,0,rp->OM_circ};
		vector v_r;
		CROSSVP(v_r,om_circ,pos); //circular velocity at particle radius
		//initial velocity is sum of shear and circular
		ADDV(vel,v_s,v_r);
    //printf("velx=%e vely=%e velz=%e \n",vel[0],vel[1],vel[2]);
		pt.vx = vel[0];
		pt.vy = vel[1];
		pt.vz = vel[2];
		reb_add(r, pt);
	}
	printf("total mass: %.18e\n",mass);
	fclose(fp);
	//r->softening 			= rad;			// m set the softening to be equivalent to particle size (UNINFLATED?)
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

void new_initial_conditions(struct reb_simulation* const r, struct run_params * rp){
  //Find new initial conditions
  //printf("\ngenerate new initial conditons with random velocities of mag %e*v_circ\n",fac);
  printf("\ngenerate new initial conditons with no random velocity\n");

	t_0=0.0; //ensure the simualtion starts recording at t_0
	double mass = rp->M_tot/rp->Ntot;
	double radius 	= pow((3*mass/(4*M_PI*rp->rho)),(1.0/3.0));
	radius=rp->f*radius; //inflate the radius by factor f
	r->softening 			= radius;			// m set the softening to be equivalent to particle size
	for(int j=0;j<Ntot;j++){
		struct reb_particle pt;

		//evenly distribute the particles in the rotating frame of radius r2
		double x1 = reb_random_uniform(0,1);
		double x2 = reb_random_uniform(0,1);
		double x3 = reb_random_uniform(0,1);
		double Rc=pow(x1,(1.0/3.0))*rp->R_c;
		pt.z = Rc*(1-(2*x2));
		pt.x = sqrt((Rc*Rc)-(pt.z*pt.z))*cos(2*M_PI*x3);
		pt.y = sqrt((Rc*Rc)-(pt.z*pt.z))*sin(2*M_PI*x3);
		vector pos={pt.x,pt.y,pt.z};
		vector vel;

		//shear velocity
		vector v_s;
		v_s[0] 		= 0;
		v_s[1] 		= -1.5*pt.x*rp->OMEGA;
		v_s[2] 		= 0;

		//uniform swarm rotation (Nesvorny)
		double om_circ[3]={0,0,rp->OM_circ};
		vector v_r;
		CROSSVP(v_r,om_circ,pos); //circular velocity at particle radius

		//initial velocity is sum of shear and circular
		ADDV(vel,v_s,v_r);

    /*
    //Calculate a random velocity vector, same as in file_input_randvel above
    double vcirc=pow(v_r[0]*v_r[0]+v_r[1]*v_r[1]+v_r[2]*v_r[2],0.5);
    double vrand=reb_random_uniform(0,fac*vcirc);//Find the magnitude of the random velocity: 0<vrand<vcirc
    double v_ranx,v_rany,v_ranz;
    rand_vec(&v_ranx,&v_rany,&v_ranz);
    vector v_ran;
    v_ran[0]=(v_ranx)*vrand;
    v_ran[1]=(v_rany)*vrand;
    v_ran[2]=(v_ranz)*vrand;
    ADDV(vel,vel,v_ran);
    */
		pt.vx = vel[0];
		pt.vy = vel[1];
		pt.vz = vel[2];

		//particle mass and size
		pt.m 		= mass; 	// kg
		pt.r 		= radius;						// m. Inflation factor f

    /*//print out values
    double mag=0.0;
    printf("%s%e\n","OM: ",OMEGA );
    printf("%s%e\n","om: ",OM_circ );
    ABSV(mag,om_circ);
    printf("%e\t%e\t%e\t%e\t\n",om_circ[0],om_circ[1],om_circ[2],mag);
    ABSV(mag,pos);
    printf("%e\t%e\t%e\t%e\t\n",pos[0],pos[1],pos[2],mag);
    printf("\n");
    ABSV(mag,v_s);
    printf("%e\t%e\t%e\t%e\t\n",v_s[0],v_s[1],v_s[2],mag);
    ABSV(mag,v_r);
    printf("%e\t%e\t%e\t%e\t\n",v_r[0],v_r[1],v_r[2],mag);
    ABSV(mag,vel);
    printf("%e\t%e\t%e\t%e\t\n",vel[0],vel[1],vel[2],mag);
    printf("\n");*/


		reb_add(r, pt);
	}
}
//------------------------------------------------------------------------------

#include <string.h>
#include <libgen.h>


int find_proc_num(int np_1, int np_2,	char *word){
  //function to set np_max to np_1 or np_2, depending on whether or not the base name contains word
  int np_max;
  char cwd1[100000];
	getcwd(cwd1, sizeof(cwd1)); //get current working directory
	char *base;
	base= basename(cwd1); //get the base of the file path
	if(strstr(base, word) != NULL) { //test if we are running in kelvin, according to base name
			printf("kelvin detected in base name\n");
	    np_max=np_1;
	}
	else{ //otherwise we are in jakita (or mac)
    printf("Running on jakita (or mac)\n");
		np_max=np_2;
	}
	printf("use a max number of cores: %d\n",np_max);
  return np_max;
}
//------------------------------------------------------------------------------

void gen_run_params(struct run_params* const rp,int iter,float run_time,int np,char* ins) {
  //Generate output
	char param_out[64];
	sprintf(param_out,"run_params_%d.txt",iter);
	FILE* out1 = fopen(param_out,"w"); //writes to file, not append!
	char *getcwd(char *buf, size_t size);
	char cwd[1024];
	if (getcwd(cwd, sizeof(cwd)) != NULL)
		fprintf(out1, "%s\n%.18e\n%.18e\n%.18e\n%.18e\n%.18e\n%.18e\n%.18e\n%.18e\n%.18e\n%.18e\n%.18e\n%.18e\n%d\n%s\n",
		cwd,rp->Ntot,rp->a,rp->R_eq,rp->rho,rp->M_tot,rp->R_c,rp->OMEGA,rp->X,rp->OM_circ,rp->f,rp->dt,run_time,np,ins);
	fclose(out1);
}
//------------------------------------------------------------------------------
void load_run_params(struct run_params* const rp,int iter){
  //Set constants (load run_params of previous iteration)
  char ignore_char;
  double ignore_float;
  int ignore_int;
  char param_in[64];
  sprintf(param_in,"run_params_%d.txt",iter-1);
  FILE *fp;
  printf("load %s\n",param_in );
  fp=fopen(param_in,"r");
  /*
  fscanf(fp,"%s",&ignore_char);
  fscanf(fp,"%lf",&Ntot); //reset Ntot from run_params
  fscanf(fp,"%lf",&a);
  fscanf(fp,"%lf",&R_eq);
  fscanf(fp,"%lf",&rho);
  fscanf(fp,"%lf",&M_tot);
  fscanf(fp,"%lf",&R_c);
  fscanf(fp,"%lf",&OMEGA);
  fscanf(fp,"%lf",&X);
  fscanf(fp,"%lf",&OM_circ);
  fscanf(fp,"%lf",&f);
  fscanf(fp,"%lf",&ignore_float);
  fscanf(fp,"%lf",&ignore_float);
  fscanf(fp,"%d",&ignore_int);
  fscanf(fp,"%s",&ignore_char);
*/
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
  //R_hill=pow(G*M_tot/(3*OMEGA*OMEGA),(1.0/3.0));
  printf("%s %lf %d\n",&ignore_char,ignore_float,ignore_int);
  printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",rp->Ntot,rp->a,rp->R_eq,rp->rho,rp->M_tot,rp->R_c,rp->OMEGA,rp->X,rp->OM_circ,rp->f,rp->dt);

}
//------------------------------------------------------------------------------
void print_constants(struct reb_simulation* const r,struct run_params* const rp){
  //Function to output the parameters used by r and rp
  double R_hill=pow(G*rp->M_tot/(3*rp->OMEGA*rp->OMEGA),(1.0/3.0));
  //Print out constants
  printf("\nConstants:\n");
  printf("G: %e m^3 kg^-1 s^-2\n",r->G);
  printf("Msun: %e kg\n",M_sun);
  printf("AU: %e m\n",AU);

  printf("\nCloud Properties:\n");
  printf("R_eq: %e m\n",rp->R_eq);
  printf("rho: %e kg m^-3\n",rp->rho);
  printf("M_tot: %e kg\n",rp->M_tot );
  printf("a: %e m\n",rp->a);
  printf("R_hill: %e m\n",R_hill );
  printf("R_c: %e m\n",rp->R_c);

  printf("\nCloud Rotation:\n");
  printf("X: %e\n",rp->X);
  printf("Omega: %e s^-1\n",rp->OMEGA );
  printf("OM_circ: %e s^-1\n",rp->OM_circ );

  struct reb_particle p = r->particles[0];
  printf("\nParticle Properties:\n");
  printf("Ntot: %e\n",rp->Ntot);
  printf("f: %e\n",rp->f);
  printf("particle mass: %e kg\n",p.m );
  printf("particle radius (simulation): %e m\n",p.r  );
  printf("dt:  %.24e s\n",r->dt);
  printf("t_max:  %e s\n\n",t_max);

}
//------------------------------------------------------------------------------
/*
void set_Ntot(struct run_params* const rp,char * ins) {
  //Function to set initial number of particles.
  //Either Ntot above, or the number of particles in the file being loaded (do this here so Ntot is correct for calculating f etc)
  char new[12]="new";
  if(!strcmp(ins,new)){ //When new conditions are required, new positions and velocities are calculated
    //printf("\n--------------\nnew initial conditions\n--------------\n");
    rp->Ntot=Ntot;
    }
  else{
    //printf("\n--------------\nload %s as initial conditions\n--------------\n",ins);
    //Input from file
    FILE *fp;
    fp=fopen(ins,"r");
    //Read number of lines
    int ch=0;
    int lines=0;
      while(!feof(fp))
    {
      ch = fgetc(fp);
      if(ch == '\n')
      {
        lines++;
      }
    }
    rewind(fp); //reset the file pointer
    //printf("%s %d\n","number of lines: ",lines);
    rp->Ntot=lines-1; //reset number of lines
    fclose(fp);
    }
}
*/
//------------------------------------------------------------------------------
