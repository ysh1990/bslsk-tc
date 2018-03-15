/**
# The 2004 Indian Ocean tsunami

The 2004 Indian Ocean tsunami was caused by a large-scale fault
rupture (> 1000 km) at the Indian–Australian and Eurasian–Andaman
plate boundaries. This example uses the fault model of [Grilli et
al, 2007](/src/references.bib#grilli2007) as initial conditions for a 
Saint-Venant solution of the subsequent tsunami. A similar setup is
discussed in [Popinet, 2011](/src/references.bib#popinet2011).

## Solver setup

The following hea ders specify that we use spherical coordinates and
the [Saint-Venant solver](/src/saint-venant.h) together with
[(dynamic) terrain reconstruction](/src/terrain.h) and the [Okada
fault model](/src/okada.h). We will use [inputs](/src/input.h) only
when restarting a simulation from a snapshot. */

//#include "spherical.h"
#include "terrain.h"
#include "input.h"
#include "TCSV-working/sv-zb.h"



/**
We then define a few useful macros and constants. */


#define ADP 1
#define OUTPUT_DAT 0
#define OUTPUT_PNG 0
#define OUTPUT_GFS 1
#define OUTPUT_BIN 0

#define MINLEVEL 5
#define ETAE     1e-1 // error on free surface elevation (1 cm)
#define HMAXE    5e-1 // error on maximum free surface elevation (5 cm)
#define CCE      1e-5
#define HEPS      3e-1
int maxlevel=10;
//#define ZB_MOD_INTERVAL 60

  double Radius = 6371220.;
#define factor_r 3.1415927/180.*Radius

//#define INIT_LR(x,y) (x<-10.5 && x>-11. && y>30.5 && y<31.)
#define INIT_LR(x,y) (x>-11.*factor_r && y<31.*factor_r)
#define INIT_ETA -2200.
#define INIT_CC 0.01

#define SAVE_INTV 0.5
#define INIT_OUTPUT_T 0 // 0 by default
#define TMAX 10


scalar hmax[];

/**
The maximum number of levels to use can be set as an argument to the
program. */

int main ()
{
  
  /**
  Here we setup the domain geometry. We choose to use metre as length
  unit, so we set the radius of the Earth (required for the [spherical
  coordinates](/src/spherical.h)) in metres. The *x* and *y*
  coordinates are longitude and latitude in degrees, so we set the
  size of the box *L0* and the coordinates of the lower-left corner
  *(X0,Y0)* in degrees. */

  // the domain is 54 degrees squared
  size (5.*factor_r);
  // centered on 94,8 longitude,latitude
  //origin (-15. - L0/2., 30. - L0/2.);
  origin (-15.*factor_r, 30.*factor_r);

  /**
  *G* is the acceleration of gravity required by the Saint-Venant
  solver. This is the only dimensional parameter. We rescale it so that
  time is in minutes. */

  // acceleration of gravity in m/min^2
  G = 9.81;//*sq(60.);
  psy = 0.15;
  d_s = 20e-6;
  Czf = 0.004;
  r_w = 0.43;
  p_0 = 0.5;
  rho_w = 1035.;


  /**
  When using a tree (i.e. adaptive) discretisation, we want to start
  with the coarsest grid, otherwise we directly refine to the maximum
  level. Note that *1 << n* is C for $2^n$. */



#if TREE
  // 32^2 grid points to start with
  init_grid (1 << maxlevel);//MINLEVEL);
#else // Cartesian
  // 1024^2 grid points
  init_grid (1 << maxlevel);
#endif

  /**
  We then call the *run()* method of the Saint-Venant solver to
  perform the integration. */

  run();
}


/**
## Boundary conditions

We set the normal velocity component on the left, right and bottom
boundaries to a "radiation condition" with a reference sealevel of
zero. The top boundary is always "dry" in this example so can be left
alone. Note that the sign is important and needs to reflect the
orientation of the boundary. */
/*
u.n[left]   = - radiation(0);
u.n[top]  = + radiation(0);
u.n[bottom] = - radiation(0);
*/
/**
## Adaptation

Here we define an auxilliary function which we will use several times
in what follows. Again we have two *#if...#else* branches selecting
whether the simulation is being run on an (adaptive) tree or a
(static) Cartesian grid.

We want to adapt according to two criteria: an estimate of the error
on the free surface position -- to track the wave in time -- and an
estimate of the error on the maximum wave height *hmax* -- to make
sure that the final maximum wave height field is properly resolved.

We first define a temporary field (in the
[automatic variable](http://en.wikipedia.org/wiki/Automatic_variable)
*η*) which we set to $h+z_b$ but only for "wet" cells. If we used
$h+z_b$ everywhere (i.e. the default $\eta$ provided by the
Saint-Venant solver) we would also refine the dry topography, which is
not useful. */

int adapt() {
#if TREE
  scalar eta[];
  foreach()
    eta[] = h[] > dry ? h[] + zb[] : zb[];
  boundary ({eta});

  /**
  We can now use wavelet adaptation on the list of scalars *{η,hmax}*
  with thresholds *{ETAE,HMAXE}*. The compiler is not clever enough yet
  and needs to be told explicitly that this is a list of *double*s,
  hence the *(double[])*
  [type casting](http://en.wikipedia.org/wiki/Type_conversion). 
  
  The function then returns the number of cells refined. */

  astats s = adapt_wavelet ({hmax, cc, h}, (double[]){HMAXE, CCE, HEPS},
			    maxlevel, MINLEVEL);
  fprintf (ferr, "# refined %d cells, coarsened %d cells\n", s.nf, s.nc);
  return s.nf;
  foreach()
		 if (cc[]<0)
		 {
				 fprintf(stderr, "!!warning!! c[]=%g at x=%g y=%g t=%g after adapt, set to 0\n",cc[], x, y, t);
		 }
#else // Cartesian
  return 0;
#endif
}

/**
## Initial conditions

We first specify the terrain database to use to reconstruct the
topography $z_b$. This KDT database needs to be built beforehand. See the
[*xyz2kdt* manual](http://gfs.sourceforge.net/wiki/index.php/Xyz2kdt)
for explanations on how to do this.

We then consider two cases, either we restart from an existing
snapshot or we start from scratch.

The next line tells the Saint-Venant solver to conserve water surface
elevation rather than volume when adapting the mesh. This is important
for tsunamis since most of the domain will be close to "lake-at-rest"
balance. */

event init (i = 0)
{
  //terrain (zb0, "/pub/home/anyi/drop/basilisk/run/TC/moroccan/terrain/mor", NULL);
  //terrain (zb0, "/home/yangsh/coding/basilisk/run/Moroccan/mor", NULL);
  //terrain (zb0, "/WORK/pp245/ysh/basilisk/run/moroccan/terrain/mor", NULL);
  terrain (zb0, "/home/ysh/basilisk/source/basrun/moroccan-terrain/mor", NULL);

  foreach()
	zb[] = zb0[];

  if (restore (file = "dump"))
    conserve_elevation();
  else {
    conserve_elevation();
    
    /**
    The initial still water surface is at $z=0$ so that the water depth
    $h$ is... */
    
    foreach() {
      h[] = INIT_LR(x,y) ? max(INIT_ETA-zb0[],0.) : 0.;
      cc[] = ((INIT_LR(x,y)) && (h[]>0)) ? INIT_CC : 0.;
    }
    boundary ({h,cc});
  }
  
				FILE * out_fp_zb0 = NULL;
				char name_zb0[10];
				sprintf (name_zb0, "zb0.dat");
				out_fp_zb0 = fopen (name_zb0, "w");
				output_field ({zb,eta,cc,h,u}, out_fp_zb0);
				fclose (out_fp_zb0);
				out_fp_zb0 = NULL;
				fprintf(stderr, "output0\n");

				FILE * out_fp_initpa = fopen ("parameter", "w");
				fprintf (out_fp_initpa, "psy = %g\n Czf = %g\n r_b = %g\n rho_s = %g\n d_s = %g\n r_w = %g\n", psy, Czf, r_b, rho_s, d_s, r_w);
				fclose (out_fp_initpa);
				out_fp_initpa = NULL;
				//fprintf(stderr, "output0\n");

}



event logfile (i++) {
  foreach() {
    if (h[] > dry && h[] + zb[] > hmax[])
      hmax[] = h[] + zb[];
	if (cc[]<0) {
			cc[] = 0;
			h[] = 0.;
	}
  }
  boundary ({hmax, u});
}




event test(t = 0; t += 1.) // note "t = t0" must be added or there will be some unknown mistake.
//event test(t = 0; t += 0.0001) // note "t = t0" must be added or there will be some unknown mistake.
{

		int count =0;
		stats check_ux = statsf(u.x);
		stats check_h = statsf(h);
		stats check_c = statsf(cc);
		stats check_zb = statsf(zb);
		fprintf(stderr, "t = %g, dt = %g, CFL = %g, minux = %g, maxux =  %g, minh = %g, maxh = %g, minc = %g, maxc = %g, minzb = %g, maxzb = %g\n", t, dt, CFL, check_ux.min, check_ux.max, check_h.min, check_h.max, check_c.min, check_c.max, check_zb.min, check_zb.max);
		//assert (abs(check_ux.max)<0.1);


#if OUTPUT_DAT
		
				int i_out=(int)(t);
				for (; (i_out%SAVE_INTV==0); i_out++)
		//if(!(t*(t-20)*(t-50)*(t-150)*(t-300)*(t-600)*(t-1200)*(t-1800))) // note the priosity of "()" and "!"
		//if(abs(t-25.)<0.0001 || abs(t-100.)<0.0001) // note the priosity of "()" and "!"
		//if((t==0) || (t==20) || (t==150) || (t==300) || (t==600) || (t==1200) ||(t-1800)) // note the priosity of "()" and "!"
				if (t>INIT_OUTPUT_T)
		{
				fprintf(stderr, "\n\n output: t=%g\n\n",t);

		foreach()
				if (h[]<dry && u.x[]!=0.)
						fprintf (stderr, "!!wrong\n");

		
				FILE * out_fp = NULL;
				char name[10];
				sprintf (name, "out_%011.2f.dat", t);
				out_fp = fopen (name, "w");
				fprintf(stderr, "output field at t = %g, i = %d, dt = %g\n", t, i, dt);
				
				scalar HU[], HC[];

#if ADP
				fprintf (stderr, "ADP is used!\n");
				output_field ({zb, h, eta, cc, u, Ri}, out_fp, linear = true);
#else
				fprintf (stderr, "ADP is not used!\n");
				foreach() {
						HU[] = h[] * u.x[];
						HC[] = h[] * cc[];
				}
				output_field ({zb, h, eta, cc, HU, HC, u, Ri}, out_fp);
#endif


				fclose (out_fp);
				out_fp = NULL;



				//count ++;
		}
#endif
}



event interface (t +=SAVE_INTV; t<=TMAX) {


		if (t>INIT_OUTPUT_T) {
  char s_gfs[20], s_mat[20];
  sprintf (s_gfs, "out-%011.2f.gfs", t);
  sprintf (s_mat, "out-%011.2f.bin", t);
#if OUTPUT_GFS
  output_gfs (file = s_gfs, t = t, list = {zb,h,eta,cc,u,Ri});
	fprintf (stderr, "\n output gfs at t=%g\n",t);
#endif
#if OUTPUT_BIN
  FILE * fp_mat =fopen (s_mat, "w");
  output_field_bin ({zb,h,eta,cc,u,Ri}, fp_mat, linear = true);
	fprintf (stderr, "\n output bin at t=%g\n",t);
  fclose (fp_mat);
  fp_mat = NULL;
#endif
  //fclose (fp);
  //system("matlab -nodesktop -nosplash -logfile `date +%Y_%m_%d-%H_%M_%S`.log -r readbin2D &");
}
}


#if OUTPUT_PNG
event pictures (t+=SAVE_INTV/2.; t<=TMAX) {


		if (t>INIT_OUTPUT_T) {
  char s1[80], s2[80];
  sprintf (s1, "outh%07.2f.png", t);
  sprintf (s2, "outl%07.2f.png", t);
  FILE * fp4 = fopen (s1, "w");
  FILE * fp5 = fopen (s2, "w");

  scalar l[];
  foreach()
		  l[] = depth();
  output_ppm (h, fp4,  spread = 10, linear = true);
  output_ppm (l, fp5,  spread = 10, linear = true);
  fclose(fp4);
  fclose(fp5);
	fprintf (stderr, "\n output png at t=%g\n",t);
}
//}
#endif

/**
## Adaptivity

And finally we apply our *adapt()* function at every timestep. */

#if ADP
/*
event do_adapt (t={64800,65700}) 
{
		adapt();
  char sadp_gfs[20], sadp_mat[20];
  sprintf (sadp_gfs, "outadp-%011.2f.gfs", t);
  sprintf (sadp_mat, "outadp-%011.2f.bin", t);
  output_gfs (file = sadp_gfs, t = t, list = {zb,h,eta,cc,u,Ri});
	fprintf (stderr, "\n !! adp output gfs at t=%g\n",t);
  FILE * fpadp_mat =fopen (sadp_mat, "w");
  output_field_bin ({zb,h,eta,cc,u,Ri}, fpadp_mat, linear = true);
	fprintf (stderr, "\n !! adp output bin at t=%g\n",t);
  fclose (fpadp_mat);
  fpadp_mat = NULL;
}
*/

event do_adapt(i++)
{
		adapt();
}
#endif


//event endtime (t = 105) {return 1;}
event endtime (t = TMAX) {return 1;}





