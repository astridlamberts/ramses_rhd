#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"



int main(int argc,char *argv[])
{
  int i;
  char filename[100];
  FILE *fd;

  /*******************************************/

  CC=     8.6;       /* halo concentration      */
  V200=   172.0;      /* circular velocity v_200 (in km/sec) */
  LAMBDA= 0.0000001;        /* spin parameter          */
  MD=     0.0000001;        /* disk mass fraction      */
  MB=     0.0;        /* bulge mass fraction     */
  JD= MD;              /* disk spin fraction      */

  GasFraction= 0.0000001;    /* relative content of gas in the disk*/ 
  DiskHeight=  0.0000001;    /* thickness of disk in units of radial scale length */
  BulgeSize=   0.0000001;    /* bulge scale length in units of disk scale length  */

  N_HALO=  200000;    /* desired number of particles in dark halo */
  N_DISK=  000;       /* desired number of collisionless particles in disk */
  N_GAS=   000;       /* number of gas particles in disk */ 
  N_BULGE= 000;       /* number of bulge particles */ 


  HI_GasMassFraction=    0.0;     /* in terms of the total gas mass */
  HI_GasDiskScaleLength= 0.0;    /* in terms of scale length of the disk */ 

  Qstabilizefactor=1.0;

  /**********************************************************/


  if(argc!=2)
    {
      fprintf(stderr,"\n\nwrong argument(s).  Specify an output filename.\n\n");
      exit(0);
    }
  strcpy(filename, argv[1]);


  init_units();        /* set system of units */
  structure();         /* determine structure of halo, disk, and bulge */
  init();              /* allocate arrays */


  set_halo_positions();
  set_disk_positions();
  set_bulge_positions();
  set_gas_positions();

  
  compute_force_field();

  compute_velocity_dispersions_disk();
  compute_velocity_dispersions_halo();  
  compute_velocity_dispersions_bulge();  
 
  compute_local_escape_speed();

  set_halo_velocities();
  set_disk_velocities();    
  set_gas_velocities();
  set_bulge_velocities();


  
  save_particles(filename);


  if(fd=fopen("curve.txt","w"))
    {
      printf("writing circular velocity curve + Toomre's Q\n");
      write_header(fd);
      plot_circular_speeds(fd);
      plot_toomre_stability(fd);
      fclose(fd);
    }
  else
    {
      fprintf(stderr,"Can't open file '%s'.\n",filename);
      exit(0);
    }
  printf("done.\n");


  printf("Disk scale length: %g\n",H);
  printf("R200: %g\n",R200);
}



int write_header(FILE *fd)
{
  fprintf(fd,"%g\n",CC);
  fprintf(fd,"%g\n",V200);
  fprintf(fd,"%g\n",LAMBDA);
  fprintf(fd,"%g\n",MD);
  fprintf(fd,"%g\n",JD);
  fprintf(fd,"%g\n",MB);
  fprintf(fd,"%g\n",DiskHeight);
  fprintf(fd,"%g\n",BulgeSize);
  fprintf(fd,"\n%g\n",R200);
  fprintf(fd,"%g\n\n",H);
}
