#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"
#include "prototypes.h"
#include "globvars.h"


#ifdef T3E
  typedef short int int4byte;   /* Note: int has 8 Bytes on the T3E ! */
#else
  typedef int int4byte;
#endif





struct io_header_1
{
  int4byte npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int4byte flag_sfr;
  int4byte flag_feedback;
  int4byte npartTotal[6];
  int4byte flag_cooling;
  int4byte num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;







void save_particles(char *fname)
{
  FILE *fd;
  int i,d;
  float xyz[3];
  double t;
  int4byte blklen;
#define BLKLEN fwrite(&blklen, sizeof(blklen), 1, fd);

  if(!(fd=fopen(fname,"w")))
    {
      printf("error opening file %s\n",fname);
      exit(0);
    }

  printf("saveing initial conditions to file `%s'\n\n",fname);

  printf("Number of particle in the gas %d\n",N_GAS);
  printf("Number of particle in the halo %d\n",N_HALO);
  printf("Number of particle in the disk %d\n",N_DISK);
  printf("Number of particle in the bulge %d\n",N_BULGE);

  for(i=1;i<=N_GAS;i++)
    {
      fprintf(fd," %g",xp_gas[i]);
      fprintf(fd," %g",yp_gas[i]);
      fprintf(fd," %g",zp_gas[i]);
      fprintf(fd," %g",vxp_gas[i]);
      fprintf(fd," %g",vyp_gas[i]);
      fprintf(fd," %g",vzp_gas[i]);
      fprintf(fd," %g\n",mp_gas[i]);

    }      
  for(i=1;i<=N_HALO;i++)
    {
      fprintf(fd," %g",xp_halo[i]);
      fprintf(fd," %g",yp_halo[i]);
      fprintf(fd," %g",zp_halo[i]);
      fprintf(fd," %g",vxp_halo[i]);
      fprintf(fd," %g",vyp_halo[i]);
      fprintf(fd," %g",vzp_halo[i]);
      fprintf(fd," %g\n",mp_halo[i]);
    }      
  for(i=1;i<=N_DISK;i++)
    {
      fprintf(fd," %g",xp_disk[i]);
      fprintf(fd," %g",yp_disk[i]);
      fprintf(fd," %g",zp_disk[i]);
      fprintf(fd," %g",vxp_disk[i]);
      fprintf(fd," %g",vyp_disk[i]);
      fprintf(fd," %g",vzp_disk[i]);
      fprintf(fd," %g\n",mp_disk[i]);
    }      
  for(i=1;i<=N_BULGE;i++)
    {
      fprintf(fd," %g",xp_bulge[i]);
      fprintf(fd," %g",yp_bulge[i]);
      fprintf(fd," %g",zp_bulge[i]);
      fprintf(fd," %g",vxp_bulge[i]);
      fprintf(fd," %g",vyp_bulge[i]);
      fprintf(fd," %g",vzp_bulge[i]);
      fprintf(fd," %g\n",mp_bulge[i]);
    }
  
  fclose(fd);
}


