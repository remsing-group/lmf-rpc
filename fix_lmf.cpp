/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Christina Payne (Vanderbilt U)
                        Stan Moore (Sandia) for dipole terms
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
  Modified by R. Remsing for Electrostatic LMF Potential and Force
   for slab systems (e.g. VR(z) ). August 1, 2020                         
------------------------------------------------------------------------- */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fix_lmf.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "modify.h"
#include "force.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "region.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixLMF::FixLMF(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
    /* ------------------------------------------------------------------------------------------------------------------------------------------- */
   /* input is "fix <label> <group> lmf <npts> <filename>"    */
  /* ------------------------------------------------------------------------------------------------------------------------------------------- */
  if (narg < 5) error->all(FLERR,"Illegal fix lmf command");

  dynamic_group_allow = 1;
  vector_flag = 1;
  scalar_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;
  extscalar = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  // number of points in lmf file and open it
  npts = force->inumeric(FLERR,arg[3]);
  lmf_fp = fopen(arg[4],"r");
  
  force_flag = 0;
  fsum[0] = fsum[1] = fsum[2] = fsum[3] = 0.0;

  maxatom = atom->nmax;
  //memory->create(efield,maxatom,4,"efield:efield");
}

/* ---------------------------------------------------------------------- */

FixLMF::~FixLMF()
{
  delete [] lmfz;
  delete [] lmfv;
  delete [] lmff;
  fclose(lmf_fp);
}

/* ---------------------------------------------------------------------- */

int FixLMF::setmask()
{
  int mask = 0;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLMF::init()
{
  qflag = 0;
  if (atom->q_flag) qflag = 1;
  if (!qflag )
    error->all(FLERR,"Fix lmf requires atom attribute q");
  
  MPI_Comm_rank(world,&me);

}

/* ---------------------------------------------------------------------- */

void FixLMF::setup(int vflag)
{

  lmfz = new double[npts];
  lmfv = new double[npts];
  lmff = new double[npts];
  read_lmf();
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixLMF::min_setup(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   apply F = qE
------------------------------------------------------------------------- */

void FixLMF::post_force(int vflag)
{
  double **f = atom->f;
  double *q = atom->q;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  // reallocate efield array if necessary

  if (varflag == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    //memory->destroy(efield);
    //memory->create(efield,maxatom,4,"efield:efield");
  }
  // fsum[0] = "potential energy" for added force
  // fsum[123] = extra force added to atoms

  fsum[0] = fsum[1] = fsum[2] = fsum[3] = 0.0;
  force_flag = 0;

  double **x = atom->x;
  double fx,fy,fz;
  double loz,dz;
  double vz;
  int iZ;
  
  loz = lmfz[0];
  dz = lmfz[1]-lmfz[0];


    if (qflag) {
      for (int i = 0; i < nlocal; i++)
       {
          iZ = int( (x[i][2]-loz)/dz );
	  
	  fz = q[i]*lmff[iZ] + (x[i][2]-lmfz[iZ])*(lmff[iZ+1]-lmff[iZ])*q[i]/(lmfz[iZ+1]-lmfz[iZ]); // interpolate the force
	  
	  vz = q[i]*lmfv[iZ] + (x[i][2]-lmfz[iZ])*(lmfv[iZ+1]-lmfv[iZ])*q[i]/(lmfz[iZ+1]-lmfz[iZ]); // interpolate the potential
          
          f[i][2] += fz;
          fsum[0] += vz;
          fsum[3] += fz;
        }
    }

}

/* ---------------------------------------------------------------------- */

void FixLMF::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixLMF::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixLMF::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = atom->nmax*4 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   return energy added by fix
------------------------------------------------------------------------- */

double FixLMF::compute_scalar(void)
{
  if (force_flag == 0) {
    MPI_Allreduce(fsum,fsum_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return fsum_all[0];
}

/* ----------------------------------------------------------------------
   return total extra force due to fix
------------------------------------------------------------------------- */

double FixLMF::compute_vector(int n)
{
  if (force_flag == 0) {
    MPI_Allreduce(fsum,fsum_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return fsum_all[n+1];
}

/* ----------------------------------------------------------------------
   read the lmf potential from file
   - RCR, Aug 1, 2020
------------------------------------------------------------------------- */
void FixLMF::read_lmf()
{
  int i = 0;
  int j = 0;
  if (me == 0) {
    int maxchar = 75;
    char line[maxchar];
    char *word;
    while(fgets(line,maxchar,lmf_fp) != NULL) {
      word = strtok(line," \t");
      while(word != NULL) 
      {
        i++;
        if (i == 1) 
	{
          lmfz[j] = atof(word);
        } 
	else if(i == 2)
	{
	  lmfv[j] = atof(word);
	}
	else if(i == 3)
	{
	  lmff[j] = -atof(word); // minus sign because the input prints dV/dz, and F = - dV/dz
	  j++;
	  i=0;
        }
        word = strtok(NULL," \t");
      }
    }
    fclose(lmf_fp);
  }
  MPI_Bcast(lmfz,npts,MPI_DOUBLE,0,world);
  MPI_Bcast(lmfv,npts,MPI_DOUBLE,0,world);
  MPI_Bcast(lmff,npts,MPI_DOUBLE,0,world);

}

