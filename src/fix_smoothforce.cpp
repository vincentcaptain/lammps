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

#include <string.h>
#include <stdlib.h>
#include "fix_smoothforce.h"
#include "atom.h"
#include "atom_masks.h"
#include "accelerator_kokkos.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixSmoothForce::FixSmoothForce(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  xstr(NULL), ystr(NULL), zstr(NULL), estr(NULL), idregion(NULL), sforce(NULL)

{
  if (narg == 7) error->all(FLERR,"Illegal fix smoothforce command");

  dynamic_group_allow = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 0;
  extvector = 0;
  respa_level_support = 0;
  ilevel_respa = 0;

  xstr = ystr = zstr = NULL;

  xlo = force->numeric(FLERR,arg[3]);
  xhi = force->numeric(FLERR,arg[4]);
  delta = force->numeric(FLERR,arg[5]);
  omega = force->numeric(FLERR,arg[6]);
  fmax = force->numeric(FLERR,arg[7]);


  
  // if (strstr(arg[4],"v_") == arg[4]) {
  //   int n = strlen(&arg[4][2]) + 1;
  //   ystr = new char[n];
  //   strcpy(ystr,&arg[4][2]);
  // } else {
  //   yvalue = force->numeric(FLERR,arg[4]);
  //   ystyle = CONSTANT;
  // }
  // if (strstr(arg[5],"v_") == arg[5]) {
  //   int n = strlen(&arg[5][2]) + 1;
  //   zstr = new char[n];
  //   strcpy(zstr,&arg[5][2]);
  // } else {
  //   zvalue = force->numeric(FLERR,arg[5]);
  //   zstyle = CONSTANT;
  // }

  // optional args

  // nevery = 1;
  // iregion = -1;

  // int iarg = 6;
  // while (iarg < narg) {
  //   if (strcmp(arg[iarg],"every") == 0) {
  //     if (iarg+2 > narg) error->all(FLERR,"Illegal fix smoothforce command");
  //     nevery = atoi(arg[iarg+1]);
  //     if (nevery <= 0) error->all(FLERR,"Illegal fix smoothforce command");
  //     iarg += 2;
  //   } else if (strcmp(arg[iarg],"region") == 0) {
  //     if (iarg+2 > narg) error->all(FLERR,"Illegal fix smoothforce command");
  //     iregion = domain->find_region(arg[iarg+1]);
  //     if (iregion == -1)
  //       error->all(FLERR,"Region ID for fix smoothforce does not exist");
  //     int n = strlen(arg[iarg+1]) + 1;
  //     idregion = new char[n];
  //     strcpy(idregion,arg[iarg+1]);
  //     iarg += 2;
  //   } else if (strcmp(arg[iarg],"energy") == 0) {
  //     if (iarg+2 > narg) error->all(FLERR,"Illegal fix smoothforce command");
  //     if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
  //       int n = strlen(&arg[iarg+1][2]) + 1;
  //       estr = new char[n];
  //       strcpy(estr,&arg[iarg+1][2]);
  //     } else error->all(FLERR,"Illegal fix smoothforce command");
  //     iarg += 2;
  //   } else error->all(FLERR,"Illegal fix smoothforce command");
  // }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

  maxatom = 1;
  memory->create(sforce,maxatom,4,"smoothforce:sforce");
}

/* ---------------------------------------------------------------------- */

FixSmoothForce::~FixSmoothForce()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] estr;
  delete [] idregion;
}

/* ---------------------------------------------------------------------- */

int FixSmoothForce::setmask()
{
  datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSmoothForce::init()
{
  // check variables

  // if (xstr) {
  //   xvar = input->variable->find(xstr);
  //   if (xvar < 0)
  //     error->all(FLERR,"Variable name for fix Smoothforce does not exist");
  //   if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
  //   else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
  //   else error->all(FLERR,"Variable for fix Smoothforce is invalid style");
  // }
  // if (ystr) {
  //   yvar = input->variable->find(ystr);
  //   if (yvar < 0)
  //     error->all(FLERR,"Variable name for fix Smoothforce does not exist");
  //   if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
  //   else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
  //   else error->all(FLERR,"Variable for fix Smoothforce is invalid style");
  // }
  // if (zstr) {
  //   zvar = input->variable->find(zstr);
  //   if (zvar < 0)
  //     error->all(FLERR,"Variable name for fix Smoothforce does not exist");
  //   if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
  //   else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
  //   else error->all(FLERR,"Variable for fix Smoothforce is invalid style");
  // }
  // if (estr) {
  //   evar = input->variable->find(estr);
  //   if (evar < 0)
  //     error->all(FLERR,"Variable name for fix Smoothforce does not exist");
  //   if (input->variable->atomstyle(evar)) estyle = ATOM;
  //   else error->all(FLERR,"Variable for fix Smoothforce is invalid style");
  // } else estyle = NONE;

  // set index and check validity of region

  // if (iregion >= 0) {
  //   iregion = domain->find_region(idregion);
  //   if (iregion == -1)
  //     error->all(FLERR,"Region ID for fix Smoothforce does not exist");
  // }

  // if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
  //   varflag = ATOM;
  // else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
  //   varflag = EQUAL;
  // else varflag = CONSTANT;

  // if (varflag == CONSTANT && estyle != NONE)
  //   error->all(FLERR,"Cannot use variable energy with "
  //              "constant force in fix Smoothforce");
  // if ((varflag == EQUAL || varflag == ATOM) &&
  //     update->whichflag == 2 && estyle == NONE)
  //   error->all(FLERR,"Must use variable energy with fix Smoothforce");

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixSmoothForce::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixSmoothForce::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSmoothForce::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  if (lmp->kokkos)
    atom->sync_modify(Host, (unsigned int) (F_MASK | MASK_MASK),
                      (unsigned int) F_MASK);

  // update region if necessary

  // Region *region = NULL;
  // if (iregion >= 0) {
  //   region = domain->regions[iregion];
  //   region->prematch();
  // }

  // reallocate sforce array if necessary

  // if ((varflag == ATOM || estyle == ATOM) && atom->nmax > maxatom) {
  //   maxatom = atom->nmax;
  //   memory->destroy(sforce);
  //   memory->create(sforce,maxatom,4,"Smoothforce:sforce");
  // }

  // foriginal[0] = "potential energy" for Smoothed force
  // foriginal[123] = force on atoms before extra force Smoothed

  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;

  // constant force
  // potential energy = - x dot f in unwrapped coords
  double ctime;
  ctime = (update->ntimestep + 1)*(update->dt);

  // if (varflag == CONSTANT) {
    // double unwrap[3];
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        // if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        // domain->unmap(x[i],image[i],unwrap);
        // foriginal[0] -= xvalue*unwrap[0] + yvalue*unwrap[1] + zvalue*unwrap[2];
        double pre_factor;
        pre_factor = 0.5*(tanh((x[i][0]-xlo)/delta)-tanh((x[i][0]-xhi)/delta));
        xvalue = pre_factor*fmax*cos(omega*ctime);
        // foriginal[1] += f[i][0];
        foriginal[2] += f[i][1];
        foriginal[3] += f[i][2];
        f[i][0] += xvalue;
        foriginal[1] += f[i][0];
        // f[i][1] += yvalue;
        // f[i][2] += zvalue;
      }

  // variable force, wrap with clear/Smooth
  // potential energy = evar if defined, else 0.0
  // wrap with clear/Smooth

  // } else {

  //   modify->clearstep_compute();

  //   if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
  //   else if (xstyle == ATOM)
  //     input->variable->compute_atom(xvar,igroup,&sforce[0][0],4,0);
  //   if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
  //   else if (ystyle == ATOM)
  //     input->variable->compute_atom(yvar,igroup,&sforce[0][1],4,0);
  //   if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
  //   else if (zstyle == ATOM)
  //     input->variable->compute_atom(zvar,igroup,&sforce[0][2],4,0);
  //   if (estyle == ATOM)
  //     input->variable->compute_atom(evar,igroup,&sforce[0][3],4,0);

    // modify->Smoothstep_compute(update->ntimestep + 1);

  //   for (int i = 0; i < nlocal; i++)
  //     if (mask[i] & groupbit) {
  //       if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
  //       if (estyle == ATOM) foriginal[0] += sforce[i][3];
  //       foriginal[1] += f[i][0];
  //       foriginal[2] += f[i][1];
  //       foriginal[3] += f[i][2];
  //       if (xstyle == ATOM) f[i][0] += sforce[i][0];
  //       else if (xstyle) f[i][0] += xvalue;
  //       if (ystyle == ATOM) f[i][1] += sforce[i][1];
  //       else if (ystyle) f[i][1] += yvalue;
  //       if (zstyle == ATOM) f[i][2] += sforce[i][2];
  //       else if (zstyle) f[i][2] += zvalue;
  //     }
  // }
}

/* ---------------------------------------------------------------------- */

void FixSmoothForce::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSmoothForce::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy of Smoothed force
------------------------------------------------------------------------- */

double FixSmoothForce::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[1];
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixSmoothForce::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n+1];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixSmoothForce::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = maxatom*4 * sizeof(double);
  return bytes;
}
