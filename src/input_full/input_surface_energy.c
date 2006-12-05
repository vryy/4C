/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/wiechert
            089/28915303
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../brick1/brick1.h"

extern struct _DESIGN *design;
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 | input of surface energetical properties                  lw 11/06    |
 *----------------------------------------------------------------------*/

void input_surface_energy()
{

#ifdef D_BRICK1
#ifdef SURFACE_ENERGY

  INT ierr, ierr1, ierr2;
  DSURF *actdsurf;
  INT dsurf_ID;
  INT number;
  INT i;

  /*--------------------------------------------------------------------*/
  #ifdef DEBUG
  dstrc_enter("input_surface_energy");
  #endif
  /*--------------------------------------------------------------------*/

  frfind("--SURFACE CONDITIONS");
  frread();
  frread();

  while(strncmp(allfiles.actplace,"------",6)!=0)
  {
    frint("E", &dsurf_ID, &ierr);
    dsassert(ierr==1,"Cannot read design-surface energy conditions");
    dsurf_ID--;

    /*-------------------------------------------------- find the dsurf */
    actdsurf=NULL;
    for (i=0; i<design->ndsurf; i++)
    {
      if (design->dsurf[i].Id ==  dsurf_ID)
      {
         actdsurf = &(design->dsurf[i]);
         break;
      }
    }
    dsassert(actdsurf!=NULL,"Cannot read design-surface energy conditions");
    actdsurf->surface=1;

    frchk("SURFACTANT", &ierr1);
    if (ierr1)
    {
      actdsurf->surface_flag=0;
      frdouble("k1", &(actdsurf->k1), &ierr);
      dsassert(ierr==1, "Cannot read adsorption coefficient");
      frdouble("k2", &(actdsurf->k2), &ierr);
      dsassert(ierr==1, "Cannot read desorption coefficient");
      frdouble("Cbulk", &(actdsurf->C), &ierr);
      dsassert(ierr==1, "Cannot read bulk surfactant concentration");
      frdouble("m1", &(actdsurf->m1), &ierr);
      dsassert(ierr==1, "Cannot read first isotherm slope");
      frdouble("m2", &(actdsurf->m2), &ierr);
      dsassert(ierr==1, "Cannot read second isotherm slope");
      frdouble("gamma_0", &(actdsurf->gamma_0), &ierr);
      dsassert(ierr==1, "Cannot read gamma_0");
      frdouble("gamma_min", &(actdsurf->gamma_min), &ierr);
      dsassert(ierr==1, "Cannot read gamma_min");
      frdouble("gamma_min_eq", &(actdsurf->gamma_min_eq), &ierr);
      dsassert(ierr==1, "Cannot read gamma_min_eq");
    }

    frchk("SURFACE TENSION", &ierr2);
    if (ierr2)
    {
      actdsurf->surface_flag=1;
      frdouble("gamma", &(actdsurf->const_gamma), &ierr);
      dsassert(ierr==1, "Cannot read gamma");
    }

    if (ierr1==0 && ierr2==0)
    {
      dserror("Cannot read design-surface energy conditions");
    }

    frread();
  }
#endif
#endif
  /*--------------------------------------------------------------------*/
  #ifdef DEBUG
  dstrc_exit();
  #endif
  return;
} /* end of input_surface_energy */
