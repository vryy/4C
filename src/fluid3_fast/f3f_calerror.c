/*-----------------------------------------------------------------------*/
/*!
\file
\brief integration loop for error calculation

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

 */
/*-----------------------------------------------------------------------*/


/*!
\addtogroup FLUID3
*//*! @{ (documentation module open)*/


#ifdef D_FLUID3_F

#include "../headers/standardtypes.h"
#include "../fluid3/fluid3_prototypes.h"
#include "../fluid3_fast/f3f_prototypes.h"
#include "../fluid3/fluid3.h"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01  |
  | vector of material laws                                            |
  | defined in global_control.c                                        |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;


/*-----------------------------------------------------------------------*/
/*!
  \brief integration loop for the error calculation


  \param ele[LOOPL]      ELEMENT    (i) the set of elements
  \param elecord        *DOUBLE     (o) vector containing nodal coordinates
  \param funct          *DOUBLE     (o) shape functions
  \param deriv          *DOUBLE     (o) deriv. of shape funcs
  \param xjm            *DOUBLE     (o) jacobian matrix
  \param evelng         *DOUBLE     (o) vels at time n+g
  \param epren          *DOUBLE     (o) pres at time n
  \param velint         *DOUBLE     (o) vel at int point
  \param xyzint         *DOUBLE     (o) coords at int point
  \param container      *CONTAINER  (o) contains variables defined in container.h
  \param sizevec[6]      INT        (i) some sizes

  \return void

  \author mn
  \date   08/05

 */
/*-----------------------------------------------------------------------*/
void f3f_int_error(
    ELEMENT         *ele[LOOPL],
    DOUBLE          *elecord,
    DOUBLE          *funct,
    DOUBLE          *deriv,
    DOUBLE          *xjm,
    DOUBLE          *evelng,
    DOUBLE          *epren,
    DOUBLE          *velint,
    DOUBLE          *xyzint,
    CONTAINER       *container,
    INT              sizevec[6]
    )

{

  INT       l;               /* a counter */
  INT       intc;            /* "integration case" for tet */
  INT       nir,nis,nit;     /* number of integration nodesin r,s,t direction */
  INT       icode=2;         /* flag for eveluation of shape functions */
  INT       lr, ls, lt;      /* counter for integration */
  INT       actmat;          /* actual material number */

  INT       ierr;
  INT       inttyp;

  DOUBLE    det[LOOPL];
  DOUBLE    preint[LOOPL];
  DOUBLE    fac[LOOPL];      /* total integration vactor */
  DOUBLE    facr=0.0, facs=0.0, fact=0.0;      /* integration weights */
  DOUBLE    e1,e2,e3;        /* natural coordinates of integr. point */
  DIS_TYP   typ;             /* element type */

  FLUID_DYNAMIC   *fdyn;
  FLUID_DATA      *data;

  /*Fortran variables - passed to fortran subroutines as parameters.*/
  DOUBLE   paravec[4];


#ifdef DEBUG
  dstrc_enter("f3f_int_error");
#endif


  /* initialisation */
  fdyn   = alldyn[genprob.numff].fdyn;
  data   = fdyn->data;
  actmat = ele[0]->mat-1;
  typ    = ele[0]->distyp;


  paravec[0] = fdyn->acttime;
  paravec[1] = mat[actmat].m.fluid->viscosity;
  paravec[2] = (DOUBLE)fdyn->init;


  intc = 0;
  nir  = 0;
  nis  = 0;
  nit  = 0;


  switch (typ)
  {
    case hex8: case hex20: case hex27:  /* hex - element */
      icode   = 2;
      /* initialise integration */
      nir = ele[0]->e.f3->nGP[0];
      nis = ele[0]->e.f3->nGP[1];
      nit = ele[0]->e.f3->nGP[2];
      intc= 0;
      break;

    case tet10: /* tet - element */
      icode   = 2;
      /* do NOT break at this point!!! */

    case tet4:    /* initialise integration */
      nir  = ele[0]->e.f3->nGP[0];
      nis  = 1;
      nit  = 1;
      intc = ele[0]->e.f3->nGP[1];
      break;

    default:
      dserror("typ unknown!");
  } /* end switch (typ) */


  /*----------------------------------------------------------------------*
    |               start loop over integration points                     |
   *----------------------------------------------------------------------*/
  for (lr=0;lr<nir;lr++)
  {
    for (ls=0;ls<nis;ls++)
    {
      for (lt=0;lt<nit;lt++)
      {


        /* get values of  shape functions and their derivatives */
        switch(typ)
        {
          case hex8: case hex20: case hex27:   /* hex - element */
            e1   = data->qxg[lr][nir-1];
            e2   = data->qxg[ls][nis-1];
            e3   = data->qxg[lt][nit-1];
            facr = data->qwgt[lr][nir-1];
            facs = data->qwgt[ls][nis-1];
            fact = data->qwgt[lt][nit-1];

            if (typ==hex8)
              inttyp=8;
            else if(typ==hex20)
              inttyp=20;
            else if(typ==hex27)
              inttyp=27;

            f3fhex(funct, deriv, NULL, &e1, &e2, &e3, &inttyp, &icode, sizevec);

            break;
          case tet4: case tet10:  /* tet - element */
            e1   = data->txgr[lr][intc];
            facr = data->twgt[lr][intc];
            e2   = data->txgs[lr][intc];
            facs = ONE;
            e3   = data->txgt[lr][intc];
            fact = ONE;

            if (typ==tet4)
              inttyp=4;
            else if(typ==tet10)
              inttyp=10;

            f3ftet(funct, deriv, NULL, &e1, &e2, &e3, &inttyp,
                &icode, sizevec);

            break;
          default:
            dserror("typ unknown!");
        } /* end switch (typ) */


        /* compute Jacobian matrix */
        f3fjaco(funct, deriv, xjm, det, elecord, sizevec);
        for(l=0;l<sizevec[4];l++)
          fac[l] = facr * facs * fact * det[l];


        /* get velocities at integraton point */
        f3fveli(velint, funct, evelng, sizevec);

        /* get coordinates at integraton point */
        f3fveli(xyzint, funct, elecord, sizevec);

        /* get pressure at integration point */
        f3fprei(preint, funct, epren, sizevec);



        f3finterr(elecord, xyzint, evelng, velint, epren, preint,
            &(container->vel_error), &(container->pre_error),
            &(container->vel_norm), &(container->pre_norm),
            &(container->error_norm), fac, paravec, sizevec, &ierr);


        if (ierr == 1)
          dserror("Error calculation for fluid3_fast only possible for beltrami problem!!");


      } /* end of loop over integration points lt*/
    } /* end of loop over integration points ls */
  } /* end of loop over integration points lr */



#ifdef DEBUG
  dstrc_exit();
#endif



  return;
} /* end of f3fint_usfem_ale */


#endif


/*! @} (documentation module close)*/



