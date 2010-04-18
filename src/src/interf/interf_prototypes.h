/*!----------------------------------------------------------------------
\file
\brief contains all prototypes for 1D interface element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>
*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_INTERF

/*!
\addtogroup INTERF
*//*! @{ (documentation module open)*/

/*-----------------------------------------------------------------------*
|  if_init.c                                                 ah 05/03  |
*-----------------------------------------------------------------------*/
void ifinit(PARTITION *actpart);
/*-----------------------------------------------------------------------*
| if_intg.c                                                 ah 05/03  |
*-----------------------------------------------------------------------*/
void ifintg(ELEMENT       *ele,
            INTERF_DATA   *data);
/*----------------------------------------------------------------------*
| if_static_ke.c                                               ah 05/03  |
*-----------------------------------------------------------------------*/
void ifstatic_ke(ELEMENT       *ele,
                 INTERF_DATA   *data,
                 MATERIAL      *mat,
                 ARRAY         *estif_global,
                 ARRAY         *emass_global,
                 DOUBLE        *force,
                 INT            init);
/*----------------------------------------------------------------------*
| if_bop.c                                                    ah 05/03  |
*-----------------------------------------------------------------------*/
void if_bop(DIS_TYP    typ,
            DOUBLE   **bop,
            DOUBLE    *funct,
            DOUBLE     co,
            DOUBLE     si,
            INT        flag);
/*----------------------------------------------------------------------*
| if_mat.c                                                    ah 05/03  |
*-----------------------------------------------------------------------*/
void if_mat(ELEMENT   *ele,        /* actual element (macro)               */
            MATERIAL  *mat,        /* actual material                      */
            DOUBLE   **bop,        /* B-operator                           */
            DOUBLE   **D,          /* material tangent                     */
            DOUBLE    *T,          /* stresses                             */
            INT        ip,         /* ID of actual GP                      */
            INT        istore,     /* is it update?                        */
            INT        newval,     /* is it stress calculation             */
            INT        smallscale, /* is it call from sm-ele in multiscale */
            ELEMENT   *actsmele,   /* if multiscale: actual sm-element     */
            DOUBLE    *jumpu_tot,  /* if multiscale: total displ. jump     */
            DOUBLE    *DELTAjumpu_tot);/* if multi: tot increm. displ.jump */
/*----------------------------------------------------------------------*
| if_cal_deltau.c                                             ah 05/03  |
*-----------------------------------------------------------------------*/
void if_jumpu(ELEMENT  *ele,
              DOUBLE  **bop,
              DOUBLE   *disjump,
              DOUBLE   *deltadisjump);
/*----------------------------------------------------------------------*
| if_stiff_fint.c                                             ah 05/03  |
*-----------------------------------------------------------------------*/
void if_ke(INT       iel,
           INT       flag,
           DOUBLE  **stiff,
           DOUBLE  **bop,
           DOUBLE  **Q,
           DOUBLE    fac);
/*----------------------------------------------------------------------*
| if_stiff_fint.c                                             ah 05/03  |
*-----------------------------------------------------------------------*/
void if_fint(INT      iel,
             DOUBLE  *T,
             DOUBLE   fac,
             DOUBLE **bop,
             DOUBLE  *fint);
/*----------------------------------------------------------------------*
| if_stress.c                                                 ah 05/03  |
*-----------------------------------------------------------------------*/
void if_stress(ELEMENT       *ele,
               INTERF_DATA   *data,
               MATERIAL      *mat,
               INT            init);
/*----------------------------------------------------------------------*
| if_mat_dyn.c                                                 ah 05/03  |
*-----------------------------------------------------------------------*/
void if_mat_dyn(ELEMENT   *ele,
                MATERIAL  *mat,
                DOUBLE   **bop,
                DOUBLE   **D,
                DOUBLE    *T,
                INT        ip,
                DOUBLE     istore,
                DOUBLE     newval);
/*----------------------------------------------------------------------*
| if_service.c                                                 ah 05/03  |
*-----------------------------------------------------------------------*/
void if_dirichnode(ELEMENT  *actele);
/*----------------------------------------------------------------------*
| if_funcderiv.c                                                ah 05/03  |
*-----------------------------------------------------------------------*/
void if_funcderiv(DOUBLE  e1,
                  DIS_TYP typ,
                  DOUBLE *x_mid,
                  DOUBLE *y_mid,
                  DOUBLE  b_parabel,
                  DOUBLE  c_parabel,
                  DOUBLE *funct,
                  DOUBLE *co,
                  DOUBLE *si,
                  DOUBLE *det);
/*----------------------------------------------------------------------*
| if_service.c                                                 ah 05/03  |
*-----------------------------------------------------------------------*/
void if_permstiff(DOUBLE **estif,
                  DOUBLE **Kdd,
                  INT      iele,
                  INT      ield);
/*----------------------------------------------------------------------*
| if_service.c                                                 ah 05/03  |
*-----------------------------------------------------------------------*/
void if_permforce(DOUBLE    *force,   /* "mixed" element int. force   */
                  DOUBLE    *fintd,    /* 2.part int. force           */
                   INT       iele,    /* num.of equiv.strain nodes   */
                   INT       ield);   /* num of displacement nodes  */
/*----------------------------------------------------------------------*
| if_restart.c                                                 ah 06/04  |
*-----------------------------------------------------------------------*/
void if_write_restart(ELEMENT   *actele,
                      MATERIAL  *mat,
                      INT        nhandle,
                      long int  *handles,
                      INT        init);
/*----------------------------------------------------------------------*
| if_restart.c                                                 ah 06/04  |
*-----------------------------------------------------------------------*/
void if_read_restart(ELEMENT  *actele,
                     MATERIAL *mat,
                     long int *handles,
                     INT       init);
/*----------------------------------------------------------------------*
| if_mat.c                                                    ah 09/04  |
*-----------------------------------------------------------------------*/
void if_mat_thermodyn(ELEMENT   *ele,        /* actual element (macro)               */
            MATERIAL  *mat,        /* actual material                      */
            DOUBLE   **bop,        /* B-operator                           */
            DOUBLE   **D,          /* material tangent                     */
            DOUBLE    *T,          /* stresses                             */
            INT        ip,         /* ID of actual GP                      */
            INT        istore,     /* is it update?                        */
            INT        newval,     /* is it stress calculation             */
            INT        smallscale, /* is it call from sm-ele in multiscale */
            ELEMENT   *actsmele,   /* if multiscale: actual sm-element     */
            DOUBLE    *jumpu_tot,  /* if multiscale: total displ. jump     */
            DOUBLE    *DELTAjumpu_tot);/* if multi: tot increm. displ.jump */
/*----------------------------------------------------------------------*/
/*! @} (documentation module close)*/

#endif /*D_INTERF*/
#endif
