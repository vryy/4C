/*!----------------------------------------------------------------------
\file
\brief contains all prototypes for 1D interface element

*----------------------------------------------------------------------*/
#ifdef D_INTERF

/*!
\addtogroup INTERF
*//*! @{ (documentation module open)*/

/*-----------------------------------------------------------------------*
|  if_init.c                                                 mn 05/03  |
*-----------------------------------------------------------------------*/
void ifinit(PARTITION *actpart);
/*-----------------------------------------------------------------------*
| if_intg.c                                                 mn 05/03  |
*-----------------------------------------------------------------------*/
void ifintg(ELEMENT       *ele,
            INTERF_DATA   *data);
/*----------------------------------------------------------------------*
| if_static_ke.c                                               mn 05/03  |
*-----------------------------------------------------------------------*/
void ifstatic_ke(ELEMENT       *ele,
                 INTERF_DATA   *data,
                 MATERIAL      *mat,
                 ARRAY         *estif_global,
                 ARRAY         *emass_global,
                 DOUBLE        *force,
                 INT            init);
/*----------------------------------------------------------------------*
| if_bop.c                                                    mn 05/03  |
*-----------------------------------------------------------------------*/
void if_bop(DIS_TYP    typ,
            DOUBLE   **bop,
            DOUBLE    *funct,
            DOUBLE     co,
            DOUBLE     si,
            INT        flag);
/*----------------------------------------------------------------------*
| if_mat.c                                                    mn 05/03  |
*-----------------------------------------------------------------------*/
void if_mat(ELEMENT   *ele,
            MATERIAL  *mat,
            DOUBLE   **bop,
            DOUBLE   **D,
            DOUBLE    *T,
            INT        ip,
            DOUBLE     istore,
            DOUBLE     newval);
/*----------------------------------------------------------------------*
| if_cal_deltau.c                                             mn 05/03  |
*-----------------------------------------------------------------------*/
void if_jumpu(ELEMENT  *ele,
              DOUBLE  **bop,
              DOUBLE   *disjump,
              DOUBLE   *deltadisjump);
/*----------------------------------------------------------------------*
| if_stiff_fint.c                                             mn 05/03  |
*-----------------------------------------------------------------------*/
void if_ke(INT       iel,
           INT       flag,
           DOUBLE  **stiff,
           DOUBLE  **bop,
           DOUBLE  **Q,
           DOUBLE    fac);
/*----------------------------------------------------------------------*
| if_stiff_fint.c                                             mn 05/03  |
*-----------------------------------------------------------------------*/
void if_fint(INT      iel,
             DOUBLE  *T,
             DOUBLE   fac,
             DOUBLE **bop,
             DOUBLE  *fint);
/*----------------------------------------------------------------------*
| if_stress.c                                                 mn 05/03  |
*-----------------------------------------------------------------------*/
void if_stress(ELEMENT       *ele,
               INTERF_DATA   *data,
               MATERIAL      *mat,
               INT            init);
/*----------------------------------------------------------------------*
| if_mat_dyn.c                                                 mn 05/03  |
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
| if_service.c                                                 mn 05/03  |
*-----------------------------------------------------------------------*/
void if_dirichnode(ELEMENT  *actele);
/*----------------------------------------------------------------------*
| if_funcderiv.c                                                 mn 05/03  |
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
| if_service.c                                                 mn 05/03  |
*-----------------------------------------------------------------------*/
void if_permstiff(DOUBLE **estif,
                  DOUBLE **Kdd,
                  INT      iele,
                  INT      ield);
/*----------------------------------------------------------------------*
| if_service.c                                                 mn 05/03  |
*-----------------------------------------------------------------------*/
void if_permforce(DOUBLE    *force,   /* "mixed" element int. force   */
                  DOUBLE    *fintd,    /* 2.part int. force           */
                   INT       iele,    /* num.of equiv.strain nodes   */
                   INT       ield);   /* num of displacement nodes  */
/*----------------------------------------------------------------------*
| if_funcderiv.c                                                 mn 05/03  |
*-----------------------------------------------------------------------*/
void if_func_bope(DOUBLE   e1,
                  DOUBLE  *x_mid,
                  DOUBLE  *y_mid,
                  DOUBLE  *functe,
                  DOUBLE  *coe,
                  DOUBLE  *sie,
                  DOUBLE  *dete,
                  INT      flag,
                  DOUBLE **bope);
/*----------------------------------------------------------------------*/
#endif /*D_INTERF*/
/*! @} (documentation module close)*/
