/*!----------------------------------------------------------------------
\file
\brief headerfile for 

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef D_MLSTRUCT

/*! 
\addtogroup MLSTRUCT
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief gradient enhanced wall element

<pre>                                                            ah 08/03
This structure contains the working array at the Gaussian points for the 
gradient enhanced wall element
</pre>

*----------------------------------------------------------------------*/

typedef struct _SM_MACRO_INFO
{
     DOUBLE       funct_ma[4];     /* Macro funct[i] at submeshnode*/
     DOUBLE       derivxy_ma[2][4];/* Macro funct[i],x funct[i],y at submeshnode*/
     DOUBLE       strain_ma[4];    /* Macro strain xx, yy, xy, zz at submeshnode*/
} SM_MACRO_INFO;

/*---------------------------------------------------------------------*/
typedef struct _SM_ELEMENT_DATA
{
  struct  _SM_GP_DATA     *sm_GPdata;
} SM_ELEMENT_DATA;

/*---------------------------------------------------------------------*/
typedef struct _SM_NODAL_DATA
{
  DOUBLE  displ_mi[2];
  DOUBLE  store_displ_mi[2];
  DOUBLE  incre_displ_mi[2];
} SM_NODAL_DATA;

/*---------------------------------------------------------------------*/
typedef struct _SM_GP_DATA
{
/*---------------------------------------------------------------- wall*/
     DOUBLE       kappa;     /* */
     DOUBLE       damage;     /* */
     DOUBLE       stress[4]; /* */
/*----------------------------------------- interface old material law */
     DOUBLE       dt;        /* */
     DOUBLE       dn;        /* */
     DOUBLE       utpl;      /* */
     DOUBLE       D[2][2];   /* */
     DOUBLE       T[2];      /* */
     INT          yip;
/*----------------------------------------- interface new material law */
     DOUBLE       kappa_n;        /* */
     DOUBLE       kappa_t;        /* */
} SM_GP_DATA;
/*! @} (documentation module close)*/

#endif /* D_MLSTRUCT */
