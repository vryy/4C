/*!---------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | neumann condition                                      m.gee 3/02    |
 |                                                                      |
 | this structure holds a neumann condition                             |
 | depend on number of dofs to a fe-node and type of element it is      |
 | is connected to, the arrays can be defined in several styles         |
 *----------------------------------------------------------------------*/
typedef struct _NEUM_CONDITION
{
     INT                       curve;        /* number of load curve associated with this conditions */        

     struct _ARRAY             neum_onoff;   /* array of on-off flags */
     struct _ARRAY             neum_val;     /* values of this condition */    
     
     enum
     {
     neum_none,
     neum_live,
     neum_dead,
     neum_live_FSI,            /* on shells surfaces it's possible to have both */
     neum_opres_FSI,            /* on shells surfaces it's possible to have both */
     neum_FSI,
     pres_domain_load,
     neum_consthydro_z,
     neum_increhydro_z,
     neum_orthopressure,
     neum_LAS
     }                         neum_type;

     enum
     {
     mid,  /*nsurf=1*/
     top,  /*nsurf=2*/
     bot   /*nsurf=3*/
     }                         neum_surf;    /* load applied on top, bottom or middle surface of shell element */

} NEUM_CONDITION;
/*----------------------------------------------------------------------*
 | dirichlet condition                                    m.gee 3/02    |
 |                                                                      |
 | this structure holds a dirichlet condition                           |
 | depend on number of dofs to a fe-node and type of element it is      |
 | is connected to, the arrays can be defined in several styles         |
 *----------------------------------------------------------------------*/
typedef struct _DIRICH_CONDITION
{
     struct _ARRAY             curve;        /* time curve */
     struct _ARRAY             funct;        /* spatial function */

     struct _ARRAY             dirich_onoff; /* array of on-off flags */
     struct _ARRAY             dirich_val;   /* values of this condition */    
     enum
     {
     dirich_none,                            /* no specification */
     dirich_FSI,                             /* fluid-structure interaction */
     dirich_FSI_pseudo,
     dirich_freesurf,
     dirich_slip
     }                         dirich_type;                       

} DIRICH_CONDITION;
/*----------------------------------------------------------------------*
 | dirichlet condition                                    m.gee 3/02    |
 |                                                                      |
 | this structure holds a dirichlet condition                           |
 | depend on number of dofs to a fe-node and type of element it is      |
 | is connected to, the arrays can be defined in several styles         |
 *----------------------------------------------------------------------*/
typedef struct _SLIPDIRICH_CONDITION
{
     INT             curve;        /* time curve */
     
     DOUBLE          dirich_val;   /* values of this condition */    
     DOUBLE          length;
     DOUBLE          alpha;       

     NODE           *firstnode;
     NODE           *lastnode;
} SLIPDIRICH_CONDITION;

/*----------------------------------------------------------------------*
 | coupling condition                                     m.gee 3/02    |
 |                                                                      |
 | this structure is assigned to nodes, which are coupled in some or    |
 | all of their dofs                                                    |
 *----------------------------------------------------------------------*/
typedef struct _COUPLE_CONDITION
{
     enum _FIELDTYP            fieldtyp;        /* type of field this structure is in */
     struct _ARRAY             couple;          /* array of the coupling conditions */
} COUPLE_CONDITION;
/*----------------------------------------------------------------------*
 | fsi coupling condition                                 genk 10/02    |
 |                                                                      |
 *----------------------------------------------------------------------*/
typedef struct _FSI_COUPLE_CONDITION
{
     enum _FIELDTYP            fieldtyp;        /* type of field this structure is in */
     INT                       fsi_coupleId;
     enum _FSI_MESH            fsi_mesh;
     enum {
           fsi_real,
	   fsi_pseudo
	  }                   fsi_typ;     
} FSI_COUPLE_CONDITION;
/*----------------------------------------------------------------------*
 | fluid freesurface condition                            genk 10/02    |
 |                                                                      |        
 *----------------------------------------------------------------------*/
typedef struct _FLUID_FREESURF_CONDITION
{
     enum _FIELDTYP            fieldtyp;        /* type of field this structure is in */
} FLUID_FREESURF_CONDITION;


/*!----------------------------------------------------------------------
\brief axishell thickness condition

<pre>                                                              mn 05/03
axishell thickness condition for linear varying thickness along a dline
</pre>

*----------------------------------------------------------------------*/
typedef struct _SAXI_THICK_CONDITION
{
  DOUBLE            value[2];          /* thickness at the two dnodes of the line */
} SAXI_THICK_CONDITION;


/*!----------------------------------------------------------------------
\brief axishell load condition

<pre>                                                              mn 05/03
axishell load condition for linear varying loads along a dline in four 
directions:
 - global vertical   pv
 - global horizontal ph
 - local tangential  px
 - local normal      pw
</pre>

*----------------------------------------------------------------------*/
typedef struct _SAXI_LOAD_CONDITION
{
  DOUBLE            pv[2];          /* vertical load at the two dnodes of the line */
  INT               interpol_pv;    /* interpolation typ: 0=along curve; 1= along vertical axis */
  DOUBLE            ph[2];          /* horizontal load at the two dnodes of the line */
  INT               interpol_ph;    /* interpolation typ: 0=along curve; 1= along vertical axis */
  DOUBLE            px[2];          /* tangential load at the two dnodes of the line */
  INT               interpol_px;    /* interpolation typ: 0=along curve; 1= along vertical axis */
  DOUBLE            pw[2];          /* normal load at the two dnodes of the line */
  INT               interpol_pw;    /* interpolation typ: 0=along curve; 1= along vertical axis */
} SAXI_LOAD_CONDITION;

/*!----------------------------------------------------------------------
\brief stabilisation parameters

<pre>                                                         chfoe 01/04
stabilisation parameters used for fluid2 or fluid3 elements in case of 
'classic' GLS stabilisation as described in the dissertation of W. Wall.
</pre>

*----------------------------------------------------------------------*/
typedef struct _STAB_PAR_GLS
{
/*---------------------------------------------- stabilisation flags ---*/
   INT	iadvec;		/*!< advection stab.: 0=no; 1=yes		*/
   INT	ipres;		/*!< pressure stab.: 0=no; 1=yes		*/
   INT	ivisc;		/*!< diffusion stab.: 0=no; 1=GLS-; 2=GLS+	*/
   INT	icont;		/*!< continuity stab.: 0=no; 1=yes		*/
   INT	istapa;		/*!< version of stab. parameter			*/
   INT	istapc;		/*!< flag for stab parameter calculation	*/
   INT	mk;		/*!< 0=mk fixed 1=min(1/3,2*C); -1 mk=1/3	*/
   INT	ihele[3];	/*!< x/y/z length-def. for vel/pres/cont stab	*/
   INT	ninths;		/*!< number of integ. points for streamlength	*/
/*----------------------------------------------- stabilisation norm ---*/
   INT	norm_p;		/*!< p-norm: p+1<=infinity; 0=Max.norm		*/

/*------------------------------------------ stabilisation constants ---*/
   DOUBLE	clamb;

/*------------------------------- statiblisation control information ---*/
   INT	istrle;		/*!< has streamlength to be computed		*/
   INT	iareavol;	/*!< calculation of area length			*/
   INT	iduring;	/*!< calculation during INT.-pt.loop		*/
   INT	itau[3];	/*!< has diagonals etc. to be computed		*/
   INT	idiaxy;		/*!< flags for tau_? calculation 
			     	  (-1: before; 1: during)		*/
} STAB_PAR_GLS;



/*!----------------------------------------------------------------------
\brief FLUID lift and drag condition

<pre>                                                              mn 03/04
FLUID lift and drag conditio
</pre>

*----------------------------------------------------------------------*/
typedef struct _FLUID_LIFTDRAG_CONDITION
{
  INT                       liftdrag;           /* flag */
  DOUBLE                    ld_center[3];	/* center of body */
  INT                       alenode;            /* id of alenode corresp. to center */
  INT                       aledline;
} FLUID_LIFTDRAG_CONDITION;

