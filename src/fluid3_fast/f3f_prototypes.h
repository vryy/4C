/*-----------------------------------------------------------------------*/
/*!
\file
\brief all prototypes for the fast 3d fluid element


<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

 */
/*-----------------------------------------------------------------------*/

#ifdef D_FLUID3_F


void f3fcalgalexf(
    DOUBLE  *eforce,
    DOUBLE  *funct,
    DOUBLE  *edeadn,
    DOUBLE  *edeadng,
    DOUBLE  *fac,
    DOUBLE  *facsl,
    DOUBLE  *facsr,
    INT     *sizevec
    );


void f3fcalstabexf(
    DOUBLE  *eforce,
    DOUBLE  *derxy,
    DOUBLE  *derxy2,
    DOUBLE  *edead,
    DOUBLE  *velint,
    DOUBLE  *fac,
    DOUBLE  *ths,
    DOUBLE  *thp,
    DOUBLE  *tau,
    DOUBLE  *paravec,
    INT     *flagvec,
    INT     *sizevec
    );


void f3fhex(
    DOUBLE  *funct,
    DOUBLE  *deriv,
    DOUBLE  *deriv2,
    DOUBLE  *r,
    DOUBLE  *s,
    DOUBLE  *t,
    INT     *typ,
    INT     *icode,
    INT     *sizevec
    );


void f3ftet(
    DOUBLE  *funct,
    DOUBLE  *deriv,
    DOUBLE  *deriv2,
    DOUBLE  *r,
    DOUBLE  *s,
    DOUBLE  *t,
    INT     *typ,
    INT     *icode,
    INT     *sizevec
    );


void f3fjaco(
    DOUBLE  *funct,
    DOUBLE  *deriv,
    DOUBLE  *xjm,
    DOUBLE  *det,
    DOUBLE  *elecord_f,
    INT     *sizevec
    );


void f3fgder(
    DOUBLE  *derxy,
    DOUBLE  *deriv,
    DOUBLE  *xjm,
    DOUBLE  *xji,
    DOUBLE  *det,
    INT     *sizevec
    );


void f3fgder2(
    DOUBLE  *elecord_f,
    DOUBLE  *xjm,
    DOUBLE  *bm_in,
    DOUBLE  *xder2_in,
    DOUBLE  *derxy,
    DOUBLE  *derxy2,
    DOUBLE  *deriv2,
    INT     *sizevec
    );


void f3fgder2loop(
    DOUBLE  *elecord_f,
    DOUBLE  *xjm,
    DOUBLE  *bm_in,
    DOUBLE  *xder2_in,
    DOUBLE  *derxy,
    DOUBLE  *derxy2,
    DOUBLE  *deriv2,
    INT     *sizevec
    );


void f3fcalgalk(
    DOUBLE  *estif,
    DOUBLE  *velint,
    DOUBLE  *gridvint,
    DOUBLE  *vderxy,
    DOUBLE  *funct,
    DOUBLE  *derxy,
    DOUBLE  *fac,
    DOUBLE  *paravec,
    INT     *flagvec,
    INT     *sizevec
    );


void f3fcalgalm(
    DOUBLE  *estif,
    DOUBLE  *funct,
    DOUBLE  *fac,
    INT     *sizevec
    );


void f3fcalif(
    DOUBLE  *eforce,
    DOUBLE  *covint,
    DOUBLE  *velint,
    DOUBLE  *funct,
    DOUBLE  *derxy,
    DOUBLE  *derxy2,
    DOUBLE  *facsl,
    DOUBLE  *tau,
    DOUBLE  *paravec,
    INT     *flagvec,
    INT     *sizevec
    );


void f3fveli(
    DOUBLE  *velint,
    DOUBLE  *funct,
    DOUBLE  *evel,
    INT     *sizevec
    );


void f3fvder(
    DOUBLE  *vderxy,
    DOUBLE  *derxy,
    DOUBLE  *evel,
    INT     *sizevec
    );


void f3fvder2(
    DOUBLE  *vderxy2,
    DOUBLE  *derxy2,
    DOUBLE  *evel,
    INT     *sizevec
    );


void f3fcovi(
    DOUBLE  *vderxy,
    DOUBLE  *velint,
    DOUBLE  *covint,
    INT     *sizevec
    );


void f3fprei(
    DOUBLE  *preint,
    DOUBLE  *funct,
    DOUBLE  *epre,
    INT     *sizevec
    );


void f3fpder(
    DOUBLE  *pderxy,
    DOUBLE  *derxy,
    DOUBLE  *epre,
    INT     *sizevec
    );


void f3fcalstabk(
    DOUBLE  *estif,
    DOUBLE  *velint,
    DOUBLE  *alecovint,
    DOUBLE  *gridvint,
    DOUBLE  *vderxy,
    DOUBLE  *funct,
    DOUBLE  *derxy,
    DOUBLE  *derxy2,
    DOUBLE  *fac,
    DOUBLE  *tau,
    DOUBLE  *paravec,
    INT     *flagvec,
    INT     *sizevec
);


void f3fcalstabm(
    DOUBLE  *estif,
    DOUBLE  *velint,
    DOUBLE  *funct,
    DOUBLE  *derxy,
    DOUBLE  *derxy2,
    DOUBLE  *fac,
    DOUBLE  *tau,
    DOUBLE  *paravec,
    INT     *flagvec,
    INT     *sizevec
    );


void f3fcaltf(
    DOUBLE  *eforce,
    DOUBLE  *velint,
    DOUBLE  *alecovint,
    DOUBLE  *covint,
    DOUBLE  *funct,
    DOUBLE  *derxy,
    DOUBLE  *derxy2,
    DOUBLE  *vderxy,
    DOUBLE  *vderxy2,
    DOUBLE  *pderxy,
    DOUBLE  *preint,
    DOUBLE  *fac,
    DOUBLE  *ths,
    DOUBLE  *thp,
    DOUBLE  *tau,
    DOUBLE  *paravec,
    INT     *flagvec,
    INT     *sizevec
    );


void f3fmassrhs(
    DOUBLE  *emass,
    DOUBLE  *eaccn,
    DOUBLE  *eiforce,
    INT     *sizevec
    );


void f3fmast(
    DOUBLE  *estif,
    DOUBLE  *emass,
    DOUBLE  *thsl,
    INT     *nis,
    INT     *sizevec
    );


void f3fcalele(
    ELEMENT        *ele[LOOPL],
    ARRAY           *estif_fast,
    ARRAY           *emass_fast,
    ARRAY           *etforce_fast,
    ARRAY           *eiforce_fast,
    ARRAY           *edforce_fast,
    INT            *hasdirich,
    INT            *hasext,
    INT             init,
    INT             loop
    );


void f3fcalint(
    ELEMENT         *ele[LOOPL],
    DOUBLE          *elecord_f,
    DOUBLE          *tau_f,
    INT             *hasext,
    DOUBLE          *estif_f,
    DOUBLE          *emass_f,
    DOUBLE          *etforce_f,
    DOUBLE          *eiforce_f,
    DOUBLE          *funct_f,
    DOUBLE          *deriv_f,
    DOUBLE          *deriv2_f,
    DOUBLE          *xjm_f,
    DOUBLE          *derxy_f,
    DOUBLE          *derxy2_f,
    DOUBLE          *eveln_f,
    DOUBLE          *evelng_f,
    DOUBLE          *epren_f,
    DOUBLE          *edeadn_f,
    DOUBLE          *edeadng_f,
    DOUBLE          *velint_f,
    DOUBLE          *vel2int_f,
    DOUBLE          *covint_f,
    DOUBLE          *vderxy_f,
    DOUBLE          *pderxy_f,
    DOUBLE          *vderxy2_f,
    DOUBLE          *wa1_f,
    DOUBLE          *wa2_f,
    INT              sizevec[5]
    );


    void f3fcalelecord(
        ELEMENT        *ele[LOOPL],
        DOUBLE         *elecord_f,
        INT             sizevec[5]
        );


    void f3fcalset(
        ELEMENT         *ele[LOOPL],
        DOUBLE          *eveln,
        DOUBLE          *evelng,
        DOUBLE          *epren,
        DOUBLE          *edeadn,
        DOUBLE          *edeadng,
        INT             *hasext,
        INT              sizevec[5]
        );


    void f3fcalelesize(
        ELEMENT         *ele[LOOPL],
        DOUBLE          *funct_f,
        DOUBLE          *deriv_f,
        DOUBLE          *deriv2_f,
        DOUBLE          *derxy_f,
        DOUBLE          *xjm_f,
        DOUBLE          *evel_f,
        DOUBLE          *velint_f,
        DOUBLE          *wa1_f,
        DOUBLE          *elecord_f,
        DOUBLE          *tau,
        INT              sizevec[5]
        );


    void f3fcalelesize2(
        ELEMENT         *ele[LOOPL],
        DOUBLE          *velint,
        DOUBLE          *derxy,
        DOUBLE          *tau,
        DOUBLE           visc,
        INT              ntyp,
        INT              sizevec[5]
        );


    void fluidf_caldirich(
        ELEMENT         *actele[LOOPL],
        DOUBLE          *dforces,
        DOUBLE          *estif,
        INT             *hasdirich,
        INT              is_relax,
        INT              sizevec[5]
        );


    void fluid3_fast(
        PARTITION   *actpart,
        INTRA       *actintra,
        ELEMENT     *ele[LOOPL],
        ARRAY       *estif_fast,
        ARRAY       *emass_fast,
        ARRAY       *etforce_fast,
        ARRAY       *eiforce_fast,
        ARRAY       *edforce_fast,
        CALC_ACTION *action,
        INT         *hasdirich,
        INT         *hasext,
        CONTAINER   *container,
        INT          counter
        );


    void f3fcalstabpar(
        ELEMENT         *ele[LOOPL],
        DOUBLE          *velint,
        DOUBLE          *tau,
        DOUBLE           visc,
        INT              ntyp,
        INT              iflag,
        INT              sizevec[5]
        );


    void f3fmake_estif(
        DOUBLE            *estif,
        DOUBLE            *emass,
        INT                sizevec[5]
        );


    void f3fcalstab(
        ELEMENT      *ele[LOOPL],
        INT           loop
        );

    void f3fcaldirich(
        ELEMENT         *actele[LOOPL],
        DOUBLE          *dforces,
        DOUBLE          *estif,
        INT             *hasdirich,
        INT              is_relax,
        INT              sizevec[5]
        );

    void f3f_write_restart(
        ELEMENT   *actele,
        INT        nhandle,
        long int  *handles
        );

    void f3f_read_restart(
        ELEMENT   *actele,
        INT        nhandle,
        long int  *handles
        );

    void f3fstress(
        FLUID_STRESS  str,
        INT           viscstr,
        ELEMENT      *ele[LOOPL],
        INT           is_relax,
        INT           aloopl
        );

void f3fsigint(
    DOUBLE   *preint,
    DOUBLE   *vderxy,
    DOUBLE   *sigint,
    DOUBLE   *visc,
    INT      *iv,
    INT      *sizevec
    );

void f3flgpl(
    INT         i,
    INT         n,
    DOUBLE     *zr,
    DOUBLE      z,
    DOUBLE     *value
    );

DOUBLE f3frsn(
    INT   node,
    INT   irs
    );


void f3fsext(
    DOUBLE   *nostr,
    DOUBLE   *sigint,
    DOUBLE   *xgr,
    DOUBLE   *xgs,
    DOUBLE   *xgt,
    INT       nir,
    INT       nis,
    INT       nit,
    INT       sizevec[6]
    );


void f3f_out_gid_sol_str(
    FILE       *out,
    FIELD *actfield,
    INT         init
    );


void f3fcalelestress(
    INT             viscstr,
    ELEMENT        *ele[LOOPL],
    DOUBLE         *evel,
    DOUBLE         *epre,
    DOUBLE         *funct,
    DOUBLE         *deriv,
    DOUBLE         *derxy,
    DOUBLE         *vderxy,
    DOUBLE         *xjm,
    DOUBLE         *wa1,
    DOUBLE         *elecord,
    DOUBLE         *sigint,
    DOUBLE         *nostr,
    INT             sizevec[6]
    );


void f3fliftdrag(
    ELEMENT       *ele,
    CONTAINER     *container
    );


void f3fhexc(
    DOUBLE      funct[MAXNOD_F3],
    DOUBLE      deriv[3][MAXNOD_F3],
    DOUBLE      r,
    DOUBLE      s,
    DOUBLE      t,
    DIS_TYP     typ
    );

void f3fjacoc(
    DOUBLE      funct[MAXNOD_F3],
    DOUBLE      deriv[3][MAXNOD_F3],
    DOUBLE      xjm[3][3],
    DOUBLE     *det,
    ELEMENT    *ele,
    INT         iel
    );

void f3fcalseta(
    ELEMENT         *ele[LOOPL],
    DOUBLE          *eveln,
    DOUBLE          *evelng,
    DOUBLE          *ealecovn,
    DOUBLE          *ealecovng,
    DOUBLE          *egridv,
    DOUBLE          *epren,
    DOUBLE          *edeadn,
    DOUBLE          *edeadng,
    INT             *hasext,
    INT              sizevec[5]
    );

void f3fcalinta(
    ELEMENT         *ele[LOOPL],
    DOUBLE          *elecord,
    DOUBLE          *tau,
    INT             *hasext,
    DOUBLE          *estif,
    DOUBLE          *emass,
    DOUBLE          *etforce,
    DOUBLE          *eiforce,
    DOUBLE          *funct,
    DOUBLE          *deriv,
    DOUBLE          *deriv2,
    DOUBLE          *xjm,
    DOUBLE          *derxy,
    DOUBLE          *derxy2,
    DOUBLE          *eveln,
    DOUBLE          *evelng,
    DOUBLE          *ealecovn,
    DOUBLE          *ealecovng,
    DOUBLE          *egridv,
    DOUBLE          *epren,
    DOUBLE          *edeadn,
    DOUBLE          *edeadng,
    DOUBLE          *velint,
    DOUBLE          *vel2int,
    DOUBLE          *covint,
    DOUBLE          *alecovint,
    DOUBLE          *gridvint,
    DOUBLE          *vderxy,
    DOUBLE          *pderxy,
    DOUBLE          *vderxy2,
    DOUBLE          *wa1,
    DOUBLE          *wa2,
    INT              sizevec[6]
    );



#endif /* ifdef D_FLUID3_F */


