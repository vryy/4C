/*----------------------------------------------------------------------*
 | materials                                              m.gee 4/01    |
 | structure to hold all types of material laws                         |
 *----------------------------------------------------------------------*/
typedef struct _MATERIAL
{
     int                       Id;           /* Id of the material */

     enum _MATERIAL_TYP        mattyp;       /* type of material */

     union
     {
     struct _LINEAR_ELASTIC   *lin_el;       /* St. Venant-Kirchhoff material */
     struct _NEO_HOOKE        *neohooke;     /* Neo-Hooke material */
     struct _FLUID            *fluid;        /* fluid material */
     }                         m;            /* union pointer to material specific structure */

} MATERIAL;





/*----------------------------------------------------------------------*
 | St. Venant-Kirchhoff material                          m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _LINEAR_ELASTIC
{
     double                    youngs;         /* Young's modulus */
     double                    possionratio;   /* Possion ratio */
     double                    spec_weight;    /* material specific weight */
} LINEAR_ELASTIC;



/*----------------------------------------------------------------------*
 | Neo Hooke material                                     m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _NEO_HOOKE
{
     double                    youngs;         /* Young's modulus */
     double                    possionratio;   /* Possion ratio */
     double                    spec_weight;    /* material specific weight */
} NEO_HOOKE;



/*----------------------------------------------------------------------*
 | fluid material                                         m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _FLUID
{
     double                    viscosity;
     double                    density;
} FLUID;
