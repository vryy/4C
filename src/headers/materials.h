/*----------------------------------------------------------------------*
 | materials                                              m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _MATERIAL
{
     int                       Id;

     enum _MATERIAL_TYP        mattyp;

     union
     {
     struct _LINEAR_ELASTIC   *lin_el;
     struct _NEO_HOOKE        *neohooke;
     struct _FLUID            *fluid;
     }                         m;

} MATERIAL;





/*----------------------------------------------------------------------*
 | linear-elastic material                                m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _LINEAR_ELASTIC
{
     double                    youngs;
     double                    possionratio;
} LINEAR_ELASTIC;



/*----------------------------------------------------------------------*
 | neo_hooke material                                     m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _NEO_HOOKE
{
     double                    youngs;
     double                    possionratio;
} NEO_HOOKE;



/*----------------------------------------------------------------------*
 | fluid material                                         m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _FLUID
{
     double                    viscosity;
     double                    density;
} FLUID;
