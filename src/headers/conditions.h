/*----------------------------------------------------------------------*
 | type definitions conditions                            m.gee 4/01    |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 |  nodal conditions                            m.gee 4/01    |
 *----------------------------------------------------------------------*/
typedef struct _COND_NODE
{
     int                       curve;

     int                       isneum;
     struct _ARRAY             neum_onoff;
     struct _ARRAY             neum_val;

     int                       isdirich;
     struct _ARRAY             dirich_onoff;
     struct _ARRAY             dirich_val;

     int                       iscoupled;                          
     struct _ARRAY             couple;

     int                       fsi_iscoupled;                          

} COND_NODE;

/*----------------------------------------------------------------------*
 | element neumann/dirichlet condition types             m.gee 10/01    |
 *----------------------------------------------------------------------*/
typedef enum _COND_ELE_TYP
{
                       ne_live,
                       ne_dead
} COND_ELE_TYP;
/*----------------------------------------------------------------------*
 |  element conditions                                    m.gee 5/01    |
 *----------------------------------------------------------------------*/
typedef struct _COND_ELEMENT
{
     enum _COND_ELE_TYP        condtyp;
     int                       curve;

     int                       isneum;
     struct _ARRAY             neum_onoff;
     struct _ARRAY             neum_val;

} COND_ELEMENT;



