/*----------------------------------------------------------------------*
 | structure for gid output data                         m.gee 12/01    |
 *----------------------------------------------------------------------*/
typedef struct _GIDSET
{
     FIELDTYP                   fieldtyp;               /* type of field */
     int                        fieldnamelenght;        /* lenght of fieldname */
     char                      *fieldname;              /* name of the field in characters */
     char                       standardrangetable[18]; /* name of standardrangetable */

     int                        is_shell8_22;           /* 4-noded shell8 2x2 GP */
     char                      *shell8_22_name;
     int                        is_shell8_33;           /* 8/9-noded shell8 3x3 GP */
     char                      *shell8_33_name;
     int                        is_wall1_22;            /* 4-noded wall1 2x2 GP */
     char                      *wall1_22_name;
     int                        is_wall1_33;            /* 8/9-noded wall1 3x3 GP */
     char                      *wall1_33_name;
     int                        is_brick1_222;          /* 8-noded brick1 2x2x2 GP */
     char                      *brick1_222_name;
     int                        is_brick1_333;          /* 20/27 noded brick1 3x3x3 GP */
     char                      *brick1_333_name;
     int                        is_fluid2_22;           /* 4-noded fluid2 2x2 GP */
     char                      *fluid2_22_name;
     int                        is_fluid2_33;           /* 8/9-noded fluid2 3x3 GP */
     char                      *fluid2_33_name;
     int                        is_fluid3_222;          /* 8-noded fluid3 2x2x2 GP */
     char                      *fluid3_222_name;
     int                        is_fluid3_333;          /* 20/27-noded fluid3 3x3x3 GP */
     char                      *fluid3_333_name;
     int                        is_ale_11;              /* 4-noded ale 1x1 GP */
     char                      *ale_11_name;
     int                        is_ale_22;              /* 4-noded ale 2x2 GP */
     char                      *ale_22_name;
     int                        is_ale_111;             /* 8-noded ale 1x1x1 GP */
     char                      *ale_111_name;
     int                        is_ale_222;             /* 8-noded ale 2x2x2 GP */
     char                      *ale_222_name;
     int                        is_ale_tri_1;           /* 3-noded tri ale 1 GP */
     char                      *ale_tri_1_name;
     int                        is_ale_tri_3;           /* 3-noded tri ale 3 GP */
     char                      *ale_tri_3_name;
     int                        is_ale_tet_1;           /* 4-noded tet ale 1 GP */
     char                      *ale_tet_1_name;
     int                        is_ale_tet_4;           /* 4-noded tet ale 4 GP */
     char                      *ale_tet_4_name;
     

} GIDSET;
/*------------------------ global variable needed by gid postprocessing */
GIDSET *gid;
/*----------------------------------------------------------------------*/
