/*----------------------------------------------------------------------*
 | structure for gid output data                         m.gee 12/01    |
 *----------------------------------------------------------------------*/
typedef struct _GIDSET
{
     FIELDTYP                   fieldtyp;               /* type of field */
     INT                        fieldnamelenght;        /* lenght of fieldname */
     char                      *fieldname;              /* name of the field in characters */
     char                       standardrangetable[18]; /* name of standardrangetable */

     INT                        is_shell8_22;           /* 4-noded shell8 2x2 GP */
     char                      *shell8_22_name;
     INT                        is_shell8_33;           /* 8/9-noded shell8 3x3 GP */
     char                      *shell8_33_name;
     INT                        is_shell9_4_22;         /* 4-noded shell9 2x2 GP */
     char                      *shell9_4_22_name;
     INT                        is_shell9_4_33;         /* 4-noded shell9 3x3 GP */
     char                      *shell9_4_33_name;
     INT                        is_shell9_8_22;         /* 8-noded shell9 2x2 GP */
     char                      *shell9_8_22_name;
     INT                        is_shell9_8_33;         /* 8-noded shell9 3x3 GP */
     char                      *shell9_8_33_name;
     INT                        is_shell9_9_22;         /* 9-noded shell9 2x2 GP */
     char                      *shell9_9_22_name;
     INT                        is_shell9_9_33;         /* 9-noded shell9 3x3 GP */
     char                      *shell9_9_33_name;
     INT                        is_wall1_22;            /* 4-noded wall1 2x2 GP */
     char                      *wall1_22_name;
     INT                        is_wall1_33;            /* 8/9-noded wall1 3x3 GP */
     char                      *wall1_33_name;
     INT                        is_brick1_222;          /* 8-noded brick1 2x2x2 GP */
     char                      *brick1_222_name;
     INT                        is_brick1_333;          /* 20/27 noded brick1 3x3x3 GP */
     char                      *brick1_333_name;
     INT                        is_fluid2_22;           /* 4-noded fluid2 2x2 GP */
     char                      *fluid2_22_name;
     INT                        is_fluid2_33;           /* 8/9-noded fluid2 3x3 GP */
     char                      *fluid2_33_name;
     INT                        is_fluid2_pro_22;       /* 4-noded fluid2_pro 2x2 GP */
     char                      *fluid2_pro_22_name;
     INT                        is_fluid2_pro_33;       /* 8/9-noded fluid2 3x3 GP */
     char                      *fluid2_pro_33_name;
     INT                        is_fluid3_222;          /* 8-noded fluid3 2x2x2 GP */
     char                      *fluid3_222_name;
     INT                        is_fluid3_333;          /* 20/27-noded fluid3 3x3x3 GP */
     char                      *fluid3_333_name;
     INT                        is_ale_11;              /* 4-noded ale 1x1 GP */
     char                      *ale_11_name;
     INT                        is_ale_22;              /* 4-noded ale 2x2 GP */
     char                      *ale_22_name;
     INT                        is_ale_111;             /* 8-noded ale 1x1x1 GP */
     char                      *ale_111_name;
     INT                        is_ale_222;             /* 8-noded ale 2x2x2 GP */
     char                      *ale_222_name;
     INT			is_beam3_21;		/* 2-noded beam3 1 GP */
     char		       *beam3_21_name;
     INT			is_beam3_22;		/* 2-noded beam3 2 GP */
     char		       *beam3_22_name;
     INT			is_beam3_32;		/* 3-noded beam3 2 GP */
     char		       *beam3_32_name;
     INT			is_beam3_33;		/* 3-noded beam3 3 GP */
     char		       *beam3_33_name;
     INT                        is_ale_tri_1;           /* 3-noded tri ale 1 GP */
     char                      *ale_tri_1_name;
     INT                        is_ale_tri_3;           /* 3-noded tri ale 3 GP */
     char                      *ale_tri_3_name;
     INT                        is_ale_tet_1;           /* 4-noded tet ale 1 GP */
     char                      *ale_tet_1_name;
     INT                        is_ale_tet_4;           /* 4-noded tet ale 4 GP */
     char                      *ale_tet_4_name;
     INT                        is_axishell;            /* 2-noded axishell */
     char                      *axishell_name;
     

} GIDSET;
/*------------------------ global variable needed by gid postprocessing */
GIDSET *gid;
/*----------------------------------------------------------------------*/
