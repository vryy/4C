#include "../headers/standardtypes.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 |  init bunkers                                         m.gee 6/01     |
 *----------------------------------------------------------------------*/
void bunker_init()
{
int        i,j;
FIELD      *actfield;
SOLVAR     *actsolv;
NODE       *actnode;
ELEMENT    *actele;
#ifdef DEBUG 
dstrc_enter("bunker_init");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------------ init db system */
db_init();
/*----------------------------------------------------------------------*/
for (i=0; i<genprob.numfld; i++)
{
   actfield = &(field[i]);
   actsolv = &(solv[i]);
/*-------------------------------------- create a bunker for each field */
/* each bunker has a bunker Id (which is type int). This Bunker Id
   is needed for writing from and to this bunker
*/   
   db_create(&(actfield->db_bunker_Id));
/*----------- loop all elements in the field and create an access_number */   
/* The elements access number and data handles are stored in the ELEMENT
   or NODE in an ARRAY named db_acess as follows:
   db_access.a.iv[0] = db_bunker_Id;
   db_access.a.iv[1] = db_access_number;
   db_access.a.iv[2] = number of actual handles;
   db_access.a.iv[3] = db_data_handle;
   .
   .
   .
   db_access.a.iv[n] = data handles;
   db_access has an initial size of 8, but can of course be redefined to 
   hold more than 8-3 = 5 data handles. 
*/
   for (j=0; j<actfield->numele; j++)
   {
      actele = &(actfield->element[j]);
      amdef("db_access",&(actele->db_access),8,1,"IV");
      actele->db_access.a.iv[0] = actfield->db_bunker_Id;
      db_create_access(actele->db_access.a.iv[0],&(actele->db_access.a.iv[1]));
      actele->db_access.a.iv[2]=0;
   }   
/*--------------------------------- now do the same thing for the NODEs */
   for (j=0; j<actfield->numnp; j++)
   {
      actnode = &(actfield->node[j]);
      amdef("db_access",&(actnode->db_access),8,1,"IV");
      actnode->db_access.a.iv[0] = actfield->db_bunker_Id;
      db_create_access(actnode->db_access.a.iv[0],&(actnode->db_access.a.iv[1]));
      actnode->db_access.a.iv[2]=0;
   }
/*----------------------------------------- create an access for SOLVAR */
      amdef("db_access",&(actsolv->db_access),8,1,"IV");
      actsolv->db_access.a.iv[0] = actfield->db_bunker_Id;
      db_create_access(actsolv->db_access.a.iv[0],&(actsolv->db_access.a.iv[1]));
      actsolv->db_access.a.iv[2]=0;
/*----------------------------------------------------------------------*/
} /* end of loop over fields */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of bunker_init */





/*----------------------------------------------------------------------*
 |  create access id for partitions                      m.gee 6/01     |
 *----------------------------------------------------------------------*/
void part_bunker_access()
{
int        i;
PARTITION  *actpart;
FIELD      *actfield;
#ifdef DEBUG 
dstrc_enter("part_bunker_access");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
for (i=0; i<genprob.numfld; i++)
{
   actpart  = &(partition[i]);
   actfield = &(field[i]);
/* db_access.a.iv[0] = db_bunker_Id;
   db_access.a.iv[1] = db_access_number;
   db_access.a.iv[2] = number of actual handles;
   db_access.a.iv[3] = db_data_handle;
   .
   .
   .
   db_access.a.iv[n] = data handles;
   db_access has an initial size of 8, but can of course be redefined to 
   hold more than 8-3 = 5 data handles. 
*/
/*---------------------------------------- create an access for actpart */
      amdef("db_access",&(actpart->db_access),8,1,"IV");
      actpart->db_access.a.iv[0] = actfield->db_bunker_Id;
      db_create_access(actpart->db_access.a.iv[0],&(actpart->db_access.a.iv[1]));
      actpart->db_access.a.iv[2]=0;
/*----------------------------------------------------------------------*/
} /* end of loop over fields */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of part_bunker_access */
