/*----------------------------------------------------------------------*
 | data bunker                                           m.gee 06/01    |
 *----------------------------------------------------------------------*/
typedef struct _BUNKER
{
int              bunker_Id;
int              used;
int              size;
/* 
style of access_handles is as follows:
access_handles[*access]    = size of shelf[*access]
*/ 
int             ahp_pss_handle;
int            *access_handles;
/* 
style of shelf is as follows:
size of shelf[*access] = access_handles[*access] 
shelf[*access][handle] = is a certain array
*/ 
int              sp_pss_handle;
struct _ARRAY  **shelf;
/* this holds the pss-file handles of all shelf entries */
} BUNKER;
