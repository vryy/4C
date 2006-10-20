/*======================================================================*/
/*!
\file
\brief Configure element THERM3

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
	    http://www.lnm.mw.tum.de/Members/bornemann
	    089-289-15237
</pre>

\author bborn
\date 09/06
*/
#ifdef D_THERM3

/*----------------------------------------------------------------------*/
/* headers */
#include "../headers/standardtypes.h"
#include "therm3.h"


/*======================================================================*/
/*!
\brief Check if pre-definable quantities (given by #define) are
       permissible

\return void

\author bborn
\date 09/06
*/
void th3_cfg_chkdef()
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_cfg_chkdef");
#endif

  /*--------------------------------------------------------------------*/
  /* check sanity of pre-definable global constants (e.g. MAXNOD) */
  if (MAXNOD_THERM3 < 4)
  {
    dserror("MAXNOD_THERM3 or MAXNOD : maximum of element nodes\n"
            "   must be at least  4 for a TET4\n"
            "   must be at least  8 for a HEX8\n"
            "   must be at least 10 for a TET10\n"
            "   must be at least 20 for a HEX20\n"
            "   must be at least 27 for a HEX27\n");
  }
  if (MAXSID_THERM3 < 4)
  {
    dserror("MAXSID_THERM3 : maximal number element sides\n"
            "   must be at least 4 for a TET4 or TET10\n"
            "   must be at least 6 for a HEX8, HEX20 or HEX27\n");
  }
  if (MAXEDG_THERM3 < 6)
  {
    dserror("MAXEDG_THERM3 : maximal number element edges\n"
            "   must be at least  6 for a TET4 or TET10\n"
            "   must be at least 12 for a HEX8, HEX20 or HEX27\n");
  }
  if (MAXNS_THERM3 < 3)
  {
    dserror("MAXNS_THERM3 : maximal number of nodes on a side\n"
            "   must be at least 3 for a TET4\n"
            "   must be at least 4 for a HEX8\n"
            "   must be at least 6 for a TET10\n"
            "   must be at least 8 for a HEX20\n"
            "   must be at least 9 for a HEX27\n");
  }  
  if (MAXNE_THERM3 < 2)
  {
    dserror("MAXNE_THERM3 : maximal number of nodes on an edge\n"
            "   must be at least 2 for a TET4 or HEX8\n"
            "   must be at least 3 for a TET10, HEX20 or HEX27\n");
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th3_cfg_chkdef */


/*======================================================================*/
/*!
\brief Initialise global variables available to all elements by accessing
       the TH3_DATA data variable.

\param   data   TH3_DATA*   (o)   element topology etc.
\return  void

\author bborn
\date 09/06
*/
void th3_cfg_init(TH3_DATA *data)
{

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_cfg_init");
#endif

  /*====================================================================*/
  /* hexhedra */
  /*--------------------------------------------------------------------*/
  /* parameter coordinates (r,s,t) of nodes 
   * of biunit cube [-1,1]x[-1,1]x[-1,1] 
   *  8-node hexahedron: node 0,1,...,7 
   * 20-node hexahedron: node 0,1,...,19
   * 27-node hexahedron: node 0,1,....26 */
  /*                                                                    
   *                      t                                             
   *                      |                                             
   *             4========|==19============7                           
   *           //|        |               /||                           
   *          // |        |              //||                            
   *         //  |        |             // ||                           
   *        16   |       25           18   ||                           
   *       //    |        |           //   ||                           
   *      //    12        |  24      //    15                         
   *     //      |        |         //     ||                           
   *     5===========17============6       ||                           
   *    ||       |        |        ||      ||                           
   *    ||  21   |       26--------||-23---------s                      
   *    ||       |       /         ||      ||                           
   *    ||       0------/----11----||------3                            
   *    ||      /      /           ||     //                            
   *    13     /     22            14    //                             
   *    ||    /      /             ||   //                              
   *    ||   8      /    20        || 10                                
   *    ||  /      /               || //                                
   *    || /      r                ||//                                 
   *    ||/                        ||/                                  
   *     1============9============2                                    
   *                                                                    
   */
  if (MAXNOD_THERM3 >= 8)
  {
    /* node 0 */
    data->nodhrst[0][0] = -1.0;  /* r-coordinate */
    data->nodhrst[0][1] = -1.0;  /* s-coordinate */
    data->nodhrst[0][2] = -1.0;  /* t-coordinate */
    /* node 1 */
    data->nodhrst[1][0] = +1.0;
    data->nodhrst[1][1] = -1.0;
    data->nodhrst[1][2] = -1.0;
    /* node 2 */
    data->nodhrst[2][0] = +1.0;
    data->nodhrst[2][1] = +1.0;
    data->nodhrst[2][2] = -1.0;
    /* node 3 */
    data->nodhrst[3][0] = -1.0;
    data->nodhrst[3][1] = +1.0;
    data->nodhrst[3][2] = -1.0;
    /* node 4 */
    data->nodhrst[4][0] = -1.0;
    data->nodhrst[4][1] = -1.0;
    data->nodhrst[4][2] = +1.0;
    /* node 5 */
    data->nodhrst[5][0] = +1.0;
    data->nodhrst[5][1] = -1.0;
    data->nodhrst[5][2] = +1.0;
    /* node 6 */
    data->nodhrst[6][0] = +1.0;
    data->nodhrst[6][1] = +1.0;
    data->nodhrst[6][2] = +1.0;
    /* node 7 */
    data->nodhrst[7][0] = -1.0;
    data->nodhrst[7][1] = +1.0;
    data->nodhrst[7][2] = +1.0;
  }
  if (MAXNOD_THERM3 >=20)
  {
    /* node 8 */
    data->nodhrst[8][0] =  0.0;
    data->nodhrst[8][1] = -1.0;
    data->nodhrst[8][2] = -1.0;
    /* node 9 */
    data->nodhrst[9][0] =  1.0;
    data->nodhrst[9][1] =  0.0;
    data->nodhrst[9][2] = -1.0;
    /* node 10 */
    data->nodhrst[10][0] =  0.0;
    data->nodhrst[10][1] =  1.0;
    data->nodhrst[10][2] = -1.0;
    /* node 11 */
    data->nodhrst[11][0] = -1.0;
    data->nodhrst[11][1] =  0.0;
    data->nodhrst[11][2] = -1.0;
    /* node 16 */
    data->nodhrst[16][0] =  0.0;
    data->nodhrst[16][1] = -1.0;
    data->nodhrst[16][2] =  1.0;
    /* node 17 */
    data->nodhrst[17][0] =  1.0;
    data->nodhrst[17][1] =  0.0;
    data->nodhrst[17][2] =  1.0;
    /* node 18 */
    data->nodhrst[18][0] =  0.0;
    data->nodhrst[18][1] =  1.0;
    data->nodhrst[18][2] =  1.0;
    /* node 19 */
    data->nodhrst[19][0] = -1.0;
    data->nodhrst[19][1] =  0.0;
    data->nodhrst[19][2] =  1.0;
    /* node 12 */
    data->nodhrst[12][0] = -1.0;
    data->nodhrst[12][1] = -1.0;
    data->nodhrst[12][2] =  0.0;
    /* node 13 */
    data->nodhrst[13][0] =  1.0;
    data->nodhrst[13][1] = -1.0;
    data->nodhrst[13][2] =  0.0;
    /* node 14 */
    data->nodhrst[14][0] =  1.0;
    data->nodhrst[14][1] =  1.0;
    data->nodhrst[14][2] =  0.0;
    /* node 15 */
    data->nodhrst[15][0] = -1.0;
    data->nodhrst[15][1] =  1.0;
    data->nodhrst[15][2] =  0.0;
  }
  if (MAXNOD_THERM3 >= 27)
  {
    /* node 20 */
    data->nodhrst[20][0] =  0.0;
    data->nodhrst[20][1] =  0.0;
    data->nodhrst[20][2] = -1.0;
    /* node 26 */
    data->nodhrst[25][0] =  0.0;
    data->nodhrst[25][1] =  0.0;
    data->nodhrst[25][2] =  1.0;
    /* node 21 */
    data->nodhrst[21][0] =  0.0;
    data->nodhrst[21][1] = -1.0;
    data->nodhrst[21][2] =  0.0;
    /* node 22 */
    data->nodhrst[22][0] =  1.0;
    data->nodhrst[22][1] =  0.0;
    data->nodhrst[22][2] =  0.0;
    /* node 23 */
    data->nodhrst[23][0] =  0.0;
    data->nodhrst[23][1] =  1.0;
    data->nodhrst[23][2] =  0.0;
    /* node 24 */
    data->nodhrst[24][0] = -1.0;
    data->nodhrst[24][1] =  0.0;
    data->nodhrst[24][2] =  0.0;
    /* node 26 */
    data->nodhrst[26][0] =  0.0;
    data->nodhrst[26][1] =  0.0;
    data->nodhrst[26][2] =  0.0;
  }
  /*--------------------------------------------------------------------*/
  /* sides
   *   - nodes on sides nodsidhl/nodsidhq
   *   - anchor point on side ancsidhl/ancsidhq
   *   - span vectors on side spasidhl/spasidhq */
  /*                                                                    
   *                      t                                             
   *                      |                                             
   *             4--------|--19------------7                           
   *            /|        |               /|                            
   *           / |        |              / |                             
   *          /  |        |             /  |                            
   *        16   |      si5           18   |                            
   *        /    |        |           /    |                            
   *       /    12        | si4      /    15                            
   *      /      |        |         /      |                            
   *     5-----------17------------6       |                            
   *     |       |        |        |       |                            
   *     | si1   |       27--------|-si3---------s                      
   *     |       |       /         |       |                            
   *     |       0------/----11----|-------3                            
   *     |      /      /           |      /                             
   *    13     /    si2           14     /                              
   *     |    /      /             |    /                               
   *     |   8      /   si0        |  10                                
   *     |  /      /               |  /                                 
   *     | /      r                | /                                  
   *     |/                        |/                                   
   *     1------------9------------2                                    
   *                                                                    
   */
  /* 8-node hexahedron */
  if ( (MAXSID_THERM3 >= 6) && (MAXNS_THERM3 >= 4) )
  {
    /* side 0 */
    data->nodsidh[0][0] = 0;
    data->nodsidh[0][1] = 1;
    data->nodsidh[0][2] = 2;
    data->nodsidh[0][3] = 3;
    /* side 5 */
    data->nodsidh[5][0] = 4;
    data->nodsidh[5][1] = 5;
    data->nodsidh[5][2] = 6;
    data->nodsidh[5][3] = 7;
    /* side 1 */
    data->nodsidh[1][0] = 0;
    data->nodsidh[1][1] = 1;
    data->nodsidh[1][2] = 5;
    data->nodsidh[1][3] = 4;
    /* side 2 */
    data->nodsidh[2][0] = 1;
    data->nodsidh[2][1] = 2;
    data->nodsidh[2][2] = 6;
    data->nodsidh[2][3] = 5;
    /* side 3 */
    data->nodsidh[3][0] = 2;
    data->nodsidh[3][1] = 3;
    data->nodsidh[3][2] = 7;
    data->nodsidh[3][3] = 6;
    /* side 4 */
    data->nodsidh[4][0] = 3;
    data->nodsidh[4][1] = 0;
    data->nodsidh[4][2] = 4;
    data->nodsidh[4][3] = 7;
  }
  /* 20-node hexahedron */
  if ( (MAXSID_THERM3 >= 6) && (MAXNS_THERM3 >=8 ) )
  {
    /* side 0 */
    data->nodsidh[0][4] = 8;
    data->nodsidh[0][5] = 9;
    data->nodsidh[0][6] = 10;
    data->nodsidh[0][7] = 11;
    /* side 5 */
    data->nodsidh[5][4] = 16;
    data->nodsidh[5][5] = 17;
    data->nodsidh[5][6] = 18;
    data->nodsidh[5][7] = 19;
    /* side 1 */
    data->nodsidh[1][4] = 8;
    data->nodsidh[1][5] = 13;
    data->nodsidh[1][6] = 16;
    data->nodsidh[1][7] = 12;
    /* side 2 */
    data->nodsidh[2][4] = 9;
    data->nodsidh[2][5] = 14;
    data->nodsidh[2][6] = 17;
    data->nodsidh[2][7] = 13;
    /* side 3 */
    data->nodsidh[3][4] = 10;
    data->nodsidh[3][5] = 15;
    data->nodsidh[3][6] = 18;
    data->nodsidh[3][7] = 14;
    /* side 4 */
    data->nodsidh[4][4] = 11;
    data->nodsidh[4][5] = 12;
    data->nodsidh[4][6] = 19;
    data->nodsidh[4][7] = 15;
  }
  /* 27-node hexahedron */
  if ( (MAXSID_THERM3 >= 6) && (MAXNS_THERM3 >=9 ) )
  {
    /* side 0 */
    data->nodsidh[0][8] = 20;
    /* side 5 */
    data->nodsidh[5][8] = 25;
    /* side 1 */
    data->nodsidh[1][8] = 21;
    /* side 2 */
    data->nodsidh[2][8] = 22;
    /* side 3 */
    data->nodsidh[3][8] = 23;
    /* side 4 */
    data->nodsidh[4][8] = 24;
  }
  /* hexahedra anchors and spans */
  if (MAXSID_THERM3 >= 6)
  {
    /* side 0 */
    data->ancsidh[0][0] = 0.0;  /* anchor r-coord */
    data->ancsidh[0][1] = 0.0;  /* anchor s-ccord */
    data->ancsidh[0][2] = -1.0;  /* anchor t-ccord */
    data->dirsidh[0][0][0] = 1.0;  /* 1st span vector in r-direction */
    data->dirsidh[0][0][1] = 0.0;  /* 1st span vector in r-direction */
    data->dirsidh[0][0][2] = 0.0;  /* 1st span vector in r-direction */
    data->dirsidh[0][1][0] = 0.0;  /* 2nd span vector in s-direction */
    data->dirsidh[0][1][1] = 1.0;  /* 2nd span vector in s-direction */
    data->dirsidh[0][1][2] = 0.0;  /* 2nd span vector in s-direction */
    /* side 5 */
    data->ancsidh[5][0] = 0.0;
    data->ancsidh[5][1] = 0.0;
    data->ancsidh[5][2] = 1.0;
    data->dirsidh[5][0][0] = 1.0;
    data->dirsidh[5][0][1] = 0.0;
    data->dirsidh[5][0][2] = 0.0;
    data->dirsidh[5][1][0] = 0.0;
    data->dirsidh[5][1][1] = 1.0;
    data->dirsidh[5][1][2] = 0.0;
    /* side 1 */
    data->ancsidh[1][0] = 0.0;
    data->ancsidh[1][1] = -1.0;
    data->ancsidh[1][2] = 0.0;
    data->dirsidh[1][0][0] = 1.0;
    data->dirsidh[1][0][1] = 0.0;
    data->dirsidh[1][0][2] = 0.0;
    data->dirsidh[1][1][0] = 0.0;
    data->dirsidh[1][1][1] = 0.0;
    data->dirsidh[1][1][2] = 1.0;
    /* side 2 */
    data->ancsidh[2][0] = 1.0;
    data->ancsidh[2][1] = 0.0;
    data->ancsidh[2][2] = 0.0;
    data->dirsidh[2][0][0] = 0.0;
    data->dirsidh[2][0][1] = 1.0;
    data->dirsidh[2][0][2] = 0.0;
    data->dirsidh[2][1][0] = 0.0;
    data->dirsidh[2][1][1] = 0.0;
    data->dirsidh[2][1][2] = 1.0;
    /* side 3 */
    data->ancsidh[3][0] = 0.0;
    data->ancsidh[3][1] = 1.0;
    data->ancsidh[3][2] = 0.0;
    data->dirsidh[3][0][0] = 1.0;
    data->dirsidh[3][0][1] = 0.0;
    data->dirsidh[3][0][2] = 0.0;
    data->dirsidh[3][1][0] = 0.0;
    data->dirsidh[3][1][1] = 0.0;
    data->dirsidh[3][1][2] = 1.0;
    /* side 4 */
    data->ancsidh[4][0] = -1.0;
    data->ancsidh[4][1] = 0.0;
    data->ancsidh[4][2] = 0.0;
    data->dirsidh[4][0][0] = 0.0;
    data->dirsidh[4][0][1] = 1.0;
    data->dirsidh[4][0][2] = 0.0;
    data->dirsidh[4][1][0] = 0.0;
    data->dirsidh[4][1][1] = 0.0;
    data->dirsidh[4][1][2] = 1.0;
  }
    

  /*--------------------------------------------------------------------*/
  /* edges
   *   -  nodes on edges nodedghl/nodedghq
   *   -  anchor on edge ancedgh
   *   -  span on edge spaedgh
  /*                                                                    
   *                      t                                             
   *                      |                                             
   *             4--------|11de------------7                           
   *            /|        |               /|                            
   *           / |        |              / |                             
   *          e  |        |            10  |                            
   *         d   4       21            d   |                            
   *        8    d        |           e    7                            
   *       /     e        |  24      /     d                            
   *      /      |        |         /      e                            
   *     5----------ed9------------6       |                            
   *     |       |        |        |       |                            
   *     |  21   |       26--------|--23---------s                      
   *     |       |       /         |       |                            
   *     |       0------/---3de----|-------3                            
   *     5      /      /           6      /                             
   *     d     /     22            d     /                              
   *     e    e      /             e    2                               
   *     |   d      /    20        |   d                                
   *     |  0      /               |  e                                 
   *     | /      r                | /                                  
   *     |/                        |/                                   
   *     1-----------ed1-----------2                                    
   *                                                                    
   */
  /* 8-node hexahedron */
  if ( (MAXEDG_THERM3 >= 12) && (MAXNE_THERM3 >= 2) )
  {
    /* edge 0 */
    data->nodedghl[0][0] = 0;
    data->nodedghl[0][1] = 1;
    /* edge 1 */
    data->nodedghl[1][0] = 1;
    data->nodedghl[1][1] = 2;
    /* edge 2 */
    data->nodedghl[2][0] = 2;
    data->nodedghl[2][1] = 3;
    /* edge 3 */
    data->nodedghl[3][0] = 3;
    data->nodedghl[3][1] = 0;
    /* edge 8 */
    data->nodedghl[8][0] = 4;
    data->nodedghl[8][1] = 5;
    /* edge 9 */
    data->nodedghl[9][0] = 5;
    data->nodedghl[9][1] = 6;
    /* edge 10 */
    data->nodedghl[10][0] = 6;
    data->nodedghl[10][1] = 7;
    /* edge 11 */
    data->nodedghl[11][0] = 7;
    data->nodedghl[11][1] = 4;
    /* edge 4 */
    data->nodedghl[4][0] = 0;
    data->nodedghl[4][1] = 4;
    /* edge 5 */
    data->nodedghl[5][0] = 1;
    data->nodedghl[5][1] = 5;
    /* edge 6 */
    data->nodedghl[6][0] = 2;
    data->nodedghl[6][1] = 6;
    /* edge 7 */
    data->nodedghl[7][0] = 3;
    data->nodedghl[7][1] = 7;
  }
  /* 20,27-node hexahedron */
  if ( (MAXEDG_THERM3 >= 12) && (MAXNE_THERM3 >= 3) )
  {
    /* edge 0 */
    data->nodedghq[0][0] = 0;
    data->nodedghq[0][1] = 8;
    data->nodedghq[0][2] = 1;
    /* edge 1 */
    data->nodedghq[1][0] = 1;
    data->nodedghq[1][1] = 9;
    data->nodedghq[1][2] = 2;
    /* edge 2 */
    data->nodedghq[2][0] = 2;
    data->nodedghq[2][1] = 10;
    data->nodedghq[2][2] = 3;
    /* edge 3 */
    data->nodedghq[3][0] = 3;
    data->nodedghq[3][1] = 11;
    data->nodedghq[3][2] = 0;
    /* edge 8 */
    data->nodedghq[8][0] = 4;
    data->nodedghq[8][1] = 16;
    data->nodedghq[8][2] = 5;
    /* edge 9 */
    data->nodedghq[9][0] = 5;
    data->nodedghq[9][1] = 17;
    data->nodedghq[9][2] = 6;
    /* edge 10 */
    data->nodedghq[10][0] = 6;
    data->nodedghq[10][1] = 18;
    data->nodedghq[10][2] = 7;
    /* edge 11 */
    data->nodedghq[11][0] = 7;
    data->nodedghq[11][1] = 19;
    data->nodedghq[11][2] = 4;
    /* edge 4 */
    data->nodedghq[4][0] = 0;
    data->nodedghq[4][1] = 12;
    data->nodedghq[4][2] = 4;
    /* edge 5 */
    data->nodedghq[5][0] = 1;
    data->nodedghq[5][1] = 13;
    data->nodedghq[5][2] = 5;
    /* edge 6 */
    data->nodedghq[6][0] = 2;
    data->nodedghq[6][1] = 14;
    data->nodedghq[6][2] = 6;
    /* edge 7 */
    data->nodedghq[7][0] = 3;
    data->nodedghq[7][1] = 15;
    data->nodedghq[7][2] = 7;
  }
  /* anchors and directions */
  if (MAXEDG_THERM3 >= 12)
  {
    /* edge 0 */
    data->ancedgh[0][0] = 0.0;
    data->ancedgh[0][1] = -1.0;
    data->ancedgh[0][2] = -1.0;
    data->diredgh[0][0] = 1.0;  /* r-compo of base vector in direct. of edge */
    data->diredgh[0][1] = 0.0;  /* s-compo of base vector in direct. of edge */
    data->diredgh[0][2] = 0.0;  /* t-compo of base vector in direct. of edge */
    /* edge 1 */
    data->ancedgh[1][0] = 1.0;
    data->ancedgh[1][1] = 0.0;
    data->ancedgh[1][2] = -1.0;
    data->diredgh[1][0] = 0.0;
    data->diredgh[1][1] = 1.0;
    data->diredgh[1][2] = 0.0;
    /* edge 2 */
    data->ancedgh[2][0] = 0.0;
    data->ancedgh[2][1] = 1.0;
    data->ancedgh[2][2] = -1.0;
    data->diredgh[2][0] = -1.0;
    data->diredgh[2][1] = 0.0;
    data->diredgh[2][2] = 0.0;
    /* edge 3 */
    data->ancedgh[3][0] = -1.0;
    data->ancedgh[3][1] = 0.0;
    data->ancedgh[3][2] = -1.0;
    data->diredgh[3][0] = 0.0;
    data->diredgh[3][1] = -1.0;
    data->diredgh[3][2] = 0.0;
    /* edge 8 */
    data->ancedgh[8][0] = 0.0;
    data->ancedgh[8][1] = -1.0;
    data->ancedgh[8][2] = 1.0;
    data->diredgh[8][0] = 1.0;
    data->diredgh[8][1] = 0.0;
    data->diredgh[8][2] = 0.0;
    /* edge 9 */
    data->ancedgh[9][0] = 1.0;
    data->ancedgh[9][1] = 0.0;
    data->ancedgh[9][2] = 1.0;
    data->diredgh[9][0] = 0.0;
    data->diredgh[9][1] = 1.0;
    data->diredgh[9][2] = 0.0;
    /* edge 10 */
    data->ancedgh[10][0] = 0.0;
    data->ancedgh[10][1] = 1.0;
    data->ancedgh[10][2] = 1.0;
    data->diredgh[10][0] = -1.0;
    data->diredgh[10][1] = 0.0;
    data->diredgh[10][2] = 0.0;
    /* edge 11 */
    data->ancedgh[11][0] = -1.0;
    data->ancedgh[11][1] = 0.0;
    data->ancedgh[11][2] = 1.0;
    data->diredgh[11][0] = 0.0;
    data->diredgh[11][1] = -1.0;
    data->diredgh[11][2] = 0.0;
    /* edge 4 */
    data->ancedgh[4][0] = -1.0;
    data->ancedgh[4][1] = -1.0;
    data->ancedgh[4][2] = 0.0;
    data->diredgh[4][0] = 0.0;
    data->diredgh[4][1] = 0.0;
    data->diredgh[4][2] = 1.0;
    /* edge 5 */
    data->ancedgh[5][0] = 1.0;
    data->ancedgh[5][1] = -1.0;
    data->ancedgh[5][2] = 0.0;
    data->diredgh[5][0] = 0.0;
    data->diredgh[5][1] = 0.0;
    data->diredgh[5][2] = 1.0;
    /* edge 6 */
    data->ancedgh[6][0] = 1.0;
    data->ancedgh[6][1] = 1.0;
    data->ancedgh[6][2] = 0.0;
    data->diredgh[6][0] = 0.0;
    data->diredgh[6][1] = 0.0;
    data->diredgh[6][2] = 1.0;
    /* edge 7 */
    data->ancedgh[7][0] = -1.0;
    data->ancedgh[7][1] = 1.0;
    data->ancedgh[7][2] = 0.0;
    data->diredgh[7][0] = 0.0;
    data->diredgh[7][1] = 0.0;
    data->diredgh[7][2] = 1.0;
  }
  

  /*====================================================================*/
  /* tetrahedra */
  /*--------------------------------------------------------------------*/
  /* parameter coordinates (r,s,t) of nodes 
   * of tetrahedron { (r,s,t) | 0<=r<=1, 0<=s<=1-r, 0<=t<=1-r-s }
   *  4-node hexahedron: node 0,1,...,3 
   * 10-node hexahedron: node 0,1,...,9 */
  /*
   * drawing
   */
  if (MAXNOD_THERM3 >= 4)
  {
    /* node 0 */
    data->nodtrst[0][0] = 1.0;
    data->nodtrst[0][1] = 0.0;
    data->nodtrst[0][2] = 0.0;
    /* node 1 */
    data->nodtrst[1][0] = 0.0;
    data->nodtrst[1][1] = 1.0;
    data->nodtrst[1][2] = 0.0;
    /* node 2 */
    data->nodtrst[2][0] = 0.0;
    data->nodtrst[2][1] = 0.0;
    data->nodtrst[2][2] = 1.0;
    /* node 4 */
    data->nodtrst[3][0] = 0.0;
    data->nodtrst[3][1] = 0.0;
    data->nodtrst[3][2] = 0.0;
  }
  if (MAXNOD_THERM3 >= 10)
  {
    /* node 4 */
    data->nodtrst[4][0] = 0.5;
    data->nodtrst[4][1] = 0.0;
    data->nodtrst[4][2] = 0.0;
    /* node 5 */
    data->nodtrst[5][0] = 0.5;
    data->nodtrst[5][1] = 0.5;
    data->nodtrst[5][2] = 0.0;
    /* node 6 */
    data->nodtrst[6][0] = 0.0;
    data->nodtrst[6][1] = 0.5;
    data->nodtrst[6][2] = 0.0;
    /* node 7 */
    data->nodtrst[7][0] = 0.5;
    data->nodtrst[7][1] = 0.0;
    data->nodtrst[7][2] = 0.5;
    /* node 8 */
    data->nodtrst[8][0] = 0.0;
    data->nodtrst[8][1] = 0.5;
    data->nodtrst[8][2] = 0.5;
    /* node 9 */
    data->nodtrst[9][0] = 0.0;
    data->nodtrst[9][1] = 0.0;
    data->nodtrst[9][2] = 0.5;
  }
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of th3_cfg_init */


/*======================================================================*/
/*!
\brief Parameter coordinates of element nodes

Obtain the (r,s,t)-coordinates (ie natural coords, or the coordinates in 
parameter space)  of the element nodes (ie hexs with 8,20,27, or tets
with 4,10)

\author bborn
\date 09/06
*/
void th3_cfg_noderst(ELEMENT *ele,
		     TH3_DATA *data,
		     INT inode,
		     DOUBLE *rst)
{

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_cfg_noderst");
#endif

  /*--------------------------------------------------------------------*/
  switch (ele->distyp)
  {
    case hex8: case hex20: case hex27:
      rst[0] = data->nodhrst[inode][0];
      rst[1] = data->nodhrst[inode][1];
      rst[2] = data->nodhrst[inode][2];
      break;
    case tet4: case tet10:
      rst[0] = data->nodtrst[inode][0];
      rst[1] = data->nodtrst[inode][1];
      rst[2] = data->nodtrst[inode][2];
      break;
    default:
      dserror("Unknown discretisation type!");
  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th3_cfg_noderst */


/*======================================================================*/
/*!
\brief Test parametric geometry and topology of THERM3 element

\author bborn
\date 10/06
*/
#ifdef TEST_THERM3 

void th3_cfg_test(TH3_DATA *data)
{
  FILE *filetest;
  INT inod, jnod, idim, isid, iedg;

    /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("th3_cfg_test");
#endif

  /*--------------------------------------------------------------------*/
  /* open test file */
  filetest = fopen("th3_cfg_test.out", "w");

  /*====================================================================*/
  /* hexahedra geometry and topology */

  /*--------------------------------------------------------------------*/
  /* loop all nodes and write to file */
  for (inod=0; inod<MAXNOD_THERM3; inod++)
  {
    fprintf(filetest, "# Node %i\n", inod);
    for (idim=0; idim<NDIM_THERM3; idim++)
    {
      fprintf(filetest, " %e", data->nodhrst[inod][idim]);
    }
    fprintf(filetest, "\n");    
  }
  fprintf(filetest, "\n\n");

  /*--------------------------------------------------------------------*/
  /* loop sides */
  for (isid=0; isid<MAXSID_THERM3; isid++)
  {
    fprintf(filetest, "# Side %i\n", isid);
    for (inod=0; inod<MAXNS_THERM3; inod++)
    {
      jnod = data->nodsidh[isid][inod];
      for (idim=0; idim<NDIM_THERM3; idim++)
      {
	fprintf(filetest, " %e", data->nodhrst[jnod][idim]);
      }
      fprintf(filetest, "\n");
    }
    fprintf(filetest, "\n\n");
  }
  fprintf(filetest, "\n\n");

  /*--------------------------------------------------------------------*/
  /* loop linear edges */
  for (iedg=0; iedg<MAXEDG_THERM3; iedg++)
  {
    fprintf(filetest, "# Edge %i\n", iedg);
    for (inod=0; inod<MAXNE_THERM3-1; inod++)
    {
      jnod = data->nodedghl[iedg][inod];
      for (idim=0; idim<NDIM_THERM3; idim++)
      {
	fprintf(filetest, " %e", data->nodhrst[jnod][idim]);
      }
      fprintf(filetest, "\n");
    }
    fprintf(filetest, "\n\n");
  }
  fprintf(filetest, "\n\n");

  /*--------------------------------------------------------------------*/
  /* loop quadratic edges */
  for (iedg=0; iedg<MAXEDG_THERM3; iedg++)
  {
    fprintf(filetest, "# Edge %i\n", iedg);
    for (inod=0; inod<MAXNE_THERM3; inod++)
    {
      jnod = data->nodedghq[iedg][inod];
      for (idim=0; idim<NDIM_THERM3; idim++)
      {
	fprintf(filetest, " %e", data->nodhrst[jnod][idim]);
      }
      fprintf(filetest, "\n");
    }
    fprintf(filetest, "\n\n");
  }
  fprintf(filetest, "\n\n");
  
  /*--------------------------------------------------------------------*/
  fclose(filetest);

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end of th3_cfg_test */

#endif  /* end of #ifdef TEST_THERM3 */

/*======================================================================*/
#endif  /* end of #ifdef D_THERM3 */
