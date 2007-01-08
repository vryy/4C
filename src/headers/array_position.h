/*======================================================================*/
/*!
\file
\brief

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>
*/


/*----------------------------------------------------------------------*/
/*!
\brief Named positions (indices) of NODE sol array in structures

\author bborn
\date 10/06
*/
typedef struct _ARRAY_POSITION_SOL
{
  /* common names */
  /* nil */
  
  /* specifically structural names */
  INT disn;  /* displacement at t_{n+1}^k */
  INT veln;  /* velocity at t_{n+1}^k */
  INT accn;
  INT disd;
  INT disdn;
  INT disdi;
  INT veldn;
  INT accdn;
  INT dis;
  INT vel;
  INT acc;

  /* specifically thermal names */
  INT tem;
} ARRAY_POSITION_SOL;


/*----------------------------------------------------------------------*/
/*!
\brief Named positions (indices) of NODE sol_increment array in structures

\author bborn
\date 10/06
*/
typedef struct _ARRAY_POSITION_SOLINC
{
  /* common names */
  /* nil */
  
  /* specifically structural names */
  INT disinc;
  INT fint;  /* internal force at t_n */
  INT fintn;  /* internal force at t_{n+1} */

  /* specifically thermal names */
  /* nil */
} ARRAY_POSITION_SOLINC;


/*----------------------------------------------------------------------*/
/*!
\brief Named positions (indices) of NODE sol_residual array in structures

\author bborn
\date 10/06
*/
typedef struct _ARRAY_POSITION_SOLRES
{
  /* common names */
  /* nil */
  
  /* specifically structural names */
  INT disres;

  /* specifically thermal names */
  /* nil */
} ARRAY_POSITION_SOLRES;


/*----------------------------------------------------------------------*/
/*!
\brief positions of physical values in node arrays

This structure contains the positions of the various fluid solutions
within the nodal array of sol_increment.a.da[pos][dim].

For the incremental acceleration gen-alpha scheme, the meaning of the
variables change. The new meaning is explained in the comments to
fluid_incr_acc_gen_alpha.                                    (gammi)

\author chfoe
\date 11/04
*/
typedef struct _ARRAY_POSITION
{
  INT num;  /*!< 1st dimension of solution field sol */
  INT numincr; /*!< number of solution fields within sol_increment (fluid)*/
  INT nummf;
  INT numres;  /*! < 1st dimension of solution field sol_residual */
#ifdef D_FLUID2_TDS
  INT accnp; /*!< most recent iteration value of new accel. at time n+1  */
#endif /* D_FLUID2_TDS */    
  INT accn;  /*!< position of sol_increment occupied by accel. at time n */
  INT accnm; /*!< position of sol_increment occup. by accel. at time n-1 */
  INT convn; /*!< position of convective velocity at n in sol vectors    */
  INT convnp; /*!< position of convective velocity at n+1 in sol vectors */
  INT dispn;
  INT dispnp;
  INT ddisp;
  INT eddy;
  INT gridv; /*!< position of grid velocity in solution vector order     */
  INT hist;  /*!< pos. of lin. comb. of hist. values needed for mass rhs */
  INT pred;  /*!< position of sol_increment occypied by predicted vels.  */
  INT relax; /*!< position of relaxation parameter solution in sol vec's */
  INT stresspro; /*! position to write projected stresses to */
  INT terr;  /*!< position of sol_increment occypied by truncation error */
  INT veln; /*!< position of sol_increment occupied by velocity at time n*/
  INT velnm; /*!< position of sol_increment occupied by vel. at time n-1 */
  INT velnp; /*!< position of sol_increment occupied by vel. at time n+1 */

  INT residuum;

  INT mf_dispn;
  INT mf_dispnp;
  INT mf_dispnm;
  INT mf_dispi;
  INT mf_posnp;
  INT mf_forcenp;
  INT mf_forcen;
  INT mf_forcecpy;
  INT mf_velnp;
  INT mf_velcpy;
  INT mf_reldisp;
  INT mf_sd_g;

  ARRAY_POSITION_SOL isol;  /* positions in sol array */
  ARRAY_POSITION_SOLINC isolinc;  /* positions in sol_increment array */
  ARRAY_POSITION_SOLRES isolres;  /* positions in sol_residual array */
} ARRAY_POSITION;
