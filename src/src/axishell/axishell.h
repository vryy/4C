/*!----------------------------------------------------------------------
\file
\brief headerfile for axisymmetric shell element,
       containing structures and prototypes

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_AXISHELL

/*!
  \addtogroup AXISHELL
  *//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief axisymmetric shell element

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

<pre>                                                            mn 05/03
This structure contains all specific information for the axisymmetric
shell element
</pre>

*----------------------------------------------------------------------*/
typedef struct _AXISHELL
{
  DOUBLE        thick[2];           /*!< thickness at node i and j       */

  DOUBLE        pv[2];              /*!< vertical load at node i and j   */
  DOUBLE        ph[2];              /*!< horizontal load at node i and j */
  DOUBLE        px[2];              /*!< tangential load at node i and j */
  DOUBLE        pw[2];              /*!< normal load at node i and j     */

  DOUBLE        statcond[6];        /*!< needed for static condensation  */
  DOUBLE        Sk;                 /*!< needed for static condensation  */

  struct _ARRAY4D  stress_GP;
  struct _ARRAY4D  stress_ND;
} AXISHELL;


/*!----------------------------------------------------------------------
  \brief axishell element data

  <pre>                                                            mn 05/03
  This structure contains the coordinates and weights for numerical
  integration of the axisymmetric shell element
  </pre>

 *----------------------------------------------------------------------*/
typedef struct _SAXI_DATA
{
  DOUBLE        xgr[5];  /*!< natural coordinates xsi of gaussian points */
  DOUBLE        wgt[5];  /*!< weights at natural coordinates xsi         */
} SAXI_DATA;


/*----------------------------------------------------------------------*
  |  PROTOTYPES OF ALL AXISHELL ROUTINES                    mn 05/03    |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
  | main axishell control routine                           mn 05/03    |
 *----------------------------------------------------------------------*/
void axishell(
    PARTITION   *actpart,
    INTRA       *actintra,
    ELEMENT     *ele,
    ARRAY       *estif_global,
    ARRAY       *intforce_global,
    CALC_ACTION *action,
    CONTAINER   *container
    );

/*----------------------------------------------------------------------*
  | read shell element                                       mn 05/03    |
 *----------------------------------------------------------------------*/
void saxi_inp(
    ELEMENT    *ele
    );

#endif /*D_AXISHELL*/
/*! @} (documentation module close)*/
#endif
