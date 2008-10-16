/*!-------------------------------------------------------------------
\file Bromotion_manager.cpp

\brief This class manages Brownian Motion including mass


<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*--------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "bromotion_manager.H"



/*-------------------------------------------------------------------*
| (public)                                                umay  09/08|
|                                                                    |
| Standard constructor                                               |
*--------------------------------------------------------------------*/
BroMotion_Manager::BroMotion_Manager(
    Teuchos::RCP<DRT::Discretization>   discretRCP,
    DRT::Discretization&                discret):
    discretRCP_(discretRCP),
    discret_(discret)
{
  
  
}



/*-------------------------------------------------------------------*
| (public)                                                umay  09/08|
|                                                                    |
| Calculate additional forces according to stochastical processes    |
*--------------------------------------------------------------------*/
void BroMotion_Manager::StochasticalForces( 
    const DRT::Element*             element,
    const DRT::UTILS::GaussRule2D&  gaussrule,
    ParameterList&                  params,
    vector<int>&                    lm,
    Epetra_SerialDenseMatrix&       K_surf,
    Epetra_SerialDenseVector&       F_int)
{

  dserror("not yet implemented");
  return;
}

#endif /*CCADISCRET*/


