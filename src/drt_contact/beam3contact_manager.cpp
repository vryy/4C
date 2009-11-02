/*!-----------------------------------------------------------------------------------------------------------
\file beam3contact_manager.cpp
\brief

<pre>
Maintainer: Alexander Popp, Christian Cyron
            {popp,cyron}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*-----------------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET

//compile only if beam3 element is complied, too, as beam3 element required for member variables of this class
#ifdef D_BEAM3

#include "beam3contact_manager.H"

/*----------------------------------------------------------------------*
 |  constructor (public)                                      popp 11/09|
 *----------------------------------------------------------------------*/
CONTACT::Beam3cmanager::Beam3cmanager(const DRT::Discretization& discret):
discret_(discret)
{
  // still empty at the moment
  return;
}

/*-------------------------------------------------------------------- -*
 |  evaluate contact (public)                                 popp 11/09|
 *----------------------------------------------------------------------*/
int CONTACT::Beam3cmanager::Evaluate()
{
  // contact search
  SearchContact();
  
  // call all element pairs to evaluate f_c and k_c
  // (...)

  return 0;
}

/*-------------------------------------------------------------------- -*
 |  search contact (public)                                   popp 11/09|
 *----------------------------------------------------------------------*/
int CONTACT::Beam3cmanager::SearchContact()
{
  // loop over all elements and find closest pairs
  // (...)

  return 0;
}

#endif  // #ifdef D_BEAM3
#endif  // #ifdef CCADISCRET 
