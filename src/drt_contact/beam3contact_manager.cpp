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
#include "beam3contact.H"

/*----------------------------------------------------------------------*
 |  constructor (public)                                      popp 11/09|
 *----------------------------------------------------------------------*/
CONTACT::Beam3cmanager::Beam3cmanager(const DRT::Discretization& discret):
discret_(discret)
{
  // initialize contact element pairs
  pairs_.resize(0);
  
#ifndef D_BEAM3
  dserror("ERROR: Beam3 contact manager called without D_BEAM3 activated");
#endif
  
  // print discretization
  //Print(cout);
  
  return;
}

/*----------------------------------------------------------------------*
 |  print beam3 contact manager (public)                      popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3cmanager::Print(ostream& os) const
{
  if (Discret().Comm().MyPID()==0)
    os << "Beam3 Contact Discretization:" << endl;
  
  os << Discret();
  return;
}

/*-------------------------------------------------------------------- -*
 |  evaluate contact (public)                                 popp 11/09|
 *----------------------------------------------------------------------*/
int CONTACT::Beam3cmanager::Evaluate(LINALG::SparseMatrix& stiffc, Epetra_Vector& fc)
{ 
  // contact search
  // loop over all elements and find closest pairs

  // (...)
  
  // pairs will be filled now
  
  // call all element pairs to evaluate f_c and k_c
  
  // (...)

  // stiffc and fc will be filled now

  return 0;
}

#endif  // #ifdef D_BEAM3
#endif  // #ifdef CCADISCRET 
