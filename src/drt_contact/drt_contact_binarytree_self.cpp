/*!----------------------------------------------------------------------
\file drt_contact_binarytree_self.cpp
\brief A class for performing contact search in 2D / 3D based
       on binary search trees

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich
              
Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed, 
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de) 
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de                   

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_contact_binarytree_self.H"
#include "drt_cnode.H"
#include "drt_celement.H"
#include "contactdefines.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/linalg_fixedsizematrix.H"

                     
/*----------------------------------------------------------------------*
 |  ctor BinaryTreeSelf (public)                              popp 11/09|
 *----------------------------------------------------------------------*/
CONTACT::BinaryTreeSelf::BinaryTreeSelf(DRT::Discretization& discret,
                                RCP<Epetra_Map> selements,
                                RCP<Epetra_Map> melements,
                                int dim, double eps) :
idiscret_(discret),
selements_(selements),
melements_(melements),
dim_(dim),
eps_(eps)
{  
  // not yet integrated into BACI
  cout << "\n\nEntering SELF CONTACT binary tree constructor!" << endl;
  dserror("ERROR: Self contact not yet available!");
  
  return;
}

/*----------------------------------------------------------------------*
 | Find min. length of master and slave elements (public)     popp 11/08|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelf::SetEnlarge(bool isinit)
{
  // not yet implemented
  return;
}
                                            
/*----------------------------------------------------------------------*
 | Search for self contact (public)                           popp 11/09|
 *----------------------------------------------------------------------*/
void CONTACT::BinaryTreeSelf::SearchContactCombined()
{
  // not yet implemented
  return;
}

#endif //#ifdef CCADISCRET
