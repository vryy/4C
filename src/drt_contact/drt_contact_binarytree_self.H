/*!----------------------------------------------------------------------
\file drt_contact_binarytree_self.H

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
#ifndef DRT_CONTACT_BINARYTREE_SELF_H
#define DRT_CONTACT_BINARYTREE_SELF_H

#include "../drt_lib/drt_discret.H"
#include "drt_celement.H"

using namespace std;

/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{

// forward declarations

/*!
\brief A class for performing self contact search in 2D / 3D based
       on a binary search tree and dual graphs
       
Refer also to the Semesterarbeit of Anh-Tu Vuong, 2009

\author popp (popp@lnm.mw.tum.de)
*/

class BinaryTreeSelf
{
public:
  
  /*!
  \brief Standard constructor
  
  Constructs an instance of this class.<br>
  For now, we only consider the serial case (1 processor) here!!!
  
  \param discret (in):    The contact interface discretization
  \param selements (in):  All slave elements (fully overlapping map)
  \param melements (in):  All master elements (fully overlapping map)
  \param dim (in):        The problem dimension 
  \param eps (in):        Parameter for enlarging dops (user-defined)
  
  */
  BinaryTreeSelf(DRT::Discretization& discret,
                 RCP<Epetra_Map> selements,
                 RCP<Epetra_Map> melements,
                 int dim, double eps);

  /*!
  \brief Destructor

  */
  virtual ~BinaryTreeSelf() {};
  
  //! @name Evaluation methods
 
  /*!
  \brief Calculate minimal element length / inflation factor

  */
  void SetEnlarge(bool isinit);

  /*!
  \brief Search for self contact

  */
  void SearchContactCombined();
  
  //@}
      
protected:
    
  // don't want = operator and cctor
  BinaryTreeSelf operator = (const BinaryTreeSelf& old);
  BinaryTreeSelf(const BinaryTreeSelf& old);

  DRT::Discretization&               idiscret_;       // contact interface discretization
  RCP<Epetra_Map>                    selements_;      // all slave elements on surface (column map)
  RCP<Epetra_Map>                    melements_;      // all master elements on surface (full map)
  const int                          dim_;            // problem dimension (2D or 3D)
  double                             eps_;            // epsilon for enlarging dops (user-defined)
     
}; // class BinaryTreeSelf
}  // namespace CONTACT


#endif  // #ifndef DRT_CONTACT_BINARYTREE_SELF_H
#endif  // #ifdef CCADISCRET
