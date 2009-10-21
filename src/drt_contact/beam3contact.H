/*!----------------------------------------------------------------------
\file beam3contact.H

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*----------------------------------------------------------------------*/
// administrative line in order to make a difference between BACI and Ccarat files
#ifdef CCADISCRET

//compile only if beam3 element is complied, too, as beam3 element required for member variables of this class
#ifdef D_BEAM3

//header file only included if not yet included!
#ifndef BEAM3C_H
#define BEAM3C_H

#include "../drt_beam3/beam3.H"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
using namespace std;
using namespace Teuchos;
#include "../drt_lib/drt_condition.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lib/linalg_fixedsizematrix.H"



/*!
\brief contact element for contact between two elements of type beam3

\author cyron (cyron @lnm.mw.tum.de)
*/
class Beam3contact 
{
public:

  //! @name Friends
  
  //no fried classes defined
  
  //@}


  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param element1 (in): reference to first element of contact pair
  \param element2 (in): reference to second element of contact pair
  */
  Beam3contact(const DRT::ELEMENTS::Beam3& element1, const DRT::ELEMENTS::Beam3& element2);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  */
  Beam3contact(const Beam3contact& old);

  /*!
  \brief Destructor
  */
  virtual ~Beam3contact();
  
  //@}

  //! @name public evaluation methods

  /*!
  \brief Evaluate an element

  this method computes internal forces and stiffness due to contact

  \param params (in/out)    : ParameterList for communication between control routine
                              and elements
  \param discretization (in): A reference to the underlying discretization
  \param lm (in)            : location vector of this element
  \param elemat1 (out)      : matrix to be filled by stiffness entries
  \param elevec1 (out)      : vector to be filled by internal force entries

  \return 0 if successful, negative otherwise
  */
  virtual int Evaluate(ParameterList&            params,
                       DRT::Discretization&      discretization,
                       vector<int>&              lm,
                       Epetra_SerialDenseMatrix& stiffmatrix,
                       Epetra_SerialDenseVector& fint);
  
  //@}

private:
  
  //! @name member variables

  //!first element of contact pair
  DRT::ELEMENTS::Beam3 element1_;
  
  //!second element of contact pair
  DRT::ELEMENTS::Beam3 element2_;

  //@}

  //! @name private evaluation methods
  
  //! calculation of nonlinear stiffness matrix and internal forces
  void nlnstiff(ParameterList& params,
                vector<double>&           disp,  //!< element displacement vector
                Epetra_SerialDenseMatrix& stiffmatrix,  //!< element stiffness matrix
                Epetra_SerialDenseVector& fint);  //!< element internal force vector


  //@}

}; // class Beam3contact

#endif  // #ifndef BEAM3C_H
#endif  // #ifdef D_BEAM3
#endif  // #ifdef CCADISCRET 
