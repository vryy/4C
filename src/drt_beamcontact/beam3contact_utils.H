/*!----------------------------------------------------------------------
\file beam3contact_utils.H

\brief A set of utility functions for beam contact

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
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*----------------------------------------------------------------------*/
#ifndef BEAM3CONTACT_UTILS_H
#define BEAM3CONTACT_UTILS_H

#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "../linalg/linalg_fixedsizematrix.H"
#include "Sacado.hpp"
#include "../drt_lib/drt_discret.H"

typedef Sacado::Fad::DFad<double> FAD;

namespace BEAMCONTACT
{

  /*!
  \brief Cast of FAD to double
  */
  double CastToDouble(FAD a);

  /*!
  \brief Cast of double to double
  */
  double CastToDouble(double a);

  /*!
  \brief Calculate Norm of a scalar FAD quantity
  */
  double Norm(double a);

  /*!
  \brief Calculate Norm of a scalar double quantity
  */
  FAD Norm(FAD a);

  /*!
  \brief Calculate Norm of a double vector
  */
  template<int length>
  double VectorNorm(LINALG::TMatrix<double,length,1> v)
  {
    return v.Norm2();
  }

  /*!
  \brief Calculate Norm of a FAD vector
  */
  template<int length>
  FAD VectorNorm(LINALG::TMatrix<FAD,length,1> v)
  {
    FAD norm = 0.0;
    for (int i=0;i<length;i++)
    {
      norm+=v(i)*v(i);
    }
    norm = pow(norm,0.5);
    return norm;
  }

  /*!
  \brief Calculate inner product of two FAD or double vectors
  */
  template<typename type>
  type ScalarProduct(LINALG::TMatrix<type,3,1> a, LINALG::TMatrix<type,3,1> b)
  {
    return a(0)*b(0)+a(1)*b(1)+a(2)*b(2);
  }

  /*!
  \brief Calculate difference of two FAD or double vectors
  */
  template<typename type>
  LINALG::TMatrix<type,3,1> DiffVector(LINALG::TMatrix<type,3,1> a, LINALG::TMatrix<type,3,1> b)
  {
    LINALG::TMatrix<type,3,1> c(true);
    for (int i=0;i<3;i++)
      c(i)=a(i)-b(i);

    return c;
  }

  /*!
  \brief Calculate signum function of FAD or double quantity
  */
  template<typename type>
  type Signum(type a)
  {
    if (CastToDouble(a) >= 0.0)
      return 1.0;
    else
      return -1.0;
  }

  /*!
  \brief Set primary displacement DoFs for automatic differentiation with Sacado
  */
  template<int numnodes, int numnodalvalues>
  void SetFADDispDofs(LINALG::TMatrix<FAD,3*numnodes*numnodalvalues,1>& ele1pos_, LINALG::TMatrix<FAD,3*numnodes*numnodalvalues,1>& ele2pos_)
  {
    // The 2*3*numnodes*numnodalvalues primary DoFs are the components of the nodal positions / tangents. The two (+2) additional
    // degrees of freedom represent the dependency on the parameter coordinates xi and eta, which is necessary in beam contact.
    for (int i=0;i<3*numnodes*numnodalvalues;i++)
      ele1pos_(i).diff(i,2*3*numnodes*numnodalvalues+2);

    for (int i=0;i<3*numnodes*numnodalvalues;i++)
      ele2pos_(i).diff(3*numnodes*numnodalvalues+i,2*3*numnodes*numnodalvalues+2);

    return;
  }

  /*!
  \brief Set primary parameter coordinate DoFs for automatic differentiation with Sacado
  */
  template<int numnodes, int numnodalvalues>
  void SetFADParCoordDofs(FAD& xi, FAD& eta)
  {
    // The 2*3*numnodes*numnodalvalues primary DoFs are the components of the nodal positions / tangents. The two (+2) additional
    // degrees of freedom represent the dependency on the parameter coordinates xi and eta, which is necessary in beam contact.
    xi.diff((2*3*numnodes*numnodalvalues+1)-1,2*3*numnodes*numnodalvalues+2);
    eta.diff((2*3*numnodes*numnodalvalues+2)-1,2*3*numnodes*numnodalvalues+2);

    return;
  }

  /*!
  \brief Check, if current node belongs to a beam element
  */
  bool BeamNode(DRT::Node& node);

  /*!
  \brief Check, if current node is a beam element
  */
  bool BeamElement(DRT::Element& element);

  /*
  \brief Check, if two elements share a node -> neighbor elements
  */
  bool ElementsShareNode(DRT::Element& element1,DRT::Element& element2);

} // namespace BEAMCONTACT

#endif  // #ifndef BEAM3CONTACT_UTILS_H
