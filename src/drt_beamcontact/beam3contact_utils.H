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
#include "beam3contact_defines.H"
#include "../drt_fem_general/largerotations.H"

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
  \brief Cast of a FAD matrix to a double matrix
  */
  template<typename type,int dim1, int dim2>
  LINALG::TMatrix<double,dim1,dim2> CastToDouble(LINALG::TMatrix<type,dim1,dim2> a)
  {
    LINALG::TMatrix<double,dim1,dim2> b(true);

    for (int i=0;i<dim1;i++)
    {
      for (int j=0;j<dim2;j++)
      {
        b(i,j)=CastToDouble(a(i,j));
      }
    }
    return b;
  }

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
  \brief Calculate vector product of two FAD or double vectors
  */
  template<typename type>
  LINALG::TMatrix<type,3,1> VectorProduct(LINALG::TMatrix<type,3,1> first_vector, LINALG::TMatrix<type,3,1> second_vector)
  {
    LINALG::TMatrix<type,3,1> result_vector;
    result_vector.Clear();
    LINALG::TMatrix<type,3,3> S_first_vector;
    S_first_vector.Clear();
    LARGEROTATIONS::computespin(S_first_vector,first_vector);

    result_vector.Multiply(S_first_vector, second_vector);

    return result_vector;
  }

  /*!
  \brief Calculate angle encompassed by two lines: returns an angle \in [0;pi/2]
  */
  double CalcAngle(LINALG::TMatrix<double,3,1> a, LINALG::TMatrix<double,3,1> b);

  /*!
  \brief Get closest distance between the endpoints of two lines
  */
  template<typename type>
  type GetClosestEndpointDist(LINALG::TMatrix<type,3,1> r1_a, LINALG::TMatrix<type,3,1> r1_b, LINALG::TMatrix<type,3,1> r2_a, LINALG::TMatrix<type,3,1> r2_b)
  {
    type minnodaldist=0.0;
    type nodaldist=0.0;

    minnodaldist=BEAMCONTACT::VectorNorm<3>(BEAMCONTACT::DiffVector(r1_a,r2_a));

    nodaldist=BEAMCONTACT::VectorNorm<3>(BEAMCONTACT::DiffVector(r1_a,r2_b));
    if(nodaldist<minnodaldist)
      minnodaldist=nodaldist;

    nodaldist=BEAMCONTACT::VectorNorm<3>(BEAMCONTACT::DiffVector(r1_b,r2_a));
    if(nodaldist<minnodaldist)
      minnodaldist=nodaldist;

    nodaldist=BEAMCONTACT::VectorNorm<3>(BEAMCONTACT::DiffVector(r1_b,r2_b));
    if(nodaldist<minnodaldist)
      minnodaldist=nodaldist;

    return minnodaldist;
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
  \brief Check, if current node belongs to a rigid sphere element
  */
  bool RigidsphereNode(DRT::Node& node);

  /*!
  \brief Check, if current node is a beam element
  */
  bool BeamElement(DRT::Element& element);

  /*!
  \brief Check, if current node is a beam element
  */
  bool RigidsphereElement(DRT::Element& element);

  /*
  \brief Check, if two elements share a node -> neighbor elements
  */
  bool ElementsShareNode(DRT::Element& element1,DRT::Element& element2);

  /*
  \brief Calculate beam radius
  */
  double CalcEleRadius(const DRT::Element* ele);

  /*
  \brief Intersect two parallel cylinders
  */
  bool IntersectParallelCylinders(LINALG::TMatrix<double,3,1>& r1_a,
                                  LINALG::TMatrix<double,3,1>& r1_b,
                                  LINALG::TMatrix<double,3,1>& r2_a,
                                  LINALG::TMatrix<double,3,1>& r2_b,
                                  double& distancelimit);

  /*
  \brief Intersect two non-parallel, arbitrary oriented cylinders
  */
  bool IntersectArbitraryCylinders(LINALG::TMatrix<double,3,1>& r1_a,
                                  LINALG::TMatrix<double,3,1>& r1_b,
                                  LINALG::TMatrix<double,3,1>& r2_a,
                                  LINALG::TMatrix<double,3,1>& r2_b,
                                  double& distancelimit,
                                  std::pair<double,double>& closestpoints, // The closest point are only set, if we have detected an intersection at a valid closest point with eta1_seg, eta2_seg \in [-1.0;1.0]
                                  bool etaset); //bool to check, if the closest point coordinates have been set or not

  /*
  \brief Calculate closest distance of a point and a line
  */
  double CalcPointLineDist( LINALG::TMatrix<double,3,1>& rline_a,
                            LINALG::TMatrix<double,3,1>& rline_b,
                            LINALG::TMatrix<double,3,1>& rp,
                            double& eta);

  /*
  \brief Determine inpute parameter representing the additive searchbox increment
  */
  double DetermineSearchboxInc(Teuchos::ParameterList& beamcontactparams);

} // namespace BEAMCONTACT

#endif  // #ifndef BEAM3CONTACT_UTILS_H
