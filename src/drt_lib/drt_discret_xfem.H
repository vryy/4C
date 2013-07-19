/*!----------------------------------------------------------------------
\file drt_discret_xfem.H

\brief a class to manage an XFEM discretization

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

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
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifndef DISCRET_XFEM_H
#define DISCRET_XFEM_H

#include "drt_discret_faces.H"
#include "drt_dserror.H"

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include <string>
#include "Teuchos_ParameterList.hpp"

namespace LINALG
{
  class MapExtractor;
  class SparseMatrix;
}

namespace DRT
{
namespace ELEMENTS
{
 class Combust3IntFace;
}
}

/*!
\brief DRT: namespace of the baci discretization module

*/
namespace DRT
{

/*!
\brief A class to manage a discretization in parallel

\author henke (henke@lnm.mw.tum.de)
*/
class DiscretizationXFEM : public DiscretizationFaces
{
public:

  /*!
  \brief Standard Constructor

  \param name (in): name of this discretization
  \param comm (in): An epetra comm object associated with this discretization
  */
  DiscretizationXFEM(const string name, RCP<Epetra_Comm> comm)
  :DiscretizationFaces(name, comm) { extension_filled_ = false; }; // use base class constructor

  /*!
  \brief Destructor

  */
  virtual ~DiscretizationXFEM()
  {
    return;
  };

  /*!
  \brief Evaluate Dirichlet boundary conditions for COMBUST problems

  this method is used instead of the base class method for XFEM discretizations

  \author henke 07/09
  */
  virtual void EvaluateDirichletCombust(Teuchos::ParameterList& params,
      Teuchos::RCP<Epetra_Vector> systemvector,
      Teuchos::RCP<Epetra_Vector> systemvectord,
      Teuchos::RCP<Epetra_Vector> systemvectordd,
      Teuchos::RCP<Epetra_Vector> toggle,
      Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor=Teuchos::null
  );


  /*!
  \brief Compute the nullspace of the discretization

  This method looks in the solver parameters whether algebraic multigrid (AMG)
  is used as preconditioner. AMG desires the nullspace of the
  system of equations which is then computed here if it does not already exist
  in the parameter list.

  \note This method is supposed to go away and live somewhere else soon....

  \param solveparams (in): List of parameters
  \param recompute (in)  : force method to recompute the nullspace
  */
  virtual void ComputeNullSpaceIfNecessary(Teuchos::ParameterList& solveparams,
                                           bool recompute = false)
  {
    // remark: the null space is not computed correctly for XFEM discretizations, since the number of
    //         degrees of freedom per node is not fixed.
    //         - it is not clear what happens with respect to the Krylov projection
    //           (having XFEM dofs seems to render the system non-singular, but it should be singular
    //           so the null space has a non-zero dimension)
    //         - the ML preconditioner also relies on a fixed number of dofs per node
    DRT::Discretization::ComputeNullSpaceIfNecessary(solveparams, recompute);
  }

  /*!
  \brief Evaluate edge-based integrals
  */
  void EvaluateEdgeBasedCombust(
       Teuchos::ParameterList&              params,
       Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
       Teuchos::RCP<Epetra_Vector>          systemvector1
  );
  void EvaluateEdgeBasedCombust(
       DRT::ELEMENTS::Combust3IntFace*      ele,             ///< internal face element
       std::vector<int>&                    nds_master,      ///< nodal dofset w.r.t. master element
       std::vector<int>&                    nds_slave,       ///< nodal dofset w.r.t. slave element
       Teuchos::ParameterList&              params,          ///< parameter list
       DRT::DiscretizationXFEM&             discretization,  ///< XFEM discretization
       Teuchos::RCP<LINALG::SparseMatrix>   systemmatrix,    ///< systemmatrix
       Teuchos::RCP<Epetra_Vector>          systemvector     ///< systemvector
       );

}; // class DiscretizationXFEM
} // namespace DRT


#endif  // #ifndef DISCRET_XFEM_H
