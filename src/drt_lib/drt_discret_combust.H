/*!----------------------------------------------------------------------
\file drt_discret_combust.H

\brief A class to manage specialized discretizations for combustion problems

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef DISCRET_COMBUST_H
#define DISCRET_COMBUST_H

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
class DiscretizationCombust : public DiscretizationFaces
{
public:

  /*!
  \brief Standard Constructor

  \param name (in): name of this discretization
  \param comm (in): An epetra comm object associated with this discretization
  */
  DiscretizationCombust(const std::string name, Teuchos::RCP<Epetra_Comm> comm)
  :DiscretizationFaces(name, comm) { extension_filled_ = false; }; // use base class constructor

  /*!
  \brief Destructor

  */
  virtual ~DiscretizationCombust()
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
       Teuchos::ParameterList&              params,          ///< parameter list
       DRT::DiscretizationCombust&             discretization,  ///< XFEM discretization
       Teuchos::RCP<LINALG::SparseMatrix>   systemmatrix,    ///< systemmatrix
       Teuchos::RCP<Epetra_Vector>          systemvector     ///< systemvector
       );

}; // class DiscretizationCombust
} // namespace DRT

#endif  // #ifndef DISCRET_XFEM_H
