/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of algebraic assemble methods for namespace LINALG

\level 0
\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/

#ifndef SRC_LINALG_LINALG_UTILS_SPARSE_ALGEBRA_ASSEMBLE_H_
#define SRC_LINALG_LINALG_UTILS_SPARSE_ALGEBRA_ASSEMBLE_H_

#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_exporter.H"
#include "linalg_blocksparsematrix.H"

namespace LINALG
{
  // forward declaration
  // class Solver;

  /*!
   \brief Assemble an Epetra_SerialDenseMatrix into an Epetra_CrsMatrix

   This is an individual call.
   Will only assemble locally and will never do any commmunication.
   All values that can not be assembled locally will be ignored.
   Will use the communicator and rowmap from matrix A to determine ownerships.
   Local matrix Aele may be \b square or \b rectangular.

   This version of 'Assemble' does not work for a matrix A that is already
   Filled()! If matrix A is not Filled(), it will be enlarged as required.

   \note The user must provide an \b additional input vector 'lmcol'
   containing the column gids for assembly seperately!

   \param A (out)         : Sparse matrix to be assembled on
   \param Aele (in)       : dense matrix to be assembled
   \param lmrow (in)      : vector with row gids
   \param lmrowowner (in) : vector with owner procs of row gids
   \param lmcol (in)      : vector with column gids
   */
  void Assemble(Epetra_CrsMatrix& A, const Epetra_SerialDenseMatrix& Aele,
      const std::vector<int>& lmrow, const std::vector<int>& lmrowowner,
      const std::vector<int>& lmcol);

  /*!
   \brief Assemble an Epetra_SerialDenseVector into an Epetra_Vector

   This is an individual call.
   Will only assemble locally and will never do any commmunication.
   All values that can not be assembled locally will be ignored.
   Will use the communicator from vector V to determine ownerships.

   \param V (out)   : Vector to be assembled on
   \param Vele (in) : dense vector to be assembled
   \param lm (in) : vector with gids
   \param lmowner (in) : vector with owner procs of gids
   */
  void Assemble(Epetra_Vector& V, const Epetra_SerialDenseVector& Vele, const std::vector<int>& lm,
      const std::vector<int>& lmowner);

  /*!
   \brief Assemble a LINALG::Matrix<3,1> into an Epetra_Vector

   This is a wrapper that converts LINALG::Matrix<3,1> to
   Epetra_SerialDenseVector (view version)

   \param V (out)   : Vector to be assembled on
   \param Vele (in) : dense vector to be assembled
   \param lm (in) : vector with gids
   \param lmowner (in) : vector with owner procs of gids
   */
  void Assemble(Epetra_Vector& V, LINALG::Matrix<3, 1>& Vele, const std::vector<int>& lm,
      const std::vector<int>& lmowner);

  /*!
   \brief Assemble a LINALG::Matrix<3,1> into an Epetra_Vector in case the
   owner is always the same

   This is a wrapper that converts LINALG::Matrix<3,1> to
   Epetra_SerialDenseVector (view version)

   \param V (out)   : Vector to be assembled on
   \param Vele (in) : dense vector to be assembled
   \param lm (in) : vector with gids
   \param lmowner (in) : vector with owner procs of gids
   */
  void Assemble(
      Epetra_Vector& V, LINALG::Matrix<3, 1>& Vele, const std::vector<int>& lm, const int& lmowner);

  /*!
   \brief Assemble a double into an Epetra_Vector

   This is an wrapper.

   \param V (out)   : Vector to be assembled on
   \param Vele (in) : dense vector to be assembled
   \param lm (in) : vector with gids
   \param lmowner (in) : vector with owner procs of gids
   */
  void Assemble(Epetra_Vector& V, double& Vele, const int& lm, const int& lmowner);

  /*!
   \brief Assemble a Epetra_SerialDenseVector into a Epetra_MultiVector

   This is an individual call.
   Will only assemble locally and will never do any commmunication.
   All values that can not be assembled locally will be ignored.
   Will use the communicator from vector V to determine ownerships.

   \param V (out)   : Vector to be assembled on
   \param n (in)   : column index of MultiVector to be assembled on
   \param Vele (in) : dense vector to be assembled
   \param lm (in) : vector with gids
   \param lmowner (in) : vector with owner procs of gids
   */
  void Assemble(Epetra_MultiVector& V, const int n, const Epetra_SerialDenseVector& Vele,
      const std::vector<int>& lm, const std::vector<int>& lmowner);

  /*! \brief Assemble a source Epetra_Vector into a target Epetra_Vector
   *
   *  The map of the source vector has to be a sub-map of the target vector and
   *  the maps must have the same processor distribution. This method does not
   *  build up any communication between different processors!
   *  The entries of the source vector are added to the target vector:
   *
   *      target_vector[GID] = source_scalar * source[GID] + target_scalar * target[GID]
   *
   *  \note The remaining GIDs, which are no part of the source map stay untouched!
   *
   *  \param scalar_target (in) : scale the target entries by this factor
   *  \param target        (out): target vector (part of the source)
   *  \param scalar_source (in) : scale the source entries by this factor
   *  \param source        (in) : source vector
   *
   *  \author hiermeier \date 03/17 */
  void AssembleMyVector(double scalar_target, Epetra_Vector& target, double scalar_source,
      const Epetra_Vector& source);

  /*!
   \brief Apply dirichlet boundary condition to a linear system of equations

   Modifies a system of equations such that dirichlet boundary conditions are enforced.
   Prescribed dirichlet BC values are supplied in dbcval and dbctoggle, where
   a prescribed value is dbcval[i] and dbctoggle[i] = 1.0. No BC is enforced in
   all places where dbctoggle[i] != 1.0.<br>
   Let us denote the \f$ A_{2 \times 2} \f$ blocks of \f$A\f$ by
   \f$A_{ff}, A_{fD}, A_{Df}, A_{DD}\f$, where \f$f\f$ stands for 'free' and
   \f$D\f$ stands for 'Dirichlet BC'. Then, after a call to this method<br>

   \f$ A_{ff} = A_{ff}, \f$<br>
   \f$ A_{fD} = A_{fD}, \f$<br>
   \f$ A_{Df} = 0_{Df}, \f$<br>
   \f$ A_{DD} = I_{DD}, \f$<br>
   \f$ x_{D} = dbcval_{D}, \f$<br>
   \f$ b_{D} = dbcval_{D} \f$<br>

   and

   \f$ A_{ff} x_f + A_{fD} x_D = b_f \f$<br>
   \f$ 0 x_f + I_{DD} x_D = x_D \f$.<br>

   \note The matrix is then nonsymmetric. When using iterative methods on this
   linear system of equations that depend on the symmetry of the matrix (such as e.g. CG),
   the initial guess supplied to the solver has to be exact at the
   Dirichlet BCs. This should be easy, as the values at the Dirichlet BCs
   are known.

   \note The mask of matrix \f$A\f$ is not modified. That is the
   entries in \f$A_{Df}\f$ and \f$A_{DD}\f$ are set to zero, not
   removed. This way the matrix can be reused in the next step.

   \param A         (in/out) : Matrix of Ax=b
   \param x         (in/out) : initial guess vector x of Ax=b
   \param b         (in/out) : rhs vector b of Ax=b
   \param dbcval    (in)     : vector holding prescribed dirichlet values
   \param dbctoggle (in)     : vector holding 1.0 where dirichlet should be applied
   and 0.0 everywhere else
     */
  void ApplyDirichlettoSystem(Teuchos::RCP<SparseOperator> A, Teuchos::RCP<Epetra_Vector>& x,
      Teuchos::RCP<Epetra_Vector>& b, const Teuchos::RCP<const Epetra_Vector> dbcval,
      const Teuchos::RCP<const Epetra_Vector> dbctoggle);

  /*!
   \brief Apply dirichlet boundary condition to a linear system of equations

   This is a flexible routine. The vectors x and dbcval might have different
   maps. The map does not need to contain all Dirichlet dofs.

   The purpose is to set Dirichlet values at a subset of all Dirichlet
   boundaries.

   \param A (in/out)         : Matrix of Ax=b
   \param x (in/out)         : vector x of Ax=b
   \param b (in/out)         : vector b of Ax=b
   \param dbcval (in)        : vector holding values that are supposed to be prescribed
   \param dbcmap (in)        : unique map of all dofs that should be constrained

   \pre The map dbcmap must be subset of the maps of the vectors.
   */
  void ApplyDirichlettoSystem(Teuchos::RCP<LINALG::SparseOperator> A,
      Teuchos::RCP<Epetra_Vector>& x, Teuchos::RCP<Epetra_Vector>& b,
      const Teuchos::RCP<const Epetra_Vector>& dbcval, const Epetra_Map& dbcmap);

  /*!
   \brief Apply dirichlet boundary condition to a linear system of equations

   This is a flexible routine. The vectors x and dbcval might have different
   maps. The map does not need to contain all Dirichlet dofs.

   The purpose is to set Dirichlet values at a subset of all Dirichlet
   boundaries.

   Special in this routine is the ability to insert rows of general rotation
   matrices (stored in #trafo) rather than simply put ones and zeros
   at the rows associated Dirichlet DOFs.

   \param A (in/out)         : Matrix of Ax=b
   \param x (in/out)         : vector x of Ax=b
   \param b (in/out)         : vector b of Ax=b
   \param trafo (in)         : global matrix holding rotation matrices to convert
   from global to local co-ordinate systems
   \param dbcval (in)        : vector holding values that are supposed to be prescribed
   \param dbcmap (in)        : unique map of all dofs that should be constrained

   \pre The map dbcmap must be subset of the maps of the vectors.
   */
  void ApplyDirichlettoSystem(Teuchos::RCP<LINALG::SparseOperator> A,
      Teuchos::RCP<Epetra_Vector>& x, Teuchos::RCP<Epetra_Vector>& b,
      Teuchos::RCP<const LINALG::SparseMatrix> trafo,
      const Teuchos::RCP<const Epetra_Vector>& dbcval, const Epetra_Map& dbcmap);

  /*!
   \brief Apply dirichlet boundary condition to a linear system of equations

   This is a flexible routine. The vectors b and dbcval might have different
   maps. The map does not need to contain all Dirichlet dofs.

   The purpose is to set Dirichlet values at a subset of all Dirichlet
   boundaries.

   Special in this routine is the ability to insert rows of general rotation
   matrices (stored in #trafo) rather than simply put ones and zeros
   at the rows associated Dirichlet DOFs.

   \param A (in/out)         : Matrix of Ax=b
   \param b (in/out)         : vector b of Ax=b
   \param trafo (in)         : global matrix holding rotation matrices to convert
   from global to local co-ordinate systems
   \param dbcval (in)        : vector holding values that are supposed to be prescribed
   \param dbcmap (in)        : unique map of all dofs that should be constrained

   \pre The map dbcmap must be subset of the maps of the vectors.
   */
  void ApplyDirichlettoSystem(Teuchos::RCP<LINALG::SparseOperator> A,
      Teuchos::RCP<Epetra_Vector>& b, Teuchos::RCP<const LINALG::SparseMatrix> trafo,
      const Teuchos::RCP<const Epetra_Vector>& dbcval, const Epetra_Map& dbcmap);

  /*!
   \brief Apply dirichlet boundary condition to a linear system of equations


   \param x (in/out)         : vector x of Ax=b
   \param b (in/out)         : vector b of Ax=b
   \param dbcval (in)        : vector holding values that are supposed to be prescribed
   \param dbctoggle (in)     : vector holding 1.0 where dirichlet should be applied
   and 0.0 everywhere else
   */
  void ApplyDirichlettoSystem(Teuchos::RCP<Epetra_Vector>& x, Teuchos::RCP<Epetra_Vector>& b,
      const Teuchos::RCP<const Epetra_Vector> dbcval,
      const Teuchos::RCP<const Epetra_Vector> dbctoggle);

  /*!
   \brief Apply dirichlet boundary condition to a linear system of equations

   This is a flexible routine. The vectors x and dbcval might have different
   maps. The dbcmap does not need to contain all Dirichlet dofs, but the vectors
   all dofs defined in it.

   The purpose is to set Dirichlet values at a subset of all Dirichlet
   boundaries.

   \param x (in/out)         : vector x of Ax=b
   \param b (in/out)         : vector b of Ax=b
   \param dbcval (in)        : vector holding values that are supposed to be prescribed
   \param dbcmap (in)        : unique map of all dofs that should be constrained

   \pre The map dbcmap must be subset of the maps of the vectors.
   */
  void ApplyDirichlettoSystem(Teuchos::RCP<Epetra_Vector>& x, Teuchos::RCP<Epetra_Vector>& b,
      Teuchos::RCP<const Epetra_Vector> dbcval, const Epetra_Map& dbcmap);

  /*!
   \brief Apply dirichlet boundary condition to a linear system of equations

   This is a flexible routine. The vectors x and dbcval might have different
   maps. The map does not need to contain all Dirichlet dofs.

   NOTE: Vector b does not need to contain all Dirichlet dofs defined in dbcmap

   The purpose is to set Dirichlet values at a subset of all Dirichlet
   boundaries.

   \param b (in/out)         : vector b of Ax=b
   \param dbcval (in)        : vector holding values that are supposed to be prescribed
   \param dbcmap (in)        : unique map of all dofs that should be constrained

   \pre The map dbcmap must be subset of the maps of the vectors.
   */
  void ApplyDirichlettoSystem(Teuchos::RCP<Epetra_Vector>& b,
      const Teuchos::RCP<const Epetra_Vector> dbcval, const Epetra_Map& dbcmap);

  /*!
   \brief Convert a Dirichlet toggle vector in a Dirichlet map

   The purpose of the routine is a smooth transition from Dirichlet toggle vectors
   to Dirichlet condition maps. Eventually, this method should be removed.

   A Dirichlet toogle vector is a real vector which holds a 1.0 at DOF subjected
   to Dirichlet boundary conditions and a 0.0 at every remaining/free DOF.

   \param dbctoggle (in)     : the Dirichlet toggle vector
   \return MapExtractor object which stores the Dirichlet condition and remaining (other) DOF map

   \author bborn
   \date 10/08
   */
  Teuchos::RCP<LINALG::MapExtractor> ConvertDirichletToggleVectorToMaps(
      const Teuchos::RCP<const Epetra_Vector>& dbctoggle);

}  // namespace LINALG

#endif /* SRC_LINALG_LINALG_UTILS_SPARSE_ALGEBRA_ASSEMBLE_H_ */
