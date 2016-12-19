/*----------------------------------------------------------------------*/
/*!
\file matpar_manager_tvsvd.H
\brief Creating reduced basis from TV approximation

<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/


/*----------------------------------------------------------------------*/
/* definitions */
#ifndef MATPAR_MANAGER_TVSVD_H
#define MATPAR_MANAGER_TVSVD_H

/*----------------------------------------------------------------------*/
/* headers */
#include "matpar_manager_elementwise.H"


// forward declarations
namespace DRT
{
  class Discretization;
}

namespace INVANA
{

/*! \class MatParManagerTVSVD
 *  \brief Reduced basis parametetrization based on eigendecomposition
 *
 *  In order to get a feasible amount of optimization parameters, basis
 *  reduction is performed on a quadratic approximation of the total
 *  variation functional, which is already used as regularization.
 *
 *  Up to now the approximation point is at the maximum a posterior
 *  solution. But this can be changed.
 *
 *  \author kehl \date 12/2016
 */
class MatParManagerTVSVD  : public MatParManagerPerElement
{

public:

  //! constructor
  MatParManagerTVSVD(Teuchos::RCP<DRT::Discretization> discret);

  //! destructor
  virtual ~MatParManagerTVSVD(){};

  //! Setup specific parametrization layout
  virtual void Setup();

  //! Account for distributed chain-rule application
  virtual void Finalize(Teuchos::RCP<Epetra_MultiVector> source,
      Teuchos::RCP<Epetra_MultiVector> target
      );

  //! Return initial covariance guess
  virtual Teuchos::RCP<Epetra_CrsMatrix> InitialCovariance();

protected:

  //! Get current material parameters
  virtual void FillParameters(Teuchos::RCP<Epetra_MultiVector> params);

  //! apply projection to INVANA::DcsMatrix and get diagonal
  virtual void ApplyParametrization(DcsMatrix& matrix, Teuchos::RCP<Epetra_MultiVector> diagonals);

  //! Initialize optimization parameters
  virtual void InitParameters();

  //! Chain rule application
  virtual void ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,
      double val, int elepos, int paraposglobal, int paraposlocal);

private:

  //! @name Setup the parametrization
  //@{

  //! Create projection to patch-wise atoms
  void CreateProjection();
  //@}

  //! @name Helper to setup the parametrization
  //@{

  //! Read elementwise parameter distribution from input
  void ReadMAPApproximation();

  //! The full approximation linear operator
  void SetupTVOperator();

  //! Factorize linearized TV operator
  void Factorize();

  //! Solve the eigenproblem for the TV operator using Anasazi
  void AnasaziEigenProblem(Teuchos::RCP<Epetra_CrsMatrix> A, ///[in] Matrix to decompose
      Teuchos::RCP<Epetra_MultiVector> L, ///[out] eigenvectors
      Teuchos::ParameterList params  ///[in] how to decompose
      );

  //! Setup restriction and prolongation operators
  void SetupRandP(int numvecs);

  //! number of eigenvectors to use
  int max_num_levels_;

  /*! \brief random initialization of eigenproblem
   *  In a nested parallel setup, one needs on each local proc
   *  the same seed across the groups, otherwise the eigenvectors
   *  come out scaled by -1.0 at random. But all groups must have
   *  the same basisvectors!
   *
   * \author kehl \date 12/2016
   */
  int seed_;

  //! Epetra_Util object for "random" drawing
  Epetra_Util util_;

  /*! \brief Random multi vector
   * Fill multivector randomly in the interval [-1,1], but
   * do consitently accros potential nested parallel groups.
   *
   * \author kehl \date 12/2016
   */
  void Random(Epetra_MultiVector& randvec);

  //! TV epsilon
  double eps_;

  // file to read the map approximation from
  std::string map_restart_file_;

  // step from which to read the map approximation
  int map_restart_step_;
  //@}

  //! @name Parameters of the initial elementwise distribution
  //@{

  //! elementwise set of optimization parameters
  Teuchos::RCP<Epetra_MultiVector> optparams_elewise_;

  //! elementwise layout
  Teuchos::RCP<Epetra_Map> elewise_map_;

  //! The linearized TV operator
  Teuchos::RCP<Epetra_CrsMatrix> lintvop_;

  //! Covariance approximation
  Teuchos::RCP<DcsMatrix> fullcovariance_;
  //@}

  //! @name Projection between sparse appoximation and full basis
  //@{

  //! Eigenvectors of the TV approximation
  Teuchos::RCP<Epetra_MultiVector> evecs_;

  //! number of eigenvectors to compute
  int nev_;

  //! number of blocks (for the block eigensolver)
  int nblocks_;

  //! block size (for the block eigensolver)
  int bsize_;
  //@}

};
}  // namespace INVANA

#endif
