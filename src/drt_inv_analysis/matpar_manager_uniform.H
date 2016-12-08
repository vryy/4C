/*----------------------------------------------------------------------*/
/*!
\file matpar_manager_uniform.H
\brief Creating patches from inputfile

<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/


/*----------------------------------------------------------------------*/
/* definitions */
#ifndef MATPAR_MANAGER_UNIFORM_H
#define MATPAR_MANAGER_UNIFORM_H

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

/*! \class MatParManagerUniform
 *  \brief Class to have patchwise constant parameters
 *
 *  Patches are created as the material blocks given
 *  in the input file.
 *
 *  \author kehl \date 08/2016
 */
class MatParManagerUniform  : public MatParManagerPerElement
{

public:

  //! constructor
  MatParManagerUniform(Teuchos::RCP<DRT::Discretization> discret);

  //! destructor
  virtual ~MatParManagerUniform(){};

  //! Setup specific parametrization layout
  virtual void Setup();

  //! Account for distributed chain-rule application
  virtual void Finalize(Teuchos::RCP<Epetra_MultiVector> source,
      Teuchos::RCP<Epetra_MultiVector> target
      );

  //! Access to the unreduced layout
  Teuchos::RCP<Epetra_Map> UnreducedMap() {return elewise_map_;}


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

  //! @name Parameters of the initial elementwise distribution
  //@{

  //! elementwise set of optimization parameters
  Teuchos::RCP<Epetra_MultiVector> optparams_elewise_;

  //! elementwise layout
  Teuchos::RCP<Epetra_Map> elewise_map_;

  //! elementwise graph
  Teuchos::RCP<Epetra_CrsMatrix> graph_;
  //@}

  //! global optimization parameter id to patch id
  std::map<int,int> pidtopatch_;

};
}  // namespace INVANA

#endif
