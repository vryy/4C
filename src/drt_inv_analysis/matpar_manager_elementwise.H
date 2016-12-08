/*----------------------------------------------------------------------*/
/*!
\file matpar_manager_elementwise.H
\brief Class to handle calls to material parameters from an optimization routine

<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifndef MATPAR_MANAGER_PERELEMENT_H
#define MATPAR_MANAGER_PERELEMENT_H

/*----------------------------------------------------------------------*/
/* headers */
#include "matpar_manager.H"


// forward declarations
namespace DRT
{
  class Discretization;
}

namespace INVANA
{

/*! \class MatParManagerPerElement
 *  \brief Class to have elementwise constant parameters
 *
 *  This class handles elementwise constant material parameters.
 *  The chainrule application in the gradient evaluation reduces
 *  to a minimum.
 *
 *  It provides itself as a base class to all the reduced basis
 *  implementation. So given the proper projections between
 *  reduced basis and full (elemental) basis, optimization
 *  parameters can be injected into the element-wise layout.
 *
 *
 *  \author kehl \date 2013
 */
class MatParManagerPerElement  : public MatParManager
{

public:

  //! constructor
  MatParManagerPerElement(Teuchos::RCP<DRT::Discretization> discret);

  //! destructor
  virtual ~MatParManagerPerElement(){};

  //! Setup specific parametrization layout
  virtual void Setup();

  //! Account for distributed chain-rule application
  virtual void Finalize(Teuchos::RCP<Epetra_MultiVector> source,
      Teuchos::RCP<Epetra_MultiVector> target
      );

  //! @name Access projection operators
  //@{

  //! Access to the restrictor
  Teuchos::RCP<Epetra_CrsMatrix> Restrictor() {return restrictor_;}

  //! Access to the prolongator (use it transposed)
  Teuchos::RCP<Epetra_CrsMatrix> Prolongator() {return prolongator_;}

  //@}

protected:

  //! Get current material parameters
  virtual void FillParameters(Teuchos::RCP<Epetra_MultiVector> params);

  //! apply projection to INVANA::DcsMatrix and get diagonal
  virtual void ApplyParametrization(DcsMatrix& matrix, Teuchos::RCP<Epetra_MultiVector> diagonals);

  //! Initialize optimization parameters
  virtual void InitParameters(int parapos, double val);

  //! Chain rule application
  virtual void ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,
      double val, int elepos, int paraposglobal, int paraposlocal);

  std::vector<int>& ParamsLIDtoeleGID() {return paramsLIDtoeleGID_;}

  std::map< int, std::vector<int> >& EleGIDtoparamsLID() {return eleGIDtoparamsLID_;}

  /* @name Projections between sparse and full
   * In here these will be just unit diagonals;
   * Its more interesting for the specializations.
   */
  //@{

  //! linear operator to project from the full signal to the sparse approximation
  Teuchos::RCP<Epetra_CrsMatrix> restrictor_;

  //! linear operator to project from the sparse approximation to the full signal
  Teuchos::RCP<Epetra_CrsMatrix> prolongator_;
  //@}

private:

  //! build a graph of elements in the optimization parameter vector
  virtual void FillAdjacencyMatrix(const Epetra_Map& paramrowmap, Teuchos::RCP<Epetra_CrsMatrix> graph);

  //! Set action for FillAdjacencyMatrix
  virtual void SetAction(Teuchos::ParameterList& p);

  //! @name Parameter access maps (processor specific)
  //@{

  /*! \brief Mapping between parameters and elements
   *
   * map keeping correspondence between params in the parameter
   * vector (parameter LID in optparams_) and eleids (i.e gids of
   * the discretizations element layout)
   *
   */
  std::vector<int> paramsLIDtoeleGID_;

  /*! \brief Mapping between elements and parameters
   *
   * map keeping correspondence between eleid (GID), parameter block id
   * (i.e. material input file parameter) and parameter LID
   *
   */
  std::map< int, std::vector<int> > eleGIDtoparamsLID_;
  //@}

};
}  // namespace INVANA

#endif
