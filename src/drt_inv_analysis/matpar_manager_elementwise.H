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
 *  The vector layout and the unique vector layout are identical.
 *
 *  \author kehl \date somewhen2013
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

protected:

  //! Get current material parameters
  virtual void FillParameters(Teuchos::RCP<Epetra_MultiVector> params);

  //! Initialize optimization parameters
  virtual void InitParameters(int parapos, double val);

  //! Chain rule application
  virtual void ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,
      double val, int elepos, int paraposglobal, int paraposlocal);

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

};
}  // namespace INVANA

#endif
