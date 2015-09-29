/*----------------------------------------------------------------------*/
/*!
\file matpar_manager_elementwise.H
\brief Class to handle calls to material parameters from an optimization routine

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            http://www.lnm.mw.tum.de
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

class MatParManagerPerElement  : public MatParManager
{

public:

  //! constructor
  MatParManagerPerElement(Teuchos::RCP<DRT::Discretization> discret);

  //! destructor
  virtual ~MatParManagerPerElement(){};

  //! Setup specific parametrization layout
  virtual void Setup();

  //! nothing to be done in case of elementwise parameters
  virtual void Finalize(Teuchos::RCP<Epetra_MultiVector> dfint){return;}

protected:

  virtual void FillParameters(Teuchos::RCP<Epetra_MultiVector> params);

  virtual void InitParameters(int parapos, double val);

  virtual void ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint, double val, int elepos, int paraposglobal, int paraposlocal);

  const std::vector<int>& ParamsLIDtoeleGID() const {return paramsLIDtoeleGID_;}

private:

  //! build a graph of elements in the optimization parameter vector
  virtual void FillAdjacencyMatrix(const Epetra_Map& paramrowmap, Teuchos::RCP<Epetra_CrsMatrix> graph);

  //! Set action for FillAdjacencyMatrix
  virtual void SetAction(Teuchos::ParameterList& p);

  //! the two following 'maps' are processor specific: each proc has its own
  //! map keeping correspondence between params in the parameter vector (parameter LID in optparams_) and
  //! eleids (i.e gids of the discretizations ele layout)
  std::vector<int> paramsLIDtoeleGID_;

  //! map keeping correspondence between eleid (GID) + parameter block id and parameter LID
  std::map< int, std::vector<int> > eleGIDtoparamsLID_;

};
}  // namespace INVANA

#endif
