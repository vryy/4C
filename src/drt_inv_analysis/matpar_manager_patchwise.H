/*----------------------------------------------------------------------*/
/*!
\file matpar_manager_patchwise.H
\brief Creating patches from an elementwise layout

<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/


/*----------------------------------------------------------------------*/
/* definitions */
#ifndef MATPAR_MANAGER_PERPATCH_H
#define MATPAR_MANAGER_PERPATCh_H

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

/*! \class MatParManagerPerPatch
 *  \brief Class to have patchwise constant parameters
 *
 *  \author kehl \date 08/2016
 */
class MatParManagerPerPatch  : public MatParManagerPerElement
{

public:

  //! constructor
  MatParManagerPerPatch(Teuchos::RCP<DRT::Discretization> discret);

  //! destructor
  virtual ~MatParManagerPerPatch(){};

  //! Setup specific parametrization layout
  virtual void Setup();

  //! Account for distributed chain-rule application
  virtual void Finalize(Teuchos::RCP<Epetra_MultiVector> source,
      Teuchos::RCP<Epetra_MultiVector> target
      );

  //! Access to the unreduced layout
  Teuchos::RCP<Epetra_Map> UnreducedMap() {return elewise_map_;}

  //! Access to the projector
  Teuchos::RCP<Epetra_CrsMatrix> Projector() {return restrictor_;}

protected:

  //! Get current material parameters
  virtual void FillParameters(Teuchos::RCP<Epetra_MultiVector> params);

  //! Initialize optimization parameters
  virtual void InitParameters();

  //! Chain rule application
  virtual void ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,
      double val, int elepos, int paraposglobal, int paraposlocal);

private:

  //! @name Setup the parametrization
  //@{

  //! Read elementwise parameter distribution from input
  void ReadMAPApproximation();

  //! Reduce to a reduced order basis
  void ReduceBasis();

  //! find connectivities within levels
  void FindLevelConnectivity(std::map<int, std::vector<int> >& neighbourhood, // [in] neighbourinformation
      std::map<int, std::vector<int> >& patches // [out] patches within this level
      );

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
  //@}

  //! @name Projection between reduced and full basis
  //@{

  Teuchos::RCP<LINALG::MultiMapExtractor> patchmap_;

  Teuchos::RCP<Epetra_CrsMatrix> restrictor_;
  //@}

  //! global optimization parameter id to patch id
  std::map<int,int> pidtopatch_;

  //! number of patch levels
  int num_levels_;

};
}  // namespace INVANA

#endif
