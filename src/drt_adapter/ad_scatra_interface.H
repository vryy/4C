/*----------------------------------------------------------------------*/
/*! \file
\brief Interface for all scatra adapters.
\level 1
\maintainer Anh-Tu Vuong
 */
/*----------------------------------------------------------------------*/


#ifndef AD_SCATRA_INTERFACE_
#define AD_SCATRA_INTERFACE_


#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>


// forward declarations
namespace DRT
{
  class Discretization;
}

namespace SCATRA
{
  class MeshtyingStrategyBase;
}


namespace ADAPTER
{
  /*! \brief General pure virtual interface for all scatra time integrators and scatra adapters.
   *
   *  The point is to keep coupled problems as far apart from our field solvers as
   *  possible. Each scatra field solver we want to use should get its own subclass
   *  of this. The coupled algorithm should be able to extract all the information
   *  from the scatra field it needs using this interface.
   *
   * \sa ScaTraTimIntImpl
   * \date 12/2016
   * \author rauch (rauch@lnm.mw.tum.de)    */
  class ScatraInterface
  {
   public:
    //! constructor
    ScatraInterface(){};

    //! virtual to get polymorph destruction
    virtual ~ScatraInterface(){};

    //! return discretization
    virtual Teuchos::RCP<DRT::Discretization> Discretization() const = 0;

    //! add parameters specific for time-integration scheme
    virtual void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver = false) = 0;

    //! return number of dofset associated with displacement dofs
    virtual int NdsDisp() = 0;

    //! return rcp ptr to neumann loads vector
    virtual Teuchos::RCP<Epetra_Vector> GetNeumannLoadsPtr() = 0;

    //! return meshtying strategy (includes standard case without meshtying)
    virtual const Teuchos::RCP<SCATRA::MeshtyingStrategyBase>& Strategy() const = 0;

    //! return scalar field phi at time n
    virtual Teuchos::RCP<Epetra_Vector> Phin() = 0;

  };  // class ScatraInterface
}  // namespace ADAPTER


#endif /* AD_SCATRA_INTERFACE_ */
