/*----------------------------------------------------------------------*/
/*! \file

\brief Structure field adapter for coupling with reduced-D airway trees


\level 3

*/

/*----------------------------------------------------------------------*/
/* macros */


#ifndef FOUR_C_ADAPTER_STR_REDAIRWAY_HPP
#define FOUR_C_ADAPTER_STR_REDAIRWAY_HPP
/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_adapter_str_wrapper.hpp"
#include "4C_fem_condition.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  class StructureRedAirway : public StructureWrapper
  {
   public:
    /// Constructor
    StructureRedAirway(Teuchos::RCP<Structure> stru);

    /// set pressure calculated from reduced-d airway tree
    void SetPressure(Teuchos::RCP<Epetra_Vector> couppres);

    /// calculate outlet fluxes for reduced-d airway tree
    void CalcFlux(
        Teuchos::RCP<Epetra_Vector> coupflux, Teuchos::RCP<Epetra_Vector> coupvol, double dt);

    /// calculate volume
    void CalcVol(std::map<int, double>& V);

    /// calculate initial volume
    void InitVol();

    //! (derived)
    void Update() override;

   private:
    /// map between coupling ID and conditions on structure
    std::map<int, Core::Conditions::Condition*> coupcond_;

    /// map of coupling IDs
    Teuchos::RCP<Epetra_Map> coupmap_;

    std::map<int, double> vn_;
    std::map<int, double> vnp_;
  };

}  // namespace Adapter
FOUR_C_NAMESPACE_CLOSE

#endif
