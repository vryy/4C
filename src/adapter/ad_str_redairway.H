/*----------------------------------------------------------------------*/
/*! \file

\brief Structure field adapter for coupling with reduced-D airway trees


\level 3

*/

/*----------------------------------------------------------------------*/
/* macros */


#ifndef AD_STR_REDAIRWAY_H
#define AD_STR_REDAIRWAY_H
/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RCP.hpp>
#include "Epetra_Vector.h"

#include "ad_str_wrapper.H"


namespace DRT
{
  class Condition;
}


namespace ADAPTER
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

    virtual void Update();

   private:
    /// map between coupling ID and conditions on structure
    std::map<int, DRT::Condition*> coupcond_;

    /// map of coupling IDs
    Teuchos::RCP<Epetra_Map> coupmap_;

    std::map<int, double> Vn_;
    std::map<int, double> Vnp_;
  };

}  // namespace ADAPTER
#endif
