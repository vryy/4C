/*----------------------------------------------------------------------------*/
/*!

\brief Solver for fluid field on a moving ALE mesh in combination with an immersed structure

\maintainer Jonas Eichinger
*/
/*----------------------------------------------------------------------------*/

#ifndef ADAPTER_FLUID_ALE_IMMERSED_H
#define ADAPTER_FLUID_ALE_IMMERSED_H


#include "ad_fld_fluid_ale.H"

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace ADAPTER
{
  /// fluid on ale in combination with immersed
  class FluidAleImmersed : public FluidAle
  {
   public:
    FluidAleImmersed(const Teuchos::ParameterList& prbdyn, std::string condname);

    void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd);

    void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove);

    void Output(const int step = -1, const double time = -1.0);
  };

}  // namespace ADAPTER

#endif
