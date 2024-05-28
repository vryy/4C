/*---------------------------------------------------------------------*/
/*! \file

\brief Control routine for coupled reduced airways and continuum tissue models


\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_RED_AIRWAYS_TISSUE_HPP
#define FOUR_C_RED_AIRWAYS_TISSUE_HPP

#include "4C_config.hpp"

#include "4C_adapter_algorithmbase.hpp"
#include "4C_inpar_bio.hpp"
#include "4C_red_airways_resulttest.hpp"
#include "4C_utils_result_test.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


// forward declarations

namespace ADAPTER
{
  class StructureRedAirway;
}

namespace AIRWAY
{
  class RedAirwayImplicitTimeInt;

  class RedAirwayTissue : public ADAPTER::AlgorithmBase
  {
   public:
    /// Standard Constructor
    RedAirwayTissue(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);


    void read_restart(const int step) override;

    void SetupRedAirways();

    /// time integration of coupled problem
    void Integrate();

    void do_structure_step();

    void RelaxPressure(int iter);

    void DoRedAirwayStep();

    /// flag whether iteration between fields should be finished
    bool NotConverged(int iter);

    void OutputIteration(Teuchos::RCP<Epetra_Vector> pres_inc,
        Teuchos::RCP<Epetra_Vector> scaled_pres_inc, Teuchos::RCP<Epetra_Vector> flux_inc,
        Teuchos::RCP<Epetra_Vector> scaled_flux_inc, int iter);

    void update_and_output();

    /// access to structural field
    Teuchos::RCP<ADAPTER::StructureRedAirway>& structure_field() { return structure_; }

    /// access to airway field
    Teuchos::RCP<RedAirwayImplicitTimeInt>& RedAirwayField() { return redairways_; }


   private:
    /// underlying structure
    Teuchos::RCP<ADAPTER::StructureRedAirway> structure_;

    Teuchos::RCP<RedAirwayImplicitTimeInt> redairways_;

    /// redundant vector of outlet pressures (new iteration step)
    Teuchos::RCP<Epetra_Vector> couppres_ip_;

    /// redundant vector of outlet pressures (old iteration step)
    Teuchos::RCP<Epetra_Vector> couppres_im_;


    // Aitken Variables:
    // Relaxation factor
    Teuchos::RCP<Epetra_Vector> omega_np_;

    /// redundant vector of outlet pressures (before old iteration step), p^{i}_{n+1}
    Teuchos::RCP<Epetra_Vector> couppres_il_;

    /// redundant vector of outlet pressures (old iteration step guess), \tilde{p}^{i+1}_{n+1}
    Teuchos::RCP<Epetra_Vector> couppres_im_tilde_;

    /// redundant vector of outlet pressures (new iteration step guess), \tilde{p}^{i+2}_{n+1}
    Teuchos::RCP<Epetra_Vector> couppres_ip_tilde_;


    /// redundant vector of outlet fluxes (new iteration step)
    Teuchos::RCP<Epetra_Vector> coupflux_ip_;

    /// redundant vector of outlet fluxes (old iteration step)
    Teuchos::RCP<Epetra_Vector> coupflux_im_;

    /// redundant vector of 3D volumes (new iteration step)
    Teuchos::RCP<Epetra_Vector> coupvol_ip_;

    /// redundant vector of 3D volumes (old iteration step)
    Teuchos::RCP<Epetra_Vector> coupvol_im_;

    /// internal iteration step
    int itermax_;

    /// restart step
    int uprestart_;

    /// internal tolerance for pressure
    double tolp_;

    /// internal tolerance for flux
    double tolq_;

    /// defined normal direction
    double normal_;

    /// fixed relaxation paramater
    double omega_;

    /// fixed relaxation paramater
    INPAR::ARTNET::Relaxtype3D0D relaxtype_;
  };
}  // namespace AIRWAY

FOUR_C_NAMESPACE_CLOSE

#endif
