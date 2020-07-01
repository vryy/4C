/*----------------------------------------------------------------------*/
/*! \file


\brief routines to calculate forcing for homogeneous isotropic turbulence simulations with
passive-scalar transport

\level 2


*----------------------------------------------------------------------*/


#ifndef TURBULENCE_HIT_SCALAR_FORCING_H
#define TURBULENCE_HIT_SCALAR_FORCING_H

#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RCP.hpp"

#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_scatra.H"

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace SCATRA
{
  // forward declarations
  class ScaTraTimIntImpl;

  class HomIsoTurbScalarForcing
  {
   public:
    //! constructor
    HomIsoTurbScalarForcing(ScaTraTimIntImpl* timeint);

    //! destructor
    ~HomIsoTurbScalarForcing() { return; };

    //    //! set time integration scheme specific vectors (required for gen-alpha)
    //    void SetTimeIntSpecificVectors(Teuchos::RCP<Epetra_Vector> phiaf,
    //                                   bool gen_alpha = false);

    //! initialize with initial spectrum
    void SetInitialSpectrum(INPAR::SCATRA::InitialField init_field_type);

    //! turn on forcing
    void ActivateForcing(const bool activate);

    //! calculate power input
    void CalculateForcing(const int step);

    //! get forcing
    void UpdateForcing(const int step);

    //! time update of energy spectrum
    void TimeUpdateForcing();

   protected:
    //! sort criterium for double values up to a tolerance of 10-9
    class LineSortCriterion
    {
     public:
      bool operator()(const double& p1, const double& p2) const { return (p1 < p2 - 1E-9); }

     protected:
     private:
    };

   private:
    //! type of forcing
    INPAR::FLUID::ForcingType forcing_type_;

    //! scatra discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! state vector of volume force to be computed
    Teuchos::RCP<Epetra_Vector> forcing_;

    //! state vectors used to compute forcing
    Teuchos::RCP<Epetra_Vector> phinp_;
    Teuchos::RCP<Epetra_Vector> phiaf_;

    //! threshold wave number for forcing
    //! i.e., forcing is applied to wave numbers <= threshold wave number
    double threshold_wavenumber_;

    //! identify gen-alpha time integration
    bool is_genalpha_;

    //! number of resolved mode
    int nummodes_;

    //! vector of coordinates in one spatial direction (same for the other two directions)
    Teuchos::RCP<std::vector<double>> coordinates_;

    //! vector of wave numbers
    Teuchos::RCP<std::vector<double>> wavenumbers_;

    //! vector scalar variance spectrum (sum over k=const) at time n
    Teuchos::RCP<std::vector<double>> scalarvariancespectrum_n_;

    //! vector scalar variance spectrum  (sum over k=const) at time n+1/n+af
    Teuchos::RCP<std::vector<double>> scalarvariancespectrum_np_;

    //! time step length
    double dt_;

    //! flag to activate forcing
    bool activate_;

    //! linear compensation factor
    Teuchos::RCP<Epetra_SerialDenseVector> force_fac_;

    //! interpolation function
    static double Interpolate(
        const double& x, const double& x_1, const double& x_2, const double& y_1, const double& y_2)
    {
      const double value = y_1 + (y_2 - y_1) / (x_2 - x_1) * (x - x_1);
      return value;
    }
  };

}  // namespace SCATRA

#endif
