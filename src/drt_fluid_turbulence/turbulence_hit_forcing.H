/*----------------------------------------------------------------------*/
/*!

\brief routines to calculate forcing for homogeneous isotropic turbulence simulations

\maintainer Martin Kronbichler

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef TURBULENCE_HIT_FORCING_H
#define TURBULENCE_HIT_FORCING_H

#include "Epetra_Vector.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RCP.hpp"

#include "../drt_inpar/inpar_fluid.H"

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace FLD
{
  // forward declarations
  class FluidImplicitTimeInt;
  class XWall;

  class ForcingInterface
  {
   public:
    //! constructor
    ForcingInterface() {}

    virtual ~ForcingInterface() {}

    //! initialize with initial spectrum
    virtual void SetInitialSpectrum(INPAR::FLUID::InitialField init_field_type) = 0;

    //! turn on forcing
    virtual void ActivateForcing(const bool activate) = 0;

    //! calculate power input
    virtual void CalculateForcing(const int step) = 0;

    //! get forcing
    virtual void UpdateForcing(const int step) = 0;

    //! time update of energy spectrum
    virtual void TimeUpdateForcing() = 0;
  };

  class HomIsoTurbForcing : public ForcingInterface
  {
   public:
    //! constructor
    HomIsoTurbForcing(FluidImplicitTimeInt& timeint);

    //! destructor
    virtual ~HomIsoTurbForcing() {}

    //! initialize with initial spectrum
    virtual void SetInitialSpectrum(INPAR::FLUID::InitialField init_field_type);

    //! turn on forcing
    void ActivateForcing(const bool activate);

    //! calculate power input
    virtual void CalculateForcing(const int step);

    //! get forcing
    virtual void UpdateForcing(const int step);

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

    //! type of forcing
    INPAR::FLUID::ForcingType forcing_type_;

    //! fluid discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! state vector of volume force to be computed
    Teuchos::RCP<Epetra_Vector> forcing_;

    //! state vectors used to compute forcing
    Teuchos::RCP<Epetra_Vector> velnp_;
    Teuchos::RCP<Epetra_Vector> velaf_;

    //! threshold wave number for forcing
    //! i.e., forcing is applied to wave numbers <= threshold wave number
    double threshold_wavenumber_;

    //! identify gen-alpha time integration
    bool is_genalpha_;

    //! number of forcing time steps for decaying case
    int num_force_steps_;

    //! specifies the special flow
    enum FlowType
    {
      decaying_homogeneous_isotropic_turbulence,
      forced_homogeneous_isotropic_turbulence
    } flow_type_;

    //! number of resolved mode
    int nummodes_;

    //! vector of coordinates in one spatial direction (same for the other two directions)
    Teuchos::RCP<std::vector<double>> coordinates_;

    //! vector of wave numbers
    Teuchos::RCP<std::vector<double>> wavenumbers_;

    //! vector energy spectrum (sum over k=const) at time n
    Teuchos::RCP<std::vector<double>> energyspectrum_n_;

    //! vector energy spectrum  (sum over k=const) at time n+1/n+af
    Teuchos::RCP<std::vector<double>> energyspectrum_np_;

    //! time step length
    double dt_;

    //! flag to activate forcing
    bool activate_;

    //! linear compensation factor
    Teuchos::RCP<Epetra_SerialDenseVector> force_fac_;

    //! fixed power input
    double Pin_;

    //! energy_contained in lowest wave numbers
    double E_kf_;

    //! fixed power input factor
    Teuchos::RCP<Epetra_SerialDenseVector> fixed_power_fac_;

    //! interpolation function
    static double Interpolate(
        const double& x, const double& x_1, const double& x_2, const double& y_1, const double& y_2)
    {
      const double value = y_1 + (y_2 - y_1) / (x_2 - x_1) * (x - x_1);
      return value;
    }
  };

  // there are quite some differences for the HDG case
  // author: bk 03/15
  class HomIsoTurbForcingHDG : public HomIsoTurbForcing
  {
   public:
    //! constructor
    HomIsoTurbForcingHDG(FluidImplicitTimeInt& timeint);

    //! destructor
    virtual ~HomIsoTurbForcingHDG() {}

    //! initialize with initial spectrum
    void SetInitialSpectrum(INPAR::FLUID::InitialField init_field_type);

    //! calculate power input
    void CalculateForcing(const int step);

    //! get forcing
    void UpdateForcing(const int step);
  };

  // this is an adaptive body force for the periodic hill benchmark such
  // that the mass flow reaches the given value.
  // author: bk 12/14
  class PeriodicHillForcing : public ForcingInterface
  {
   public:
    //! constructor
    PeriodicHillForcing(FluidImplicitTimeInt& timeint);

    //! destructor
    virtual ~PeriodicHillForcing() {}

    //! initialize with initial spectrum
    void SetInitialSpectrum(INPAR::FLUID::InitialField init_field_type) { return; }

    //! turn on forcing
    void ActivateForcing(const bool activate) { return; }

    //! calculate power input
    void CalculateForcing(const int step) { return; }

    //! get forcing
    void UpdateForcing(const int step);

    //! time update of energy spectrum
    void TimeUpdateForcing();


   private:
    //! fluid discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! state vector of volume force to be computed
    Teuchos::RCP<Epetra_Vector> forcing_;

    //! state vectors used to compute forcing
    Teuchos::RCP<Epetra_Vector> velnp_;
    Teuchos::RCP<Epetra_Vector> velaf_;

    //! xwall object is required for evaluating inner element planes of channel
    Teuchos::RCP<FLD::XWall> myxwall_;

    //! values of previous step
    double oldforce_;
    double oldflow_;

    //! reference value for optimal control
    double idealmassflow_;

    // length of overall flow domain
    double length_;

    //! step and statistical data
    int step_;
    int count_;
    double sum_;
  };

}  // namespace FLD

#endif
