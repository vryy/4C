/*----------------------------------------------------------------------*/
/*! \file

\brief routines to initialize homogeneous isotropic turbulence simulations


\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_TURBULENCE_HIT_INITIAL_FIELD_HPP
#define FOUR_C_FLUID_TURBULENCE_HIT_INITIAL_FIELD_HPP

#include "baci_config.hpp"

#include "baci_inpar_fluid.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace FLD
{
  // forward declarations
  class FluidImplicitTimeInt;

  // inital condition for homogeneous isotropic turbulence
  // based on the Comte-Bellot - Corrsion experiment
  class HomIsoTurbInitialField
  {
   public:
    //! constructor
    HomIsoTurbInitialField(
        FluidImplicitTimeInt& timeint, const INPAR::FLUID::InitialField initfield);

    //! destructor
    virtual ~HomIsoTurbInitialField() = default;

    //! calculate initial field
    virtual void CalculateInitialField();

   protected:
    //! sort criterium for double values up to a tolerance of 10-9
    class LineSortCriterion
    {
     public:
      bool operator()(const double& p1, const double& p2) const { return (p1 < p2 - 1E-9); }

     protected:
     private:
    };

    //! non-dimensionalize and store experimental data
    void PrepareExparimentalData();

    //! estimate energy form given energy spectrum (experimental data)
    double InterpolateEnergyFromSpectrum(double k);

    //! estimate energy form given energy spectrum (function for E)
    double CalculateEnergyFromSpectrum(double k);

    //! fluid discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! state vectors to be initialized
    Teuchos::RCP<Epetra_Vector> velnp_;
    Teuchos::RCP<Epetra_Vector> veln_;
    Teuchos::RCP<Epetra_Vector> velnm_;
    //! type of energy spectrum for initialization
    INPAR::FLUID::InitialField type_;

    //! number of resolved mode
    int nummodes_;

    //! vector of coordinates in one spatial direction (same for the other two directions)
    Teuchos::RCP<std::vector<double>> coordinates_;

    //! vector containing wave numbers of experiment
    std::vector<double> k_exp_;
    //! vector containing corresponding energy
    std::vector<double> E_exp_;
  };

  // inital condition for homogeneous isotropic turbulence
  // based on the Comte-Bellot - Corrsion experiment
  class HomIsoTurbInitialFieldHDG : public HomIsoTurbInitialField
  {
   public:
    //! constructor
    HomIsoTurbInitialFieldHDG(
        FluidImplicitTimeInt& timeint, const INPAR::FLUID::InitialField initfield);


    //! calculate initial field
    void CalculateInitialField() override;

   protected:
    Teuchos::RCP<Epetra_Vector> intveln_;
    Teuchos::RCP<Epetra_Vector> intvelnm_;
    Teuchos::RCP<Epetra_Vector> intvelnp_;
  };

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
