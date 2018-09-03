/*!

\file turbulence_hit_initial_field.H

\brief routines to initialize homogeneous isotropic turbulence simulations

References are

<pre>

</pre>

\level 2
<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*/


#ifndef TURBULENCE_HIT_INITIAL_FIELD_H
#define TURBULENCE_HIT_INITIAL_FIELD_H

#include "Epetra_Vector.h"
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

  // inital condition for homogeneous isotropic turbulence
  // based on the Comte-Bellot - Corrsion experiment
  class HomIsoTurbInitialField
  {
   public:
    //! constructor
    HomIsoTurbInitialField(
        FluidImplicitTimeInt& timeint, const INPAR::FLUID::InitialField initfield);

    //! destructor
    virtual ~HomIsoTurbInitialField() { return; };

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

    //! destructor
    ~HomIsoTurbInitialFieldHDG() { return; };

    //! calculate initial field
    void CalculateInitialField();

   protected:
    Teuchos::RCP<Epetra_Vector> intveln_;
    Teuchos::RCP<Epetra_Vector> intvelnm_;
    Teuchos::RCP<Epetra_Vector> intvelnp_;
  };

}  // namespace FLD

#endif
