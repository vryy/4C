/*!

\file turbulence_hit_initial_scalar_field.H

\brief routines to initialize homogeneous isotropic turbulence simulations with passive scalar transport

\level 2

<pre>
\maintainer Benjamin Krank
            krank@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*/


#ifndef TURBULENCE_HIT_INITIAL_SCALAR_FIELD_H
#define TURBULENCE_HIT_INITIAL_SCALAR_FIELD_H

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

#include "../drt_inpar/inpar_scatra.H"

//forward declarations
namespace DRT
{
  class Discretization;
}

namespace SCATRA
{

// forward declarations
class ScaTraTimIntImpl;

// inital condition for homogeneous isotropic turbulence
// based on the Comte-Bellot - Corrsion experiment
class HomIsoTurbInitialScalarField
{
  public:

    //! constructor
    HomIsoTurbInitialScalarField(ScaTraTimIntImpl& timeint,
                                 const INPAR::SCATRA::InitialField initfield);

    //! destructor
    ~HomIsoTurbInitialScalarField(){return;};

    //! calculate initial field
    void CalculateInitialField();

  protected:

    //! sort criterium for double values up to a tolerance of 10-9
    class LineSortCriterion
    {
      public:
        bool operator() (const double& p1,
                         const double& p2) const { return (p1 < p2 - 1E-9);}

      protected:
      private:
    };

  private:

    //! estimate energy form given scalar variance spectrum (function for E_phi)
    double CalculateEnergyFromSpectrum(double k);

    //! scatra discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! state vectors to be initialized
    Teuchos::RCP<Epetra_Vector> phinp_;
    Teuchos::RCP<Epetra_Vector> phin_;

    //! type of energy spectrum for initialization
    INPAR::SCATRA::InitialField type_;

    //! number of resolved mode
    int nummodes_;

    //! vector of coordinates in one spatial direction (same for the other two directions)
    Teuchos::RCP<std::vector<double> > coordinates_;

};

} // end FLUID

#endif
