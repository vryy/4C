/*---------------------------------------------------------------------*/
/*!
\file spherebeamlinking_params.H

\brief data container holding all contractile cells input parameters

\level 3

\maintainer Jonas Eichinger
*/
/*---------------------------------------------------------------------*/

#ifndef INTEGRINS_H
#define INTEGRINS_H

#include "../drt_lib/drt_dserror.H"
#include<Teuchos_RCP.hpp>


// forward declaration

namespace STR {
namespace TIMINT {
  class BaseDataGlobalState;
}
}
namespace MAT
{
  class CrosslinkerMat;
}
namespace BEAMINTERACTION {

  /*!
   * data container for input file parameters for submodel crosslinking in beam interaction */
  class SphereBeamLinkingParams
  {
  public:
    //! constructor
    SphereBeamLinkingParams();

    //! destructor
    virtual ~SphereBeamLinkingParams() {};

    //! initialize with the stuff coming from input file
    void Init( STR::TIMINT::BaseDataGlobalState const & gstate );

    //! setup member variables
    void Setup();

    //! returns the isinit_ flag
    inline const bool& IsInit() const
    { return isinit_; };

    //! returns the issetup_ flag
    inline const bool& IsSetup() const
    { return issetup_; };

    //! Checks the init and setup status
    inline void CheckInitSetup() const
    {
      if (!IsInit() or !IsSetup())
        dserror("Call Init() and Setup() first!");
    }

    //! Checks the init status
    inline void CheckInit() const
    {
      if (!IsInit())
        dserror("Init() has not been called, yet!");
    }

    /// max number of integrins per cell
    int MaxNumberOfIntegrinsPerCell() const
    {
      CheckInitSetup();
      return maxnumintegrins_;
    };

    /// integrin material id
    Teuchos::RCP< MAT::CrosslinkerMat const > GetLinkerMaterial() const
    {
      CheckInitSetup();
      return mat_;
    };

    /// time step for stochastic events concerning crosslinking
    double const & DeltaTime() const
    {
      CheckInitSetup();
      return deltatime_;
    };

    /// time step for stochastic events concerning crosslinking
    double const & ContractionRate() const
    {
      CheckInitSetup();
      return contractionrate_;
    };

  private:

    bool isinit_;

    bool issetup_;

    // maximum number integrins per cell
    int maxnumintegrins_;
    /// time step for stochastic events concerning integrins, e.g. catch-slip-bond behavior
    double deltatime_;
    /// contraction rate of cell (integrin linker) in [μm/s]
    double contractionrate_;
    /// crosslinker material
    Teuchos::RCP< MAT::CrosslinkerMat > mat_;

  };

}

#endif
