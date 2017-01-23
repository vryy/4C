/*---------------------------------------------------------------------*/
/*!
\file contractilecells_params.H

\brief data container holding all contractile cells input parameters

\level 3

\maintainer Jonas Eichinger
*/
/*---------------------------------------------------------------------*/

#ifndef CONTRACTILECELLS_H
#define CONTRACTILECELLS_H

#include "../drt_lib/drt_dserror.H"


// forward declaration
namespace BEAMINTERACTION {

  /*!
   * data container for input file parameters for submodel crosslinking in beam interaction */
  class ContractileCellsParams
  {
  public:
    //! constructor
    ContractileCellsParams();

    //! destructor
    virtual ~ContractileCellsParams() {};

    //! initialize with the stuff coming from input file
    void Init();

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

    /// number of contractile cells
    int const& NumCells() const
    {
      CheckInitSetup();
      return numcells_;
    };



  private:

    bool isinit_;

    bool issetup_;

    // number of cells in the simulated volume
    int numcells_;

  };

}

#endif
