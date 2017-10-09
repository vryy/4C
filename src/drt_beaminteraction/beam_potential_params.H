/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_potential_params.H

\brief data container holding all input parameters relevant for potential based beam interactions

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef BEAMPOTENTIALPARAMS_H
#define BEAMPOTENTIALPARAMS_H

#include "../drt_inpar/inpar_beampotential.H"

// forward declaration
namespace BEAMINTERACTION {

  /*!
   *  */
  class BeamPotentialParams
  {
  public:
    //! constructor
    BeamPotentialParams();

    //! destructor
    virtual ~BeamPotentialParams() {};

    //! initialize with the stuff coming from input file
    void Init();

    //! setup member variables
    void Setup();

    //! returns the isinit_ flag
    inline bool IsInit() const
    {
      return isinit_;
    }

    //! returns the issetup_ flag
    inline bool IsSetup() const
    {
      return issetup_;
    }

    //! asserts the init and setup status
    void ThrowErrorIfNotInitAndSetup() const;

    //! asserts the init status
    void ThrowErrorIfNotInit() const;

    inline std::vector<double> const& PotentialLawExponents() const
    {
      ThrowErrorIfNotInit();
      return *pot_law_exponents_;
    }

    inline std::vector<double> const& PotentialLawPrefactors() const
    {
      ThrowErrorIfNotInit();
      return *pot_law_prefactors_;
    }

    inline enum INPAR::BEAMPOTENTIAL::BeamPotentialType PotentialType() const
    {
      return potential_type_;
    }

    inline enum INPAR::BEAMPOTENTIAL::BeamPotentialStrategy Strategy() const
    {
      return strategy_;
    }

    inline int NumberGaussPoints() const
    {
      return num_GPs_;
    }

    inline bool UseFAD() const
    {
      return useFAD_;
    }


  private:

    bool isinit_;

    bool issetup_;

    //! exponents of the summands of a potential law in form of a power law
    // Todo maybe change to integer?
    Teuchos::RCP<std::vector<double> > pot_law_exponents_;

    //! prefactors of the summands of a potential law in form of a power law
    Teuchos::RCP<std::vector<double> > pot_law_prefactors_;

    //! type of applied potential (volume, surface)
    enum INPAR::BEAMPOTENTIAL::BeamPotentialType potential_type_;

    //! strategy to evaluate interaction potential
    enum INPAR::BEAMPOTENTIAL::BeamPotentialStrategy strategy_;

    //! number of Gauss points to be used per beam element
    int num_GPs_;

    //! use automatic differentiation via FAD
    bool useFAD_;

  };

}

#endif
