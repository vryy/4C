/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief data container holding all input parameters relevant for potential based beam interactions

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef BEAMPOTENTIALPARAMS_H
#define BEAMPOTENTIALPARAMS_H

#include "../drt_inpar/inpar_beampotential.H"

// forward declaration
namespace BEAMINTERACTION
{
  class BeamToBeamPotentialRuntimeVtkParams;

  /*!
   *  */
  class BeamPotentialParams
  {
   public:
    //! constructor
    BeamPotentialParams();

    //! destructor
    virtual ~BeamPotentialParams(){};

    //! initialize with the stuff coming from input file
    void Init();

    //! setup member variables
    void Setup();

    //! returns the isinit_ flag
    inline bool IsInit() const { return isinit_; }

    //! returns the issetup_ flag
    inline bool IsSetup() const { return issetup_; }

    //! asserts the init and setup status
    void ThrowErrorIfNotInitAndSetup() const;

    //! asserts the init status
    void ThrowErrorIfNotInit() const;

    inline std::vector<double> const& PotentialLawExponents() const
    {
      ThrowErrorIfNotInitAndSetup();
      return *pot_law_exponents_;
    }

    inline std::vector<double> const& PotentialLawPrefactors() const
    {
      ThrowErrorIfNotInitAndSetup();
      return *pot_law_prefactors_;
    }

    inline enum INPAR::BEAMPOTENTIAL::BeamPotentialType PotentialType() const
    {
      ThrowErrorIfNotInitAndSetup();
      return potential_type_;
    }

    inline enum INPAR::BEAMPOTENTIAL::BeamPotentialStrategy Strategy() const
    {
      ThrowErrorIfNotInitAndSetup();
      return strategy_;
    }

    inline double CutoffRadius() const
    {
      ThrowErrorIfNotInitAndSetup();
      return cutoff_radius_;
    }

    inline enum INPAR::BEAMPOTENTIAL::BeamPotentialRegularizationType RegularizationType() const
    {
      ThrowErrorIfNotInitAndSetup();
      return regularization_type_;
    }

    inline double RegularizationSeparation() const
    {
      ThrowErrorIfNotInitAndSetup();
      return regularization_separation_;
    }

    inline int NumberIntegrationSegments() const
    {
      ThrowErrorIfNotInitAndSetup();
      return num_integration_segments_;
    }

    inline int NumberGaussPoints() const
    {
      ThrowErrorIfNotInitAndSetup();
      return num_GPs_;
    }

    inline bool UseFAD() const
    {
      ThrowErrorIfNotInitAndSetup();
      return useFAD_;
    }

    //! whether to write vtk output for beam contact
    inline bool VtkRuntimeOutput() const
    {
      ThrowErrorIfNotInitAndSetup();
      return vtk_output_;
    }

    //! get the data container for parameters regarding vtk output
    inline Teuchos::RCP<const BEAMINTERACTION::BeamToBeamPotentialRuntimeVtkParams>
    GetBeamPotentialVtkParams() const
    {
      ThrowErrorIfNotInitAndSetup();
      return params_runtime_vtk_BTB_potential_;
    }

   private:
    bool isinit_;

    bool issetup_;

    //! exponents of the summands of a potential law in form of a power law
    // Todo maybe change to integer?
    Teuchos::RCP<std::vector<double>> pot_law_exponents_;

    //! prefactors of the summands of a potential law in form of a power law
    Teuchos::RCP<std::vector<double>> pot_law_prefactors_;

    //! type of applied potential (volume, surface)
    enum INPAR::BEAMPOTENTIAL::BeamPotentialType potential_type_;

    //! strategy to evaluate interaction potential
    enum INPAR::BEAMPOTENTIAL::BeamPotentialStrategy strategy_;

    //! neglect all contributions at separation larger than this cutoff radius
    double cutoff_radius_;

    //! type of regularization to use for force law at separations below specified separation
    enum INPAR::BEAMPOTENTIAL::BeamPotentialRegularizationType regularization_type_;

    //! use specified regularization type for separations smaller than this value
    double regularization_separation_;

    //! number of integration segments to be used per beam element
    int num_integration_segments_;

    //! number of Gauss points to be used per integration segment
    int num_GPs_;

    //! use automatic differentiation via FAD
    bool useFAD_;

    //! whether to write vtk output at runtime
    bool vtk_output_;

    //! data container for input parameters related to VTK output of beam contact at runtime
    Teuchos::RCP<BEAMINTERACTION::BeamToBeamPotentialRuntimeVtkParams>
        params_runtime_vtk_BTB_potential_;
  };

}  // namespace BEAMINTERACTION

#endif
