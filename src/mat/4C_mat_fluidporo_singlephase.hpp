/*----------------------------------------------------------------------*/
/*! \file
 \brief single phase material for multiphase porous flow

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_FLUIDPORO_SINGLEPHASE_HPP
#define FOUR_C_MAT_FLUIDPORO_SINGLEPHASE_HPP



/*---------------------------------------------------------------------*
 | headers                                                              |
 *---------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_fluidporo_relpermeability_law.hpp"
#include "4C_mat_fluidporo_viscosity_law.hpp"
#include "4C_mat_material.hpp"
#include "4C_mat_par_parameter.hpp"
#include "4C_mat_poro_density_law.hpp"

FOUR_C_NAMESPACE_OPEN


/*---------------------------------------------------------------------*
 *---------------------------------------------------------------------*/
namespace MAT
{
  namespace PAR
  {
    // forward declaration
    class FluidPoroPhaseDof;

    /*----------------------------------------------------------------------*/
    /// material parameters for a single phase of porous multiphase fluid
    ///
    /// This object exists only once for each read fluid.
    class FluidPoroSinglePhase : public Parameter
    {
     public:
      /// standard constructor
      FluidPoroSinglePhase(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      /// initialize
      void Initialize();

      /// @name material parameters
      //@{

      /// density
      const double density_;
      /// type of degree of freedom
      FluidPoroPhaseDof* phasedof_;
      /// density law
      PoroDensityLaw* densitylaw_;
      /// permeability law
      FluidPoroRelPermeabilityLaw* relpermeabilitylaw_;
      /// viscosity law
      FluidPoroViscosityLaw* viscositylaw_;

      //@}

     private:
      bool isinit_;

    };  // class FluidPoroSinglePhase

    /*----------------------------------------------------------------------*/
    /// material parameters for a single volfrac of porous multiphase fluid
    ///
    /// This object exists only once for each read fluid.
    class FluidPoroSingleVolFrac : public Parameter
    {
     public:
      /// standard constructor
      FluidPoroSingleVolFrac(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      /// initialize
      void Initialize();

      /// @name material parameters
      //@{

      /// density
      const double density_;
      /// bulk diffusivity
      const double diffusivity_;
      /// do we have additional scalar dependent flux
      const bool scalardependentflux_;
      /// number of scalars
      const int numscal_;
      /// their diffusivities
      const std::vector<double> scalardiffs_;
      /// constant omega_half for receptor kinetic law
      const std::vector<double> omega_half_;

      //@}

     private:
      bool isinit_;

    };  // class FluidPoroSingleVolFrac

    /*----------------------------------------------------------------------*/
    /// material parameters for a single volfrac pressure of porous multiphase fluid
    ///
    /// This object exists only once for each read fluid.
    class FluidPoroVolFracPressure : public Parameter
    {
     public:
      /// standard constructor
      FluidPoroVolFracPressure(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      /// initialize
      void Initialize();

      /// @name material parameters
      //@{

      /// permeability
      const double permeability_;

      /// minimum volume fraction
      const double min_volfrac_;

      /// viscosity law
      FluidPoroViscosityLaw* viscositylaw_;

      //@}

     private:
      bool isinit_;

    };  // class FluidPoroVolFracPressure



  }  // namespace PAR

  /*----------------------------------------------------------------------*
   | instance access method                                   vuong 08/16 |
   *----------------------------------------------------------------------*/
  class FluidPoroSinglePhaseType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "FluidPoroSinglePhaseType"; }

    static FluidPoroSinglePhaseType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static FluidPoroSinglePhaseType instance_;
  };

  /*----------------------------------------------------------------------*
   | instance access method                                   vuong 08/16 |
   *----------------------------------------------------------------------*/
  class FluidPoroSingleVolFracType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "FluidPoroSingleVolFracType"; }

    static FluidPoroSingleVolFracType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static FluidPoroSingleVolFracType instance_;
  };

  /*----------------------------------------------------------------------*
   | instance access method                              kremheller 02/18 |
   *----------------------------------------------------------------------*/
  class FluidPoroVolFracPressureType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "FluidPoroVolFracPressureType"; }

    static FluidPoroVolFracPressureType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static FluidPoroVolFracPressureType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Base class for a single porous fluid phase within multiphase porous flow
  ///
  /// This object exists (several times) at every element
  class FluidPoroSinglePhaseBase : public Material
  {
   public:
    /// construct empty material object
    FluidPoroSinglePhaseBase(){};

    /// initialize
    virtual void Initialize() = 0;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for a single porous fluid phase within multiphase porous flow
  ///
  /// This object exists (several times) at every element
  class FluidPoroSinglePhase : public FluidPoroSinglePhaseBase
  {
   public:
    /// construct empty material object
    FluidPoroSinglePhase();

    /// construct the material object given material parameters
    explicit FluidPoroSinglePhase(MAT::PAR::FluidPoroSinglePhase* params);

    //! @name Packing and Unpacking

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of parobject.H (this file) and should return it in this method.
     */
    int UniqueParObjectId() const override
    {
      return FluidPoroSinglePhaseType::Instance().UniqueParObjectId();
    }

    /*!
     \brief Pack this class so it can be communicated

     Resizes the vector data and stores all information of a class in it.
     The first information to be stored in data has to be the
     unique parobject id delivered by UniqueParObjectId() which will then
     identify the exact class on the receiving processor.

     \param data (in/out): char vector to store class information
     */
    void Pack(CORE::COMM::PackBuffer& data) const override;

    /*!
     \brief Unpack data from a char vector into this class

     The vector data contains all information to rebuild the
     exact copy of an instance of a class on a different processor.
     The first entry in data has to be an integer which is the unique
     parobject id defined at the top of this file and delivered by
     UniqueParObjectId().

     \param data (in) : vector storing all data to be unpacked into this
     instance.
     */
    void Unpack(const std::vector<char>& data) override;

    //@}

    /// initialize
    void Initialize() override;

    /// material type
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_fluidporo_singlephase;
    }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new FluidPoroSinglePhase(*this));
    }

    /// return density
    double Density() const override { return params_->density_; }

    /// return relative permeability
    double RelPermeability(const double saturation) const
    {
      return params_->relpermeabilitylaw_->GetRelPermeability(saturation);
    }

    // check for constant relative permeability
    bool HasConstantRelPermeability() const
    {
      return params_->relpermeabilitylaw_->HasConstantRelPermeability();
    }

    /// return derivative of relative permeabilty w.r.t. saturation
    double EvaluateDerivOfRelPermeabilityWrtSaturation(const double saturation) const
    {
      return params_->relpermeabilitylaw_->GetDerivOfRelPermeabilityWrtSaturation(saturation);
    }

    // check for constant viscosity
    bool HasConstantViscosity() const { return params_->viscositylaw_->HasConstantViscosity(); }

    /// return viscosity
    double Viscosity(const double abspressgrad) const
    {
      return params_->viscositylaw_->GetViscosity(abspressgrad);
    }

    /// return derivative of viscosity w.r.t. to absolute value of pressure gradient
    double ViscosityDeriv(const double abspressgrad) const
    {
      return params_->viscositylaw_->GetDerivOfViscosityWrtAbsPressGrad(abspressgrad);
    }

    /// return inverse bulk modulus (compressibility)
    double InvBulkmodulus() const { return params_->densitylaw_->InvBulkmodulus(); }

    /// return type of degree of freedom
    CORE::Materials::MaterialType PoroDofType() const;

    /// return type of phase law
    CORE::Materials::MaterialType PoroPhaseLawType() const;

    /// mark dofs associated with this phase in a given row (=numphase) in a matrix
    void FillDoFMatrix(CORE::LINALG::SerialDenseMatrix& dofmat, int numphase) const;

    /// evaluate saturation of the phase
    double EvaluateSaturation(
        int phasenum, const std::vector<double>& state, const std::vector<double>& pressure) const;

    /// evaluate the generalized(!) pressure of this phase
    double EvaluateGenPressure(int phasenum, const std::vector<double>& state) const;

    //! evaluate derivative of saturation with respect to pressure
    double EvaluateDerivOfSaturationWrtPressure(
        int phasenum, int doftoderive, const std::vector<double>& pressure) const;

    //! evaluate 2nd derivative of saturation with respect to pressure
    double EvaluateSecondDerivOfSaturationWrtPressure(int phasenum, int firstdoftoderive,
        int seconddoftoderive, const std::vector<double>& pressure) const;

    //! evaluate derivative of degree of freedom with respect to pressure
    double EvaluateDerivOfDofWrtPressure(
        int phasenum, int doftoderive, const std::vector<double>& state) const;

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

   private:
    /// my material parameters
    MAT::PAR::FluidPoroSinglePhase* params_;
  };


  /*----------------------------------------------------------------------*/
  /// Wrapper for a single volfrac within multiphase porous flow
  ///
  /// This object exists (several times) at every element
  class FluidPoroSingleVolFrac : public FluidPoroSinglePhaseBase
  {
   public:
    /// construct empty material object
    FluidPoroSingleVolFrac();

    /// construct the material object given material parameters
    explicit FluidPoroSingleVolFrac(MAT::PAR::FluidPoroSingleVolFrac* params);

    //! @name Packing and Unpacking

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of parobject.H (this file) and should return it in this method.
     */
    int UniqueParObjectId() const override
    {
      return FluidPoroSingleVolFracType::Instance().UniqueParObjectId();
    }

    /*!
     \brief Pack this class so it can be communicated

     Resizes the vector data and stores all information of a class in it.
     The first information to be stored in data has to be the
     unique parobject id delivered by UniqueParObjectId() which will then
     identify the exact class on the receiving processor.

     \param data (in/out): char vector to store class information
     */
    void Pack(CORE::COMM::PackBuffer& data) const override;

    /*!
     \brief Unpack data from a char vector into this class

     The vector data contains all information to rebuild the
     exact copy of an instance of a class on a different processor.
     The first entry in data has to be an integer which is the unique
     parobject id defined at the top of this file and delivered by
     UniqueParObjectId().

     \param data (in) : vector storing all data to be unpacked into this
     instance.
     */
    void Unpack(const std::vector<char>& data) override;

    //@}

    /// initialize
    void Initialize() override;

    /// material type
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_fluidporo_singlevolfrac;
    }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new FluidPoroSingleVolFrac(*this));
    }

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

    /// return density
    double Density() const override { return params_->density_; }

    /// return diffusivity
    double Diffusivity() const { return params_->diffusivity_; }

    /// return number of scalars
    int NumScal() const { return params_->numscal_; }

    /// return scalardependentflux_
    bool HasAddScalarDependentFlux() const { return params_->scalardependentflux_; }

    /// return diffusivities for scalar-dependent flux
    std::vector<double> ScalarDiffs() const { return params_->scalardiffs_; }

    /// return omega_half for scalar-dependent flux
    std::vector<double> OmegaHalf() const { return params_->omega_half_; }

   private:
    /// my material parameters
    MAT::PAR::FluidPoroSingleVolFrac* params_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for a single volfrac pressure within multiphase porous flow
  ///
  /// This object exists (several times) at every element
  class FluidPoroVolFracPressure : public FluidPoroSinglePhaseBase
  {
   public:
    /// construct empty material object
    FluidPoroVolFracPressure();

    /// construct the material object given material parameters
    explicit FluidPoroVolFracPressure(MAT::PAR::FluidPoroVolFracPressure* params);

    //! @name Packing and Unpacking

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of parobject.H (this file) and should return it in this method.
     */
    int UniqueParObjectId() const override
    {
      return FluidPoroVolFracPressureType::Instance().UniqueParObjectId();
    }

    /*!
     \brief Pack this class so it can be communicated

     Resizes the vector data and stores all information of a class in it.
     The first information to be stored in data has to be the
     unique parobject id delivered by UniqueParObjectId() which will then
     identify the exact class on the receiving processor.

     \param data (in/out): char vector to store class information
     */
    void Pack(CORE::COMM::PackBuffer& data) const override;

    /*!
     \brief Unpack data from a char vector into this class

     The vector data contains all information to rebuild the
     exact copy of an instance of a class on a different processor.
     The first entry in data has to be an integer which is the unique
     parobject id defined at the top of this file and delivered by
     UniqueParObjectId().

     \param data (in) : vector storing all data to be unpacked into this
     instance.
     */
    void Unpack(const std::vector<char>& data) override;

    //@}

    /// initialize
    void Initialize() override;

    /// material type
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_fluidporo_volfracpressure;
    }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new FluidPoroVolFracPressure(*this));
    }

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

    /// return permeability
    double Permeability() const { return params_->permeability_; }

    /// return minimum volume fraction
    double MinVolFrac() const { return params_->min_volfrac_; }

    // check for constant viscosity
    bool HasConstantViscosity() const { return params_->viscositylaw_->HasConstantViscosity(); }

    /// return viscosity
    double Viscosity(const double abspressgrad) const
    {
      return params_->viscositylaw_->GetViscosity(abspressgrad);
    }

    /// return derivative of viscosity w.r.t. to absolute value of pressure gradient
    double ViscosityDeriv(const double abspressgrad) const
    {
      return params_->viscositylaw_->GetDerivOfViscosityWrtAbsPressGrad(abspressgrad);
    }

   private:
    /// my material parameters
    MAT::PAR::FluidPoroVolFracPressure* params_;
  };

}  // namespace MAT


FOUR_C_NAMESPACE_CLOSE

#endif
