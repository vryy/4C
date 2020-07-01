/*! \file
\brief A thermo-mechnical pseudo-solid material used for additive manufacturing applications

\level 3

*/

#ifndef AM_THREEPHASE_THERMOMECH_H
#define AM_THREEPHASE_THERMOMECH_H

#include "material_thermomechanical.H"
#include "matpar_material.H"

namespace MAT
{
  class Consolidation;

  class FourierVar;
}  // namespace MAT

namespace MAT
{
  namespace PAR
  {
    //! @brief material parameters for the three phase pseudo solid material
    //!
    //! <h3>Input line</h3>
    //! MAT 1 MAT_Struct_ThrStVenantK YOUNG 400 NUE 0.3 DENS 1 THEXPANS 1 INITTEMP 20
    class ThermoMechThreePhase : public Parameter
    {
     public:
      //! construct parameter instance from input data
      explicit ThermoMechThreePhase(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! destructor
      ~ThermoMechThreePhase() override { ; }

      //! @name material parameters
      //! @{

      //! functions for Youngs' modulus
      const std::vector<int> youngsfunct_;
      //! Possion's ratio \f$ \nu \f$
      const double poissonratio_;
      //! mass density \f$ \rho \f$
      const double density_;
      //! functions for thermal expansion
      const std::vector<int> thermexpansfunct_;
      //! reference temperature for thermal expansion
      const double thetaref_;
      //! id of thermal fourier material wiht three phases
      const int fouriervarmat_;
      //! id of material to use for consolidation tracking
      const int consolmat_;

      //! @}

      //! create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;
    };
  }  // namespace PAR

  class ThermoMechThreePhaseType : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "ThermoMechThreePhaseType"; }

    static ThermoMechThreePhaseType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ThermoMechThreePhaseType instance_;
  };

  /*!
   * Thermo-mechanical Material with three phases (powder, melt and solid)
   *
   * All three phases are modeled as solids with varying material parameters. The material can
   * determine the correct parameters by tracking the temperature history.
   */
  class ThermoMechThreePhase : public ThermoMechanicalMaterial
  {
   public:
    //! construct empty material object
    ThermoMechThreePhase();

    //! default copy constructor
    ThermoMechThreePhase(const ThermoMechThreePhase& other) = default;

    //! construct material object from parameters
    explicit ThermoMechThreePhase(PAR::ThermoMechThreePhase* params);


    int UniqueParObjectId() const override
    {
      return ThermoMechThreePhaseType::Instance().UniqueParObjectId();
    }

    INPAR::MAT::MaterialType MaterialType() const override
    {
      return INPAR::MAT::m_thermomechthreephase;
    }

    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new ThermoMechThreePhase(*this));
    }

    void Reinit(const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* glstrain,
        double temperature, unsigned gp) override;


    void Reinit(double temperature, unsigned gp) override;

    void Setup(int numgp, DRT::INPUT::LineDefinition* linedef) override;

    void Evaluate(const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* glstrain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat,
        int gp, int eleGID) override;

    void Evaluate(const LINALG::Matrix<3, 1>& gradtemp, LINALG::Matrix<3, 3>& cmat,
        LINALG::Matrix<3, 1>& heatflux) const override;
    void Evaluate(const LINALG::Matrix<2, 1>& gradtemp, LINALG::Matrix<2, 2>& cmat,
        LINALG::Matrix<2, 1>& heatflux) const override;
    void Evaluate(const LINALG::Matrix<1, 1>& gradtemp, LINALG::Matrix<1, 1>& cmat,
        LINALG::Matrix<1, 1>& heatflux) const override;

    void ConductivityDerivT(LINALG::Matrix<3, 3>& dCondDT) const override;

    void ConductivityDerivT(LINALG::Matrix<2, 2>& dCondDT) const override;

    void ConductivityDerivT(LINALG::Matrix<1, 1>& dCondDT) const override;

    double CapacityDerivT() const override;

    void GetdSdT(LINALG::Matrix<6, 1>* dS_dT) override;

    double Capacity() const override;

    double Density() const override;

    void StressTemperatureModulusAndDeriv(
        LINALG::Matrix<6, 1>& stm, LINALG::Matrix<6, 1>& stm_dT) override;

    PAR::Parameter* Parameter() const override { return params_; }

    void Pack(DRT::PackBuffer& data) const override;
    void Unpack(const std::vector<char>& data) override;

    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (!(kinem == INPAR::STR::kinem_linear || kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    void ResetCurrentState() override;

    void CommitCurrentState() override;

   private:
    //! my material parameters
    MAT::PAR::ThermoMechThreePhase* params_;

    //! pointer to the consolidation material managing the evaluation of functions
    Teuchos::RCP<MAT::Consolidation> consolidation_;

    //! pointer to the FourierVar material for purely thermal stuff
    Teuchos::RCP<MAT::FourierVar> fouriervar_;


    void SetupCmat(LINALG::Matrix<6, 6>& cmat);
    void SetupCmat_dT(LINALG::Matrix<6, 6>& derivcmat);

    double STModulus();
    double STModulus_dT();

    inline double GetMaterialParameter(const std::vector<int>& functions) const;
    inline double GetMaterialParameterThermalDerivative(const std::vector<int>& functions) const;

    //! Create thermo material member from  input parameters
    void CreateThermoMaterial();

    //! current temperature (set by Reinit())
    double currentTemperature_{};

    //! current Gauss point (set by Reinit())
    unsigned currentGp_{};

    //! current Green-Lagrange strain
    const LINALG::Matrix<6, 1>* currentGlstrain_{};
  };

}  // namespace MAT



#endif
