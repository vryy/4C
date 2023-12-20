/*----------------------------------------------------------------------------*/
/*! \file
\brief material stores parameters for ion species in electrolyte solution. The scl material is
derived for a binary electrolyte assuming a mobile ionic species in a fixed anion lattice.
Local electroneutrality is dismissed, which enables formation of Space-Charge-Layers (SCLs).

\level 2


*/
/*----------------------------------------------------------------------------*/

#ifndef BACI_MAT_SCL_H
#define BACI_MAT_SCL_H

#include "baci_config.H"

#include "baci_comm_parobjectfactory.H"
#include "baci_mat_elchsinglemat.H"

BACI_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for electrolytes including space-charge-layer formation
    class Scl : public ElchSingleMat
    {
     public:
      Scl(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{
      /// valence (= charge number)
      const double valence_;

      /// definition of transference number
      /// (by function number or implemented concentration dependence)
      const int transnrcurve_;

      /// number of parameter needed for implemented concentration dependence
      const int transnrparanum_;

      /// parameter needed for implemented concentration dependence
      const std::vector<double> transnr_;

      //! maximum concentration of species
      const double cmax_;

      //! strategy for extrapolation of diffusion coefficient
      const int extrapolation_diffussion_coeff_strategy_;

      //! limit concentration for extrapolation strategy
      const double clim_;

      //! bulk concentration i.e. anion concentration for equal transference numbers
      const double cbulk_;

      //! dieelectric susceptibility of electrolyte material
      const double susceptibility_;

      //! difference in partial molar volumes (vacancy <=> interstitial)
      const double delta_nu_;

      //! Faraday constant
      const double faraday_;

      //! vacuum Permittivity
      const double epsilon_0_;
      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;
    };

  }  // namespace PAR

  class SclType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "SclType"; }

    static SclType& Instance() { return instance_; }

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static SclType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for the material properties of an ion species in an electrolyte solution
  class Scl : public ElchSingleMat
  {
   public:
    /// construct empty material object
    Scl();

    /// construct the material object given material parameters
    explicit Scl(MAT::PAR::Scl* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override { return SclType::Instance().UniqueParObjectId(); }

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

    /// material type
    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_scl; }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override { return Teuchos::rcp(new Scl(*this)); }

    /// valence (= charge number)
    double Valence() const { return params_->valence_; }

    /// computation of the transference number based on the defined curve
    double ComputeTransferenceNumber(const double cint) const;

    /// computation of the first derivative of the transference number based on the defined curve
    double ComputeFirstDerivTrans(const double cint) const;

    double ComputeDiffusionCoefficient(
        const double concentration, const double temperature) const override;

    double ComputeConcentrationDerivativeOfDiffusionCoefficient(
        const double concentration, const double temperature) const override;

    /// computation of dielectric susceptibility (currently a constant)
    double ComputeSusceptibility() const { return params_->susceptibility_; }

    /// computation of 1/(z^2F^2) with valence of cations
    double InvValValenceFaradaySquared() const;

    /// Computation of dielectric permittivity based on dieelectric susceptibility
    double ComputePermittivity() const;

    /// Returns Value of cation concentration in the neutral bulk (= anion concentration)
    double BulkConcentration() const { return params_->cbulk_; }

    /// computation of mobility factor in linear onsager ansatz
    double ComputeOnsagerCoefficient(const double concentration, const double temperature) const;

    /// computation of the derivative of the mobility factor w.r.t to cation concentration
    double ComputeConcentrationDerivativeOfOnsagerCoefficient(
        const double concentration, const double temperature) const;

   private:
    /// return curve defining the transference number
    int TransNrCurve() const { return params_->transnrcurve_; }

    /// parameter needed for implemented concentration dependence
    const std::vector<double>& TransNrParams() const { return params_->transnr_; }

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

    /// my material parameters
    MAT::PAR::Scl* params_;
  };
}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif
