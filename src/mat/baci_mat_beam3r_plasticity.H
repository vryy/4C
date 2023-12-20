/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief constitutive relations for beam cross-section resultants (hyperelastic stored energy
function)


\level 3
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef BACI_MAT_BEAM3R_PLASTICITY_H
#define BACI_MAT_BEAM3R_PLASTICITY_H

#include "baci_config.H"

#include "baci_comm_parobjectfactory.H"
#include "baci_inpar_material.H"
#include "baci_linalg_fixedsizematrix.H"
#include "baci_mat_beam_elasthyper.H"
#include "baci_mat_beam_elasthyper_parameter.H"
#include "baci_mat_material.H"

#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN


// forward declaration
namespace DRT
{
  class ParObject;
}

namespace MAT
{
  namespace PAR
  {
    /*-------------------------------------------------------------------------------------------*/
    /** constitutive parameters for a Reissner beam formulation (hyperelastic stored energy
     * function)
     */
    class BeamReissnerElastPlasticMaterialParams : public BeamReissnerElastHyperMaterialParams
    {
     public:
      //! standard constructor
      BeamReissnerElastPlasticMaterialParams(Teuchos::RCP<MAT::PAR::Material> matdata);

      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      //! @name Access to plasticity parameters
      //@{

      //! yield stress for forces
      double GetYieldStressN() const override { return yield_stress_n_; }

      //! yield stress momentum
      double GetYieldStressM() const override { return yield_stress_m_; }

      //! hardening rigidity axial direction
      double GetHardeningAxialRigidity() const override
      {
        return isohard_modulus_n_ * GetCrossSectionArea();
      };

      //! hardening rigidity shear one direction
      double GetHardeningShearRigidity2() const override
      {
        return GetShearModulus() * GetCrossSectionArea() * GetShearCorrectionFactor();
      };

      //! hardening rigidity shear other direction
      double GetHardeningShearRigidity3() const override
      {
        return GetShearModulus() * GetCrossSectionArea() * GetShearCorrectionFactor();
      };

      //! hardening rigidity for momentum
      double GetHardeningMomentalRigidity() const override
      {
        return isohard_modulus_m_ * GetMomentInertia2();
      }

      //! consider torsion plasticity
      bool GetTorsionPlasticity() const override { return torsion_plasticity_; }
      //@}

     private:
      //! @name plasticity parameters
      //@{
      //! Yield stress of forces
      double yield_stress_n_;
      //! Yield stress of moments
      double yield_stress_m_;
      //! Isotropic hardening modulus of forces
      const double isohard_modulus_n_;
      //! Isotropic hardening modulus of moments
      const double isohard_modulus_m_;
      //! defines whether torsional moment contributes to plasticity
      const bool torsion_plasticity_;
      //@}
    };
  }  // namespace PAR

  //! singleton for constitutive law of a beam formulation (hyperelastic stored energy function)
  template <typename T>
  class BeamElastPlasticMaterialType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return typeid(this).name(); }

    //! get instance for beam material
    static BeamElastPlasticMaterialType& Instance() { return instance_; };

    //! create material object
    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static BeamElastPlasticMaterialType instance_;
  };

  /*---------------------------------------------------------------------------------------------*/
  //! constitutive relations for beam cross-section resultants (hyperelastic stored energy function)
  template <typename T>
  class BeamPlasticMaterial : public BeamElastHyperMaterial<T>
  {
   public:
    //! construct empty material object
    BeamPlasticMaterial() = default;

    //! construct the material object from given material parameters
    explicit BeamPlasticMaterial(MAT::PAR::BeamReissnerElastPlasticMaterialParams* params);

    /**
     * \brief Initialize and setup element specific variables
     *
     */
    void Setup(int numgp_force, int numgp_moment) override;

    //! @name Packing and Unpacking
    //@{

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return BeamElastPlasticMaterialType<T>::Instance().UniqueParObjectId();
    }

    /*!
     * \brief Pack this class so it can be communicated
     *
     * Resizes the vector data and stores all information of a class in it. The first information
     * to be stored in data has to be the unique parobject id delivered by UniqueParObjectId() which
     * will then identify the exact class on the receiving processor.
     *
     * @param data (in/out): char vector to store class information
     */
    void Pack(CORE::COMM::PackBuffer& data) const override;

    /*!
     * \brief Unpack data from a char vector into this class
     *
     * The vector data contains all information to rebuild the exact copy of an instance of a class
     * on a different processor. The first entry in data hast to be an integer which is the unique
     * parobject id defined at the top of this file and delivered by UniqueParObjetId().
     *
     * @param data (in) : vector storing all data to be unpacked into this instance
     */
    void Unpack(const std::vector<char>& data) override;

    //@}

    //! @name Stress contributions
    //@{

    /*!
     * \brief Compute axial stress contributions
     *
     *\param[out] stressN axial stress
     *\param[in] CN constitutive matrix
     *\param[in] Gamma strain
     */
    void EvaluateForceContributionsToStress(CORE::LINALG::Matrix<3, 1, T>& stressN,
        const CORE::LINALG::Matrix<3, 3, T>& CN, const CORE::LINALG::Matrix<3, 1, T>& Gamma,
        const unsigned int gp) override;

    /*!
     * \brief Compute moment stress contributions
     *
     *\param[out] stressM moment stress
     *\param[in] CM constitutive matrix
     *\param[in] Cur curvature
     */
    void EvaluateMomentContributionsToStress(CORE::LINALG::Matrix<3, 1, T>& stressM,
        const CORE::LINALG::Matrix<3, 3, T>& CM, const CORE::LINALG::Matrix<3, 1, T>& Cur,
        const unsigned int gp) override;

    /** \brief return copy of this material object
     */
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new BeamPlasticMaterial(*this));
    }

    //@}

    //! @name Constitutive relations
    //@{

    /** \brief get constitutive matrix relating stress force resultants and translational strain
     *         measures, expressed w.r.t. material frame
     */
    void GetConstitutiveMatrixOfForcesMaterialFrame(
        CORE::LINALG::Matrix<3, 3, T>& C_N) const override;

    /** \brief get constitutive matrix relating stress moment resultants and rotational strain
     *         measures, expressed w.r.t. material frame
     *
     */
    void GetConstitutiveMatrixOfMomentsMaterialFrame(
        CORE::LINALG::Matrix<3, 3, T>& C_M) const override;

    /** \brief get mass inertia factor with respect to translational accelerations
     *         (usually: density * cross-section area)
     *
     */
    double GetTranslationalMassInertiaFactor() const override;

    /** \brief get the radius of a circular cross-section that is ONLY to be used for evaluation of
     *         any kinds of beam interactions (contact, potentials, viscous drag forces ...)
     *
     */
    double GetInteractionRadius() const override;

    /** \brief compute stiffness matrix of moments for plastic regime
     */
    void GetStiffnessMatrixOfMoments(CORE::LINALG::Matrix<3, 3, T>& stiffM,
        const CORE::LINALG::Matrix<3, 3, T>& C_M, const int gp) override;

    /** \brief compute stiffness matrix of forces for plastic regime
     */
    void GetStiffnessMatrixOfForces(CORE::LINALG::Matrix<3, 3, T>& stiffN,
        const CORE::LINALG::Matrix<3, 3, T>& C_N, const int gp) override;

    /** \brief update the plastic strain and curvature vectors
     */
    void Update() override;

    /** \brief reset the values for current plastic strain and curvature
     */
    void Reset() override;

    /** \brief get hardening constitutive parameters depending on the type of plasticity
     */
    void ComputeConstitutiveParameter(
        CORE::LINALG::Matrix<3, 3, T>& C_N, CORE::LINALG::Matrix<3, 3, T>& C_M) override;

   protected:
    /** \brief get the constitutive matrix of forces during kinematic hardening
     *
     */
    void GetHardeningConstitutiveMatrixOfForcesMaterialFrame(
        CORE::LINALG::Matrix<3, 3, T>& CN_eff) const;

    /** \brief get the constitutive matrix of moments during kinematic hardening
     */
    void GetHardeningConstitutiveMatrixOfMomentsMaterialFrame(
        CORE::LINALG::Matrix<3, 3, T>& CM_eff) const;

    /** \brief returns current effective yield stress of forces depending on plastic deformation
     */
    void GetEffectiveYieldStressN(
        T& eff_yieldN, T init_yieldN, T CN_0, T CN_eff_0, const unsigned int gp) const;

    /** \brief returns current effective yield stress of moments depending on plastic deformation
     */
    void GetEffectiveYieldStressM(
        T& eff_yieldM, T init_yieldM, T CM_1, T CM_eff_1, const unsigned int gp) const;

   private:
    //! effective constitutive matrices forces
    std::vector<CORE::LINALG::Matrix<3, 3, T>> cN_eff_;
    //! effective constitutive matrices moments
    std::vector<CORE::LINALG::Matrix<3, 3, T>> cM_eff_;

    //! converged plastic strain vectors at GPs
    std::vector<CORE::LINALG::Matrix<3, 1, T>> gammaplastconv_;
    //! new plastic strain vectors at GPs
    std::vector<CORE::LINALG::Matrix<3, 1, T>> gammaplastnew_;
    //! accumulated plastic strain vectors at GPs
    std::vector<T> gammaplastaccum_;

    //! converged plastic curvature vectors at GPs
    std::vector<CORE::LINALG::Matrix<3, 1, T>> kappaplastconv_;
    //! new plastic curvature vectors at GPs
    std::vector<CORE::LINALG::Matrix<3, 1, T>> kappaplastnew_;
    //! accumulated plastic curvature vectors at GPs
    std::vector<T> kappaplastaccum_;

    //! effective yield force depending on accumulated plastic strain
    std::vector<T> effyieldstressN_;
    //! effective yield moment depending on accumulated plastic curvature
    std::vector<T> effyieldstressM_;

    //! norm of material plastic curvature increment
    std::vector<T> deltaKappaplast_;
    //! norm of the moment vector
    std::vector<T> normstressM_;
    //! fraction of the norm of the moment vector exceeding the current yield moment
    std::vector<T> deltastressM_;

    //! material elastic curvature
    std::vector<CORE::LINALG::Matrix<3, 1, T>> kappaelast_;
    /** copy of material elastic curvature needed to determine flow direction when computing the
     * stiffness matrix (first entry is 0 if torsional plasticity is turned off)
     */
    std::vector<CORE::LINALG::Matrix<3, 1, T>> kappaelastflow_;
    //! unit vector in the direction of the elastic curvature
    std::vector<CORE::LINALG::Matrix<3, 1, T>> elastic_curvature_;
    //! material plastic strain increment
    std::vector<CORE::LINALG::Matrix<3, 1, T>> deltaGammaplast_;
    //! fraction of the stress exceeding the current yield force
    std::vector<CORE::LINALG::Matrix<3, 1, T>> deltastressN_;

    //! axial stress
    std::vector<T> stressN_;

    /// Number of integration points for forces
    unsigned int numgp_force_;

    /// Number of integration points for moments
    unsigned int numgp_moment_;
    //@}
  };
}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif
