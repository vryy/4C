/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief constitutive relations for beam cross-section resultants (hyperelastic stored energy
function)


\level 3
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef MAT_BEAM_PLASTICITY_H_
#define MAT_BEAM_PLASTICITY_H_

#include "material.H"
#include "../drt_inpar/inpar_material.H"
#include "beam_elasthyper_parameter.H"
#include "beam_elasthyper.H"

#include "../drt_lib/drt_parobjectfactory.H"

#include "./beam_elasthyper.H"

#include "../linalg/linalg_fixedsizematrix.H"

#include <Teuchos_RCP.hpp>


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
      virtual double GetYieldStressN() const { return yield_stress_n_; }

      //! yield stress momentum
      virtual double GetYieldStressM() const { return yield_stress_m_; }

      //! hardening rigidity axial direction
      virtual double GetHardeningAxialRigidity() const
      {
        return isohard_modulus_n_ * GetCrossSectionArea();
      };

      //! hardening rigidity shear one direction
      virtual double GetHardeningShearRigidity2() const
      {
        return GetShearModulus() * GetCrossSectionArea() * GetShearCorrectionFactor();
      };

      //! hardening rigidity shear other direction
      virtual double GetHardeningShearRigidity3() const
      {
        return GetShearModulus() * GetCrossSectionArea() * GetShearCorrectionFactor();
      };

      //! hardening rigidity for momentum
      virtual double GetHardeningMomentalRigidity() const
      {
        return isohard_modulus_m_ * GetMomentInertia2();
      }

      //! consider torsion plasticity
      virtual bool GetTorsionPlasticity() const { return torsion_plasticity_; }
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
  class BeamElastPlasticMaterialType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return typeid(this).name(); }

    //! get instance for beam material
    static BeamElastPlasticMaterialType& Instance() { return instance_; };

    //! create material object
    virtual DRT::ParObject* Create(const std::vector<char>& data);

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

    //! @name Packing and Unpacking
    //@{

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H and should return it in this method.
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
    void Pack(DRT::PackBuffer& data) const override;

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
    void EvaluateForceContributionsToStress(LINALG::Matrix<3, 1, T>& stressN,
        const LINALG::Matrix<3, 3, T>& CN, const LINALG::Matrix<3, 1, T>& Gamma) override;

    /*!
     * \brief Compute moment stress contributions
     *
     *\param[out] stressM moment stress
     *\param[in] CM constitutive matrix
     *\param[in] Cur curvature
     */
    void EvaluateMomentContributionsToStress(LINALG::Matrix<3, 1, T>& stressM,
        const LINALG::Matrix<3, 3, T>& CM, const LINALG::Matrix<3, 1, T>& Cur) override;

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
    void GetConstitutiveMatrixOfForcesMaterialFrame(LINALG::Matrix<3, 3, T>& C_N) const override;

    /** \brief get constitutive matrix relating stress moment resultants and rotational strain
     *         measures, expressed w.r.t. material frame
     *
     */
    void GetConstitutiveMatrixOfMomentsMaterialFrame(LINALG::Matrix<3, 3, T>& C_M) const override;

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
    void GetStiffnessMatrixOfMoments(
        LINALG::Matrix<3, 3, T>& stiffM, const LINALG::Matrix<3, 3, T>& C_M) override;

    /** \brief compute stiffness matrix of forces for plastic regime
     */
    void GetStiffnessMatrixOfForces(
        LINALG::Matrix<3, 3, T>& stiffN, const LINALG::Matrix<3, 3, T>& C_N) override;

    /** \brief update the plastic strain and curvature vectors
     */
    void Update() override;

    /** \brief reset the values for current plastic strain and curvature
     */
    void Reset() override;

    /** \brief get hardening constitutive parameters depending on the type of plasticity
     */
    void ComputeConstitutiveParameter(
        LINALG::Matrix<3, 3, T>& C_N, LINALG::Matrix<3, 3, T>& C_M) override;

   protected:
    /** \brief get the constitutive matrix of forces during kinematic hardening
     *
     */
    void GetHardeningConstitutiveMatrixOfForcesMaterialFrame(LINALG::Matrix<3, 3, T>& CN_eff) const;

    /** \brief get the constitutive matrix of moments during kinematic hardening
     */
    void GetHardeningConstitutiveMatrixOfMomentsMaterialFrame(
        LINALG::Matrix<3, 3, T>& CM_eff) const;

    /** \brief returns current effective yield stress of forces depending on plastic deformation
     */
    void GetEffectiveYieldStressN(T& eff_yieldN, T init_yieldN, T CN_0, T CN_eff_0) const;

    /** \brief returns current effective yield stress of moments depending on plastic deformation
     */
    void GetEffectiveYieldStressM(T& eff_yieldM, T init_yieldM, T CM_1, T CM_eff_1) const;

   private:
    //! effective constitutive matrices forces
    LINALG::Matrix<3, 3, T> cN_eff_{true};
    //! effective constitutive matrices moments
    LINALG::Matrix<3, 3, T> cM_eff_{true};

    //! converged plastic strain vectors at GPs
    LINALG::Matrix<3, 1, T> gammaplastconv_{true};
    //! new plastic strain vectors at GPs
    LINALG::Matrix<3, 1, T> gammaplastnew_{true};
    //! accumulated plastic strain vectors at GPs
    T gammaplastaccum_{0};

    //! converged plastic curvature vectors at GPs
    LINALG::Matrix<3, 1, T> kappaplastconv_{true};
    //! new plastic curvature vectors at GPs
    LINALG::Matrix<3, 1, T> kappaplastnew_{true};
    //! accumulated plastic curvature vectors at GPs
    T kappaplastaccum_{0};

    //! effective yield force depending on accumulated plastic strain
    T effyieldstressN_{0};
    //! effective yield moment depending on accumulated plastic curvature
    T effyieldstressM_{0};

    //! norm of material plastic curvature increment
    T deltaKappaplast_{0};
    //! norm of the moment vector
    T normstressM_{0};
    //! fraction of the norm of the moment vector exceeding the current yield moment
    T deltastressM_{0};

    //! copy of material curvature K (but first entry is 0 if torsional plasticity is turned off)
    LINALG::Matrix<3, 1, T> kappa_{true};
    //! material elastic curvature
    LINALG::Matrix<3, 1, T> kappaelast_{true};
    /** copy of material elastic curvature needed to determine flow direction when computing the
     * stiffness matrix (first entry is 0 if torsional plasticity is turned off)
     */
    LINALG::Matrix<3, 1, T> kappaelastflow_{true};
    //! unit vector in the direction of the elastic curvature
    LINALG::Matrix<3, 1, T> elastic_curvature_{true};
    //! material plastic strain increment
    LINALG::Matrix<3, 1, T> deltaGammaplast_{true};
    //! fraction of the stress exceeding the current yield force
    LINALG::Matrix<3, 1, T> deltastressN_{true};

    //! axial stress
    T stressN_{0};

    //@}
  };
}  // namespace MAT

#endif
