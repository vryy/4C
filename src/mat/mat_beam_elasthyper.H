/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief constitutive relations for beam cross-section resultants (hyperelastic stored energy
function)


\level 3
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef MAT_BEAM_ELASTHYPER_H
#define MAT_BEAM_ELASTHYPER_H

#include "inpar_material.H"

#include "lib_parobjectfactory.H"

#include "linalg_fixedsizematrix.H"

#include <Teuchos_RCP.hpp>
#include "mat_beam_templated_material_generic.H"
#include <typeinfo>


// forward declaration
namespace DRT
{
  class ParObject;
}

namespace MAT
{
  // forward declaration
  namespace PAR
  {
    class BeamElastHyperMaterialParameterGeneric;
  }

  /// singleton for constitutive law of a beam formulation (hyperelastic stored energy function)
  template <typename T>
  class BeamElastHyperMaterialType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return typeid(this).name(); }

    static BeamElastHyperMaterialType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static BeamElastHyperMaterialType instance_;
  };



  /*---------------------------------------------------------------------------------------------*/
  /// constitutive relations for beam cross-section resultants (hyperelastic stored energy function)
  template <typename T>
  class BeamElastHyperMaterial : public BeamMaterialTemplated<T>
  {
   public:
    /// construct empty material object
    BeamElastHyperMaterial();

    /// construct the material object from given material parameters
    explicit BeamElastHyperMaterial(MAT::PAR::BeamElastHyperMaterialParameterGeneric* params);

    //! @name Packing and Unpacking
    //@{

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of parobject.H and should return it in this method.
    */
    int UniqueParObjectId() const override
    {
      return BeamElastHyperMaterialType<T>::Instance().UniqueParObjectId();
    }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    void Pack(DRT::PackBuffer& data) const override;

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

    //! @name derived access methods
    //@{

    /** \brief get type of this material
     *
     */
    INPAR::MAT::MaterialType MaterialType() const override
    {
      /* the idea is that we have a generic type of material (this class), but two
       * possible types of material parameter definitions and hence types
       * (BeamReissnerElastHyperMaterialParams and BeamReissnerElastHyperMaterialParamsByMode) */
      return INPAR::MAT::m_beam_elast_hyper_generic;
    }

    /** \brief return copy of this material object
     *
     */
    Teuchos::RCP<Material> Clone() const override
    {
      return Teuchos::rcp(new BeamElastHyperMaterial(*this));
    }

    /** \brief get quick accessible material parameter data
     *
     */
    MAT::PAR::Parameter* Parameter() const override;

    //@}


    //! @name Access methods
    //@{

    /** \brief get quick accessible material parameter data
     *
     */
    const MAT::PAR::BeamElastHyperMaterialParameterGeneric& Params() const;

    /*
     * \brief Compute axial stress contributions
     *
     *\param[out] stressM axial stress
     *
     *\param[in] CM constitutive matrix
     *
     *\param[in] Cur curvature
     */
    void EvaluateMomentContributionsToStress(LINALG::Matrix<3, 1, T>& stressM,
        const LINALG::Matrix<3, 3, T>& CM, const LINALG::Matrix<3, 1, T>& Cur) override;

    /*
     * \brief Compute axial stress contributions
     *
     *\param[out] stressN axial stress
     *
     *\param[in] CN constitutive matrix
     *
     *\param[in] Gamma triad
     */

    void EvaluateForceContributionsToStress(LINALG::Matrix<3, 1, T>& stressN,
        const LINALG::Matrix<3, 3, T>& CN, const LINALG::Matrix<3, 1, T>& Gamma) override;

    /*
     * \brief Update material-dependent variables
     */
    void ComputeConstitutiveParameter(
        LINALG::Matrix<3, 3, T>& C_N, LINALG::Matrix<3, 3, T>& C_M) override;

    /** \brief get constitutive matrix relating stress force resultants and translational strain
     *         measures, expressed w.r.t. material frame
     *
     */
    void GetConstitutiveMatrixOfForcesMaterialFrame(LINALG::Matrix<3, 3, T>& C_N) const override;

    /** \brief get constitutive matrix relating stress moment resultants and rotational strain
     *         measures, expressed w.r.t. material frame
     *
     */
    void GetConstitutiveMatrixOfMomentsMaterialFrame(LINALG::Matrix<3, 3, T>& C_M) const override;

    /** \brief get mass inertia factor with respect to translational accelerations
     *         (usually: density * cross-section area)
     */
    double GetTranslationalMassInertiaFactor() const override;

    /** \brief get mass moment of inertia tensor, expressed w.r.t. material frame
     *
     */
    void GetMassMomentOfInertiaTensorMaterialFrame(LINALG::Matrix<3, 3>& J) const override;

    /** \brief get mass moment of inertia tensor, expressed w.r.t. material frame
     */
    void GetMassMomentOfInertiaTensorMaterialFrame(
        LINALG::Matrix<3, 3, Sacado::Fad::DFad<double>>& J) const override;



    /** \brief get the radius of a circular cross-section that is ONLY to be used for evaluation of
     *         any kinds of beam interactions (contact, potentials, viscous drag forces ...)
     *
     */
    double GetInteractionRadius() const;

    void GetStiffnessMatrixOfMoments(
        LINALG::Matrix<3, 3, T>& stiffness_matrix, const LINALG::Matrix<3, 3, T>& C_M) override;

    void GetStiffnessMatrixOfForces(
        LINALG::Matrix<3, 3, T>& stiffness_matrix, const LINALG::Matrix<3, 3, T>& C_N) override;

    void Update() override{};

    void Reset() override{};

   protected:
    void SetParameter(MAT::PAR::BeamElastHyperMaterialParameterGeneric* parameter)
    {
      params_ = parameter;
    }

   private:
    /// my material parameters
    MAT::PAR::BeamElastHyperMaterialParameterGeneric* params_;
  };
}  // namespace MAT

#endif
