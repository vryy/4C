/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief constitutive relations for beam cross-section resultants (hyperelastic stored energy
function)

\maintainer Maximilian Grill

\level 3
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef MAT_BEAM_ELASTHYPER_H_
#define MAT_BEAM_ELASTHYPER_H_

#include "material.H"
#include "../drt_inpar/inpar_material.H"

#include "../drt_lib/drt_parobjectfactory.H"

#include "../linalg/linalg_fixedsizematrix.H"

#include <Teuchos_RCP.hpp>


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
  class BeamElastHyperMaterialType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "BeamElastHyperMaterialType"; }

    static BeamElastHyperMaterialType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static BeamElastHyperMaterialType instance_;
  };



  /*---------------------------------------------------------------------------------------------*/
  /// constitutive relations for beam cross-section resultants (hyperelastic stored energy function)
  class BeamElastHyperMaterial : public Material
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
      top of drt_parobject.H and should return it in this method.
    */
    virtual int UniqueParObjectId() const
    {
      return BeamElastHyperMaterialType::Instance().UniqueParObjectId();
    }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    virtual void Pack(DRT::PackBuffer& data) const;

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
    virtual void Unpack(const std::vector<char>& data);

    //@}


    //! @name derived access methods
    //@{

    /** \brief get type of this material
     *
     *  \author grill
     *  \date 02/17 */
    virtual INPAR::MAT::MaterialType MaterialType() const
    {
      /* the idea is that we have a generic type of material (this class), but two
       * possible types of material parameter definitions and hence types
       * (BeamReissnerElastHyperMaterialParams and BeamReissnerElastHyperMaterialParamsByMode) */
      return INPAR::MAT::m_beam_elast_hyper_generic;
    }

    /** \brief return copy of this material object
     *
     *  \author grill
     *  \date 02/17 */
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new BeamElastHyperMaterial(*this));
    }

    /** \brief get quick accessible material parameter data
     *
     *  \author grill
     *  \date 02/17 */
    virtual MAT::PAR::Parameter* Parameter() const;

    //@}


    //! @name Access methods
    //@{

    /** \brief get quick accessible material parameter data
     *
     *  \author grill
     *  \date 02/17 */
    const MAT::PAR::BeamElastHyperMaterialParameterGeneric& Params() const;

    /** \brief throw error if pointer to parameter class is not set
     *
     *  \author grill
     *  \date 02/17 */
    void ThrowErrorIfParamsPointerIsNull() const;

    /** \brief get constitutive matrix relating stress force resultants and translational strain
     *         measures, expressed w.r.t. material frame
     *
     *  \author grill
     *  \date 02/17 */
    template <typename T>
    void GetConstitutiveMatrixOfForcesMaterialFrame(LINALG::Matrix<3, 3, T>& C_N) const;

    /** \brief get constitutive matrix relating stress moment resultants and rotational strain
     *         measures, expressed w.r.t. material frame
     *
     *  \author grill
     *  \date 02/17 */
    template <typename T>
    void GetConstitutiveMatrixOfMomentsMaterialFrame(LINALG::Matrix<3, 3, T>& C_M) const;

    /** \brief get mass inertia factor with respect to translational accelerations
     *         (usually: density * cross-section area)
     *
     *  \author grill
     *  \date 02/17 */
    double GetTranslationalMassInertiaFactor() const;

    /** \brief get mass moment of inertia tensor, expressed w.r.t. material frame
     *
     *  \author grill
     *  \date 02/17 */
    template <typename T>
    void GetMassMomentOfInertiaTensorMaterialFrame(LINALG::Matrix<3, 3, T>& J) const;


    /** \brief get the radius of a circular cross-section that is ONLY to be used for evaluation of
     *         any kinds of beam interactions (contact, potentials, viscous drag forces ...)
     *
     *  \author grill
     *  \date 02/17 */
    double GetInteractionRadius() const;

    //@}


    //! @name Evaluation methods
    //@{

    // Todo @grill implement and call this instead of calculating stress resultant in element itself
    //  /** \brief evaluate constitutive relations, i.e. compute stress resultants from given
    //  strains
    //   *
    //   *  \author grill
    //   *  \date 02/17 */
    //  void Evaluate() const;

    //@}

   private:
    /// my material parameters
    MAT::PAR::BeamElastHyperMaterialParameterGeneric* params_;
  };
}  // namespace MAT

#endif
