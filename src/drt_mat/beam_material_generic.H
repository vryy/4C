/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief interface for constitutive relations for beam cross-section resultants


\level 1
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef MAT_BEAM_H_
#define MAT_BEAM_H_

#include "material.H"
#include "../drt_inpar/inpar_material.H"

#include "../drt_lib/drt_parobjectfactory.H"

#include "../linalg/linalg_fixedsizematrix.H"

#include <Teuchos_RCP.hpp>
#include <Sacado.hpp>


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

  /*---------------------------------------------------------------------------------------------*/
  /// constitutive relations for beam cross-section resultants (hyperelastic stored energy function)
  class BeamMaterial : public Material
  {
   public:
    //! @name Packing and Unpacking
    //@{

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H and should return it in this method.
    */
    virtual int UniqueParObjectId() const = 0;

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    virtual void Pack(DRT::PackBuffer& data) const = 0;

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
    virtual void Unpack(const std::vector<char>& data) = 0;

    //@}

    //! @name derived access methods
    //@{

    /** \brief get type of this material
     *
     */
    virtual INPAR::MAT::MaterialType MaterialType() const = 0;

    /** \brief return copy of this material object
     *
     */
    virtual Teuchos::RCP<Material> Clone() const = 0;

    /** \brief get quick accessible material parameter data
     *
     */
    virtual MAT::PAR::Parameter* Parameter() const = 0;

    //@}


    //! @name Access methods
    //@{

    /** \brief get quick accessible material parameter data
     *
     */
    virtual const MAT::PAR::BeamElastHyperMaterialParameterGeneric& Params() const = 0;

    /** \brief throw error if pointer to parameter class is not set
     *
     */
    virtual void ThrowErrorIfParamsPointerIsNull() const = 0;

    /** \brief get the radius of a circular cross-section that is ONLY to be used for evaluation of
     *         any kinds of beam interactions (contact, potentials, viscous drag forces ...)
     *
     */
    virtual double GetInteractionRadius() const = 0;

    /** \brief get mass inertia factor with respect to translational accelerations
     *         (usually: density * cross-section area)
     *
     */
    virtual double GetTranslationalMassInertiaFactor() const = 0;

    /** \brief get mass moment of inertia tensor, expressed w.r.t. material frame
     *
     */
    virtual void GetMassMomentOfInertiaTensorMaterialFrame(LINALG::Matrix<3, 3>& J) const = 0;

    /** \brief get mass moment of inertia tensor, expressed w.r.t. material frame
     *
     */
    virtual void GetMassMomentOfInertiaTensorMaterialFrame(
        LINALG::Matrix<3, 3, Sacado::Fad::DFad<double>>& J) const = 0;

    /** \brief Update all material related variables at the end of a time step
     *
     */
    virtual void Update() const = 0;

    /** \brief Resets all material related variables i.e. in case of adaptive time stepping
     *
     */
    virtual void Reset() const = 0;
  };
}  // namespace MAT

#endif
