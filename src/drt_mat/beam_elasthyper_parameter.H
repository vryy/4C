/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_elasthyper_parameter.H

\brief constitutive parameters for a beam material based on hyperelastic stored energy function

\maintainer Maximilian Grill

\level 3
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef MAT_BEAM_ELASTHYPER_PARAMETER_H_
#define MAT_BEAM_ELASTHYPER_PARAMETER_H_

#include "matpar_parameter.H"


namespace MAT
{
namespace PAR
{
// forward declaration
class Material;


/** \brief determine shear modulus which is either given directly or via Young's modulus and
 *         Poisson's ratio
 *
 *  \author grill
 *  \date 02/17 */
double DetermineShearModulus(const Teuchos::RCP<MAT::PAR::Material>& matdata);

/** \brief determine default value for interaction radius from area moment of inertia just in case
 *         that no value was explicitly specified
 *
 *  \author grill
 *  \date 02/17 */
double DetermineDefaultInteractionRadius(const Teuchos::RCP<MAT::PAR::Material>& matdata);

/** \brief determine default value for interaction radius from area moment of inertia just in case
 *         that no value was explicitly specified: isotropic case, i.e. only one moment of inertia
 *
 *  \author grill
 *  \date 02/17 */
double DetermineDefaultInteractionRadiusIsotropic(const Teuchos::RCP<MAT::PAR::Material>& matdata);

/*-------------------------------------------------------------------------------------------*/
/// (generic) constitutive parameters for a beam material based on hyperelastic stored energy function
class BeamElastHyperMaterialParameterGeneric : public Parameter
{
public:

  /// standard constructor
  BeamElastHyperMaterialParameterGeneric(
      Teuchos::RCP<MAT::PAR::Material> matdata
      );

  /// destructor
  virtual ~BeamElastHyperMaterialParameterGeneric() { }

  /// @name accessors to 'modal' constitutive parameters
  //@{
  virtual double GetAxialRigidity() const = 0;

  virtual double GetShearRigidity2() const = 0;

  virtual double GetShearRigidity3() const = 0;


  virtual double GetTorsionalRigidity() const = 0;

  virtual double GetBendingRigidity2() const = 0;

  virtual double GetBendingRigidity3() const = 0;


  virtual double GetTranslationalMassInertia() const = 0;

  virtual double GetPolarMassMomentOfInertia() const = 0;

  virtual double GetMassMomentOfInertia2() const = 0;

  virtual double GetMassMomentOfInertia3() const = 0;
  //@}

  virtual double GetInteractionRadius() const = 0;

  /// create material instance of matching type with my parameters
  virtual Teuchos::RCP<MAT::Material> CreateMaterial();

};

/*-------------------------------------------------------------------------------------------*/
/// constitutive parameters for a Reissner beam formulation (hyperelastic stored energy function)
class BeamReissnerElastHyperMaterialParams : public BeamElastHyperMaterialParameterGeneric
{
public:

  /// standard constructor
  BeamReissnerElastHyperMaterialParams(
    Teuchos::RCP<MAT::PAR::Material> matdata
    );

  /// destructor
  virtual ~BeamReissnerElastHyperMaterialParams() { }


  /// @name derived: accessors to 'modal' constitutive parameters
  //@{
  virtual double GetAxialRigidity() const
  {
    return youngs_modulus_ * cross_section_area_;
  }

  virtual double GetShearRigidity2() const
  {
    return shear_modulus_ * cross_section_area_ * shear_correction_factor_;
  }

  virtual double GetShearRigidity3() const
  {
    return shear_modulus_ * cross_section_area_ * shear_correction_factor_;
  }


  virtual double GetTorsionalRigidity() const
  {
    return shear_modulus_ * area_moment_inertia_polar_;
  }

  virtual double GetBendingRigidity2() const
  {
    return youngs_modulus_ * area_moment_inertia_2_;
  }

  virtual double GetBendingRigidity3() const
  {
    return youngs_modulus_ * area_moment_inertia_3_;
  }


  virtual double GetTranslationalMassInertia() const
  {
    return density_ * cross_section_area_;
  }


  virtual double GetPolarMassMomentOfInertia() const
  {
    return density_ * (area_moment_inertia_2_ + area_moment_inertia_3_);
  }

  virtual double GetMassMomentOfInertia2() const
  {
    return density_ * area_moment_inertia_2_;
  }

  virtual double GetMassMomentOfInertia3() const
  {
    return density_ * area_moment_inertia_3_;
  }


  virtual double GetInteractionRadius() const
  {
    if (radius_interaction_ == -1.0)
      dserror("the radius of a beam which is to be used for interactions (contact, potentials, "
          "viscous drag in background fluid ...) has not been specified in the material definition!");

    return radius_interaction_;
  }
  //@}

private:
  /// @name constitutive parameters
  //@{

  /// Young's modulus
  const double youngs_modulus_;
  /// shear modulus
  const double shear_modulus_;
  /// mass density
  const double density_;

  /// cross-section area
  const double cross_section_area_;
  /// shear correction factor
  const double shear_correction_factor_;

  /// polar/axial area moment of inertia
  const double area_moment_inertia_polar_;
  /// area moment of inertia w.r.t. first principal axis of inertia (i.e. second base vector)
  const double area_moment_inertia_2_;
  /// area moment of inertia w.r.t. second principal axis of inertia (i.e. third base vector)
  const double area_moment_inertia_3_;

  /// radius of a circular cross-section which is ONLY used to evaluate interactions
  //  such as contact, potentials, viscous drag in background fluid ...
  //  it is an optional input parameter and defaults to -1.0 if not required (no interactions at all)
  //  of course, this should be generalized as soon as we allow for other cross-section
  //  shapes in our models for beam-to-X interactions
  const double radius_interaction_;
  //@}

};

/*-------------------------------------------------------------------------------------------*/
/// constitutive parameters for a Reissner beam formulation (hyperelastic stored energy function),
// specified individually 'by mode', i.e. axial tension, torsion, bending (2x)
class BeamReissnerElastHyperMaterialParamsByMode : public BeamElastHyperMaterialParameterGeneric
{
public:

  /// standard constructor
  BeamReissnerElastHyperMaterialParamsByMode(
    Teuchos::RCP<MAT::PAR::Material> matdata
    );

  /// destructor
  virtual ~BeamReissnerElastHyperMaterialParamsByMode() { }


  /// @name derived: accessors to 'modal' constitutive parameters
  //@{
  virtual double GetAxialRigidity() const
  {
    return axial_rigidity_;
  }

  virtual double GetShearRigidity2() const
  {
    return shear_rigidity_2_;
  }

  virtual double GetShearRigidity3() const
  {
    return shear_rigidity_3_;
  }


  virtual double GetTorsionalRigidity() const
  {
    return torsional_rigidity_;
  }

  virtual double GetBendingRigidity2() const
  {
    return bending_rigidity_2_;
  }

  virtual double GetBendingRigidity3() const
  {
    return bending_rigidity_3_;
  }


  virtual double GetTranslationalMassInertia() const
  {
    return translational_mass_inertia_;
  }


  virtual double GetPolarMassMomentOfInertia() const
  {
    return mass_moment_inertia_polar_;
  }

  virtual double GetMassMomentOfInertia2() const
  {
    return mass_moment_inertia_2_;
  }

  virtual double GetMassMomentOfInertia3() const
  {
    return mass_moment_inertia_3_;
  }

  virtual double GetInteractionRadius() const
  {
    if (radius_interaction_ == -1.0)
      dserror("the radius of a beam which is to be used for interactions (contact, potentials, "
          "viscous drag in background fluid ...) has not been specified in the material definition!");

    return radius_interaction_;
  }
  //@}

private:
  /// @name constitutive parameters
  //@{

  /// axial rigidity
  const double axial_rigidity_;
  /// shear rigidity w.r.t first principal axis of inertia (i.e. second base vector)
  const double shear_rigidity_2_;
  /// shear rigidity w.r.t second principal axis of inertia (i.e. third base vector)
  const double shear_rigidity_3_;

  /// torsional rigidity
  const double torsional_rigidity_;
  /// flexural/bending rigidity w.r.t. first principal axis of inertia (i.e. second base vector)
  const double bending_rigidity_2_;
  /// flexural/bending rigidity w.r.t. second principal axis of inertia (i.e. third base vector)
  const double bending_rigidity_3_;

  /// translational mass inertia: mass density * cross-section area
  const double translational_mass_inertia_;

  /// polar mass moment of inertia, i.e. w.r.t. rotation around beam axis
  const double mass_moment_inertia_polar_;
  /// mass moment of inertia w.r.t. first principal axis of inertia (i.e. second base vector)
  const double mass_moment_inertia_2_;
  /// mass moment of inertia w.r.t. second principal axis of inertia (i.e. third base vector)
  const double mass_moment_inertia_3_;

  /// radius of a circular cross-section which is ONLY used to evaluate interactions
  //  such as contact, potentials, viscous drag in background fluid ...
  //  it is an optional input parameter and defaults to -1.0 if not required (no interactions at all)
  //  of course, this should be generalized as soon as we allow for other cross-section
  //  shapes in our models for beam-to-X interactions
  const double radius_interaction_;
  //@}
};



/*-------------------------------------------------------------------------------------------*/
/// constitutive parameters for a Kirchhoff beam formulation (hyperelastic stored energy function)
class BeamKirchhoffElastHyperMaterialParams : public BeamElastHyperMaterialParameterGeneric
{
public:

  /// standard constructor
  BeamKirchhoffElastHyperMaterialParams(
    Teuchos::RCP<MAT::PAR::Material> matdata
    );

  /// destructor
  virtual ~BeamKirchhoffElastHyperMaterialParams() { }


  /// @name derived: accessors to 'modal' constitutive parameters
  //@{
  virtual double GetAxialRigidity() const
  {
    return youngs_modulus_ * cross_section_area_;
  }

  virtual double GetShearRigidity2() const
  {
    return 0.0;
  }

  virtual double GetShearRigidity3() const
  {
    return 0.0;
  }


  virtual double GetTorsionalRigidity() const
  {
    return shear_modulus_ * area_moment_inertia_polar_;
  }

  virtual double GetBendingRigidity2() const
  {
    return youngs_modulus_ * area_moment_inertia_2_;
  }

  virtual double GetBendingRigidity3() const
  {
    return youngs_modulus_ * area_moment_inertia_3_;
  }


  virtual double GetTranslationalMassInertia() const
  {
    return density_ * cross_section_area_;
  }


  virtual double GetPolarMassMomentOfInertia() const
  {
    return density_ * (area_moment_inertia_2_ + area_moment_inertia_3_);
  }

  virtual double GetMassMomentOfInertia2() const
  {
    return density_ * area_moment_inertia_2_;
  }

  virtual double GetMassMomentOfInertia3() const
  {
    return density_ * area_moment_inertia_3_;
  }


  virtual double GetInteractionRadius() const
  {
    if (radius_interaction_ == -1.0)
      dserror("the radius of a beam which is to be used for interactions (contact, potentials, "
          "viscous drag in background fluid ...) has not been specified in the material definition!");

    return radius_interaction_;
  }
  //@}

private:
  /// @name constitutive parameters
  //@{

  /// Young's modulus
  const double youngs_modulus_;
  /// shear modulus
  const double shear_modulus_;
  /// mass density
  const double density_;

  /// cross-section area
  const double cross_section_area_;

  /// polar/axial area moment of inertia
  const double area_moment_inertia_polar_;
  /// area moment of inertia w.r.t. first principal axis of inertia (i.e. second base vector)
  const double area_moment_inertia_2_;
  /// area moment of inertia w.r.t. second principal axis of inertia (i.e. third base vector)
  const double area_moment_inertia_3_;

  /// radius of a circular cross-section which is ONLY used to evaluate interactions
  //  such as contact, potentials, viscous drag in background fluid ...
  //  it is an optional input parameter and defaults to -1.0 if not required (no interactions at all)
  //  of course, this should be generalized as soon as we allow for other cross-section
  //  shapes in our models for beam-to-X interactions
  const double radius_interaction_;
  //@}

};

/*-------------------------------------------------------------------------------------------*/
/// constitutive parameters for a Kirchhoff beam formulation (hyperelastic stored energy function),
// specified individually 'by mode', i.e. axial tension, torsion, bending (2x)
class BeamKirchhoffElastHyperMaterialParamsByMode : public BeamElastHyperMaterialParameterGeneric
{
public:

  /// standard constructor
  BeamKirchhoffElastHyperMaterialParamsByMode(
    Teuchos::RCP<MAT::PAR::Material> matdata
    );

  /// destructor
  virtual ~BeamKirchhoffElastHyperMaterialParamsByMode() { }


  /// @name derived: accessors to 'modal' constitutive parameters
  //@{
  virtual double GetAxialRigidity() const
  {
    return axial_rigidity_;
  }

  virtual double GetShearRigidity2() const
  {
    return 0.0;
  }

  virtual double GetShearRigidity3() const
  {
    return 0.0;
  }


  virtual double GetTorsionalRigidity() const
  {
    return torsional_rigidity_;
  }

  virtual double GetBendingRigidity2() const
  {
    return bending_rigidity_2_;
  }

  virtual double GetBendingRigidity3() const
  {
    return bending_rigidity_3_;
  }


  virtual double GetTranslationalMassInertia() const
  {
    return translational_mass_inertia_;
  }


  virtual double GetPolarMassMomentOfInertia() const
  {
    return mass_moment_inertia_polar_;
  }

  virtual double GetMassMomentOfInertia2() const
  {
    return mass_moment_inertia_2_;
  }

  virtual double GetMassMomentOfInertia3() const
  {
    return mass_moment_inertia_3_;
  }

  virtual double GetInteractionRadius() const
  {
    if (radius_interaction_ == -1.0)
      dserror("the radius of a beam which is to be used for interactions (contact, potentials, "
          "viscous drag in background fluid ...) has not been specified in the material definition!");

    return radius_interaction_;
  }
  //@}

private:
  /// @name constitutive parameters
  //@{

  /// axial rigidity
  const double axial_rigidity_;

  /// torsional rigidity
  const double torsional_rigidity_;
  /// flexural/bending rigidity w.r.t. first principal axis of inertia (i.e. second base vector)
  const double bending_rigidity_2_;
  /// flexural/bending rigidity w.r.t. second principal axis of inertia (i.e. third base vector)
  const double bending_rigidity_3_;

  /// translational mass inertia: mass density * cross-section area
  const double translational_mass_inertia_;

  /// polar mass moment of inertia, i.e. w.r.t. rotation around beam axis
  const double mass_moment_inertia_polar_;
  /// mass moment of inertia w.r.t. first principal axis of inertia (i.e. second base vector)
  const double mass_moment_inertia_2_;
  /// mass moment of inertia w.r.t. second principal axis of inertia (i.e. third base vector)
  const double mass_moment_inertia_3_;

  /// radius of a circular cross-section which is ONLY used to evaluate interactions
  //  such as contact, potentials, viscous drag in background fluid ...
  //  it is an optional input parameter and defaults to -1.0 if not required (no interactions at all)
  //  of course, this should be generalized as soon as we allow for other cross-section
  //  shapes in our models for beam-to-X interactions
  const double radius_interaction_;
  //@}
};



/*-------------------------------------------------------------------------------------------*/
/// constitutive parameters for a torsion-free, isotropic Kirchhoff beam formulation
// (hyperelastic stored energy function)
class BeamKirchhoffTorsionFreeElastHyperMaterialParams : public BeamElastHyperMaterialParameterGeneric
{
public:

  /// standard constructor
  BeamKirchhoffTorsionFreeElastHyperMaterialParams(
    Teuchos::RCP<MAT::PAR::Material> matdata
    );

  /// destructor
  virtual ~BeamKirchhoffTorsionFreeElastHyperMaterialParams() { }


  /// @name derived: accessors to 'modal' constitutive parameters
  //@{
  virtual double GetAxialRigidity() const
  {
    return youngs_modulus_ * cross_section_area_;
  }

  virtual double GetShearRigidity2() const
  {
    return 0.0;
  }

  virtual double GetShearRigidity3() const
  {
    return 0.0;
  }


  virtual double GetTorsionalRigidity() const
  {
    return 0.0;
  }

  virtual double GetBendingRigidity2() const
  {
    return youngs_modulus_ * area_moment_inertia_;
  }

  virtual double GetBendingRigidity3() const
  {
    return youngs_modulus_ * area_moment_inertia_;
  }


  virtual double GetTranslationalMassInertia() const
  {
    return density_ * cross_section_area_;
  }


  virtual double GetPolarMassMomentOfInertia() const
  {
    return 0.0;
  }

  virtual double GetMassMomentOfInertia2() const
  {
    return 0.0;
  }

  virtual double GetMassMomentOfInertia3() const
  {
    return 0.0;
  }


  virtual double GetInteractionRadius() const
  {
    if (radius_interaction_ == -1.0)
      dserror("the radius of a beam which is to be used for interactions (contact, potentials, "
          "viscous drag in background fluid ...) has not been specified in the material definition!");

    return radius_interaction_;
  }
  //@}

private:
  /// @name constitutive parameters
  //@{

  /// Young's modulus
  const double youngs_modulus_;
  /// mass density
  const double density_;

  /// cross-section area
  const double cross_section_area_;
  /// area moment of inertia
  const double area_moment_inertia_;

  /// radius of a circular cross-section which is ONLY used to evaluate interactions
  //  such as contact, potentials, viscous drag in background fluid ...
  //  it is an optional input parameter and defaults to -1.0 if not required (no interactions at all)
  //  of course, this should be generalized as soon as we allow for other cross-section
  //  shapes in our models for beam-to-X interactions
  const double radius_interaction_;
  //@}

};

/*-------------------------------------------------------------------------------------------*/
/// constitutive parameters for a torsion-free, isotropic Kirchhoff beam formulation
// (hyperelastic stored energy function),
// specified individually 'by mode', i.e. axial tension and bending
class BeamKirchhoffTorsionFreeElastHyperMaterialParamsByMode : public BeamElastHyperMaterialParameterGeneric
{
public:

  /// standard constructor
  BeamKirchhoffTorsionFreeElastHyperMaterialParamsByMode(
    Teuchos::RCP<MAT::PAR::Material> matdata
    );

  /// destructor
  virtual ~BeamKirchhoffTorsionFreeElastHyperMaterialParamsByMode() { }


  /// @name derived: accessors to 'modal' constitutive parameters
  //@{
  virtual double GetAxialRigidity() const
  {
    return axial_rigidity_;
  }

  virtual double GetShearRigidity2() const
  {
    return 0.0;
  }

  virtual double GetShearRigidity3() const
  {
    return 0.0;
  }


  virtual double GetTorsionalRigidity() const
  {
    return 0.0;
  }

  virtual double GetBendingRigidity2() const
  {
    return bending_rigidity_;
  }

  virtual double GetBendingRigidity3() const
  {
    return bending_rigidity_;
  }


  virtual double GetTranslationalMassInertia() const
  {
    return translational_mass_inertia_;
  }


  virtual double GetPolarMassMomentOfInertia() const
  {
    return 0.0;
  }

  virtual double GetMassMomentOfInertia2() const
  {
    return 0.0;
  }

  virtual double GetMassMomentOfInertia3() const
  {
    return 0.0;
  }


  virtual double GetInteractionRadius() const
  {
    if (radius_interaction_ == -1.0)
      dserror("the radius of a beam which is to be used for interactions (contact, potentials, "
          "viscous drag in background fluid ...) has not been specified in the material definition!");

    return radius_interaction_;
  }
  //@}

private:
  /// @name constitutive parameters
  //@{

  /// axial rigidity
  const double axial_rigidity_;
  /// flexural/bending rigidity
  const double bending_rigidity_;

  /// translational mass inertia: mass density * cross-section area
  const double translational_mass_inertia_;

  /// radius of a circular cross-section which is ONLY used to evaluate interactions
  //  such as contact, potentials, viscous drag in background fluid ...
  //  it is an optional input parameter and defaults to -1.0 if not required (no interactions at all)
  //  of course, this should be generalized as soon as we allow for other cross-section
  //  shapes in our models for beam-to-X interactions
  const double radius_interaction_;
  //@}
};


}

}

#endif
