/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief constitutive parameters for a beam material based on hyperelastic stored energy function


\level 3
*/
/*-----------------------------------------------------------------------------------------------*/

#include "4C_mat_beam_elasthyper_parameter.hpp"

#include "4C_mat_beam_elasthyper.hpp"
#include "4C_mat_par_material.hpp"

#include <Sacado.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
double MAT::PAR::DetermineShearModulus(const Teuchos::RCP<MAT::PAR::Material>& matdata)
{
  double shearmodulus = 0.0;
  double poissonratio = 0.0;

  // We want the flexibility to either specify the shear modulus or the Poisson's ratio.
  // Therefore, both parameters are defined as optional in the definition of the input file line

  shearmodulus = *matdata->Get<double>("SHEARMOD");
  poissonratio = *matdata->Get<double>("POISSONRATIO");

  if (shearmodulus != -1.0 and poissonratio == -1.0)
  {
    // all good, only a value for shear modulus was given directly
  }
  else if (shearmodulus == -1.0 and poissonratio != -1.0)
  {
    // compute shear modulus from Young's modulus and given Poisson's ratio
    shearmodulus = *matdata->Get<double>("YOUNG") / (2.0 * (1.0 + poissonratio));
  }
  else if (shearmodulus != -1.0 and poissonratio != -1.0)
  {
    FOUR_C_THROW(
        "You specified both of the redundant material parameters SHEARMOD and POISSONRATIO! "
        "Specify exactly one of them in the material definition of your input file!");
  }
  else
  {
    FOUR_C_THROW(
        "You specified none of the material parameters SHEARMOD and POISSONRATIO! "
        "Specify exactly one of them in the material definition of your input file!");
  }

  return shearmodulus;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
double MAT::PAR::DetermineDefaultInteractionRadius(const Teuchos::RCP<MAT::PAR::Material>& matdata)
{
  double radius = *matdata->Get<double>("INTERACTIONRADIUS");

  double Iyy = *matdata->Get<double>("MOMIN2");
  double Izz = *matdata->Get<double>("MOMIN3");

  // determine default value for interaction radius if no value was given:
  // assume circular cross-section and compute from the area moment of inertia
  if (radius == -1.0 and Iyy == Izz) radius = std::pow(4.0 * Iyy / M_PI, 0.25);

  return radius;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
double MAT::PAR::DetermineDefaultInteractionRadiusIsotropic(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
{
  double radius = *matdata->Get<double>("INTERACTIONRADIUS");

  double Iyy = *matdata->Get<double>("MOMIN");

  // determine default value for interaction radius if no value was given:
  // assume circular cross-section and compute from the area moment of inertia
  if (radius == -1.0) radius = std::pow(4.0 * Iyy / M_PI, 0.25);

  return radius;
}



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::PAR::BeamElastHyperMaterialParameterGeneric::BeamElastHyperMaterialParameterGeneric(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata), use_fad_(*matdata->Get<bool>("FAD"))
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::BeamElastHyperMaterialParameterGeneric::CreateMaterial()
{
  /* all the different parameter sets (Reissner/Kirchhoff/..., 'classic'/'by modes') are used to
   * parameterize the same constitutive relations based on a hyperelastic stored energy function
   * formulated for cross-section resultants which are implemented in BeamElastHyperMaterial */
  Teuchos::RCP<MAT::Material> matobject;

  if (Uses_FAD())
  {
    matobject = Teuchos::rcp(new MAT::BeamElastHyperMaterial<Sacado::Fad::DFad<double>>(this));
  }
  else
    matobject = Teuchos::rcp(new MAT::BeamElastHyperMaterial<double>(this));
  return matobject;
}



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::PAR::BeamReissnerElastHyperMaterialParams::BeamReissnerElastHyperMaterialParams(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : BeamElastHyperMaterialParameterGeneric(matdata),
      youngs_modulus_(*matdata->Get<double>("YOUNG")),
      shear_modulus_(DetermineShearModulus(matdata)),
      density_(*matdata->Get<double>("DENS")),
      cross_section_area_(*matdata->Get<double>("CROSSAREA")),
      shear_correction_factor_(*matdata->Get<double>("SHEARCORR")),
      area_moment_inertia_polar_(*matdata->Get<double>("MOMINPOL")),
      area_moment_inertia_2_(*matdata->Get<double>("MOMIN2")),
      area_moment_inertia_3_(*matdata->Get<double>("MOMIN3")),
      radius_interaction_(DetermineDefaultInteractionRadius(matdata))
{
  if (youngs_modulus_ <= 0.0) FOUR_C_THROW("Young's modulus must be positive value");

  if (shear_modulus_ <= 0.0) FOUR_C_THROW("shear modulus must be positive value");

  if (density_ < 0.0) FOUR_C_THROW("density must not be negative value");

  if (cross_section_area_ <= 0.0) FOUR_C_THROW("cross-section area must be positive value");

  if (shear_correction_factor_ <= 0.0)
    FOUR_C_THROW("shear correction factor must be positive value");

  if (area_moment_inertia_polar_ <= 0.0)
    FOUR_C_THROW("polar/axial area moment of inertia must be positive value");

  if (area_moment_inertia_2_ <= 0.0) FOUR_C_THROW("area moment of inertia must be positive value");

  if (area_moment_inertia_3_ <= 0.0) FOUR_C_THROW("area moment of inertia must be positive value");


  /* the radius of an assumed circular cross-section is only used for the evaluation
   * of all kinds of interactions. it can hence be ignored if no interaction are considered. */
  if (radius_interaction_ != -1.0 and radius_interaction_ <= 0.0)
    FOUR_C_THROW(
        "if specified (only required if any kind of beam interactions are considered and you "
        "don't want to use the default radius computed from the area moment of inertia), the "
        "given interaction radius must be a positive value");
}



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::PAR::BeamReissnerElastHyperMaterialParamsByMode::BeamReissnerElastHyperMaterialParamsByMode(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : BeamElastHyperMaterialParameterGeneric(matdata),
      axial_rigidity_(*matdata->Get<double>("EA")),
      shear_rigidity_2_(*matdata->Get<double>("GA2")),
      shear_rigidity_3_(*matdata->Get<double>("GA3")),
      torsional_rigidity_(*matdata->Get<double>("GI_T")),
      bending_rigidity_2_(*matdata->Get<double>("EI2")),
      bending_rigidity_3_(*matdata->Get<double>("EI3")),
      translational_mass_inertia_(*matdata->Get<double>("RhoA")),
      mass_moment_inertia_polar_(*matdata->Get<double>("MASSMOMINPOL")),
      mass_moment_inertia_2_(*matdata->Get<double>("MASSMOMIN2")),
      mass_moment_inertia_3_(*matdata->Get<double>("MASSMOMIN3")),
      radius_interaction_(*matdata->Get<double>("INTERACTIONRADIUS"))
{
  if (axial_rigidity_ <= 0.0) FOUR_C_THROW("axial rigidity must be positive value");

  if (shear_rigidity_2_ <= 0.0 or shear_rigidity_3_ <= 0.0)
    FOUR_C_THROW("shear rigidity must be positive value");

  if (torsional_rigidity_ <= 0.0) FOUR_C_THROW("torsional rigidity must be positive value");

  if (bending_rigidity_2_ <= 0.0 or bending_rigidity_3_ <= 0.0)
    FOUR_C_THROW("bending rigidity must be positive value");

  if (translational_mass_inertia_ < 0.0)
    FOUR_C_THROW("translational mass inertia must not be negative value");

  if (mass_moment_inertia_polar_ < 0.0)
    FOUR_C_THROW("polar mass moment of inertia must not be negative value");

  if (mass_moment_inertia_2_ < 0.0 or mass_moment_inertia_3_ < 0.0)
    FOUR_C_THROW("mass moment of inertia must not be negative value");


  /* the radius of an assumed circular cross-section is only used for the evaluation
   * of all kinds of interactions. it can hence be ignored if no interaction are considered. */
  if (radius_interaction_ != -1.0 and radius_interaction_ <= 0.0)
    FOUR_C_THROW(
        "if specified (only required if any kind of beam interactions are considered), the "
        "given interaction radius must be a positive value");
}



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::PAR::BeamKirchhoffElastHyperMaterialParams::BeamKirchhoffElastHyperMaterialParams(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : BeamElastHyperMaterialParameterGeneric(matdata),
      youngs_modulus_(*matdata->Get<double>("YOUNG")),
      shear_modulus_(DetermineShearModulus(matdata)),
      density_(*matdata->Get<double>("DENS")),
      cross_section_area_(*matdata->Get<double>("CROSSAREA")),
      area_moment_inertia_polar_(*matdata->Get<double>("MOMINPOL")),
      area_moment_inertia_2_(*matdata->Get<double>("MOMIN2")),
      area_moment_inertia_3_(*matdata->Get<double>("MOMIN3")),
      radius_interaction_(DetermineDefaultInteractionRadius(matdata))
{
  if (youngs_modulus_ <= 0.0) FOUR_C_THROW("Young's modulus must be positive value");

  if (shear_modulus_ <= 0.0) FOUR_C_THROW("shear modulus must be positive value");

  if (density_ < 0.0) FOUR_C_THROW("density must not be negative value");

  if (cross_section_area_ <= 0.0) FOUR_C_THROW("cross-section area must be positive value");

  if (area_moment_inertia_polar_ <= 0.0)
    FOUR_C_THROW("polar/axial area moment of inertia must be positive value");

  if (area_moment_inertia_2_ <= 0.0) FOUR_C_THROW("area moment of inertia must be positive value");

  if (area_moment_inertia_3_ <= 0.0) FOUR_C_THROW("area moment of inertia must be positive value");


  /* the radius of an assumed circular cross-section is only used for the evaluation
   * of all kinds of interactions. it can hence be ignored if no interaction are considered. */
  if (radius_interaction_ != -1.0 and radius_interaction_ <= 0.0)
    FOUR_C_THROW(
        "if specified (only required if any kind of beam interactions are considered and you "
        "don't want to use the default radius computed from the area moment of inertia), the "
        "given interaction radius must be a positive value");
}



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::PAR::BeamKirchhoffElastHyperMaterialParamsByMode::BeamKirchhoffElastHyperMaterialParamsByMode(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : BeamElastHyperMaterialParameterGeneric(matdata),
      axial_rigidity_(*matdata->Get<double>("EA")),
      torsional_rigidity_(*matdata->Get<double>("GI_T")),
      bending_rigidity_2_(*matdata->Get<double>("EI2")),
      bending_rigidity_3_(*matdata->Get<double>("EI3")),
      translational_mass_inertia_(*matdata->Get<double>("RhoA")),
      mass_moment_inertia_polar_(*matdata->Get<double>("MASSMOMINPOL")),
      mass_moment_inertia_2_(*matdata->Get<double>("MASSMOMIN2")),
      mass_moment_inertia_3_(*matdata->Get<double>("MASSMOMIN3")),
      radius_interaction_(*matdata->Get<double>("INTERACTIONRADIUS"))
{
  if (axial_rigidity_ <= 0.0) FOUR_C_THROW("axial rigidity must be positive value");

  if (torsional_rigidity_ <= 0.0) FOUR_C_THROW("torsional rigidity must be positive value");

  if (bending_rigidity_2_ <= 0.0 or bending_rigidity_3_ <= 0.0)
    FOUR_C_THROW("bending rigidity must be positive value");

  if (translational_mass_inertia_ < 0.0)
    FOUR_C_THROW("translational mass inertia must not be negative value");

  if (mass_moment_inertia_polar_ < 0.0)
    FOUR_C_THROW("polar mass moment of inertia must not be negative value");

  if (mass_moment_inertia_2_ < 0.0 or mass_moment_inertia_3_ < 0.0)
    FOUR_C_THROW("mass moment of inertia must not be negative value");


  /* the radius of an assumed circular cross-section is only used for the evaluation
   * of all kinds of interactions. it can hence be ignored if no interaction are considered. */
  if (radius_interaction_ != -1.0 and radius_interaction_ <= 0.0)
    FOUR_C_THROW(
        "if specified (only required if any kind of beam interactions are considered), the "
        "given interaction radius must be a positive value");
}



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParams::
    BeamKirchhoffTorsionFreeElastHyperMaterialParams(Teuchos::RCP<MAT::PAR::Material> matdata)
    : BeamElastHyperMaterialParameterGeneric(matdata),
      youngs_modulus_(*matdata->Get<double>("YOUNG")),
      density_(*matdata->Get<double>("DENS")),
      cross_section_area_(*matdata->Get<double>("CROSSAREA")),
      area_moment_inertia_(*matdata->Get<double>("MOMIN")),
      radius_interaction_(DetermineDefaultInteractionRadiusIsotropic(matdata))
{
  if (youngs_modulus_ <= 0.0) FOUR_C_THROW("Young's modulus must be positive value");

  if (density_ < 0.0) FOUR_C_THROW("density must not be negative value");

  if (cross_section_area_ <= 0.0) FOUR_C_THROW("cross-section area must be positive value");

  if (area_moment_inertia_ <= 0.0) FOUR_C_THROW("area moment of inertia must be positive value");


  /* the radius of an assumed circular cross-section is only used for the evaluation
   * of all kinds of interactions. it can hence be ignored if no interaction are considered. */
  if (radius_interaction_ != -1.0 and radius_interaction_ <= 0.0)
    FOUR_C_THROW(
        "if specified (only required if any kind of beam interactions are considered and you "
        "don't want to use the default radius computed from the area moment of inertia), the "
        "given interaction radius must be a positive value");
}



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParamsByMode::
    BeamKirchhoffTorsionFreeElastHyperMaterialParamsByMode(Teuchos::RCP<MAT::PAR::Material> matdata)
    : BeamElastHyperMaterialParameterGeneric(matdata),
      axial_rigidity_(*matdata->Get<double>("EA")),
      bending_rigidity_(*matdata->Get<double>("EI")),
      translational_mass_inertia_(*matdata->Get<double>("RhoA")),
      radius_interaction_(*matdata->Get<double>("INTERACTIONRADIUS"))
{
  if (axial_rigidity_ <= 0.0) FOUR_C_THROW("axial rigidity must be positive value");

  if (bending_rigidity_ <= 0.0) FOUR_C_THROW("bending rigidity must be positive value");

  if (translational_mass_inertia_ < 0.0)
    FOUR_C_THROW("translational mass inertia must not be negative value");


  /* the radius of an assumed circular cross-section is only used for the evaluation
   * of all kinds of interactions. it can hence be ignored if no interaction are considered. */
  if (radius_interaction_ != -1.0 and radius_interaction_ <= 0.0)
    FOUR_C_THROW(
        "if specified (only required if any kind of beam interactions are considered), the "
        "given interaction radius must be a positive value");
}

FOUR_C_NAMESPACE_CLOSE
