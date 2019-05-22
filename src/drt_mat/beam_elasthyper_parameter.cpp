/*-----------------------------------------------------------------------------------------------*/
/*!
\brief constitutive parameters for a beam material based on hyperelastic stored energy function

\maintainer Maximilian Grill

\level 3
*/
/*-----------------------------------------------------------------------------------------------*/

#include "beam_elasthyper_parameter.H"
#include "beam_elasthyper.H"

#include "../drt_mat/matpar_material.H"

#include <Sacado.hpp>
#include <Teuchos_RCP.hpp>


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
double MAT::PAR::DetermineShearModulus(const Teuchos::RCP<MAT::PAR::Material>& matdata)
{
  double shearmodulus = 0.0;
  double poissonratio = 0.0;

  // We want the flexibility to either specify the shear modulus or the Poisson's ratio.
  // Therefore, both parameters are defined as optional in the definition of the input file line

  shearmodulus = matdata->GetDouble("SHEARMOD");
  poissonratio = matdata->GetDouble("POISSONRATIO");

  if (shearmodulus != -1.0 and poissonratio == -1.0)
  {
    // all good, only a value for shear modulus was given directly
  }
  else if (shearmodulus == -1.0 and poissonratio != -1.0)
  {
    // compute shear modulus from Young's modulus and given Poisson's ratio
    shearmodulus = matdata->GetDouble("YOUNG") / (2.0 * (1.0 + poissonratio));
  }
  else if (shearmodulus != -1.0 and poissonratio != -1.0)
  {
    dserror(
        "You specified both of the redundant material parameters SHEARMOD and POISSONRATIO! "
        "Specify exactly one of them in the material definition of your input file!");
  }
  else
  {
    dserror(
        "You specified none of the material parameters SHEARMOD and POISSONRATIO! "
        "Specify exactly one of them in the material definition of your input file!");
  }

  return shearmodulus;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
double MAT::PAR::DetermineDefaultInteractionRadius(const Teuchos::RCP<MAT::PAR::Material>& matdata)
{
  double radius = matdata->GetDouble("INTERACTIONRADIUS");

  double Iyy = matdata->GetDouble("MOMIN2");
  double Izz = matdata->GetDouble("MOMIN3");

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
  double radius = matdata->GetDouble("INTERACTIONRADIUS");

  double Iyy = matdata->GetDouble("MOMIN");

  // determine default value for interaction radius if no value was given:
  // assume circular cross-section and compute from the area moment of inertia
  if (radius == -1.0) radius = std::pow(4.0 * Iyy / M_PI, 0.25);

  return radius;
}



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::PAR::BeamElastHyperMaterialParameterGeneric::BeamElastHyperMaterialParameterGeneric(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata)
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
  return Teuchos::rcp(new MAT::BeamElastHyperMaterial(this));
}



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::PAR::BeamReissnerElastHyperMaterialParams::BeamReissnerElastHyperMaterialParams(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : BeamElastHyperMaterialParameterGeneric(matdata),
      youngs_modulus_(matdata->GetDouble("YOUNG")),
      shear_modulus_(DetermineShearModulus(matdata)),
      density_(matdata->GetDouble("DENS")),
      cross_section_area_(matdata->GetDouble("CROSSAREA")),
      shear_correction_factor_(matdata->GetDouble("SHEARCORR")),
      area_moment_inertia_polar_(matdata->GetDouble("MOMINPOL")),
      area_moment_inertia_2_(matdata->GetDouble("MOMIN2")),
      area_moment_inertia_3_(matdata->GetDouble("MOMIN3")),
      radius_interaction_(DetermineDefaultInteractionRadius(matdata))
{
  if (youngs_modulus_ <= 0.0) dserror("Young's modulus must be positive value");

  if (shear_modulus_ <= 0.0) dserror("shear modulus must be positive value");

  if (density_ < 0.0) dserror("density must not be negative value");

  if (cross_section_area_ <= 0.0) dserror("cross-section area must be positive value");

  if (shear_correction_factor_ <= 0.0) dserror("shear correction factor must be positive value");

  if (area_moment_inertia_polar_ <= 0.0)
    dserror("polar/axial area moment of inertia must be positive value");

  if (area_moment_inertia_2_ <= 0.0) dserror("area moment of inertia must be positive value");

  if (area_moment_inertia_3_ <= 0.0) dserror("area moment of inertia must be positive value");


  /* the radius of an assumed circular cross-section is only used for the evaluation
   * of all kinds of interactions. it can hence be ignored if no interaction are considered. */
  if (radius_interaction_ != -1.0 and radius_interaction_ <= 0.0)
    dserror(
        "if specified (only required if any kind of beam interactions are considered and you "
        "don't want to use the default radius computed from the area moment of inertia), the "
        "given interaction radius must be a positive value");
}



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::PAR::BeamReissnerElastHyperMaterialParamsByMode::BeamReissnerElastHyperMaterialParamsByMode(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : BeamElastHyperMaterialParameterGeneric(matdata),
      axial_rigidity_(matdata->GetDouble("EA")),
      shear_rigidity_2_(matdata->GetDouble("GA2")),
      shear_rigidity_3_(matdata->GetDouble("GA3")),
      torsional_rigidity_(matdata->GetDouble("GI_T")),
      bending_rigidity_2_(matdata->GetDouble("EI2")),
      bending_rigidity_3_(matdata->GetDouble("EI3")),
      translational_mass_inertia_(matdata->GetDouble("RhoA")),
      mass_moment_inertia_polar_(matdata->GetDouble("MASSMOMINPOL")),
      mass_moment_inertia_2_(matdata->GetDouble("MASSMOMIN2")),
      mass_moment_inertia_3_(matdata->GetDouble("MASSMOMIN3")),
      radius_interaction_(matdata->GetDouble("INTERACTIONRADIUS"))
{
  if (axial_rigidity_ <= 0.0) dserror("axial rigidity must be positive value");

  if (shear_rigidity_2_ <= 0.0 or shear_rigidity_3_ <= 0.0)
    dserror("shear rigidity must be positive value");

  if (torsional_rigidity_ <= 0.0) dserror("torsional rigidity must be positive value");

  if (bending_rigidity_2_ <= 0.0 or bending_rigidity_3_ <= 0.0)
    dserror("bending rigidity must be positive value");

  if (translational_mass_inertia_ < 0.0)
    dserror("translational mass inertia must not be negative value");

  if (mass_moment_inertia_polar_ < 0.0)
    dserror("polar mass moment of inertia must not be negative value");

  if (mass_moment_inertia_2_ < 0.0 or mass_moment_inertia_3_ < 0.0)
    dserror("mass moment of inertia must not be negative value");


  /* the radius of an assumed circular cross-section is only used for the evaluation
   * of all kinds of interactions. it can hence be ignored if no interaction are considered. */
  if (radius_interaction_ != -1.0 and radius_interaction_ <= 0.0)
    dserror(
        "if specified (only required if any kind of beam interactions are considered), the "
        "given interaction radius must be a positive value");
}



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::PAR::BeamKirchhoffElastHyperMaterialParams::BeamKirchhoffElastHyperMaterialParams(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : BeamElastHyperMaterialParameterGeneric(matdata),
      youngs_modulus_(matdata->GetDouble("YOUNG")),
      shear_modulus_(DetermineShearModulus(matdata)),
      density_(matdata->GetDouble("DENS")),
      cross_section_area_(matdata->GetDouble("CROSSAREA")),
      area_moment_inertia_polar_(matdata->GetDouble("MOMINPOL")),
      area_moment_inertia_2_(matdata->GetDouble("MOMIN2")),
      area_moment_inertia_3_(matdata->GetDouble("MOMIN3")),
      radius_interaction_(DetermineDefaultInteractionRadius(matdata))
{
  if (youngs_modulus_ <= 0.0) dserror("Young's modulus must be positive value");

  if (shear_modulus_ <= 0.0) dserror("shear modulus must be positive value");

  if (density_ < 0.0) dserror("density must not be negative value");

  if (cross_section_area_ <= 0.0) dserror("cross-section area must be positive value");

  if (area_moment_inertia_polar_ <= 0.0)
    dserror("polar/axial area moment of inertia must be positive value");

  if (area_moment_inertia_2_ <= 0.0) dserror("area moment of inertia must be positive value");

  if (area_moment_inertia_3_ <= 0.0) dserror("area moment of inertia must be positive value");


  /* the radius of an assumed circular cross-section is only used for the evaluation
   * of all kinds of interactions. it can hence be ignored if no interaction are considered. */
  if (radius_interaction_ != -1.0 and radius_interaction_ <= 0.0)
    dserror(
        "if specified (only required if any kind of beam interactions are considered and you "
        "don't want to use the default radius computed from the area moment of inertia), the "
        "given interaction radius must be a positive value");
}



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::PAR::BeamKirchhoffElastHyperMaterialParamsByMode::BeamKirchhoffElastHyperMaterialParamsByMode(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : BeamElastHyperMaterialParameterGeneric(matdata),
      axial_rigidity_(matdata->GetDouble("EA")),
      torsional_rigidity_(matdata->GetDouble("GI_T")),
      bending_rigidity_2_(matdata->GetDouble("EI2")),
      bending_rigidity_3_(matdata->GetDouble("EI3")),
      translational_mass_inertia_(matdata->GetDouble("RhoA")),
      mass_moment_inertia_polar_(matdata->GetDouble("MASSMOMINPOL")),
      mass_moment_inertia_2_(matdata->GetDouble("MASSMOMIN2")),
      mass_moment_inertia_3_(matdata->GetDouble("MASSMOMIN3")),
      radius_interaction_(matdata->GetDouble("INTERACTIONRADIUS"))
{
  if (axial_rigidity_ <= 0.0) dserror("axial rigidity must be positive value");

  if (torsional_rigidity_ <= 0.0) dserror("torsional rigidity must be positive value");

  if (bending_rigidity_2_ <= 0.0 or bending_rigidity_3_ <= 0.0)
    dserror("bending rigidity must be positive value");

  if (translational_mass_inertia_ < 0.0)
    dserror("translational mass inertia must not be negative value");

  if (mass_moment_inertia_polar_ < 0.0)
    dserror("polar mass moment of inertia must not be negative value");

  if (mass_moment_inertia_2_ < 0.0 or mass_moment_inertia_3_ < 0.0)
    dserror("mass moment of inertia must not be negative value");


  /* the radius of an assumed circular cross-section is only used for the evaluation
   * of all kinds of interactions. it can hence be ignored if no interaction are considered. */
  if (radius_interaction_ != -1.0 and radius_interaction_ <= 0.0)
    dserror(
        "if specified (only required if any kind of beam interactions are considered), the "
        "given interaction radius must be a positive value");
}



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParams::
    BeamKirchhoffTorsionFreeElastHyperMaterialParams(Teuchos::RCP<MAT::PAR::Material> matdata)
    : BeamElastHyperMaterialParameterGeneric(matdata),
      youngs_modulus_(matdata->GetDouble("YOUNG")),
      density_(matdata->GetDouble("DENS")),
      cross_section_area_(matdata->GetDouble("CROSSAREA")),
      area_moment_inertia_(matdata->GetDouble("MOMIN")),
      radius_interaction_(DetermineDefaultInteractionRadiusIsotropic(matdata))
{
  if (youngs_modulus_ <= 0.0) dserror("Young's modulus must be positive value");

  if (density_ < 0.0) dserror("density must not be negative value");

  if (cross_section_area_ <= 0.0) dserror("cross-section area must be positive value");

  if (area_moment_inertia_ <= 0.0) dserror("area moment of inertia must be positive value");


  /* the radius of an assumed circular cross-section is only used for the evaluation
   * of all kinds of interactions. it can hence be ignored if no interaction are considered. */
  if (radius_interaction_ != -1.0 and radius_interaction_ <= 0.0)
    dserror(
        "if specified (only required if any kind of beam interactions are considered and you "
        "don't want to use the default radius computed from the area moment of inertia), the "
        "given interaction radius must be a positive value");
}



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
MAT::PAR::BeamKirchhoffTorsionFreeElastHyperMaterialParamsByMode::
    BeamKirchhoffTorsionFreeElastHyperMaterialParamsByMode(Teuchos::RCP<MAT::PAR::Material> matdata)
    : BeamElastHyperMaterialParameterGeneric(matdata),
      axial_rigidity_(matdata->GetDouble("EA")),
      bending_rigidity_(matdata->GetDouble("EI")),
      translational_mass_inertia_(matdata->GetDouble("RhoA")),
      radius_interaction_(matdata->GetDouble("INTERACTIONRADIUS"))
{
  if (axial_rigidity_ <= 0.0) dserror("axial rigidity must be positive value");

  if (bending_rigidity_ <= 0.0) dserror("bending rigidity must be positive value");

  if (translational_mass_inertia_ < 0.0)
    dserror("translational mass inertia must not be negative value");


  /* the radius of an assumed circular cross-section is only used for the evaluation
   * of all kinds of interactions. it can hence be ignored if no interaction are considered. */
  if (radius_interaction_ != -1.0 and radius_interaction_ <= 0.0)
    dserror(
        "if specified (only required if any kind of beam interactions are considered), the "
        "given interaction radius must be a positive value");
}
