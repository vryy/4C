/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief constitutive parameters for a beam material based on hyperelastic stored energy function


\level 3
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_BEAM_ELASTHYPER_PARAMETER_HPP
#define FOUR_C_MAT_BEAM_ELASTHYPER_PARAMETER_HPP

#include "baci_config.hpp"

#include "baci_mat_par_parameter.hpp"
#include "baci_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


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

    /** \brief determine default value for interaction radius from area moment of inertia just in
     * case that no value was explicitly specified
     *
     *  \author grill
     *  \date 02/17 */
    double DetermineDefaultInteractionRadius(const Teuchos::RCP<MAT::PAR::Material>& matdata);

    /** \brief determine default value for interaction radius from area moment of inertia just in
     * case that no value was explicitly specified: isotropic case, i.e. only one moment of inertia
     *
     *  \author grill
     *  \date 02/17 */
    double DetermineDefaultInteractionRadiusIsotropic(
        const Teuchos::RCP<MAT::PAR::Material>& matdata);

    /*-------------------------------------------------------------------------------------------*/
    /// (generic) constitutive parameters for a beam material based on hyperelastic stored energy
    /// function
    class BeamElastHyperMaterialParameterGeneric : public Parameter
    {
     public:
      /// standard constructor
      BeamElastHyperMaterialParameterGeneric(Teuchos::RCP<MAT::PAR::Material> matdata);

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


      virtual double GetYieldStressN() const { return -1.0; };

      virtual double GetYieldStressM() const { return -1.0; };


      virtual double GetHardeningAxialRigidity() const { return -1.0; };

      virtual double GetHardeningShearRigidity2() const { return -1.0; };

      virtual double GetHardeningShearRigidity3() const { return -1.0; };


      virtual double GetHardeningMomentalRigidity() const { return -1.0; };

      virtual bool GetTorsionPlasticity() const { return false; };

      bool Uses_FAD() { return use_fad_; };

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

     private:
      /// Flag to store if the material is used with a beam material and automatic differentiation
      bool use_fad_;
    };

    /*-------------------------------------------------------------------------------------------*/
    /// constitutive parameters for a Reissner beam formulation (hyperelastic stored energy
    /// function)
    class BeamReissnerElastHyperMaterialParams : public BeamElastHyperMaterialParameterGeneric
    {
     public:
      /// standard constructor
      BeamReissnerElastHyperMaterialParams(Teuchos::RCP<MAT::PAR::Material> matdata);


      /// @name derived: accessors to 'modal' constitutive parameters
      //@{
      double GetAxialRigidity() const override { return youngs_modulus_ * cross_section_area_; }

      double GetShearRigidity2() const override
      {
        return shear_modulus_ * cross_section_area_ * shear_correction_factor_;
      }

      double GetShearRigidity3() const override
      {
        return shear_modulus_ * cross_section_area_ * shear_correction_factor_;
      }


      double GetTorsionalRigidity() const override
      {
        return shear_modulus_ * area_moment_inertia_polar_;
      }

      double GetBendingRigidity2() const override
      {
        return youngs_modulus_ * area_moment_inertia_2_;
      }

      double GetBendingRigidity3() const override
      {
        return youngs_modulus_ * area_moment_inertia_3_;
      }


      double GetTranslationalMassInertia() const override { return density_ * cross_section_area_; }


      double GetPolarMassMomentOfInertia() const override
      {
        return density_ * (area_moment_inertia_2_ + area_moment_inertia_3_);
      }

      double GetMassMomentOfInertia2() const override { return density_ * area_moment_inertia_2_; }

      double GetMassMomentOfInertia3() const override { return density_ * area_moment_inertia_3_; }


      double GetInteractionRadius() const override
      {
        if (radius_interaction_ == -1.0)
          FOUR_C_THROW(
              "the radius of a beam which is to be used for interactions (contact, potentials, "
              "viscous drag in background fluid ...) has not been specified in the material "
              "definition!");

        return radius_interaction_;
      }
      //@}
     protected:
      double GetCrossSectionArea() const { return cross_section_area_; }
      double GetShearCorrectionFactor() const { return shear_correction_factor_; }
      double GetShearModulus() const { return shear_modulus_; }
      double GetMomentInertia2() const { return area_moment_inertia_2_; }
      double GetMomentInertia3() const { return area_moment_inertia_3_; }
      double GetMomentInertiaPolar() const { return area_moment_inertia_polar_; }
      double GetYoungsModulus() const { return youngs_modulus_; }

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
      //  it is an optional input parameter and defaults to -1.0 if not required (no interactions at
      //  all) of course, this should be generalized as soon as we allow for other cross-section
      //  shapes in our models for beam-to-X interactions
      const double radius_interaction_;
      //@}
    };

    /*-------------------------------------------------------------------------------------------*/
    /// constitutive parameters for a Reissner beam formulation (hyperelastic stored energy
    /// function),
    // specified individually 'by mode', i.e. axial tension, torsion, bending (2x)
    class BeamReissnerElastHyperMaterialParamsByMode : public BeamElastHyperMaterialParameterGeneric
    {
     public:
      /// standard constructor
      BeamReissnerElastHyperMaterialParamsByMode(Teuchos::RCP<MAT::PAR::Material> matdata);


      /// @name derived: accessors to 'modal' constitutive parameters
      //@{
      double GetAxialRigidity() const override { return axial_rigidity_; }

      double GetShearRigidity2() const override { return shear_rigidity_2_; }

      double GetShearRigidity3() const override { return shear_rigidity_3_; }


      double GetTorsionalRigidity() const override { return torsional_rigidity_; }

      double GetBendingRigidity2() const override { return bending_rigidity_2_; }

      double GetBendingRigidity3() const override { return bending_rigidity_3_; }


      double GetTranslationalMassInertia() const override { return translational_mass_inertia_; }


      double GetPolarMassMomentOfInertia() const override { return mass_moment_inertia_polar_; }

      double GetMassMomentOfInertia2() const override { return mass_moment_inertia_2_; }

      double GetMassMomentOfInertia3() const override { return mass_moment_inertia_3_; }

      double GetInteractionRadius() const override
      {
        if (radius_interaction_ == -1.0)
          FOUR_C_THROW(
              "the radius of a beam which is to be used for interactions (contact, potentials, "
              "viscous drag in background fluid ...) has not been specified in the material "
              "definition!");

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
      //  it is an optional input parameter and defaults to -1.0 if not required (no interactions at
      //  all) of course, this should be generalized as soon as we allow for other cross-section
      //  shapes in our models for beam-to-X interactions
      const double radius_interaction_;
      //@}
    };



    /*-------------------------------------------------------------------------------------------*/
    /// constitutive parameters for a Kirchhoff beam formulation (hyperelastic stored energy
    /// function)
    class BeamKirchhoffElastHyperMaterialParams : public BeamElastHyperMaterialParameterGeneric
    {
     public:
      /// standard constructor
      BeamKirchhoffElastHyperMaterialParams(Teuchos::RCP<MAT::PAR::Material> matdata);


      /// @name derived: accessors to 'modal' constitutive parameters
      //@{
      double GetAxialRigidity() const override { return youngs_modulus_ * cross_section_area_; }

      double GetShearRigidity2() const override { return 0.0; }

      double GetShearRigidity3() const override { return 0.0; }


      double GetTorsionalRigidity() const override
      {
        return shear_modulus_ * area_moment_inertia_polar_;
      }

      double GetBendingRigidity2() const override
      {
        return youngs_modulus_ * area_moment_inertia_2_;
      }

      double GetBendingRigidity3() const override
      {
        return youngs_modulus_ * area_moment_inertia_3_;
      }


      double GetTranslationalMassInertia() const override { return density_ * cross_section_area_; }


      double GetPolarMassMomentOfInertia() const override
      {
        return density_ * (area_moment_inertia_2_ + area_moment_inertia_3_);
      }

      double GetMassMomentOfInertia2() const override { return density_ * area_moment_inertia_2_; }

      double GetMassMomentOfInertia3() const override { return density_ * area_moment_inertia_3_; }


      double GetInteractionRadius() const override
      {
        if (radius_interaction_ == -1.0)
          FOUR_C_THROW(
              "the radius of a beam which is to be used for interactions (contact, potentials, "
              "viscous drag in background fluid ...) has not been specified in the material "
              "definition!");

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
      //  it is an optional input parameter and defaults to -1.0 if not required (no interactions at
      //  all) of course, this should be generalized as soon as we allow for other cross-section
      //  shapes in our models for beam-to-X interactions
      const double radius_interaction_;
      //@}
    };

    /*-------------------------------------------------------------------------------------------*/
    /// constitutive parameters for a Kirchhoff beam formulation (hyperelastic stored energy
    /// function),
    // specified individually 'by mode', i.e. axial tension, torsion, bending (2x)
    class BeamKirchhoffElastHyperMaterialParamsByMode
        : public BeamElastHyperMaterialParameterGeneric
    {
     public:
      /// standard constructor
      BeamKirchhoffElastHyperMaterialParamsByMode(Teuchos::RCP<MAT::PAR::Material> matdata);


      /// @name derived: accessors to 'modal' constitutive parameters
      //@{
      double GetAxialRigidity() const override { return axial_rigidity_; }

      double GetShearRigidity2() const override { return 0.0; }

      double GetShearRigidity3() const override { return 0.0; }


      double GetTorsionalRigidity() const override { return torsional_rigidity_; }

      double GetBendingRigidity2() const override { return bending_rigidity_2_; }

      double GetBendingRigidity3() const override { return bending_rigidity_3_; }


      double GetTranslationalMassInertia() const override { return translational_mass_inertia_; }


      double GetPolarMassMomentOfInertia() const override { return mass_moment_inertia_polar_; }

      double GetMassMomentOfInertia2() const override { return mass_moment_inertia_2_; }

      double GetMassMomentOfInertia3() const override { return mass_moment_inertia_3_; }

      double GetInteractionRadius() const override
      {
        if (radius_interaction_ == -1.0)
          FOUR_C_THROW(
              "the radius of a beam which is to be used for interactions (contact, potentials, "
              "viscous drag in background fluid ...) has not been specified in the material "
              "definition!");

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
      //  it is an optional input parameter and defaults to -1.0 if not required (no interactions at
      //  all) of course, this should be generalized as soon as we allow for other cross-section
      //  shapes in our models for beam-to-X interactions
      const double radius_interaction_;
      //@}
    };



    /*-------------------------------------------------------------------------------------------*/
    /// constitutive parameters for a torsion-free, isotropic Kirchhoff beam formulation
    // (hyperelastic stored energy function)
    class BeamKirchhoffTorsionFreeElastHyperMaterialParams
        : public BeamElastHyperMaterialParameterGeneric
    {
     public:
      /// standard constructor
      BeamKirchhoffTorsionFreeElastHyperMaterialParams(Teuchos::RCP<MAT::PAR::Material> matdata);


      /// @name derived: accessors to 'modal' constitutive parameters
      //@{
      double GetAxialRigidity() const override { return youngs_modulus_ * cross_section_area_; }

      double GetShearRigidity2() const override { return 0.0; }

      double GetShearRigidity3() const override { return 0.0; }


      double GetTorsionalRigidity() const override { return 0.0; }

      double GetBendingRigidity2() const override { return youngs_modulus_ * area_moment_inertia_; }

      double GetBendingRigidity3() const override { return youngs_modulus_ * area_moment_inertia_; }


      double GetTranslationalMassInertia() const override { return density_ * cross_section_area_; }


      double GetPolarMassMomentOfInertia() const override { return 0.0; }

      double GetMassMomentOfInertia2() const override { return 0.0; }

      double GetMassMomentOfInertia3() const override { return 0.0; }


      double GetInteractionRadius() const override
      {
        if (radius_interaction_ == -1.0)
          FOUR_C_THROW(
              "the radius of a beam which is to be used for interactions (contact, potentials, "
              "viscous drag in background fluid ...) has not been specified in the material "
              "definition!");

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
      //  it is an optional input parameter and defaults to -1.0 if not required (no interactions at
      //  all) of course, this should be generalized as soon as we allow for other cross-section
      //  shapes in our models for beam-to-X interactions
      const double radius_interaction_;
      //@}
    };

    /*-------------------------------------------------------------------------------------------*/
    /// constitutive parameters for a torsion-free, isotropic Kirchhoff beam formulation
    // (hyperelastic stored energy function),
    // specified individually 'by mode', i.e. axial tension and bending
    class BeamKirchhoffTorsionFreeElastHyperMaterialParamsByMode
        : public BeamElastHyperMaterialParameterGeneric
    {
     public:
      /// standard constructor
      BeamKirchhoffTorsionFreeElastHyperMaterialParamsByMode(
          Teuchos::RCP<MAT::PAR::Material> matdata);


      /// @name derived: accessors to 'modal' constitutive parameters
      //@{
      double GetAxialRigidity() const override { return axial_rigidity_; }

      double GetShearRigidity2() const override { return 0.0; }

      double GetShearRigidity3() const override { return 0.0; }


      double GetTorsionalRigidity() const override { return 0.0; }

      double GetBendingRigidity2() const override { return bending_rigidity_; }

      double GetBendingRigidity3() const override { return bending_rigidity_; }


      double GetTranslationalMassInertia() const override { return translational_mass_inertia_; }


      double GetPolarMassMomentOfInertia() const override { return 0.0; }

      double GetMassMomentOfInertia2() const override { return 0.0; }

      double GetMassMomentOfInertia3() const override { return 0.0; }


      double GetInteractionRadius() const override
      {
        if (radius_interaction_ == -1.0)
          FOUR_C_THROW(
              "the radius of a beam which is to be used for interactions (contact, potentials, "
              "viscous drag in background fluid ...) has not been specified in the material "
              "definition!");

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
      //  it is an optional input parameter and defaults to -1.0 if not required (no interactions at
      //  all) of course, this should be generalized as soon as we allow for other cross-section
      //  shapes in our models for beam-to-X interactions
      const double radius_interaction_;
      //@}
    };


  }  // namespace PAR

}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
