/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for the viscous contribution of a generalized maxwell model

\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MATELAST_VISCO_GENERALIZEDGENMAX_HPP
#define FOUR_C_MATELAST_VISCO_GENERALIZEDGENMAX_HPP

#include "4C_config.hpp"

#include "4C_matelast_summand.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*!
       * @brief material parameters for viscous contribution to a viscoelastic branch of a
       * generalized Maxwell model
       *
       * <h3>Input line</h3>
       *  MAT 1 VISCO_GeneralizedGenMax NUMBRANCH 3 MATIDS 4 5 6 SOLVE CONVOL
       *  MAT 1 VISCO_GeneralizedGenMax NUMBRANCH 3 MATIDS 4 5 6 SOLVE OST
       */
      class GeneralizedGenMax : public CORE::MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        GeneralizedGenMax(const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// material parameters
        int numbranch_;
        const std::vector<int> matids_;
        std::string solve_;
        //@}

        /// create material instance of matching type with my parameters

        Teuchos::RCP<CORE::MAT::Material> create_material() override { return Teuchos::null; };
      };  // class GeneralizedGenMax


      /*!
       * @brief material parameters for viscous contribution to a viscoelastic branch of a
       * generalized Maxwell modell
       *
       * <h3>Input line</h3>
       * MAT 1 VISCO_BRANCH NUMMAT 2 MATIDS 2 3
       */
      class ViscoBranch : public CORE::MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        ViscoBranch(const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// material parameters

        double nummat_;
        const std::vector<int> matids_;


        //@}

        /// Override this method and throw error, as the material should be created in within the
        /// Factory method of the elastic summand
        Teuchos::RCP<CORE::MAT::Material> create_material() override
        {
          FOUR_C_THROW(
              "Cannot create a material from this method, as it should be created in "
              "MAT::ELASTIC::Summand::Factory.");
          return Teuchos::null;
        };
      };  // class ViscoBranch

      /*!
       * @brief material parameters for viscous contribution to a viscoelastic branch of a
       * generalized Maxwell modell
       *
       * <h3>Input line</h3>
       * MAT 1 VISCO_PART TAU 1.5
       */
      class ViscoPart : public CORE::MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        ViscoPart(const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata);

        /// @name material parameters
        //@{

        /// material parameters
        double tau_;

        //@}

        /// create material instance of matching type with my parameters
        Teuchos::RCP<CORE::MAT::Material> create_material() override { return Teuchos::null; };
      };  // class ViscoPart

    }  // namespace PAR



    class GeneralizedGenMax : public Summand
    {
     public:
      /// constructor with given material parameters
      GeneralizedGenMax(MAT::ELASTIC::PAR::GeneralizedGenMax* params);

      /// @name Access material constants
      //@{

      /// material type
      CORE::Materials::MaterialType MaterialType() const override
      {
        return CORE::Materials::mes_generalizedgenmax;
      }
      //@}

      /// Read material parameters
      void read_material_parameters(int& numbranch,  ///< number of viscoelastic branches
          const std::vector<int>*& matids,           ///< material IDs of the viscoelastic branches
          std::string& solve  /// variant of the solution of the evolution integral
          ) override;

      /// @name Access methods
      //@{
      std::vector<std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>> GetBranchespotsum() const
      {
        return branchespotsum_;
      }
      //@}

      /// Indicator for formulation
      void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< general indicator, if one viscoelastic formulation is used
          ) override
      {
        viscogeneral = true;
        return;
      };

      /// Indicator for the chosen viscoelastic formulations
      void specify_visco_formulation(
          bool& isovisco,     ///< global indicator for isotropic, splitted and viscous formulation
          bool& viscogenmax,  ///< global indicator for viscous contribution according the SLS-Model
          bool& viscogeneralizedgenmax,  ///< global indicator for viscoelastic contribution
                                         ///< according to the generalized Maxwell Model
          bool& viscofract  ///< global indicator for viscous contribution according the FSLS-Model
          ) override
      {
        viscogeneralizedgenmax = true;
        return;
      };


     private:
      /// my material parameters
      MAT::ELASTIC::PAR::GeneralizedGenMax* params_;

     protected:
      /// summands of the GeneralizedGenMax material or each branch
      std::vector<std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>> branchespotsum_;
      /// summands in one particular branch
      std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>> internalpotsum_;
    };

    class ViscoBranch : public Summand
    {
     public:
      /// constructor with given material parameters
      ViscoBranch(MAT::ELASTIC::PAR::ViscoBranch* params);

      /// @name Access material constants
      //@{

      /// material type
      CORE::Materials::MaterialType MaterialType() const override
      {
        return CORE::Materials::mes_viscobranch;
      }

      //@}

      /// Read material parameters
      void read_material_parameters(double& nummat,  ///< number of materials in one branch
          const std::vector<int>*& matids            ///< matierial IDs of each part of the branch
          ) override;

      /// Indicator for formulation
      void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< general indicator, if one viscoelastic formulation is used
          ) override
      {
        viscogeneral = true;
        return;
      };


     private:
      /// my material parameters
      MAT::ELASTIC::PAR::ViscoBranch* params_;

    };  // class ViscoBranch


    class ViscoPart : public Summand
    {
     public:
      /// constructor with given material parameters
      ViscoPart(MAT::ELASTIC::PAR::ViscoPart* params);

      /// @name Access material constants
      //@{

      /// material type
      CORE::Materials::MaterialType MaterialType() const override
      {
        return CORE::Materials::mes_viscopart;
      }

      //@}

      /// Read material parameters
      virtual void read_material_parameters(
          double& tau  ///< viscous contribution to viscoelastic part
      );

      /// Indicator for formulation
      void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< general indicator, if one viscoelastic formulation is used
          ) override
      {
        viscogeneral = true;
        return;
      };

      /// Indicator for the chosen viscoelastic formulations
      void specify_visco_formulation(
          bool& isovisco,     ///< global indicator for isotropic, splitted and viscous formulation
          bool& viscogenmax,  ///< global indicator for viscous contribution according the SLS-Model
          bool& viscogeneralizedgenmax,  ///< global indicator for viscoelastic contribution
                                         ///< according to the generalized Maxwell Model
          bool& viscofract  ///< global indicator for viscous contribution according the FSLS-Model
          ) override
      {
        return;
      };

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::ViscoPart* params_;

    };  // class ViscoPart

  }  // namespace ELASTIC
}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
