/*----------------------------------------------------------------------*/
/*!
\file visco_generalizedgenmax.H

\brief generalized maxwell model
\level 2
\maintainer Andreas Rauch

*/
/*----------------------------------------------------------------------*/
#ifndef VISCO_GENERALIZEDGENMAX_H_
#define VISCO_GENERALIZEDGENMAX_H_

#include "../drt_mat/matpar_parameter.H"
#include "../drt_matelast/elast_summand.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*----------------------------------------------------------------------*/
      /// material parameters for viscous contribution to a viscoelastic branch of a generalized
      /// Maxwell model
      ///
      /// <h3>Input line</h3>
      ///  MAT 1 VISCO_GeneralizedGenMax NUMBRANCH 3 MATIDS 4 5 6 SOLVE CONVOL
      ///  MAT 1 VISCO_GeneralizedGenMax NUMBRANCH 3 MATIDS 4 5 6 SOLVE OST
      class GeneralizedGenMax : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        GeneralizedGenMax(Teuchos::RCP<MAT::PAR::Material> matdata);

        /// destructor
        virtual ~GeneralizedGenMax() { ; }

        /// @name material parameters
        //@{

        /// material parameters
        int numbranch_;
        const std::vector<int>* matids_;
        std::string solve_;
        //@}

        /// create material instance of matching type with my parameters

        virtual Teuchos::RCP<MAT::Material> CreateMaterial() { return Teuchos::null; };

        /// @name Write methods
        //@{
        void SetNumbranch(const int& numbranch) { numbranch_ = numbranch; }
        void SetMatids(const std::vector<int>* matids) { matids_ = matids; }
        void SetSolve(const std::string& solve) { solve_ = solve; }
        //@}

      };  // class GeneralizedGenMax


      /*----------------------------------------------------------------------*/
      /// material parameters for viscous contribution to a viscoelastic branch of a generalized
      /// Maxwell modell
      ///
      /// <h3>Input line</h3>
      /// MAT 1 VISCO_BRANCH NUMMAT 2 MATIDS 2 3
      class ViscoBranch : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        ViscoBranch(Teuchos::RCP<MAT::PAR::Material> matdata);

        /// destructor
        virtual ~ViscoBranch() { ; }

        /// @name material parameters
        //@{

        /// material parameters

        double nummat_;
        const std::vector<int>* matids_;


        //@}

        /// create material instance of matching type with my parameters
        virtual Teuchos::RCP<MAT::Material> CreateMaterial() { return Teuchos::null; };

        /// @name Write methods
        //@{
        void SetNummat(const double nummat) { nummat_ = nummat; }
        void SetMatids(const std::vector<int>* matids) { matids_ = matids; }
        //@}

      };  // class ViscoBranch

      /*----------------------------------------------------------------------*/
      /// material parameters for viscous contribution to a viscoelastic branch of a generalized
      /// Maxwell modell
      ///
      /// <h3>Input line</h3>
      /// MAT 1 VISCO_PART TAU 1.5
      class ViscoPart : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        ViscoPart(Teuchos::RCP<MAT::PAR::Material> matdata);

        /// destructor
        virtual ~ViscoPart() { ; }

        /// @name material parameters
        //@{

        /// material parameters
        double tau_;

        //@}

        /// create material instance of matching type with my parameters
        virtual Teuchos::RCP<MAT::Material> CreateMaterial() { return Teuchos::null; };

        /// @name Write methods
        //@{

        void SetTau(const double tau) { tau_ = tau; }
        //@}

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
      virtual INPAR::MAT::MaterialType MaterialType() const
      {
        return INPAR::MAT::mes_generalizedgenmax;
      }


      //@}

      /// Read material parameters
      virtual void ReadMaterialParameters(int& numbranch,  ///< number of viscoelastic branches
          const std::vector<int>*& matids,  ///< material IDs of the viscoelastic branches
          std::string& solve                /// variant of the solution of the evolution integral
      );

      /// @name Access methods
      //@{
      double Numbranch() const { return params_->numbranch_; }
      const std::vector<int>* Matids() const { return params_->matids_; }
      const std::vector<std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>> GetBranchespotsum() const
      {
        return branchespotsum_;
      }
      std::string Solve() const { return params_->solve_; }
      //@}


      /// Indicator for formulation
      virtual void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< general indicator, if one viscoelastic formulation is used
      )
      {
        viscogeneral = true;
        return;
      };

      /// Indicator for the chosen viscoelastic formulations
      virtual void SpecifyViscoFormulation(
          bool& isovisco,     ///< global indicator for isotropic, splitted and viscous formulation
          bool& viscogenmax,  ///< global indicator for viscous contribution according the SLS-Model
          bool& viscogeneralizedgenmax,  ///< global indicator for viscoelastic contribution
                                         ///< according to the generalized Maxwell Model
          bool& viscofract  ///< global indicator for viscous contribution according the FSLS-Model
      )
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
      virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_viscobranch; }


      //@}

      /// Read material parameters
      virtual void ReadMaterialParameters(double& nummat,  ///< number of materials in one branch
          const std::vector<int>*& matids  ///< matierial IDs of each part of the branch
      );

      /// @name Access methods
      //@{
      double Nummat() const { return params_->nummat_; }
      const std::vector<int>* Matids() const { return params_->matids_; }

      //@}

      /// Indicator for formulation
      virtual void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< general indicator, if one viscoelastic formulation is used
      )
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
      virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_viscopart; }


      //@}

      /// Read material parameters
      virtual void ReadMaterialParameters(
          double& tau  ///< viscous contribution to viscoelastic part
      );

      /// @name Access methods
      //@{
      double Tau() const { return params_->tau_; }

      //@}

      /// Indicator for formulation
      virtual void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< general indicator, if one viscoelastic formulation is used
      )
      {
        viscogeneral = true;
        return;
      };

      /// Indicator for the chosen viscoelastic formulations
      virtual void SpecifyViscoFormulation(
          bool& isovisco,     ///< global indicator for isotropic, splitted and viscous formulation
          bool& viscogenmax,  ///< global indicator for viscous contribution according the SLS-Model
          bool& viscogeneralizedgenmax,  ///< global indicator for viscoelastic contribution
                                         ///< according to the generalized Maxwell Model
          bool& viscofract  ///< global indicator for viscous contribution according the FSLS-Model
      )
      {
        return;
      };


     private:
      /// my material parameters
      MAT::ELASTIC::PAR::ViscoPart* params_;

    };  // class ViscoPart

  }  // namespace ELASTIC
}  // namespace MAT



#endif /* VISCO_GENERALIZEDGENMAX_H_ */
