/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains the viscohyperelastic material toolbox.
The availabele viscous models can be applied to any hyperelastic law
of the Elasthyper toolbox.
The viscous part is summed up with the hyperelastic laws
to build a viscohyperelastic material model.

The input line should read
MAT 0   MAT_ViscoElastHyper   NUMMAT 2 MATIDS 1 2 DENS 0

\level 1

*/

/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_VISCOELASTHYPER_HPP
#define FOUR_C_MAT_VISCOELASTHYPER_HPP



#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_elasthyper.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration due to avoid header definition
namespace MAT
{
  namespace ELASTIC
  {
    class Summand;
  }

  // forward declaration
  class ViscoElastHyper;

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// Collection of viscohyperelastic materials
    ///
    /// Storage map of hyperelastic summands.


    class ViscoElastHyper : public MAT::PAR::ElastHyper
    //    class ViscoElastHyper : public CORE::MAT::PAR::Parameter
    {
      friend class MAT::ViscoElastHyper;

     public:
      /// standard constructor
      ///
      /// This constructor recursively calls the constructors of the
      /// parameter sets of the hyperelastic summands.
      ViscoElastHyper(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<CORE::MAT::Material> CreateMaterial() override;

      //@}

    };  // class ViscoElastHyper

  }  // namespace PAR

  class ViscoElastHyperType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "ViscoElastHyperType"; }

    static ViscoElastHyperType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ViscoElastHyperType instance_;
  };


  /*----------------------------------------------------------------------*/
  /// Collection of hyperelastic materials
  ///
  /// This collection offers to additively compose a stress response
  /// based on summands defined separately.  This is possible, because
  /// we deal with hyperelastic materials, which are composed
  /// of (Helmholtz free energy density) potentials.  Effectively, we want
  ///\f[
  ///  \Psi(\boldsymbol{C}) = \sum_i \Psi_i(\boldsymbol{C})
  ///\f]
  /// in which the individual \f$\Psi_i\f$ is implemented as #MAT::ELASTIC::Summand.
  ///
  /// Quite often the right Cauchy-Green 2-tensor \f$\boldsymbol{C}\f$
  /// is replaced by its various invariant forms as argument.
  ///
  /// The task of ElastHyper is the evaluation of the
  /// potential energies and their derivatives to obtain the actual
  /// stress response and the elasticity tensor. The storage is located
  /// at the associated member #params_.
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
  /// </ul>
  ///
  /// \author rausch,tk,bborn
  /// \date 05/09
  class Material;

  class ViscoElastHyper : public MAT::ElastHyper
  {
   public:
    /// construct empty material object
    ViscoElastHyper();

    /// construct the material object given material parameters
    explicit ViscoElastHyper(MAT::PAR::ViscoElastHyper* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of parobject.H (this file) and should return it in this method.
    int UniqueParObjectId() const override
    {
      return ViscoElastHyperType::Instance().UniqueParObjectId();
    }

    /// \brief Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by UniqueParObjectId() which will then
    /// identify the exact class on the receiving processor.
    ///
    /// \param data (in/out): char vector to store class information
    void Pack(CORE::COMM::PackBuffer& data) const override;

    /// \brief Unpack data from a char vector into this class
    ///
    /// The vector data contains all information to rebuild the
    /// exact copy of an instance of a class on a different processor.
    /// The first entry in data has to be an integer which is the unique
    /// parobject id defined at the top of this file and delivered by
    /// UniqueParObjectId().
    ///
    /// \param data (in) : vector storing all data to be unpacked into this
    ///                    instance.
    void Unpack(const std::vector<char>& data) override;

    //@}

    /// material type
    CORE::Materials::MaterialType MaterialType() const override
    {
      return CORE::Materials::m_viscoelasthyper;
    }

    /// check if element kinematics and material kinematics are compatible
    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (!(kinem == INPAR::STR::KinemType::nonlinearTotLag))
        FOUR_C_THROW("element and material kinematics are not compatible");
    }

    /// return copy of this material object
    Teuchos::RCP<CORE::MAT::Material> Clone() const override
    {
      return Teuchos::rcp(new ViscoElastHyper(*this));
    }

    /// Check if history variables are already initialized
    virtual bool Initialized() const
    {
      return isinitvis_ && (histscgcurr_ != Teuchos::null);
      return isinitvis_ && (histstresscurr_ != Teuchos::null);
      return isinitvis_ && (histartstresscurr_ != Teuchos::null);
    }

    /// hyperelastic stress response plus elasticity tensor
    void Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const CORE::LINALG::Matrix<6, 1>* glstrain,          ///< Green-Lagrange strain
        Teuchos::ParameterList& params,      ///< Container for additional information
        CORE::LINALG::Matrix<6, 1>* stress,  ///< 2nd Piola-Kirchhoff stresses
        CORE::LINALG::Matrix<6, 6>* cmat,
        int gp,                      ///< Gauss point
        const int eleGID) override;  ///< Constitutive matrix

    /// setup material description
    void Setup(int numgp, INPUT::LineDefinition* linedef) override;

    /// update history variables
    void Update() override;

   protected:
    /// calculates the kinematic quantities and tensors used afterwards for viscous part
    virtual void EvaluateKinQuantVis(CORE::LINALG::Matrix<6, 1>& rcg,
        CORE::LINALG::Matrix<6, 1>& scg, CORE::LINALG::Matrix<6, 1>& icg,
        CORE::LINALG::Matrix<3, 1>& prinv, CORE::LINALG::Matrix<7, 1>& rateinv,
        CORE::LINALG::Matrix<6, 1>& modrcg, Teuchos::ParameterList& params,
        CORE::LINALG::Matrix<6, 1>& scgrate, CORE::LINALG::Matrix<6, 1>& modrcgrate,
        CORE::LINALG::Matrix<7, 1>& modrateinv, int gp);

    /// calculates the factors associated to the viscous laws
    virtual void EvaluateMuXi(CORE::LINALG::Matrix<3, 1>& inv, CORE::LINALG::Matrix<3, 1>& modinv,
        CORE::LINALG::Matrix<8, 1>& mu, CORE::LINALG::Matrix<8, 1>& modmu,
        CORE::LINALG::Matrix<33, 1>& xi, CORE::LINALG::Matrix<33, 1>& modxi,
        CORE::LINALG::Matrix<7, 1>& rateinv, CORE::LINALG::Matrix<7, 1>& modrateinv,
        Teuchos::ParameterList& params, int gp, int eleGID);

    /// calculates the isotropic stress and elasticity tensor for viscous principal configuration
    virtual void EvaluateIsoViscoPrincipal(CORE::LINALG::Matrix<6, 1>& stress,
        CORE::LINALG::Matrix<6, 6>& cmat, CORE::LINALG::Matrix<8, 1>& mu,
        CORE::LINALG::Matrix<33, 1>& xi, CORE::LINALG::Matrix<6, 6>& id4sharp,
        CORE::LINALG::Matrix<6, 1>& scgrate);

    /// calculates the isotropic stress and elasticity tensor for viscous decoupled configuration
    virtual void EvaluateIsoViscoModified(CORE::LINALG::Matrix<6, 1>& stressisomodisovisco,
        CORE::LINALG::Matrix<6, 1>& stressisomodvolvisco,
        CORE::LINALG::Matrix<6, 6>& cmatisomodisovisco,
        CORE::LINALG::Matrix<6, 6>& cmatisomodvolvisco, CORE::LINALG::Matrix<3, 1>& prinv,
        CORE::LINALG::Matrix<3, 1>& modinv, CORE::LINALG::Matrix<8, 1>& modmu,
        CORE::LINALG::Matrix<33, 1>& modxi, CORE::LINALG::Matrix<6, 1>& rcg,
        CORE::LINALG::Matrix<6, 1>& id2, CORE::LINALG::Matrix<6, 1>& icg,
        CORE::LINALG::Matrix<6, 6>& id4, CORE::LINALG::Matrix<6, 1>& modrcgrate);

    /// calculates the stress and elasticitiy tensor for the viscous Genmax-material
    /// depending on the viscoelastic material isochoric-principal, isochoric-modified
    /// (volumetric and isochoric equal treaded) or anisotropic stresses and elasticity
    /// tensors are added
    virtual void EvaluateViscoGenMax(CORE::LINALG::Matrix<6, 1>* stress,
        CORE::LINALG::Matrix<6, 6>* cmat, CORE::LINALG::Matrix<6, 1>& Q,
        CORE::LINALG::Matrix<6, 6>& cmatq, Teuchos::ParameterList& params, int gp);

    /// calculates the stress and elasticitiy tensor for the GeneneralizedMax-material
    virtual void EvaluateViscoGeneralizedGenMax(CORE::LINALG::Matrix<6, 1>& Q,
        CORE::LINALG::Matrix<6, 6>& cmatq, Teuchos::ParameterList& params,
        const CORE::LINALG::Matrix<6, 1>* glstrain, int gp, int eleGID);

    /// calculates the stress and elasticity tensor for the viscos Fract-material
    /// depending on the viscoelastic material isochoric-principal, isochoric-modified
    /// (volumetric and isochoric equal treaded) or anisotropic stresses and elasticity
    /// tensors are added
    virtual void EvaluateViscoFract(CORE::LINALG::Matrix<6, 1> stress,
        CORE::LINALG::Matrix<6, 6> cmat, CORE::LINALG::Matrix<6, 1>& Q,
        CORE::LINALG::Matrix<6, 6>& cmatq, Teuchos::ParameterList& params, int gp);

    /// @name Flags to specify the viscous formulations
    //@{
    bool isovisco_;     ///< global indicator for isotropic splitted viscous formulation
    bool viscogenmax_;  ///< global indicator for viscous contribution according the SLS-Model
    bool viscogeneralizedgenmax_;  ///< global indicator for viscous contribution of the branches
                                   ///< according to the generalized Maxwell model
    bool viscofract_;  ///< global indicator for viscous contribution according the FSLS-Model
                       //@}

   private:
    /// viscous history Cauchy-Green-Tensor
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<6, 1>>>
        histscgcurr_;  ///< current Cauchy-Green-Tensor
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<6, 1>>>
        histscglast_;  ///< Cauchy-Green-Tensor of last converged state
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<6, 1>>>
        histmodrcgcurr_;  ///< current decoupled Cauchy-Green-Tensor
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<6, 1>>>
        histmodrcglast_;  ///< decoupled Cauchy-Green-Tensor of last converged state
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>>
        histstresscurr_;  ///< current stress
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>>
        histstresslast_;  ///< stress of last converged state
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>>
        histartstresscurr_;  ///< current stress
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>>
        histartstresslast_;  ///< stress of last converged state
    Teuchos::RCP<std::vector<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>>>
        histbranchstresscurr_;  ///< current stress of the branches
    Teuchos::RCP<std::vector<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>>>
        histbranchstresslast_;  ///< stress of the branches of last converged state
    Teuchos::RCP<std::vector<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>>>
        histbranchelaststresscurr_;  ///< current elastic stress of the branches
    Teuchos::RCP<std::vector<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>>>
        histbranchelaststresslast_;  ///< elastic stress of the branches of last converged state
    Teuchos::RCP<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>>
        histfractartstresscurr_;  ///< current artificial stress of fractional SLS-model
    Teuchos::RCP<std::vector<std::vector<CORE::LINALG::Matrix<NUM_STRESS_3D, 1>>>>
        histfractartstresslastall_;  ///< artificial fractional stress of all last converged states

    bool isinitvis_;  ///< indicates if #Initialized routine has been called
  };                  // class ViscoElastHyper

}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
