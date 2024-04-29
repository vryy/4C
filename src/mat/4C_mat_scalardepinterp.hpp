/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains a framework, which allows for the linear interpolation between two different
strain energy functions, where the interpolation ratio come from 'outside'. For now this is supposed
to come from the idea, that there is a scalar quantity (e.g. foam cells) which induces growth. The
grown volume fraction thereby consists of another material (e.g. softer or stiffer material
parameters or totally different strain energy function).

\level 2



*/

/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_SCALARDEPINTERP_HPP
#define FOUR_C_MAT_SCALARDEPINTERP_HPP



#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_input_base.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration due to avoid header definition
namespace MAT
{
  // forward declaration
  class ScalarDepInterp;

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /*! \class ScalarDepInterp
     *  \brief Common parameters for scalar dependent interpolation between two material laws
     *
     *  \author thon
     *  \date 12/2015
     */
    class ScalarDepInterp : public CORE::MAT::PAR::Parameter
    {
     public:
      /// standard constructor
      ScalarDepInterp(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{
      /// elastic material number
      const int id_lambda_zero_;
      /// growthlaw material number
      const int id_lambda_unit_;
      //@}

      // create material instance of matching type with my parameters
      Teuchos::RCP<CORE::MAT::Material> CreateMaterial() override;

    };  // class ScalarDepInterp

  }  // namespace PAR

  class ScalarDepInterpType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "ScalarDepInterpType"; }

    static ScalarDepInterpType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ScalarDepInterpType instance_;
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

  class ScalarDepInterp : public So3Material
  {
   public:
    /// construct empty material object
    ScalarDepInterp();

    /// construct the material object given material parameters
    explicit ScalarDepInterp(MAT::PAR::ScalarDepInterp* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of parobject.H (this file) and should return it in this method.
    int UniqueParObjectId() const override
    {
      return ScalarDepInterpType::Instance().UniqueParObjectId();
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
      return CORE::Materials::m_sc_dep_interp;
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
      return Teuchos::rcp(new ScalarDepInterp(*this));
    }

    /// number of materials
    //    virtual int NumMat() const { return params_->nummat_; }

    /// Update
    void Update() override;

    /// Reset time step
    void ResetStep() override;

    double Density() const override
    {
      // Note: we have already check that both densities are equal in Setuo()
      return lambda_zero_mat_->Density();
    };

    /// provide access to material by its ID
    //     virtual Teuchos::RCP<const MAT::ELASTIC::Summand> MaterialById(const int id) const {
    //     return params_->MaterialById(id); }


    /// hyperelastic stress response plus elasticity tensor
    void Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const CORE::LINALG::Matrix<6, 1>* glstrain,          ///< Green-Lagrange strain
        Teuchos::ParameterList& params,      ///< Container for additional information
        CORE::LINALG::Matrix<6, 1>* stress,  ///< 2nd Piola-Kirchhoff stresses
        CORE::LINALG::Matrix<6, 6>* cmat,    ///< Constitutive matrix
        int gp,                              ///< Gauss point
        int eleGID                           ///< Element GID
        ) override;

    /// setup
    void Setup(int numgp, INPUT::LineDefinition* linedef) override;

    /// Return quick accessible material parameter data
    CORE::MAT::PAR::Parameter* Parameter() const override { return params_; }

    /// Return names of visualization data
    void VisNames(std::map<std::string, int>& names) override;

    /// Return visualization data
    bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID) override;

    /// Reading material weights (lambda) from input file
    void ReadLambda(INPUT::LineDefinition* linedef, std::string specifier, double& lambda);

    /// evaluate strain energy function
    void StrainEnergy(const CORE::LINALG::Matrix<6, 1>& glstrain, double& psi, const int gp,
        const int eleGID) override;

   private:
    /// my material parameters
    MAT::PAR::ScalarDepInterp* params_;

    /// indicates if material is initialized
    bool isinit_;

    /// elastic material for zero lambda
    Teuchos::RCP<MAT::So3Material> lambda_zero_mat_;

    /// elastic material for unit lambda
    Teuchos::RCP<MAT::So3Material> lambda_unit_mat_;

    /// interpolation parameter (0-->IDMATZEROSC,1-->IDMATUNITSC)
    std::vector<double> lambda_;
  };


}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
