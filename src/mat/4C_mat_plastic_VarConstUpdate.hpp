/*----------------------------------------------------------------------*/
/*! \file
\level 3
\brief finite deformation plasticity algorithm based on
       variational constitutive update
*/

/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_PLASTIC_VARCONSTUPDATE_HPP
#define FOUR_C_MAT_PLASTIC_VARCONSTUPDATE_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_plasticelasthyper.hpp"
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
  class PlasticElastHyperVCU;

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// Collection of hyperelastic materials
    ///
    /// Storage map of hyperelastic summands.
    class PlasticElastHyperVCU : public MAT::PAR::PlasticElastHyper
    {
      friend class MAT::PlasticElastHyperVCU;

     public:
      /// standard constructor
      ///
      /// This constructor recursively calls the constructors of the
      /// parameter sets of the hyperelastic summands.
      PlasticElastHyperVCU(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      /// @name material parameters
      //@{

      /// provide access to material/summand by its ID
      Teuchos::RCP<const MAT::ELASTIC::Summand> MaterialById(
          const int id  ///< ID to look for in collection of summands
      ) const;

      /// create material instance of matching type with my parameters
      Teuchos::RCP<CORE::MAT::Material> CreateMaterial() override;
      //@}

    };  // class PlasticElastHyper

  }  // namespace PAR

  class PlasticElastHyperVCUType : public CORE::COMM::ParObjectType
  {
   public:
    std::string Name() const override { return "PlasticElastHyperVCUType"; }

    static PlasticElastHyperVCUType& Instance() { return instance_; };

    CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

   private:
    static PlasticElastHyperVCUType instance_;
  };


  /*----------------------------------------------------------------------*/
  /// Collection of hyperelastic materials
  class PlasticElastHyperVCU : public MAT::PlasticElastHyper
  {
   public:
    /// construct empty material object
    PlasticElastHyperVCU();

    /// construct the material object given material parameters
    explicit PlasticElastHyperVCU(MAT::PAR::PlasticElastHyperVCU* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of parobject.H (this file) and should return it in this method.
    int UniqueParObjectId() const override
    {
      return PlasticElastHyperVCUType::Instance().UniqueParObjectId();
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
      return CORE::Materials::m_plelasthyperVCU;
    }

    /// return copy of this material object
    Teuchos::RCP<CORE::MAT::Material> Clone() const override
    {
      return Teuchos::rcp(new PlasticElastHyperVCU(*this));
    }


    /// hyperelastic stress response plus elasticity tensor
    /// (pure virtual in material base class. Not allowed here)
    void Evaluate(const CORE::LINALG::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const CORE::LINALG::Matrix<6, 1>* glstrain,          ///< Green-Lagrange strain
        Teuchos::ParameterList& params,      ///< Container for additional information
        CORE::LINALG::Matrix<6, 1>* stress,  ///< 2nd Piola-Kirchhoff stresses
        CORE::LINALG::Matrix<6, 6>* cmat,    ///< Constitutive matrix
        int gp,                              ///< Gauss point
        int eleGID) override;                ///< Element GID

    /// setup material data
    void Setup(int numgp, INPUT::LineDefinition* linedef) override;

    /// update sumands
    void Update() override;

    //! return names of visualization data
    void VisNames(std::map<std::string, int>& names) override;

    //! return visualization data
    bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID) override;

    // ********************************************************* //
    // ** here are some functions we don't need variationally ** //
    // ** consistent constitutive updates. They all give an   ** //
    // ** error when called. They are only here because of the * //
    // ** inheritance                                         ** //
    // ********************************************************* //

    /// evaluate quantities for elastic stiffness matrix
    /// in consideration of plastic history/deformation
    virtual void EvaluateElast(const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<3, 3>* deltaLp, Teuchos::ParameterList& params,
        CORE::LINALG::Matrix<6, 1>* pk2, CORE::LINALG::Matrix<6, 6>* cmat, const int gp,
        const int eleGID)
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
    }

    /// evaluate stresses and stiffness contribution
    /// due to thermal expansion
    virtual void evaluate_thermal_stress(const CORE::LINALG::Matrix<3, 3>* defgrd,
        const double temp, Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>* pk2,
        CORE::LINALG::Matrix<6, 6>* cmat, const int gp, const int eleGID)
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
    }

    /// evaluate stresses and stiffness contribution
    /// due to thermal expansion
    virtual void EvaluateCTvol(const CORE::LINALG::Matrix<3, 3>* defgrd,
        Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>* cTvol,
        CORE::LINALG::Matrix<6, 6>* dCTvoldE, const int gp, const int eleGID)
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
    }

    /// evaluate everything needed for the condensation of the plastic deformation
    /// at element level. (with zero plastic spin)
    virtual void EvaluatePlast(const CORE::LINALG::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const CORE::LINALG::Matrix<3, 3>* deltaDp,  ///< symmetric part of plastic flow increment
        const double temp,                          ///< current temperature
        Teuchos::ParameterList& params,             ///< Container for additional information
        CORE::LINALG::Matrix<6, 6>* dPK2dDp,        ///< derivative of PK2 w.r.t. F_p^{-1}
        CORE::LINALG::Matrix<6, 1>* NCP,            ///< NCP function
        CORE::LINALG::Matrix<6, 6>* dNCPdC,         ///< derivative of NCP function w.r.t. RCG
        CORE::LINALG::Matrix<6, 6>* dNCPdDp,        ///< derivative of NCP function w.r.t. deltaLp
        bool* active,                               ///< gauss point is active
        bool* elast,         ///< gauss point needs condensation if it is not elast
        bool* as_converged,  ///< convergence of active set (false, if as has changed)
        const int gp,        ///< gauss point
        CORE::LINALG::Matrix<6, 1>*
            dNCPdT,  ///< derivative of NCP function w.r.t. temperature (only in TSI case)
        CORE::LINALG::Matrix<6, 1>* dHdC,  ///< derivative of Heating w.r.t. RCG (only in TSI case)
        CORE::LINALG::Matrix<6, 1>*
            dHdDp,        ///< derivative of Heating w.r.t. deltaLp (only in TSI case)
        const double dt,  ///< time step size
        const int eleGID  ///< global ID of element
    )
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
    }

    /// evaluate everything needed for the condensation of the plastic deformation
    /// at element level. (with plastic spin)
    virtual void EvaluatePlast(const CORE::LINALG::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const CORE::LINALG::Matrix<3, 3>*
            deltaLp,                          ///< plastic deformation gradient (non-symmetric)
        const double temp,                    ///< current temperature
        Teuchos::ParameterList& params,       ///< Container for additional information
        CORE::LINALG::Matrix<6, 9>* dPK2dLp,  ///< derivative of PK2 w.r.t. F_p^{-1}
        CORE::LINALG::Matrix<9, 1>* NCP,      ///< NCP function
        CORE::LINALG::Matrix<9, 6>* dNCPdC,   ///< derivative of NCP function w.r.t. RCG
        CORE::LINALG::Matrix<9, 9>* dNCPdLp,  ///< derivative of NCP function w.r.t. deltaLp
        bool* active,                         ///< gauss point is active
        bool* elast,                          ///< gauss point needs condensation if it is not elast
        bool* as_converged,  ///< convergence of active set (false, if as has changed)
        const int gp,        ///< gauss point
        CORE::LINALG::Matrix<9, 1>*
            dNCPdT,  ///< derivative of NCP function w.r.t. temperature (only in TSI case)
        CORE::LINALG::Matrix<6, 1>* dHdC,  ///< derivative of Heating w.r.t. RCG (only in TSI case)
        CORE::LINALG::Matrix<9, 1>*
            dHdLp,        ///< derivative of Heating w.r.t. deltaLp (only in TSI case)
        const double dt,  ///< time step size
        const int eleGID  ///< global ID of element
    )
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
    }

    /// update plastic history variables
    void UpdateGP(const int gp, const CORE::LINALG::Matrix<3, 3>* deltaDp) override
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
    }

    /// get plastic algorithm parameters
    void GetParams(double s, double cpl) override
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
    }

    /// is this GP active
    virtual bool Active(int gp)
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
      return false;
    }

    /// heating at this gp
    double& HepDiss(int gp) override
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
      static double a = 0.;
      return a;
    }

    /// derivative of heating at this gp
    CORE::LINALG::SerialDenseVector& dHepDissDd(int gp) override
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
      static CORE::LINALG::SerialDenseVector tmp(0);
      return tmp;
    }

    // derivative of heating w.r.t. temperature
    double& dHepDT(int gp) override
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
      static double a = 0.;
      return a;
    }

    // derivative of heating at each gp w.r.t. nodal temperature vector
    // (only EAS contribution)
    Teuchos::RCP<std::vector<CORE::LINALG::SerialDenseVector>> dHepDTeas() override
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
      return Teuchos::null;
    }

   protected:
    virtual void EvaluateNCP(const CORE::LINALG::Matrix<3, 3>* mStr,
        const CORE::LINALG::Matrix<6, 6>* dMdC, const CORE::LINALG::Matrix<6, 9>* dMdFpinv,
        const CORE::LINALG::Matrix<6, 9>* dPK2dFpinv, const CORE::LINALG::Matrix<3, 3>* deltaDp,
        const int gp, const double temp, CORE::LINALG::Matrix<6, 1>* NCP,
        CORE::LINALG::Matrix<6, 6>* dNCPdC, CORE::LINALG::Matrix<6, 6>* dNCPdDp,
        CORE::LINALG::Matrix<6, 1>* dNCPdT, CORE::LINALG::Matrix<6, 6>* dPK2dDp, bool* active,
        bool* elast, bool* as_converged, CORE::LINALG::Matrix<6, 1>* dHdC,
        CORE::LINALG::Matrix<6, 1>* dHdDp, Teuchos::ParameterList& params, const double dt)
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
    }


    virtual void EvaluateNCPandSpin(const CORE::LINALG::Matrix<3, 3>* mStr,
        const CORE::LINALG::Matrix<6, 6>* dMdC, const CORE::LINALG::Matrix<6, 9>* dMdFpinv,
        const CORE::LINALG::Matrix<6, 9>* dPK2dFpinv, const CORE::LINALG::Matrix<3, 3>* deltaLp,
        const int gp, CORE::LINALG::Matrix<9, 1>* NCP, CORE::LINALG::Matrix<9, 6>* dNCPdC,
        CORE::LINALG::Matrix<9, 9>* dNCPdLp, CORE::LINALG::Matrix<6, 9>* dPK2dLp, bool* active,
        bool* elast, bool* as_converged, const double dt)
    {
      FOUR_C_THROW("Don't need this for Variationally consistent constitutive update");
    }



    virtual void EvalDceDlp(const CORE::LINALG::Matrix<3, 3> fpi,
        const CORE::LINALG::Matrix<3, 3>* defgrd, const CORE::LINALG::Matrix<6, 6> Dexp,
        const CORE::LINALG::Matrix<3, 3> cetrial, const CORE::LINALG::Matrix<3, 3> explp,
        CORE::LINALG::Matrix<6, 6>& dceDdeltalp, CORE::LINALG::Matrix<9, 6>& dFpiDdeltaDp);



    virtual void EvaluateRHS(const int gp, const CORE::LINALG::Matrix<3, 3> dLp,
        const CORE::LINALG::Matrix<3, 3> defgrd, CORE::LINALG::Matrix<6, 1>& eeOut,
        CORE::LINALG::Matrix<5, 1>& rhs, CORE::LINALG::Matrix<5, 1>& rhsElast,
        CORE::LINALG::Matrix<6, 6>& dcedlp, CORE::LINALG::Matrix<9, 6>& dFpiDdeltaDp,
        Teuchos::ParameterList& params, const int eleGID);


    virtual void YieldFunction(const double last_ai, const double norm_dLp,
        const CORE::LINALG::Matrix<3, 3> ExpEqui, const CORE::LINALG::Matrix<3, 3> cetr,
        const CORE::LINALG::Matrix<6, 1> str, double* yieldFunc,
        CORE::LINALG::Matrix<3, 3>& devMandelStr, CORE::LINALG::Matrix<3, 3>& MandelStr);

    virtual void CompElastQuant(const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<3, 3> fpi, const CORE::LINALG::Matrix<3, 3> MatExp,
        CORE::LINALG::Matrix<3, 3>* cetrial, CORE::LINALG::Matrix<6, 1>* Ee);

    virtual void matrix_exponential_second_derivative_sym3x3x6(
        const CORE::LINALG::Matrix<3, 3> MatrixIn, CORE::LINALG::Matrix<3, 3>& exp,
        CORE::LINALG::Matrix<6, 6>& dexp_mat, CORE::LINALG::Matrix<6, 6>* MatrixExp2ndDerivVoigt);

    virtual void matrix_exponential_second_derivative_sym3x3(
        const CORE::LINALG::Matrix<3, 3> MatrixIn, CORE::LINALG::Matrix<3, 3>& exp,
        std::vector<CORE::LINALG::Matrix<3, 3>>& MatrixExp1stDeriv,
        std::vector<std::vector<CORE::LINALG::Matrix<3, 3>>>& MatrixExp2ndDeriv);

    virtual void EvaluatePlast(CORE::LINALG::Matrix<6, 9>& dPK2dFpinvIsoprinc,
        const CORE::LINALG::Matrix<3, 1>& gamma, const CORE::LINALG::Matrix<8, 1>& delta,
        const CORE::LINALG::Matrix<3, 3>& id2, const CORE::LINALG::Matrix<6, 1>& Cpi,
        const CORE::LINALG::Matrix<3, 3>& Fpi, const CORE::LINALG::Matrix<3, 3>& CpiC,
        const CORE::LINALG::Matrix<9, 1>& CFpi, const CORE::LINALG::Matrix<9, 1>& CFpiCei,
        const CORE::LINALG::Matrix<6, 1>& ircg, const CORE::LINALG::Matrix<3, 3>& FpiCe,
        const CORE::LINALG::Matrix<9, 1>& CFpiCe, const CORE::LINALG::Matrix<6, 1>& CpiCCpi);

    virtual void evaluate_kin_quant_plast(int gp, int eleGID,
        const CORE::LINALG::Matrix<3, 3>* defgrd, const CORE::LINALG::Matrix<3, 3>* fpi,
        CORE::LINALG::Matrix<3, 1>& gamma, CORE::LINALG::Matrix<8, 1>& delta,
        CORE::LINALG::Matrix<3, 3>& id2, CORE::LINALG::Matrix<6, 1>& Cpi,
        CORE::LINALG::Matrix<3, 3>& CpiC, CORE::LINALG::Matrix<9, 1>& CFpi,
        CORE::LINALG::Matrix<9, 1>& CFpiCei, CORE::LINALG::Matrix<6, 1>& ircg,
        CORE::LINALG::Matrix<3, 3>& FpiCe, CORE::LINALG::Matrix<9, 1>& CFpiCe,
        CORE::LINALG::Matrix<6, 1>& CpiCCpi);

    virtual void Dpk2dFpi(int gp, int eleGID, const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<3, 3>* fpi, CORE::LINALG::Matrix<6, 9>& dPK2dFpinvIsoprinc);

    virtual void DpsiplastDalphaiso(const double norm_dLp, const double last_alphaiso,
        const double isoHardMod, const double initYield, const double infYield,
        const double expIsoHard, double* dpsiplastdalphaiso);

    virtual void Ce2ndDeriv(const CORE::LINALG::Matrix<3, 3>* defgrd,
        const CORE::LINALG::Matrix<3, 3> fpi, const CORE::LINALG::Matrix<3, 3> dLp,
        CORE::LINALG::Matrix<6, 6>* DDceDdLpDdLpVoigt);


    /// Access to material params
    MAT::PAR::PlasticElastHyperVCU* MatParams() const override { return params_; }

    /// get dissipation mode
    INPAR::TSI::DissipationMode DisMode() const override { return INPAR::TSI::pl_flow; }

    /// inverse plastic deformation gradient for each Gauss point at current state
    std::vector<CORE::LINALG::Matrix<3, 3>> plastic_defgrd_inverse_;

    /// my material parameters
    MAT::PAR::PlasticElastHyperVCU* params_;
  };

}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
