/*----------------------------------------------------------------------*/
/*!
\level 3
\brief finite deformation plasticity algorithm based on
       variational constitutive update
\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------*/
#ifndef PLASTIC_VARCONSTUPDATE_H
#define PLASTIC_VARCONSTUPDATE_H

#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "plasticelasthyper.H"

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
      PlasticElastHyperVCU(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~PlasticElastHyperVCU() { ; }

      /// @name material parameters
      //@{

      /// provide access to material/summand by its ID
      Teuchos::RCP<const MAT::ELASTIC::Summand> MaterialById(
          const int id  ///< ID to look for in collection of summands
          ) const;

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();
      //@}

    };  // class PlasticElastHyper

  }  // namespace PAR

  class PlasticElastHyperVCUType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "PlasticElastHyperVCUType"; }

    static PlasticElastHyperVCUType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

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
    /// top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const
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
    virtual void Pack(DRT::PackBuffer& data) const;

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
    virtual void Unpack(const std::vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_plelasthyperVCU; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new PlasticElastHyperVCU(*this));
    }


    /// hyperelastic stress response plus elasticity tensor
    /// (pure virtual in material base class. Not allowed here)
    virtual void Evaluate(const LINALG::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const LINALG::Matrix<6, 1>* glstrain,                  ///< Green-Lagrange strain
        Teuchos::ParameterList& params,  ///< Container for additional information
        LINALG::Matrix<6, 1>* stress,    ///< 2nd Piola-Kirchhoff stresses
        LINALG::Matrix<6, 6>* cmat,      ///< Constitutive matrix
        const int eleGID);               ///< Element GID

    /// setup material data
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    /// update sumands
    virtual void Update();

    //! return names of visualization data
    virtual void VisNames(std::map<std::string, int>& names);

    //! return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

    // ********************************************************* //
    // ** here are some functions we don't need variationally ** //
    // ** consistent constitutive updates. They all give an   ** //
    // ** error when called. They are only here because of the * //
    // ** inheritance                                         ** //
    // ********************************************************* //

    /// evaluate quantities for elastic stiffness matrix
    /// in consideration of plastic history/deformation
    virtual void EvaluateElast(const LINALG::Matrix<3, 3>* defgrd,
        const LINALG::Matrix<3, 3>* deltaLp, Teuchos::ParameterList& params,
        LINALG::Matrix<6, 1>* pk2, LINALG::Matrix<6, 6>* cmat, const int gp, const int eleGID)
    {
      dserror("Don't need this for Variationally consistent constitutive update");
    }

    /// evaluate stresses and stiffness contribution
    /// due to thermal expansion
    virtual void EvaluateThermalStress(const LINALG::Matrix<3, 3>* defgrd, const double temp,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>* pk2, LINALG::Matrix<6, 6>* cmat,
        const int gp, const int eleGID)
    {
      dserror("Don't need this for Variationally consistent constitutive update");
    }

    /// evaluate stresses and stiffness contribution
    /// due to thermal expansion
    virtual void EvaluateCTvol(const LINALG::Matrix<3, 3>* defgrd, Teuchos::ParameterList& params,
        LINALG::Matrix<6, 1>* cTvol, LINALG::Matrix<6, 6>* dCTvoldE, const int gp, const int eleGID)
    {
      dserror("Don't need this for Variationally consistent constitutive update");
    }

    /// evaluate everything needed for the condensation of the plastic deformation
    /// at element level. (with zero plastic spin)
    virtual void EvaluatePlast(const LINALG::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const LINALG::Matrix<3, 3>* deltaDp,  ///< symmetric part of plastic flow increment
        const double temp,                    ///< current temperature
        Teuchos::ParameterList& params,       ///< Container for additional information
        LINALG::Matrix<6, 6>* dPK2dDp,        ///< derivative of PK2 w.r.t. F_p^{-1}
        LINALG::Matrix<6, 1>* NCP,            ///< NCP function
        LINALG::Matrix<6, 6>* dNCPdC,         ///< derivative of NCP function w.r.t. RCG
        LINALG::Matrix<6, 6>* dNCPdDp,        ///< derivative of NCP function w.r.t. deltaLp
        bool* active,                         ///< gauss point is active
        bool* elast,                          ///< gauss point needs condensation if it is not elast
        bool* as_converged,  ///< convergence of active set (false, if as has changed)
        const int gp,        ///< gauss point
        LINALG::Matrix<6, 1>*
            dNCPdT,  ///< derivative of NCP function w.r.t. temperature (only in TSI case)
        LINALG::Matrix<6, 1>* dHdC,   ///< derivative of Heating w.r.t. RCG (only in TSI case)
        LINALG::Matrix<6, 1>* dHdDp,  ///< derivative of Heating w.r.t. deltaLp (only in TSI case)
        const double dt,              ///< time step size
        const int eleGID              ///< global ID of element
    )
    {
      dserror("Don't need this for Variationally consistent constitutive update");
    }

    /// evaluate everything needed for the condensation of the plastic deformation
    /// at element level. (with plastic spin)
    virtual void EvaluatePlast(const LINALG::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const LINALG::Matrix<3, 3>* deltaLp,  ///< plastic deformation gradient (non-symmetric)
        const double temp,                    ///< current temperature
        Teuchos::ParameterList& params,       ///< Container for additional information
        LINALG::Matrix<6, 9>* dPK2dLp,        ///< derivative of PK2 w.r.t. F_p^{-1}
        LINALG::Matrix<9, 1>* NCP,            ///< NCP function
        LINALG::Matrix<9, 6>* dNCPdC,         ///< derivative of NCP function w.r.t. RCG
        LINALG::Matrix<9, 9>* dNCPdLp,        ///< derivative of NCP function w.r.t. deltaLp
        bool* active,                         ///< gauss point is active
        bool* elast,                          ///< gauss point needs condensation if it is not elast
        bool* as_converged,  ///< convergence of active set (false, if as has changed)
        const int gp,        ///< gauss point
        LINALG::Matrix<9, 1>*
            dNCPdT,  ///< derivative of NCP function w.r.t. temperature (only in TSI case)
        LINALG::Matrix<6, 1>* dHdC,   ///< derivative of Heating w.r.t. RCG (only in TSI case)
        LINALG::Matrix<9, 1>* dHdLp,  ///< derivative of Heating w.r.t. deltaLp (only in TSI case)
        const double dt,              ///< time step size
        const int eleGID              ///< global ID of element
    )
    {
      dserror("Don't need this for Variationally consistent constitutive update");
    }

    /// update plastic history variables
    virtual void UpdateGP(const int gp, const LINALG::Matrix<3, 3>* deltaDp)
    {
      dserror("Don't need this for Variationally consistent constitutive update");
    }

    /// get plastic algorithm parameters
    virtual void GetParams(double s, double cpl)
    {
      dserror("Don't need this for Variationally consistent constitutive update");
    }

    /// is this GP active
    virtual bool Active(int gp)
    {
      dserror("Don't need this for Variationally consistent constitutive update");
      return false;
    }

    /// heating at this gp
    virtual double& HepDiss(int gp)
    {
      dserror("Don't need this for Variationally consistent constitutive update");
      static double a = 0.;
      return a;
    }

    /// derivative of heating at this gp
    virtual Epetra_SerialDenseVector& dHepDissDd(int gp)
    {
      dserror("Don't need this for Variationally consistent constitutive update");
      static Epetra_SerialDenseVector tmp(0);
      return tmp;
    }

    // derivative of heating w.r.t. temperature
    virtual double& dHepDT(int gp)
    {
      dserror("Don't need this for Variationally consistent constitutive update");
      static double a = 0.;
      return a;
    }

    // derivative of heating at each gp w.r.t. nodal temperature vector
    // (only EAS contribution)
    virtual Teuchos::RCP<std::vector<Epetra_SerialDenseVector>> dHepDTeas()
    {
      dserror("Don't need this for Variationally consistent constitutive update");
      return Teuchos::null;
    }

   protected:
    virtual void EvaluateNCP(const LINALG::Matrix<3, 3>* mStr, const LINALG::Matrix<6, 6>* dMdC,
        const LINALG::Matrix<6, 9>* dMdFpinv, const LINALG::Matrix<6, 9>* dPK2dFpinv,
        const LINALG::Matrix<3, 3>* deltaDp, const int gp, const double temp,
        LINALG::Matrix<6, 1>* NCP, LINALG::Matrix<6, 6>* dNCPdC, LINALG::Matrix<6, 6>* dNCPdDp,
        LINALG::Matrix<6, 1>* dNCPdT, LINALG::Matrix<6, 6>* dPK2dDp, bool* active, bool* elast,
        bool* as_converged, LINALG::Matrix<6, 1>* dHdC, LINALG::Matrix<6, 1>* dHdDp,
        Teuchos::ParameterList& params, const double dt)
    {
      dserror("Don't need this for Variationally consistent constitutive update");
    }


    virtual void EvaluateNCPandSpin(const LINALG::Matrix<3, 3>* mStr,
        const LINALG::Matrix<6, 6>* dMdC, const LINALG::Matrix<6, 9>* dMdFpinv,
        const LINALG::Matrix<6, 9>* dPK2dFpinv, const LINALG::Matrix<3, 3>* deltaLp, const int gp,
        LINALG::Matrix<9, 1>* NCP, LINALG::Matrix<9, 6>* dNCPdC, LINALG::Matrix<9, 9>* dNCPdLp,
        LINALG::Matrix<6, 9>* dPK2dLp, bool* active, bool* elast, bool* as_converged,
        const double dt)
    {
      dserror("Don't need this for Variationally consistent constitutive update");
    }



    virtual void EvalDceDlp(const LINALG::Matrix<3, 3> fpi, const LINALG::Matrix<3, 3>* defgrd,
        const LINALG::Matrix<6, 6> Dexp, const LINALG::Matrix<3, 3> cetrial,
        const LINALG::Matrix<3, 3> explp, LINALG::Matrix<6, 6>& dceDdeltalp,
        LINALG::Matrix<9, 6>& dFpiDdeltaDp);



    virtual void EvaluateRHS(const int gp, const LINALG::Matrix<3, 3> dLp,
        const LINALG::Matrix<3, 3> defgrd, LINALG::Matrix<6, 1>& eeOut, LINALG::Matrix<5, 1>& rhs,
        LINALG::Matrix<5, 1>& rhsElast, LINALG::Matrix<6, 6>& dcedlp,
        LINALG::Matrix<9, 6>& dFpiDdeltaDp, Teuchos::ParameterList& params, const int eleGID);


    virtual void YieldFunction(const double last_ai, const double norm_dLp,
        const LINALG::Matrix<3, 3> ExpEqui, const LINALG::Matrix<3, 3> cetr,
        const LINALG::Matrix<6, 1> str, double* yieldFunc, LINALG::Matrix<3, 3>& devMandelStr,
        LINALG::Matrix<3, 3>& MandelStr);

    virtual void CompElastQuant(const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<3, 3> fpi,
        const LINALG::Matrix<3, 3> MatExp, LINALG::Matrix<3, 3>* cetrial, LINALG::Matrix<6, 1>* Ee);

    virtual void MatrixExponentialSecondDerivativeSym3x3x6(const LINALG::Matrix<3, 3> MatrixIn,
        LINALG::Matrix<3, 3>& exp, LINALG::Matrix<6, 6>& dexp_mat,
        LINALG::Matrix<6, 6>* MatrixExp2ndDerivVoigt);

    virtual void MatrixExponentialSecondDerivativeSym3x3(const LINALG::Matrix<3, 3> MatrixIn,
        LINALG::Matrix<3, 3>& exp, std::vector<LINALG::Matrix<3, 3>>& MatrixExp1stDeriv,
        std::vector<std::vector<LINALG::Matrix<3, 3>>>& MatrixExp2ndDeriv);

    virtual void EvaluatePlast(LINALG::Matrix<6, 9>& dPK2dFpinvIsoprinc,
        const LINALG::Matrix<3, 1>& gamma, const LINALG::Matrix<8, 1>& delta,
        const LINALG::Matrix<3, 3>& id2, const LINALG::Matrix<6, 1>& Cpi,
        const LINALG::Matrix<3, 3>& Fpi, const LINALG::Matrix<3, 3>& CpiC,
        const LINALG::Matrix<9, 1>& CFpi, const LINALG::Matrix<9, 1>& CFpiCei,
        const LINALG::Matrix<6, 1>& ircg, const LINALG::Matrix<3, 3>& FpiCe,
        const LINALG::Matrix<9, 1>& CFpiCe, const LINALG::Matrix<6, 1>& CpiCCpi);

    virtual void EvaluateKinQuantPlast(const int eleGID, const LINALG::Matrix<3, 3>* defgrd,
        const LINALG::Matrix<3, 3>* fpi, LINALG::Matrix<3, 1>& gamma, LINALG::Matrix<8, 1>& delta,
        LINALG::Matrix<3, 3>& id2, LINALG::Matrix<6, 1>& Cpi, LINALG::Matrix<3, 3>& CpiC,
        LINALG::Matrix<9, 1>& CFpi, LINALG::Matrix<9, 1>& CFpiCei, LINALG::Matrix<6, 1>& ircg,
        LINALG::Matrix<3, 3>& FpiCe, LINALG::Matrix<9, 1>& CFpiCe, LINALG::Matrix<6, 1>& CpiCCpi);

    virtual void Dpk2dFpi(const int eleGID, const LINALG::Matrix<3, 3>* defgrd,
        const LINALG::Matrix<3, 3>* fpi, LINALG::Matrix<6, 9>& dPK2dFpinvIsoprinc);

    virtual void DpsiplastDalphaiso(const double norm_dLp, const double last_alphaiso,
        const double isoHardMod, const double initYield, const double infYield,
        const double expIsoHard, double* dpsiplastdalphaiso);

    virtual void Ce2ndDeriv(const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<3, 3> fpi,
        const LINALG::Matrix<3, 3> dLp, LINALG::Matrix<6, 6>* DDceDdLpDdLpVoigt);


    /// Access to material params
    virtual MAT::PAR::PlasticElastHyperVCU* MatParams() const { return params_; }

    /// get dissipation mode
    virtual INPAR::TSI::DissipationMode DisMode() const { return INPAR::TSI::pl_flow; }

    /// inverse plastic deformation gradient for each Gauss point at current state
    std::vector<LINALG::Matrix<3, 3>> plastic_defgrd_inverse_;

    /// my material parameters
    MAT::PAR::PlasticElastHyperVCU* params_;
  };

}  // namespace MAT

#endif  // PLASTIC_VARCONSTUPDATE_H
