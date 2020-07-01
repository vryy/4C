/*-----------------------------------------------------------*/
/*! \file
\brief
This file contains the hyperelastic toolbox. It allows summing up several summands
of several types (isotropic or anisotropic, splitted or not) to build a hyperelastic
strain energy function.

The input line should read
MAT 0   MAT_ElastHyper   NUMMAT 2 MATIDS 1 2 DENS 0

\level 1


*/
/*-----------------------------------------------------------*/
#ifndef ELASTHYPER_H
#define ELASTHYPER_H



#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "elasthyper_service.H"
#include "../drt_mat/anisotropy.H"

// forward declaration due to avoid header definition
namespace MAT
{
  namespace ELASTIC
  {
    class Summand;
  }

  // forward declaration
  class ElastHyper;

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// Collection of hyperelastic materials
    ///
    /// Storage map of hyperelastic summands.
    ///
    /// \author rausch,tk,bborn
    /// \date 05/09

    class ElastHyper : public Parameter
    {
      friend class MAT::ElastHyper;

     public:
      /// standard constructor
      ///
      /// This constructor recursively calls the constructors of the
      /// parameter sets of the hyperelastic summands.
      explicit ElastHyper(const Teuchos::RCP<MAT::PAR::Material>& matdata);

      /// destructor
      ~ElastHyper() override = default;

      /// @name material parameters
      //@{

      //       /// provide access to material/summand by its ID
      //       Teuchos::RCP<const MAT::ELASTIC::Summand> MaterialById(
      //         const int id  ///< ID to look for in collection of summands
      //         ) const;

      /// length of material list
      const int nummat_;

      /// the list of material IDs
      const std::vector<int>* matids_;

      /// material mass density
      const double density_;

      /// 1.0 if polyconvexity of system is checked (0 = no = default)
      const int polyconvex_;

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override;

      /// for stat inverse analysis: add paramters for optimization
      void OptParams(std::map<std::string, int>* pnames) override;

      /// for stat inverse analysis: set pointer to elasthyper
      void SetMaterialPtrSIA(MAT::ElastHyper* in) { statiaelasthyper_ = in; }

     protected:
      /// pointer to elasthyper for stat inverse analysis
      MAT::ElastHyper* statiaelasthyper_;

      //@}

    };  // class ElastHyper

  }  // namespace PAR

  class ElastHyperType : public DRT::ParObjectType
  {
   public:
    std::string Name() const override { return "ElastHyperType"; }

    static ElastHyperType& Instance() { return instance_; };

    DRT::ParObject* Create(const std::vector<char>& data) override;

   private:
    static ElastHyperType instance_;
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

  // forward declaration of unit test
  class ElastHyper_TestSuite;
  class ElastHyperService_TestSuite;

  class ElastHyper : public So3Material
  {
   public:
    friend class MAT::PAR::ElastHyper;

    /// construct empty material object
    ElastHyper();

    /// construct the material object given material parameters
    explicit ElastHyper(MAT::PAR::ElastHyper* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of drt_parobject.H (this file) and should return it in this method.
    int UniqueParObjectId() const override
    {
      return ElastHyperType::Instance().UniqueParObjectId();
    }

    /// \brief Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by UniqueParObjectId() which will then
    /// identify the exact class on the receiving processor.
    ///
    /// \param data (in/out): char vector to store class information
    void Pack(DRT::PackBuffer& data) const override;

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
    INPAR::MAT::MaterialType MaterialType() const override { return INPAR::MAT::m_elasthyper; }

    /// check if element kinematics and material kinematics are compatible
    void ValidKinematics(INPAR::STR::KinemType kinem) override
    {
      if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const override { return Teuchos::rcp(new ElastHyper(*this)); }

    /// number of materials
    virtual int NumMat() const { return params_->nummat_; }

    /// deliver material ID of index i'th potential summand in collection
    virtual int MatID(unsigned index  ///< index i
        ) const;

    /// provide access to material by its ID
    //     virtual Teuchos::RCP<const MAT::ELASTIC::Summand> MaterialById(const int id) const {
    //     return params_->MaterialById(id); }

    /// material mass density
    double Density() const override { return params_->density_; }

    /// a shear modulus equivalent
    virtual double ShearMod() const;

    /// a young's modulus equivalent
    virtual double GetYoung();

    /// evaluate strain energy function
    void StrainEnergy(const LINALG::Matrix<6, 1>& glstrain,  ///< Green-Lagrange strain
        double& psi,                                         ///< Strain energy function
        int gp,                                              ///< Gauss point
        int eleGID                                           ///< Element GID
        ) override;

    /// hyperelastic stress response plus elasticity tensor
    void Evaluate(const LINALG::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const LINALG::Matrix<6, 1>* glstrain,          ///< Green-Lagrange strain
        Teuchos::ParameterList& params,                ///< Container for additional information
        LINALG::Matrix<6, 1>* stress,                  ///< 2nd Piola-Kirchhoff stresses
        LINALG::Matrix<6, 6>* cmat,                    ///< Constitutive matrix
        int gp,                                        ///< Gauss point
        int eleGID                                     ///< Element GID
        ) override;

    /// evaluate for GEMM
    void EvaluateGEMM(LINALG::Matrix<6, 1>* stress,  ///< 2nd Piola-Kirchhoff stresses
        LINALG::Matrix<6, 6>* cmat,                  ///< Constitutive matrix
        double* density,                             ///< Density
        LINALG::Matrix<6, 1>* glstrain_m,            ///< Green-Lagrange strain at t_{n+1/2}
        LINALG::Matrix<6, 1>* glstrain_new,          ///< Green-Lagrange strain at t_{n+1}
        LINALG::Matrix<6, 1>* glstrain_old,          ///< Green-Lagrange strain at t_{n}
        LINALG::Matrix<3, 3>* rcg_new,               ///< Right Cauchy-Green tensor at t_{n+1}
        LINALG::Matrix<3, 3>* rcg_old,               ///< Right Cauchy-Green tensor at t_{n}
        int gp,                                      ///< Gauss point
        int eleGID) override;                        ///< Element GID

    /**
     * \brief Evaluate the Cauchy stress and linearizations with given deformation gradient.
     *
     * Cauchy stress is evaluated within this function call. If requested, the required
     * linearizations are calculated. A potential thermal dependency is handled if the temperature
     * is handed in.
     *
     * \param[in]     defgrd    deformation gradient (\f[\bold{F}\f])
     * \param[in]     n         vector n (\f[\bold{n}\f])
     * \param[in]     t         vector t (\f[\bold{t}\f])
     * \param[out]    snt       cauchy stress tensor contracted using the vectors n and t (\f[
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t} \f])
     * \param[out]    DsntDn    derivative of snt w.r.t. vector n (\f[ \frac{\mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{n}}
                                \f])
     * \param[out]    DsntDt    derivative of snt w.r.t. vector t (\f[ \frac{\mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{t}}
                                \f])
     * \param[out]    DsntDF    derivative of snt w.r.t. deformation gradient (\f[\frac{\mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d}\bold{F}}\f])
     * \param[out]    D2sntDF2  second derivative of snt w.r.t. deformation gradient (\f[ \frac{
                                \mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}
                                {\mathrm{d} \bold{F}^2} \f])
     * \param[out]    D2sntDFDn second derivative of snt w.r.t. deformation gradient and vector n
                                (\f[ \frac{ \mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot
                                \bold{t}} {\mathrm{d} \bold{F} \mathrm{d} \bold{n} } \f])
     * \param[out]    D2sntDFDt second derivative of snt w.r.t. deformation gradient and vector n
                                (\f[ \frac{ \mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot
                                \bold{t}} {\mathrm{d} \bold{F} \mathrm{d} \bold{n} } \f])
     * \oaram[in]     gp        Gauss point
     * \param[in]     eleGID    global ID of element
     * \param[in]     temp      temperature (if thermal dependency shall be handeled)
     * \param[out]    DsntDT    derivative of snt w.r.t. temperature (\f[ \frac{ \mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t} }{\mathrm{d} T} \f])
     * \param[out]    D2sntDFDT second derivative of snt w.r.t. deformation gradient and temperature
                                (\f[ \frac{ \mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}
                                } {\mathrm{d} \bold{F} \mathrm{d} T } \f])
     */
    void EvaluateCauchy(const LINALG::Matrix<3, 3>& defgrd, const LINALG::Matrix<3, 1>& n,
        const LINALG::Matrix<3, 1>& t, double& snt, LINALG::Matrix<3, 1>* DsntDn,
        LINALG::Matrix<3, 1>* DsntDt, LINALG::Matrix<9, 1>* DsntDF, LINALG::Matrix<9, 9>* D2sntDF2,
        LINALG::Matrix<9, 3>* D2sntDFDn, LINALG::Matrix<9, 3>* D2sntDFDt, const int gp,
        const int eleGID, const double* temp = nullptr, double* DsntDT = nullptr,
        LINALG::Matrix<9, 1>* D2sntDFDT = nullptr) override;

    /**
     * \brief The derivatives of the strain energy function w.r.t. the principle invariants are
     * filled up to the third derivative, can be overloaded if temperature is given.
     *
     * \param prinv[in]       principle invariants of the tensor
     * \param gp[in]          Gauss point
     * \param eleGID[in]      global ID of the element
     * \param dPI[in,out]     first derivative of strain energy function w.r.t. principle invariants
     * \param ddPII[in,out]   second derivative of strain energy function w.r.t. principle
                              invariants
     * \param dddPIII[in,out] third derivative of strain energy function w.r.t. principle invariants
     * \param temp[in]        temperature
     */
    virtual void EvaluateCauchyDerivs(const LINALG::Matrix<3, 1>& prinv, int gp, int eleGID,
        LINALG::Matrix<3, 1>& dPI, LINALG::Matrix<6, 1>& ddPII, LINALG::Matrix<10, 1>& dddPIII,
        const double* temp);

    /**
     * \brief Evaluate the thermal dependency of the linearizations of the cauchy stress, to be
     *        overloaded
     *
     * \param prinv[in]      principle invariants of the tensor
     * \param ndt[in]        dot product of vectors n and t (\f[\bold{t} \cdot \bold{t}\f])
     * \param bdndt[in]      left Cauchy-Green tensor contracted using the vector n and t (\f[
                             \bold{b} \cdot \bold{n} \cdot \bold{t} \f])
     * \param ibdndt[in]     inverse of left Cauchy-Green tensor contracted using the vector n and t
                             (\f[ \bold{b}^{-1} \cdot \bold{n} \cdot \bold{t} \f])
     * \param temp[in]       temperature
     * \param DsntDT[out]    derivative of cauchy stress contracted with vectors n and t w.r.t.
                             temperature (\f[\frac{\mathrm{d} \bold{\sigma} \cdot \bold{n} \cdot
                             \bold{t} }{ \mathrm{d} T }\f])
     * \param iFTV[in]       inverse transposed of the deformation gradient
     * \param DbdndtDFV[in]  derivative of bdndt w.r.t. deformation gradient (\f[\frac{ \mathrm{d}
                             \bold{b} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{F} }\f])
     * \param DibdndtDFV[in] derivative of ibdndt w.r.t. deformation gradient (\f[\frac{ \mathrm{d}
                             \bold{b}^{-1} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{F} }\f])
     * \param DI1DF[in]      derivative of first invariant w.r.t. deformation gradient
     * \param DI2DF[in]      derivative of second invariant w.r.t. deformation gradient
     * \param DI3DF[in]      derivative of third invariant w.r.t. deformation gradient
     * \param D2sntDFDT[out] second derivative of cauchy stress contracted with vectors n and t
                             w.r.t. temperature and deformation gradient (\f[ \frac{\mathrm{d}^2
                             \bold{\sigma} \cdot \bold{n} \cdot \bold{t}} {\mathrm{d} \bold{F}
                             \mathrm{d} T } \f])
     */
    virtual void EvaluateCauchyTempDeriv(const LINALG::Matrix<3, 1>& prinv, const double ndt,
        const double bdndt, const double ibdndt, const double* temp, double* DsntDT,
        const LINALG::Matrix<9, 1>& iFTV, const LINALG::Matrix<9, 1>& DbdndtDFV,
        const LINALG::Matrix<9, 1>& DibdndtDFV, const LINALG::Matrix<9, 1>& DI1DF,
        const LINALG::Matrix<9, 1>& DI2DF, const LINALG::Matrix<9, 1>& DI3DF,
        LINALG::Matrix<9, 1>* D2sntDFDT)
    {
    }

    /// setup
    void Setup(int numgp, DRT::INPUT::LineDefinition* linedef) override;

    /*!
     * Post setup routine that is executed after the input. This class will forward the call to all
     * summands
     *
     * @param params Container for additional information
     * @param eleGID Global element id
     */
    void PostSetup(Teuchos::ParameterList& params, int eleGID) override;

    /// update
    void Update() override;

    /// setup patient-specific AAA stuff
    virtual void SetupAAA(Teuchos::ParameterList& params, int eleGID);

    /// return if anisotropic not splitted formulation
    virtual bool AnisotropicPrincipal() const { return summandProperties_.anisoprinc; }

    /// return if anisotropic and splitted formulation
    virtual bool AnisotropicModified() const { return summandProperties_.anisomod; }

    /// get fiber vectors
    virtual void GetFiberVecs(std::vector<LINALG::Matrix<3, 1>>& fibervecs);

    /// evaluate fiber directions from locsys and alignment angle, pull back
    virtual void EvaluateFiberVecs(double newgamma,  ///< new angle
        const LINALG::Matrix<3, 3>& locsys,          ///< local coordinate system
        const LINALG::Matrix<3, 3>& defgrd           ///< deformation gradient
    );

    /// Return potential summand pointer for the given material type
    Teuchos::RCP<const MAT::ELASTIC::Summand> GetPotSummandPtr(
        const INPAR::MAT::MaterialType& materialtype) const;

    /// Return quick accessible material parameter data
    MAT::PAR::Parameter* Parameter() const override { return params_; }

    /// Return names of visualization data
    void VisNames(std::map<std::string, int>& names) override;

    /// Return visualization data
    bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID) override;

    /// Return whether the material requires the deformation gradient for its evaluation
    bool NeedsDefgrd() override
    {
      // only the polyconvexity check needs the deformation gradient. Regular materials don't need
      // it
      return params_->polyconvex_ != 0;
    };

   protected:
    /// Sat Inverse Analysis: add params of summands
    virtual void ElastOptParams(std::map<std::string, int>* pnames);

    /// @name Flags to specify the elastic formulations (initialize with false)
    //@{
    /// Class indicating the formulation used by the summands
    SummandProperties summandProperties_;
    //@}

    /// my material parameters
    MAT::PAR::ElastHyper* params_;

    /// map to materials/potential summands
    std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>> potsum_;

    /// Holder of anisotropy
    MAT::Anisotropy anisotropy_;
  };

}  // namespace MAT

#endif
