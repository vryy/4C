/*-----------------------------------------------------------*/
/*!
\brief
This file contains the hyperelastic toolbox. It allows summing up several summands
of several types (isotropic or anisotropic, splitted or not) to build a hyperelastic
strain energy function.

The input line should read
MAT 0   MAT_ElastHyper   NUMMAT 2 MATIDS 1 2 DENS 0

\level 1

\maintainer Fabian Braeu

*/
/*-----------------------------------------------------------*/
#ifndef ELASTHYPER_H
#define ELASTHYPER_H



#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"

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
      ElastHyper(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~ElastHyper() { ; }

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
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// for stat inverse analysis: add paramters for optimization
      virtual void OptParams(std::map<std::string, int>* pnames);

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
    std::string Name() const { return "ElastHyperType"; }

    static ElastHyperType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

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
    virtual int UniqueParObjectId() const { return ElastHyperType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_elasthyper; }

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem)
    {
      if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new ElastHyper(*this)); }

    /// number of materials
    virtual int NumMat() const { return params_->nummat_; }

    /// deliver material ID of index i'th potential summand in collection
    virtual int MatID(const unsigned index  ///< index i
        ) const;

    /// provide access to material by its ID
    //     virtual Teuchos::RCP<const MAT::ELASTIC::Summand> MaterialById(const int id) const {
    //     return params_->MaterialById(id); }

    /// material mass density
    virtual double Density() const { return params_->density_; }

    /// a shear modulus equivalent
    virtual double ShearMod() const;

    /// a young's modulus equivalent
    virtual double GetYoung();

    /// evaluate strain energy function
    virtual void StrainEnergy(const LINALG::Matrix<6, 1>& glstrain,  ///< Green-Lagrange strain
        double& psi,                                                 ///< Strain energy function
        const int eleGID                                             ///< Element GID
    );

    /// hyperelastic stress response plus elasticity tensor
    virtual void Evaluate(const LINALG::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const LINALG::Matrix<6, 1>* glstrain,                  ///< Green-Lagrange strain
        Teuchos::ParameterList& params,  ///< Container for additional information
        LINALG::Matrix<6, 1>* stress,    ///< 2nd Piola-Kirchhoff stresses
        LINALG::Matrix<6, 6>* cmat,      ///< Constitutive matrix
        const int eleGID                 ///< Element GID
    );

    /*!
     * \brief Evaluate the stress response and the elasticity tensor
     *
     * This static method evaluates the 2nd Piola-Kirchhoff stress tensor and its linearization with
     * respect to the Green-Lagrange strain.
     *
     * @param defgrd      (in)      : deformation gradient
     * @param glstrain    (in)      : Green lagrange strain in strain-like Voigt notation
     * @param params      (in/out)  : Container for additional information
     * @param stress      (out)     : 2nd Piola Kirchhoff stress
     * @param cmat        (out)     : Elasticity tensor
     * @param eleGID      (in)      : Element id
     * @param potsum      (in)      : Summands of the Free-energy function
     * @param isoprinc    (in)      : Flag, whether at least one summand is formulated in the
     * isotropic principal invariants
     * @param isomod      (in)      : Flag, whether at least one summand is formulated in the
     * isotropic modified principal invariants
     * @param anisoprinc  (in)      : Flag, whether at least one summand is formulated in the
     * anisotropic principal invariants
     * @param anisomod    (in)      : Flag, whether at least one summand is formulated in the
     * anisotropic modified principal invariants
     * @param checkpolyconvexity (in) : Flag, whether to check the polyconvexity
     */
    static void EvaluateStressCmatStatic(
        const LINALG::Matrix<3, 3>* defgrd,    ///< Deformation gradient
        const LINALG::Matrix<6, 1>* glstrain,  ///< Green-Lagrange strain
        Teuchos::ParameterList& params,        ///< Holder for additional parameters
        LINALG::Matrix<6, 1>* stress,          ///< 2nd Piola-Kirchhoff stresses
        LINALG::Matrix<6, 6>* cmat,            ///< Constitutive matrix
        int eleGID,                            ///< Element GID
        const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>&
            potsum,     ///< Map to material/potential summands
        bool isoprinc,  ///< flag, whether material uses principal invariants (isotropic case)
        bool
            isomod,  ///< flag, whether material uses modified principal invariants (isotropic case)
        bool anisoprinc,  ///< flag, whether material uses principal invariants (anisotropic case)
        bool anisomod,  ///< flag, whether material uses modified principal invariants (anisotropic
                        ///< case)
        bool checkpolyconvexity  ///< flag, whether polyconvexity should be checked
    );

    /// evaluate for GEMM
    virtual void EvaluateGEMM(LINALG::Matrix<6, 1>* stress,  ///< 2nd Piola-Kirchhoff stresses
        LINALG::Matrix<6, 6>* cmat,                          ///< Constitutive matrix
        double* density,                                     ///< Density
        LINALG::Matrix<6, 1>* glstrain_m,                    ///< Green-Lagrange strain at t_{n+1/2}
        LINALG::Matrix<6, 1>* glstrain_new,                  ///< Green-Lagrange strain at t_{n+1}
        LINALG::Matrix<6, 1>* glstrain_old,                  ///< Green-Lagrange strain at t_{n}
        LINALG::Matrix<3, 3>* rcg_new,  ///< Right Cauchy-Green tensor at t_{n+1}
        LINALG::Matrix<3, 3>* rcg_old,  ///< Right Cauchy-Green tensor at t_{n}
        const int eleGID);              ///< Element GID

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
        LINALG::Matrix<9, 3>* D2sntDFDn, LINALG::Matrix<9, 3>* D2sntDFDt, const int eleGID,
        const double* temp = NULL, double* DsntDT = NULL,
        LINALG::Matrix<9, 1>* D2sntDFDT = NULL) override;

    /**
     * \brief The derivatives of the strain energy function w.r.t. the principle invariants are
     * filled up to the third derivative, can be overloaded if temperature is given.
     *
     * \param prinv[in]       principle invariants of the tensor
     * \param eleGID[in]      global ID of the element
     * \param dPI[in,out]     first derivative of strain energy function w.r.t. principle invariants
     * \param ddPII[in,out]   second derivative of strain energy function w.r.t. principle
                              invariants
     * \param dddPIII[in,out] third derivative of strain energy function w.r.t. principle invariants
     * \param temp[in]        temperature
     */
    virtual void EvaluateCauchyDerivs(const LINALG::Matrix<3, 1>& prinv, const int eleGID,
        LINALG::Matrix<3, 1>& dPI, LINALG::Matrix<6, 1>& ddPII, LINALG::Matrix<10, 1>& dddPIII,
        const double* temp = NULL);

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
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    /// update
    virtual void Update();

    /// setup patient-specific AAA stuff
    virtual void SetupAAA(Teuchos::ParameterList& params, const int eleGID);

    /// return if isotropic not splitted formulation
    virtual bool IsotropicPrincipal() const { return isoprinc_; }

    /// return if isotropic and splitted formulation
    virtual bool IsotropicModified() const { return isomod_; }

    /// return if anisotropic not splitted formulation
    virtual bool AnisotropicPrincipal() const { return anisoprinc_; }

    /// return if anisotropic and splitted formulation
    virtual bool AnisotropicModified() const { return anisomod_; }

    /// get fiber vectors
    virtual void GetFiberVecs(std::vector<LINALG::Matrix<3, 1>>& fibervecs);

    /// evaluate fiber directions from locsys and alignment angle, pull back
    virtual void EvaluateFiberVecs(const double newgamma,  ///< new angle
        const LINALG::Matrix<3, 3>& locsys,                ///< local coordinate system
        const LINALG::Matrix<3, 3>& defgrd                 ///< deformation gradient
    );

    /// Return potential summand pointer for the given material type
    Teuchos::RCP<const MAT::ELASTIC::Summand> GetPotSummandPtr(
        const INPAR::MAT::MaterialType& materialtype) const;

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string, int>& names);

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

    /// Return whether the material requires the deformation gradient for its evaluation
    virtual bool NeedsDefgrd()
    {
      if (params_->polyconvex_) return true;
      /* usual materials are based on GL-strains */
      else
        return false;
    };

   protected:
    /*!
     * \brief Calculates the kinematic quantities and tensors
     *
     * @param glstrain (in) : Green-Lagrange strain tensor
     * @param id2 (out) : Cartesian identity 2-tensor in stress/strain-like Voigt notation
     * @param scg (out) : Right Cauchy-Green strain tensor in stress-like Voigt notation
     * @param rcg (out) : Right Cauchy-Green strain tensor in strain-like Voigt notation
     * @param icg (out) : Inverse of the Right Cauchy Green strain tensor in stress-like Voigt
     * notation
     * @param id4 (out) : Cartesian identity 4-tensor in Voigt notation (rows are stress-like,
     * columns are strain-like) This tensor is 'mixed co- and contra-variant'.
     * @param id4sharp (out) : Cartesian identity 4-tensor in Voigt notation (rows and columns are
     * stress-like). This tensor is 'contra-variant'.
     * @param prinv (out) : Principal invariants of the Right Cauchy-Green strain tensor
     */
    static void EvaluateKinQuant(const LINALG::Matrix<6, 1>& glstrain, LINALG::Matrix<6, 1>& id2,
        LINALG::Matrix<6, 1>& scg, LINALG::Matrix<6, 1>& rcg, LINALG::Matrix<6, 1>& icg,
        LINALG::Matrix<6, 6>& id4, LINALG::Matrix<6, 6>& id4sharp, LINALG::Matrix<3, 1>& prinv);


    /*!
     * \brief Calculates the derivatives of the hyperelastic laws with respect to the (modified)
     * principal invariants
     *
     * @param prinv (in) : Principal invariants of the Right Cauchy-Green tensor
     * @param dPI (out) : First derivative of the Free-energy function with respect to the principal
     * invariants
     * @param ddPII (out) : Second derivative of the Free-energy function with respect to the
     * principal invariants
     * @param eleGID (in) : Element id
     * @param potsum (in) : List of summands
     * @param isoprinc (in) : Flag, whether at least one summand is formulated in the isotropic
     * principal invariants
     * @param isomod (in) : Flag, whether at least one summand is formulated in the isotropic
     * modified principal invariants
     */
    static void EvaluateInvariantDerivatives(const LINALG::Matrix<3, 1>& prinv,
        LINALG::Matrix<3, 1>& dPI, LINALG::Matrix<6, 1>& ddPII, int eleGID,
        const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum, bool isoprinc, bool isomod);

    /*!
     * \brief Converts the derivatives with respect to the modified principal invariants in
     * derivatives with respect to the principal invariants.
     *
     * The used conversions are:
     * \f[
     *   \overline{I}_{\boldsymbol{C}} = J^{-2/3} I_{\boldsymbol{C}},
     * \f]
     * \f[
     *   \overline{II}_{\boldsymbol{C}} = J^{-4/3} II_{\boldsymbol{C}},
     * \f]
     * \f[
     *   J = \sqrt{III_{\boldsymbol{C}}}
     * \f]
     *
     * @param prinv (in) : Principal invariants of the Right Cauchy-Green strain tensor
     * @param dPmodI (in) : First derivatives of the Free-energy function with respect to the
     * modified principal invariants
     * @param ddPmodII (in) : Second derivatives of the Free-energy function with respect to the
     * modified principal invariants
     * @param dPI (out) : First derivatives of the Free-energy function with respect to the
     * principal invariants
     * @param ddPII (out) : Second derivatives of the Free-energy function with respect to the
     * principal invariants
     */
    static void ConvertModToPrinc(const LINALG::Matrix<3, 1>& prinv,
        const LINALG::Matrix<3, 1>& dPmodI, const LINALG::Matrix<6, 1>& ddPmodII,
        LINALG::Matrix<3, 1>& dPI, LINALG::Matrix<6, 1>& ddPII);

    /*!
     * \brief Calculate coefficients gamma and delta from partial derivatives w.r.t. invariants.
     *
     * The coefficients \f$\gamma_i\f$ and \f$\delta_j\f$ are based
     *  on the summand \f$\Psi(I_{\boldsymbol{C}},II_{\boldsymbol{C}},III_{\boldsymbol{C}})\f$
     *  in which the principal invariants of the right Cauchy-Green tensor \f$\boldsymbol{C}\f$
     *  are the arguments, cf. Holzapfel [1],
     * \f[
     *   I_{\boldsymbol{C}} = C_{AA},
     *   \quad
     *   II_{\boldsymbol{C}} = 1/2 \big( \mathrm{trace}^2(\boldsymbol{C}) -
     *   \mathrm{trace}(\boldsymbol{C}^2) \big), \quad III_{\boldsymbol{C}} = \det(\boldsymbol{C})
     * \f]
     *
     *  cf. Holzapfel [1], p. 216 (6.32) and p. 248
     * \f[
     *   \bold{S} = \gamma_1 \ \bold{Id} + \gamma_2 \ \bold{C} + \gamma_3 \ \bold{C}^{-1}
     * \f]
     * \f[
     *  \gamma_1 = 2\left( \frac{\partial \Psi}{\partial I_{\boldsymbol{C}}}
     *           + I_{\boldsymbol{C}}\frac{\partial \Psi}{\partial II_{\boldsymbol{C}}} \right);
     * \f]
     * \f[
     *  \gamma_2 = -2\frac{\partial \Psi}{\partial II_{\boldsymbol{C}}};
     * \f]
     * \f[
     *  \gamma_3 = 2III_{\boldsymbol{C}} \frac{\partial \Psi}{\partial III_{\boldsymbol{C}}};
     * \f]
     *
     *  material constitutive tensor coefficients
     *  cf. Holzapfel [1], p. 261
     * \f[
     *   \mathbb{C} = \delta_1 \left( \bold{Id} \otimes \bold{Id} \right) + \delta_2 \left(
     *   \bold{Id} \otimes \bold{C} + \bold{C} \otimes \bold{Id} \right)
     *   + \delta_3 \left( \bold{Id} \otimes \bold{C}^{-1} + \bold{C}^{-1} \otimes \bold{Id} \right)
     *   + \delta_4 \left( \bold{C} \otimes \bold{C} \right)
     *   + \delta_5 \left( \bold{C} \otimes \bold{C}^{-1} + \bold{C}^{-1} \otimes \bold{C} \right) +
     *   \delta_6 \left( \bold{C}^{-1} \otimes \bold{C}^{-1} \right)
     *   + \delta_7 \mathbb{P}
     * \f]
     * \f[
     *  \delta_1 = 4\left( \frac{\partial^2 \Psi}{\partial I_{\boldsymbol{C}}^2}
     *           + 2 Ic \frac{\partial^2 \Psi}{\partial I_{\boldsymbol{C}} \partial
     *           II_{\boldsymbol{C}}}
     *           + \frac{\partial \Psi}{\partial II_{\boldsymbol{C}}}
     *           + I_{\boldsymbol{C}}^2 \frac{\partial^2 \Psi}{\partial II_{\boldsymbol{C}}^2}
     *           \right)
     * \f]
     * \f[
     *  \delta_2 = -4\left( \frac{\partial^2 \Psi}{\partial I_{\boldsymbol{C}} \partial
     *  II_{\boldsymbol{C}}}
     *           + I_{\boldsymbol{C}}\frac{\partial^2 \Psi}{\partial II_{\boldsymbol{C}}^2} \right)
     * \f]
     * \f[
     *  \delta_3 = 4\left( III_{\boldsymbol{C}} \frac{\partial^2 \Psi}{\partial I_{\boldsymbol{C}}
     *  \partial III_{\boldsymbol{C}}}
     *           + I_{\boldsymbol{C}} III_{\boldsymbol{C}} \frac{\partial^2 \Psi}{\partial IIc
     *           \partial III_{\boldsymbol{C}}} \right)
     * \f]
     * \f[
     *  \delta_4 = 4\frac{\partial^2 \Psi}{\partial II_{\boldsymbol{C}}^2}
     * \f]
     * \f[
     *  \delta_5 = -4 III_{\boldsymbol{C}}\frac{\partial^2 \Psi}{\partial II_{\boldsymbol{C}}
     *  \partial III_{\boldsymbol{C}}}
     * \f]
     * \f[
     *  \delta_6 = 4\left( III_{\boldsymbol{C}} \frac{\partial \Psi}{\partial III_{\boldsymbol{C}}}
     *           + III_{\boldsymbol{C}}^2 \frac{\partial^2 \Psi}{\partial III_{\boldsymbol{C}}^2}
     *           \right)
     * \f]
     * \f[
     *  \delta_7 = -4 III_{\boldsymbol{C}} \frac{\partial \Psi}{\partial III_{\boldsymbol{C}}}
     * \f]
     * \f[
     *  \delta_8 = -4 \frac{\partial \Psi}{\partial II_{\boldsymbol{C}}}
     * \f]
     *
     * @param gamma (out) : Coefficients gamma, cf. Holzapfel [1]
     * @param delta (out) : Coefficients delta, cf. Holzapfel [1]
     * @param prinv (in) : Principal invariants of the Right Cauchy-Green strain tensor
     * @param dPI (in) : First derivatives of the Free energy function w.r.t the principal
     * invariants of the Right Cauchy-Green strain tensor
     * @param ddPII (in) : Second derivatives of the Free energy function w.r.t the principal
     * invariants of the Right Cauchy-Green strain tensor
     */
    static void CalculateGammaDelta(LINALG::Matrix<3, 1>& gamma, LINALG::Matrix<8, 1>& delta,
        const LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<3, 1>& dPI,
        const LINALG::Matrix<6, 1>& ddPII);

    /// calculates the isotropic stress and elasticity tensor

    /*!
     * \brief Evaluation of the isotropic contribution to the stress tensor and the elasticity
     * tensor
     *
     * @param stress (out) : 2nd Piola-Kirchhoff stress tensor in stress-like Voigt notation
     * @param cmat (out) : Elasticity tensor
     * @param scg (in) : Right Cauchy-Green strain tensor in stress-like Voigt notation
     * @param id2 (in) : Cartesian identity 2-tensor in stress/strain-like Voigt notation
     * @param icg (in) : Inverse of the Right Cauchy Green strain tensor in stress-like Voigt
     * notation
     * @param id4sharp (in) : Cartesian identity 4-tensor in Voigt notation (rows and columns are
     * stress-like). This tensor is 'contra-variant'.
     * @param prinv (in) : Principal invariants of the Right Cauchy-Green strain tensor
     * @param dPI (in) : First derivative of the Free-energy function with respect to the principal
     * invariants
     * @param ddPII (in) : Second derivative of the Free-energy function with respect to the
     * principal invariants
     */
    static void EvaluateIsotropicStressCmat(LINALG::Matrix<6, 1>& stress,
        LINALG::Matrix<6, 6>& cmat, const LINALG::Matrix<6, 1>& scg,
        const LINALG::Matrix<6, 1>& id2, const LINALG::Matrix<6, 1>& icg,
        const LINALG::Matrix<6, 6>& id4sharp, const LINALG::Matrix<3, 1>& prinv,
        const LINALG::Matrix<3, 1>& dPI, const LINALG::Matrix<6, 1>& ddPII);

    /*!
     * \brief Evaluation of the anisotropic stress and elasticity tensor for coupled configuration
     *
     * @param stressanisoprinc (out) : 2nd Piola-Kirchhoff stress tensor in stress-like Voigt
     * notation
     * @param cmatanisoprinc (out) : Elasticity tensor
     * @param rcg (in) : Right Cauchy-Green tensor in strain-like Voigt notation
     * @param params (in/out) : Container for additional information
     * @param eleGID (in) : Element id
     * @param potsum (in) : List of summands in the Free-energy function
     */
    static void EvaluateAnisotropicPrinc(LINALG::Matrix<6, 1>& stressanisoprinc,
        LINALG::Matrix<6, 6>& cmatanisoprinc, const LINALG::Matrix<6, 1>& rcg,
        Teuchos::ParameterList& params, int eleGID,
        const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum);

    /*!
     * \brief Evaluates the anisotropic stress and elasticity tensor for decoupled configuration
     * @param stressanisoprinc (out) : 2nd Piola-Kirchhoff stress tensor in stress-like Voigt
     * notation
     * @param cmatanisoprinc (out) : Elasticity tensor
     * @param rcg (in) : Right Cauchy-Green tensor in strain-like Voigt notation
     * @param icg (in) : Inverse of the Right Cauchy Green strain tensor in stress-like Voigt
     * notation
     * @param prinv (in) : Principal invariants of the Right Cauchy-Green strain tensor
     * @param eleGID (in) : Element id
     * @param potsum (in) : List of summands in the Free-energy function
     */
    static void EvaluateAnisotropicMod(LINALG::Matrix<6, 1>& stressanisoprinc,
        LINALG::Matrix<6, 6>& cmatanisoprinc, const LINALG::Matrix<6, 1>& rcg,
        const LINALG::Matrix<6, 1>& icg, const LINALG::Matrix<3, 1>& prinv, int eleGID,
        const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum);

    /*!
     * \brief Evaluates the principal invariants of the Right Cauchy-Green tensor
     *
     * @param prinv (out) : Principal invariants of the Right Cauchy-Green strain tensor
     * @param rcg (in) : Right Cauchy-Green tensor in strain-like Voigt notation
     */
    static void InvariantsPrincipal(LINALG::Matrix<3, 1>& prinv,  ///< principal invariants
        const LINALG::Matrix<6, 1>&
            rcg  ///< symmetric Cartesian 2-tensor in strain-like 6-Voigt notation
    );

    /*!
     * \brief Converts the principal invariants to the modified principal invariants.
     *
     * @param modinv (out) : Modified principal invariants of the Right Cauchy-Green strain tensor
     * @param prinv (in) : Principal invariants of the Right Cauchy-Green strain tensor
     */
    static void InvariantsModified(LINALG::Matrix<3, 1>& modinv,  ///< modified principal invariants
        const LINALG::Matrix<3, 1>& prinv                         ///< principal invariants
    );

    /*!
     * \brief Evaluates principal stretches \f$\lambda_\alpha\f$ and material principal directions
     *
     * @param prstr (out) : Principal stretches
     * @param prdir (out) : Principal directions
     * @param rcg (in) : Right Cauchy-Green strain tensor in strain-like Voigt notation
     */
    static void StretchesPrincipal(LINALG::Matrix<3, 1>& prstr,  ///< principal stretches
        LINALG::Matrix<3, 3>& prdir,  ///< principal directions, stored col-by-col
        const LINALG::Matrix<6, 1>&
            rcg);  ///< symmetric Cartesian 2-tensor in strain-like 6-Voigt notation


    /*!
     * \brief Evaluates the modified principal stretches \f$\bar{\lambda}_\alpha\f$
     *
     * @param modstr (out) : Modified principal stretches
     * @param prstr (in) : Principal stretches
     */
    static void StretchesModified(LINALG::Matrix<3, 1>& modstr,  ///< modified principal stretches
        const LINALG::Matrix<3, 1>& prstr                        ///< principal stretches
    );

    /// determine if at least one summand uses principal stretches

    /*!
     * \brief determine if at least one summand uses principal stretches
     * @param potsum (in) : Summands of the Free-energy function
     * @return true if at least one summand uses principal stretches, otherwise false
     */
    static bool HaveCoefficientsStretchesPrincipal(
        const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum);

    /*!
     * \brief determine if at least one summand uses modified principal stretches
     * @param potsum (in) : Summands of the Free-energy function
     * @return true if at least one summand uses modified principal stretches, otherwise false
     */
    static bool HaveCoefficientsStretchesModified(
        const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum);

    /// Determine PK2 stress response and material elasticity tensor
    /// due to energy densities (MAT::ELASTIC::Summand) described
    /// in (modified) principal stretches.
    ///
    /// The stress response is achieved by collecting the coefficients
    /// \f$\gamma_\alpha\f$ and \f$\delta_{\alpha\beta}\f$ due to
    /// MAT::ELASTIC::Summand::AddCoefficientsStretchesPrincipal()
    /// (and/or \f$\bar{\gamma}_\alpha\f$ and \f$\bar{\delta}_{\alpha\beta}\f$
    /// due to MAT::ELASTIC::Summand::AddCoefficientsStretchesModified()).
    /// The collected coefficients build the principal 2nd Piola--Kirchhoff
    /// stresses which are transformed into ordinary Cartesian co-ordinates
    /// applying the principal directions. Similarly, the elasticity
    /// 4-tensor is obtained.
    ///
    /// Please note, unlike suggested by Holzapfel, p 263-264, the modified
    /// coefficients are transformed to unmodified coefficients and than added
    /// onto the respective quantities. The transformation goes along the following
    /// lines. The first derivatives or \f$\gamma_\alpha\f$ coefficients are related
    /// to the modified coefficients \f$\bar{\gamma}_\alpha\f$ via the chain rule
    ///\f[
    ///   \gamma_\alpha
    ///   = \frac{\partial \Psi(\boldsymbol{\lambda})}{\partial\lambda_\alpha}
    ///   = \frac{\partial \Psi(\bar{\boldsymbol{\lambda}})}{\partial \bar{\lambda}_\eta}
    ///   \, \frac{\partial \bar{\lambda}_\eta}{\partial\lambda_\alpha}
    ///   = \bar{\gamma}_\eta \, \frac{\partial \bar{\lambda}_\eta}{\partial\lambda_\alpha}
    ///\f]
    /// utilising Holzapfel Eq (6.142):
    ///\f[
    ///   \frac{\partial\bar{\lambda}_\alpha}{\partial\lambda_\beta}
    ///   = J^{-1/3} \Big( 1_{\alpha\beta} - \frac{1}{3} \lambda_\alpha \lambda_\beta^{-1} \Big)
    ///\f]
    /// in which Kronecker's delta is denoted \f$1_{\alpha\beta}\f$.
    /// Likewise (and once without Holzapfel), the second derivatives \f$\delta_{\alpha\beta}\f$ can
    /// be recovered directly by knowledge of the modified coefficients
    /// \f$\bar{\gamma}_\alpha\f$ and \f$\bar{\delta}_{\alpha\beta}\f$, i.e.
    ///\f[
    ///   \delta_{\alpha\beta}
    ///   = \frac{\partial^2\Psi(\boldsymbol{\lambda})}{\partial\lambda_\alpha
    ///   \partial\lambda_\beta} =
    ///   \frac{\partial^2\Psi(\bar{\boldsymbol{\lambda}})}{\partial\bar{\lambda}_\eta\partial\bar{\lambda}_\epsilon}
    ///     \, \frac{\partial \bar{\lambda}_\eta}{\partial\lambda_\alpha}
    ///     \, \frac{\partial \bar{\lambda}_\epsilon}{\partial\lambda_\beta}
    ///   + \frac{\partial \Psi(\bar{\boldsymbol{\lambda}})}{\partial \bar{\lambda}_\eta}
    ///     \, \frac{\partial^2 \bar{\lambda}_\eta}{\partial\lambda_\alpha\partial\lambda_\beta}
    ///   = \bar{\delta}_{\eta\epsilon}
    ///     \, \frac{\partial \bar{\lambda}_\eta}{\partial\lambda_\alpha}
    ///     \, \frac{\partial \bar{\lambda}_\epsilon}{\partial\lambda_\beta}
    ///   + \bar{\gamma}_\eta
    ///     \, \frac{\partial^2 \bar{\lambda}_\eta}{\partial\lambda_\alpha\partial\lambda_\beta}
    ///\f]
    /// with
    ///\f[
    ///   \frac{\partial^2 \bar{\lambda}_\alpha}{\partial\lambda_\beta\partial\lambda_\eta}
    ///   = -\frac{1}{3} J^{-1/3} \Big(
    ///     1_{\alpha\beta} \lambda_\eta^{-1}
    ///     + 1_{\alpha\eta}\lambda_\beta^{-1}
    ///     - \big( 1_{\beta\eta} + \frac{1}{3} \big) \lambda_\alpha \lambda_\beta^{-1}
    ///     \lambda_\eta^{-1}
    ///   \Big)
    ///\f]
    ///
    /// <h3>References</h3>
    /// See Holzapfel, p 245-246, p 257-259, p 263-264
    ///
    /// \author bborn
    /// \date 08/09

    /*!
     * \brief Determine PK2 stress response and material elasticity tensor
     *  due to energy densities (MAT::ELASTIC::Summand) described
     *  in (modified) principal stretches.
     *  The stress response is achieved by collecting the coefficients
     *  \f$\gamma_\alpha\f$ and \f$\delta_{\alpha\beta}\f$ due to
     *  MAT::ELASTIC::Summand::AddCoefficientsStretchesPrincipal()
     *  (and/or \f$\bar{\gamma}_\alpha\f$ and \f$\bar{\delta}_{\alpha\beta}\f$
     *  due to MAT::ELASTIC::Summand::AddCoefficientsStretchesModified()).
     *  The collected coefficients build the principal 2nd Piola--Kirchhoff
     *  stresses which are transformed into ordinary Cartesian co-ordinates
     *  applying the principal directions. Similarly, the elasticity
     *  4-tensor is obtained.
     *
     *  Please note, unlike suggested by Holzapfel, p 263-264, the modified
     *  coefficients are transformed to unmodified coefficients and than added
     *  onto the respective quantities. The transformation goes along the following
     *  lines. The first derivatives or \f$\gamma_\alpha\f$ coefficients are related
     *  to the modified coefficients \f$\bar{\gamma}_\alpha\f$ via the chain rule
     * \f[
     *    \gamma_\alpha
     *    = \frac{\partial \Psi(\boldsymbol{\lambda})}{\partial\lambda_\alpha}
     *    = \frac{\partial \Psi(\bar{\boldsymbol{\lambda}})}{\partial \bar{\lambda}_\eta}
     *    \, \frac{\partial \bar{\lambda}_\eta}{\partial\lambda_\alpha}
     *    = \bar{\gamma}_\eta \, \frac{\partial \bar{\lambda}_\eta}{\partial\lambda_\alpha}
     * \f]
     *  utilising Holzapfel Eq (6.142):
     * \f[
     *    \frac{\partial\bar{\lambda}_\alpha}{\partial\lambda_\beta}
     *    = J^{-1/3} \Big( 1_{\alpha\beta} - \frac{1}{3} \lambda_\alpha \lambda_\beta^{-1} \Big)
     * \f]
     *  in which Kronecker's delta is denoted \f$1_{\alpha\beta}\f$.
     *  Likewise (and once without Holzapfel), the second derivatives \f$\delta_{\alpha\beta}\f$ can
     *  be recovered directly by knowledge of the modified coefficients
     *  \f$\bar{\gamma}_\alpha\f$ and \f$\bar{\delta}_{\alpha\beta}\f$, i.e.
     * \f[
     *    \delta_{\alpha\beta}
     *    = \frac{\partial^2\Psi(\boldsymbol{\lambda})}{\partial\lambda_\alpha
     *    \partial\lambda_\beta} =
     *    \frac{\partial^2\Psi(\bar{\boldsymbol{\lambda}})}{\partial\bar{\lambda}_\eta\partial\bar{\lambda}_\epsilon}
     *      \, \frac{\partial \bar{\lambda}_\eta}{\partial\lambda_\alpha}
     *      \, \frac{\partial \bar{\lambda}_\epsilon}{\partial\lambda_\beta}
     *    + \frac{\partial \Psi(\bar{\boldsymbol{\lambda}})}{\partial \bar{\lambda}_\eta}
     *      \, \frac{\partial^2 \bar{\lambda}_\eta}{\partial\lambda_\alpha\partial\lambda_\beta}
     *    = \bar{\delta}_{\eta\epsilon}
     *      \, \frac{\partial \bar{\lambda}_\eta}{\partial\lambda_\alpha}
     *      \, \frac{\partial \bar{\lambda}_\epsilon}{\partial\lambda_\beta}
     *    + \bar{\gamma}_\eta
     *      \, \frac{\partial^2 \bar{\lambda}_\eta}{\partial\lambda_\alpha\partial\lambda_\beta}
     * \f]
     *  with
     * \f[
     *    \frac{\partial^2 \bar{\lambda}_\alpha}{\partial\lambda_\beta\partial\lambda_\eta}
     *    = -\frac{1}{3} J^{-1/3} \Big(
     *      1_{\alpha\beta} \lambda_\eta^{-1}
     *      + 1_{\alpha\eta}\lambda_\beta^{-1}
     *      - \big( 1_{\beta\eta} + \frac{1}{3} \big) \lambda_\alpha \lambda_\beta^{-1}
     *      \lambda_\eta^{-1}
     *    \Big)
     * \f]
     *
     *  <h3>References</h3>
     *  See Holzapfel, p 245-246, p 257-259, p 263-264
     *
     *  \author bborn
     *  \date 08/09
     * @param cmat (out) : Material elasticity tensor in Voigt notation
     * @param stress (out) : 2n Piola-Kirchhoff stress tensor in stress-like Voigt notation
     * @param rcg (in) : Right Cauchy-Green strain tensor in strain-like Voigt notation
     * @param havecoeffstrpr (in) : Flag, whether at least one summand is based in principal
     * stretches
     * @param havecoeffstrmod (in) : Flag, whether at least one summand is based in modified
     * principal stretches
     * @param eleGID (in) : Element id
     * @param potsum (in) : Summands of the Free-energy function
     */
    static void ResponseStretches(LINALG::Matrix<6, 6>& cmat,  ///< material elasticity matrix
        LINALG::Matrix<6, 1>& stress,     ///< 2nd Piola-Kirchhoff stress vector
        const LINALG::Matrix<6, 1>& rcg,  ///< Right Cauchy-Green tensor
        const bool& havecoeffstrpr,  ///< true, if at least one summand is based on pr.stretches
        const bool&
            havecoeffstrmod,  ///< true, if at least one summand is based on mod.pr.stretches
        const int eleGID,     ///< element GID
        const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum);

    /*!
     * \brief Check if material state is polyconvex
     *
     *  Polyconvexity of isotropic hyperelastic material
     *  dependent on principal or modified invariants)
     *  is tested with eq. (5.31) of Vera Ebbing - PHD-thesis (p. 79).
     *  Partial derivatives of SEF are used.
     *  \author birzle
     *  \date 04/2016
     *
     * @param defgrd (in) : Deformation gradient
     * @param prinv (in) : Principal invariants of the Right Cauchy-Green tensor
     * @param dPI (in) : First derivative of the Free-energy function with respect to the principal
     * invariants
     * @param ddPII (in) : Second derivative of the Free-energy function with respect to the
     * principal invariants
     * @param params (in/out) : Container for additional information
     * @param eleGID (in) : Element id
     * @param isoprinc (in) : Flag, whether at least one summand is formulated in the isotropic
     * principal invariants
     * @param isomod : Flag, whether at least one summand is formulated in the isotropic modified
     * principal invariants
     */
    static void CheckPolyconvexity(const LINALG::Matrix<3, 3>& defgrd,
        const LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<3, 1>& dPI,
        const LINALG::Matrix<6, 1>& ddPII, Teuchos::ParameterList& params, const int eleGID,
        bool isoprinc, bool isomod);

    /// Sat Inverse Analysis: add params of summands
    virtual void ElastOptParams(std::map<std::string, int>* pnames);

    /// @name Flags to specify the elastic formulations (initialize with false)
    //@{
    bool isoprinc_;    ///< global indicator for isotropic principal formulation
    bool isomod_;      ///< global indicator for isotropic splitted formulation
    bool anisoprinc_;  ///< global indicator for anisotropic principal formulation
    bool anisomod_;    ///< global indicator for anisotropic splitted formulation
    //@}

    /// @name Voigt index arrays
    //@{
    static const int VOIGT6ROW_[];    ///< 6-Voigt row index of corresponding 2-tensor
    static const int VOIGT6COL_[];    ///< 6-Voigt column index of corresponding 2-tensor
    static const int VOIGT3X3SYM_[];  ///< go from symmetric 2-tensor index pair to 6-Voigt index
                                      ///< by [3*i+j] for any i,j=0,1,2
    static const int
        VOIGT3X3NONSYM_[3][3];  ///< from non-symmetric 2-tensor index pair to 9-Voigt index
    //@}

    /// my material parameters
    MAT::PAR::ElastHyper* params_;

    /// map to materials/potential summands
    std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>> potsum_;
  };

}  // namespace MAT

#endif
