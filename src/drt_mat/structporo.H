/*----------------------------------------------------------------------*/
/*! \file
 \brief wrapper for structure material of porous media

\maintainer Johannes Kremheller

\level 2
 *-----------------------------------------------------------------------*/

#ifndef STRUCTPORO_H_
#define STRUCTPORO_H_

/*---------------------------------------------------------------------*
 | headers                                                              |
 *---------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  // forward declaration
  class StructPoro;

  namespace PAR
  {
    class PoroLaw;

    class StructPoro : public Parameter
    {
      friend class MAT::StructPoro;

     public:
      /// standard constructor
      StructPoro(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~StructPoro() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// @name material parameters
      //@{

      /// material ID of sub-material
      int matid_;

      /// initial porosity
      int porolawID_;

      /// initial porosity
      double initporosity_;

      //@}

      // implementation of porosity law
      PoroLaw* porolaw_;

    };  // class StructPoro

  }  // namespace PAR

  class StructPoroType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "StructPoroType"; }

    static StructPoroType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static StructPoroType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for StructPoro material
  ///
  /// This object exists (several times) at every element

  /*!
    The idea is to use any material formulation within the poro framework.
    Therefore, a poro material wraps the 'real' material and holds it as
    a private member. For most evaluation routines it will just call this material.
    In addition it provides poro specific functions, as giving the constitutive law
    for the porosity.

    Main methods of this material are the ComputePorosity(...) methods, providing
    the porosity and its derivatives. If the constitutive law is a governing equation
    itself (for poro P1 elements, for instance), the material evaluates the
    consitutive law itself and its derivatives in the ConstituitiveDerivatives(...)
    methods.
    All other Evaluate() methods are basically passed through to the underlying
    structure material.

    The poro material can save the porosity gauss point wise. Therefore it
    has an additional setup method, giving the number of gauss points. This is
    only (!) meant for post processing/visualization processes! The gauss point
    wise saved porosity must not be used during simulation as it is not
    guaranteed (and actually not the case) that the gauss point numbering
    is the same for every element (especially for e.g. fluid and solid elements).

   \author vuong 10/14
   */
  class StructPoro : public So3Material
  {
   public:
    /// construct empty material object
    StructPoro();

    /// construct the material object given material parameters
    explicit StructPoro(MAT::PAR::StructPoro* params);

    //! @name Packing and Unpacking

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of drt_parobject.H (this file) and should return it in this method.
     */
    virtual int UniqueParObjectId() const { return StructPoroType::Instance().UniqueParObjectId(); }

    /*!
     \brief Pack this class so it can be communicated

     Resizes the vector data and stores all information of a class in it.
     The first information to be stored in data has to be the
     unique parobject id delivered by UniqueParObjectId() which will then
     identify the exact class on the receiving processor.

     \param data (in/out): char vector to store class information
     */
    virtual void Pack(DRT::PackBuffer& data) const;

    /*!
     \brief Unpack data from a char vector into this class

     The vector data contains all information to rebuild the
     exact copy of an instance of a class on a different processor.
     The first entry in data has to be an integer which is the unique
     parobject id defined at the top of this file and delivered by
     UniqueParObjectId().

     \param data (in) : vector storing all data to be unpacked into this
     instance.
     */
    virtual void Unpack(const std::vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_structporo; }

    /// poro law type
    virtual INPAR::MAT::MaterialType PoroLawType() const;

    /// return inverse bulkmodulus (=compressibility)
    double InvBulkmodulus() const;

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem) { mat_->ValidKinematics(kinem); }

    /// return material
    Teuchos::RCP<MAT::Material> GetMaterial() const { return mat_; }

    /// return material ID
    int MatID() const { return params_->matid_; }

    /// return porosity average (for post processing only!)
    double PorosityAv() const;

    /// return initial porosity
    double Initporosity() const { return params_->initporosity_; }

    /// return time derivative of reference porosity (only nonzero with reaction)
    virtual double RefPorosityTimeDeriv() const { return 0.0; }

    /// compute current porosity and save it
    virtual void ComputePorosity(Teuchos::ParameterList& params,  ///< (i) element parameter list
        double press,                                             ///< (i) pressure at gauss point
        double J,          ///< (i) determinant of jacobian at gauss point
        int gp,            ///< (i) number of current gauss point
        double& porosity,  ///< (o) porosity at gauss point
        double* dphi_dp,   ///< (o) first derivative of porosity w.r.t. pressure at gauss point
        double* dphi_dJ,   ///< (o) first derivative of porosity w.r.t. jacobian at gauss point
        double*
            dphi_dJdp,  ///< (o) derivative of porosity w.r.t. pressure and jacobian at gauss point
        double* dphi_dJJ,  ///< (o) second derivative of porosity w.r.t. jacobian at gauss point
        double* dphi_dpp,  ///< (o) second derivative of porosity w.r.t. pressure at gauss point
        bool save = true);

    /// compute current porosity and save it
    void ComputePorosity(Teuchos::ParameterList& params,  ///< (i) element parameter list
        double press,                                     ///< (i) pressure at gauss point
        double J,          ///< (i) determinant of jacobian at gauss point
        int gp,            ///< (i) number of current gauss point
        double& porosity,  ///< (o) porosity at gauss point
        bool save = true);

    /// compute current surface porosity and save it
    void ComputeSurfPorosity(Teuchos::ParameterList& params,  ///< (i) element parameter list
        double press,                                         ///< (i) pressure at gauss point
        double J,           ///< (i) determinant of jacobian at gauss point
        const int surfnum,  ///< (i) number of surface
        int gp,             ///< (i) number of current gauss point
        double& porosity,   ///< (o) porosity at gauss point
        double* dphi_dp,    ///< (o) first derivative of porosity w.r.t. pressure at gauss point
        double* dphi_dJ,    ///< (o) first derivative of porosity w.r.t. jacobian at gauss point
        double*
            dphi_dJdp,  ///< (o) derivative of porosity w.r.t. pressure and jacobian at gauss point
        double* dphi_dJJ,  ///< (o) second derivative of porosity w.r.t. jacobian at gauss point
        double* dphi_dpp,  ///< (o) second derivative of porosity w.r.t. pressure at gauss point
        bool save = true);

    /// compute current surface porosity and save it
    void ComputeSurfPorosity(Teuchos::ParameterList& params,  ///< (i) element parameter list
        double press,                                         ///< (i) pressure at gauss point
        double J,           ///< (i) determinant of jacobian at gauss point
        const int surfnum,  ///< (i) number of surface
        int gp,             ///< (i) number of current gauss point
        double& porosity,   ///< (o) porosity at gauss point
        bool save = true);

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new StructPoro(*this)); }

    /// Initialize internal variables
    virtual void PoroSetup(int numgp,  ///< number of Gauss points
        DRT::INPUT::LineDefinition* linedef);

    //! Calculate coupling part of homogenized 2 Piola-Kirchhoff stress (3D)
    void CouplStress(const LINALG::Matrix<3, 3>& defgrad,  ///< (i) deformation gradient
        const LINALG::Matrix<3, 1>& fluidvel,              ///< (i) fluid velocity in gausspoint
        const double& press,                               ///< (i) pressure in gauss point
        LINALG::Matrix<6, 1>& couplstress                  ///< (o) coupling stress in gauss point
        ) const;

    //! Calculate coupling part of homogenized 2 Piola-Kirchhoff stress (2D)
    void CouplStress(const LINALG::Matrix<2, 2>& defgrad,  ///< (i) deformation gradient
        const LINALG::Matrix<2, 1>& fluidvel,              ///< (i) fluid velocity in gausspoint
        const double& press,                               ///< (i) pressure in gauss point
        LINALG::Matrix<3, 1>& couplstress                  ///< (o) coupling stress in gauss point
        ) const;

    //! evaluate constitutive relation for porosity and compute derivatives
    virtual void ConstitutiveDerivatives(Teuchos::ParameterList& params,  ///< (i) parameter list
        double press,        ///< (i) fluid pressure at gauss point
        double J,            ///< (i) Jacobian determinant at gauss point
        double porosity,     ///< (i) porosity at gauss point
        double* dW_dp,       ///< (o) derivative of potential w.r.t. pressure
        double* dW_dphi,     ///< (o) derivative of potential w.r.t. porosity
        double* dW_dJ,       ///< (o) derivative of potential w.r.t. jacobian
        double* dW_dphiref,  ///< (o) derivative of potential w.r.t. reference porosity
        double* W            ///< (o) inner potential
    );

    //! evaluate constitutive relation for porosity and compute derivatives using reference porosity
    void ConstitutiveDerivatives(Teuchos::ParameterList& params,  ///< (i) parameter list
        double press,        ///< (i) fluid pressure at gauss point
        double J,            ///< (i) Jacobian determinant at gauss point
        double porosity,     ///< (i) porosity at gauss point
        double refporosity,  ///< (i) porosity at gauss point
        double* dW_dp,       ///< (o) derivative of potential w.r.t. pressure
        double* dW_dphi,     ///< (o) derivative of potential w.r.t. porosity
        double* dW_dJ,       ///< (o) derivative of potential w.r.t. jacobian
        double* dW_dphiref,  ///< (o) derivative of potential w.r.t. reference porosity
        double* W            ///< (o) inner potential
    );

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    //! @name Evaluation methods

    /// evaluate material law
    virtual void Evaluate(const LINALG::Matrix<3, 3>* defgrd,  ///< (i) deformation gradient
        const LINALG::Matrix<6, 1>* glstrain,                  ///< (i) green lagrange strain
        Teuchos::ParameterList& params,                        ///< (i) parameter list
        LINALG::Matrix<6, 1>* stress,  ///< (o) second piola kirchhoff stress
        LINALG::Matrix<6, 6>* cmat,    ///< (o) constitutive matrix
        int EleID)                     ///< (i) element GID
    {
      mat_->Evaluate(defgrd, glstrain, params, stress, cmat, EleID);
    };

    /// Evaluate strain energy function (for hyperelastic materials only)
    virtual void StrainEnergy(const LINALG::Matrix<6, 1>& glstrain,  ///< (i) green lagrange strain
        double& psi,                                                 ///< (o) strain energy
        int EleID                                                    ///< (i) element GID
    )
    {
      mat_->StrainEnergy(glstrain, psi, EleID);
    };

    /// Evaluate material law for GEMM
    virtual void EvaluateGEMM(LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* stress,
        LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>* cmat, double* density,
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_m,
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_new,
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_old, LINALG::Matrix<3, 3>* rcg_new,
        LINALG::Matrix<3, 3>* rcg_old, const int eleGID)
    {
      mat_->EvaluateGEMM(
          stress, cmat, density, glstrain_m, glstrain_new, glstrain_old, rcg_new, rcg_old, eleGID);
    };

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
        LINALG::Matrix<9, 1>* D2sntDFDT = NULL) override
    {
      mat_->EvaluateCauchy(defgrd, n, t, snt, DsntDn, DsntDt, DsntDF, D2sntDF2, D2sntDFDn,
          D2sntDFDt, eleGID, temp, DsntDT, D2sntDFDT);
    }

    //@}

    /// Return material density (if provided by the specific material)
    virtual double Density() const;
    virtual double DensitySolidPhase() const;

    //! @name Handling of Gauss point data. Here, the poro material just calls the underlying
    //! material

    /// Set up for materials with GP data (e.g., history variables)
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
    {
      // setup the underlying material
      // Note: poro material itself is setup when calling PoroSetup()
      mat_->Setup(numgp, linedef);
    };

    /// Update of GP data (e.g., history variables)
    virtual void Update() { mat_->Update(); };

    /// Reset time step (for time adaptivity)
    virtual void ResetStep() { mat_->ResetStep(); };

    /// Reset internal variables to state in the beginning of the computation (needed for inverse
    /// analysis)
    virtual void ResetAll(int numgp) { mat_->ResetAll(numgp); };

    //@}

    //! @name Visualization methods

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string, int>& names);

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

    //@}

   protected:
    /// compute current porosity and save it
    void ComputePorosity(
        const double& refporosity,  ///< (i) initial/reference porosity at gauss point
        const double& press,        ///< (i) pressure at gauss point
        const double& J,            ///< (i) determinant of jacobian at gauss point
        const int& gp,              ///< (i) number of current gauss point
        double& porosity,           ///< (o) porosity at gauss point
        double* dphi_dp,  ///< (o) first derivative of porosity w.r.t. pressure at gauss point
        double* dphi_dJ,  ///< (o) first derivative of porosity w.r.t. jacobian at gauss point
        double*
            dphi_dJdp,  ///< (o) derivative of porosity w.r.t. pressure and jacobian at gauss point
        double* dphi_dJJ,      ///< (o) second derivative of porosity w.r.t. jacobian at gauss point
        double* dphi_dpp,      ///< (o) second derivative of porosity w.r.t. pressure at gauss point
        double* dphi_dphiref,  ///< (o) derivative of porosity w.r.t. reference porosity (only
                               ///< nonzero with reaction)
        bool save = true);

    /// my material parameters
    MAT::PAR::StructPoro* params_;

    /// actual material
    Teuchos::RCP<MAT::So3Material> mat_;

    /// porosity at gauss points
    Teuchos::RCP<std::vector<double>> porosity_;

    /// porosity at gauss points of surface element
    Teuchos::RCP<std::map<int, std::vector<double>>> surfporosity_;

    /// flag indicating initialization of attributes
    bool isinitialized_;
  };

}  // namespace MAT

#endif /* STRUCTPORO_H_ */
