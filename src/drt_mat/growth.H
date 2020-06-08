/*----------------------------------------------------------------------*/
/*! \file
\brief This file contains routines for integration point based isotropic and anisotropic volumetric
growth laws.

\level 2

   \maintainer Amadeus Gebauer
 *----------------------------------------------------------------------*/
#ifndef GROWTH_H
#define GROWTH_H


#include "so3_material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_container.H"

/*------------------------------------------------------------------------*/
/* forward declarations */

namespace MAT
{
  class GrowthLaw;

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /*! \class Growth
     *  \brief Common parameters for volumetric growth
     *
     *  \author kehl
     *  \date 6/2015
     */
    class Growth : public Parameter
    {
     public:
      /// standard constructor
      Growth(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~Growth() { ; }

      /// @name material parameters
      //@{
      /// elastic material number
      const int idmatelastic_;
      /// growthlaw material number
      const int idgrowthlaw_;
      /// start growth after starttime
      const double starttime_;
      /// stop growth after endtime
      const double endtime_;
      /// growth law
      Teuchos::RCP<MAT::GrowthLaw> growthlaw_;
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class Growth

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /*! \class Growth
       \brief Wrapper material for volumetric growth

       Base class for isovolumetric growth materials.

       It keeps a pointer to the elastic material and to a specific
       growth law (via its Parameters) and some variables common to
       all growth materials.

       Further it basically passes through the interfaces functions
       Evaluate(...) and EvaluateNonLinMass(...).

       \author kehl
       \date 6/2015
  */
  class Growth : public So3Material
  {
   public:
    /// construct empty material object
    Growth();

    /// construct the material object given material parameters
    explicit Growth(MAT::PAR::Growth* params);

    //! @name Packing and Unpacking
    //@{

    /*!
      \brief Return unique ParObject id

      \sa DRT::ParObject
    */
    virtual int UniqueParObjectId() const = 0;

    /*!
      \brief Pack this class so it can be communicated

      \sa DRT::ParObject
      \param data (in/out): char vector to store class information
    */
    void Pack(DRT::PackBuffer& data) const;

    /*!
      \brief Unpack data from a char vector into this class

      \sa DRT::ParObject
      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    void Unpack(const std::vector<char>& data);

    //@}

    /// material type
    INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_volumetric; }

    /// check if element kinematics and material kinematics are compatible
    void ValidKinematics(INPAR::STR::KinemType kinem)
    {
      if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const = 0;

    /// Setup
    void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    /// Reset internal variables to state in the beginning of the computation (needed for inverse
    /// analysis)
    void ResetAll(int numgp);

    /// Update
    void Update();

    /// Reset time step
    void ResetStep();

    /// Store history/internal variables
    void StoreHistory(int timestep);

    /// Set history/internal variables
    void SetHistory(int timestep);


    //! @name Evaluation methods
    //@{

    /*! \brief Evaluate material stresses and cmat
     *
     * \param In
     * defgrd - deformation gradient
     * \param In
     * glstrain - green lagrange strain
     * \param In
     * params - a parameter list as handed in from the element
     * \param Out
     * stress - stresses from material evaluation
     * \param Out
     * cmat - material stiffness matrix
     * \param In
     * gp - Gauss point
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 06/2015
     */
    void Evaluate(const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* glstrain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat,
        const int gp, const int eleGID) = 0;

    /*! \brief Evaluate mass change
     *
     * \param In
     * defgrd - deformation gradient
     * \param In
     * glstrain - green lagrange strain
     * \param In
     * params - a parameter list as handed in from the element
     * \param Out
     * linmass_disp - linearization wrt to displacements
     * \param Out
     * linmass_vel - linearization wrt to velocities
     * gp - Gauss point
     * \param In
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 06/2015
     */
    void EvaluateNonLinMass(const LINALG::Matrix<3, 3>* defgrd,
        const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        LINALG::Matrix<6, 1>* linmass_disp, LINALG::Matrix<6, 1>* linmass_vel, const int gp,
        const int eleGID) = 0;

    /*! \brief Evaluate elastic material stresses and cmat
     *
     * \param In
     * defgrd - deformation gradient
     * \param In
     * glstrain - green lagrange strain
     * \param Out
     * stress - stresses from material evaluation
     * \param Out
     * cmat - material stiffness matrix
     * \param In
     * params - a parameter list as handed in from the element
     * \param In gp - Gauss point
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 06/2015
     */
    void EvaluateElastic(const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* glstrain,
        LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat, Teuchos::ParameterList& params,
        int gp, int eleGID);
    //@}


    //! @name Query  methods
    //@{

    /// Return density
    double Density() const
    {
      dserror("growth material needs gauss point data for density!");
      return -1.0;
    }

    // growth law-specific!
    double Density(int gp) const;

    /// Return whether material has a varying material density
    bool VaryingDensity() const;

    /// Return quick accessible material parameter data
    MAT::PAR::Growth* Parameter() const { return params_; }

    /// access to the elastic material
    const Teuchos::RCP<MAT::So3Material> Matelastic() const { return matelastic_; }

    // read access to thetaold_
    Teuchos::RCP<const std::vector<double>> ThetaOld() const { return thetaold_; }

    // read access to theta_
    Teuchos::RCP<const std::vector<double>> Theta() const { return theta_; }

    /// Return theta at Gauss point
    double Theta_atgp(int gp) const { return (*theta_)[gp]; }

    // Return thetaold at Gauss-point
    double ThetaOld_atgp(int gp) const { return ThetaOld()->at(gp); }

    // read access to isinit_
    bool IsInit() const { return isinit_; }

    //@}

   protected:
    /// growth stretch
    Teuchos::RCP<std::vector<double>> theta_;

    /// indicates if material is initialized
    bool isinit_;

   private:
    /// my material parameters
    MAT::PAR::Growth* params_;

    /// elastic material
    Teuchos::RCP<MAT::So3Material> matelastic_;

    /// growth stretch at the time step before
    Teuchos::RCP<std::vector<double>> thetaold_;

    /// data container to store and recall history data
    DRT::Container histdata_;

  };  // class Growth


  /*! \class GrowthVolumetricType
       \brief ParObjectType instance

       \sa DRT::ParObjectType

    \author kehl
    \date 6/2015
  */
  class GrowthVolumetricType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "GrowthVolumetricType"; }

    static GrowthVolumetricType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static GrowthVolumetricType instance_;
  };


  /*----------------------------------------------------------------------*/
  /*! \class GrowthVolumetric
      \brief Split of growth and elastic part of the deformation

      Implementation of the multiplicative split into isovolumetric growth
      and elastic part of the deformation.

      Evaluation of the growth factor \f$\theta\f$ and its derivative
      \f$\frac{\partial \vartheta}{\partial C}\f$ is passed to a generic growth law.

      \author kehl
      \date 6/2015
  */
  class GrowthVolumetric : public Growth
  {
   public:
    /// construct empty material object
    GrowthVolumetric();

    /// construct the material object given material parameters
    explicit GrowthVolumetric(MAT::PAR::Growth* params);

    //! @name Packing and Unpacking
    //@{

    /*!
      \brief Return unique ParObject id

      \sa DRT::ParObject
    */
    int UniqueParObjectId() const { return GrowthVolumetricType::Instance().UniqueParObjectId(); }

    /*!
      \brief Pack this class so it can be communicated

      \sa DRT::ParObject
      \param data (in/out): char vector to store class information
    */
    void Pack(DRT::PackBuffer& data) const;

    /*!
      \brief Unpack data from a char vector into this class

      \sa DRT::ParObject
      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    void Unpack(const std::vector<char>& data);

    //@}

    /// return copy of this material object
    Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new GrowthVolumetric(*this)); }

    /// Setup
    void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    /// Update
    void Update();

    /// Reset internal variables to state in the beginning of the computation (needed for inverse
    /// analysis)
    void ResetAll(int numgp);

    /// Return names of visualization data
    void VisNames(std::map<std::string, int>& names);

    /// Return visualization data
    bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

    //! @name Evaluation methods
    //@{

    /*! \brief Evaluate material stresses and cmat
     *
     * \param In
     * defgrd - deformation gradient
     * \param In
     * glstrain - green lagrange strain
     * \param In
     * params - a parameter list as handed in from the element
     * \param Out
     * stress - stresses from material evaluation
     * \param Out
     * cmat - material stiffness matrix
     * \param In
     * gp - Gauss point
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 06/2015
     */
    void Evaluate(const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* glstrain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat,
        const int gp, const int eleGID);

    /// calculate stresses and elastic material tangent (both in Voigt notation)
    void GetSAndCmatdach(const double theta, const LINALG::Matrix<3, 3>* defgrd,
        LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmatdach,
        Teuchos::ParameterList& params, int gp, int eleGID);

    /// Function which reads in the given fiber value due to the FIBER1 nomenclature
    void ReadFiber(DRT::INPUT::LineDefinition* linedef, std::string specifier,
        LINALG::Matrix<3, 1>& fiber_vector);

    /// pull back of a symmetric elastic 4th order tensor (in matrix/voigt notation)
    /// via the 2th order deformation gradient (also in matrix notation)
    static LINALG::Matrix<6, 6> PullBack4Tensor(
        const LINALG::Matrix<3, 3>& defgr, const LINALG::Matrix<6, 6>& Cmat);

    /// pull back the ijkl-th entry of a symmetric elastic 4th order tensor (in matrix/voigt
    /// notation) via the 2th order deformation gradient (also in matrix notation)
    static double PullBack4Tensorijkl(const LINALG::Matrix<3, 3>& defgr,
        const double (&FourTensor)[3][3][3][3], const double& i, const double& j, const double& k,
        const double& l);

    // Setup 4-Tensor from 6x6 Voigt notation
    static void Setup4Tensor(
        double (&FourTensor)[3][3][3][3], const LINALG::Matrix<6, 6>& VoigtMatrix);

    // Print Four Tensor
    void PrintFourTensor(double (&FourTensor)[3][3][3][3]);



    /*! \brief Evaluate mass change
     *
     * \param In
     * defgrd - deformation gradient
     * \param In
     * glstrain - green lagrange strain
     * \param In
     * params - a parameter list as handed in from the element
     * \param Out
     * linmass_disp - linearization wrt to displacements
     * \param Out
     * linmass_vel - linearization wrt to velocities
     * \param In
     * gp - Gauss point
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 06/2015
     */
    void EvaluateNonLinMass(const LINALG::Matrix<3, 3>* defgrd,
        const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        LINALG::Matrix<6, 1>* linmass_disp, LINALG::Matrix<6, 1>* linmass_vel, const int gp,
        const int eleGID);

    /*! \brief evaluate the volumetric growth factor
     *
     * This function just hands over the comuptation of the growth factor
     * theta to the attached growthlaw.
     *
     * \param Out
     * theta - the volumetric growth factor
     * \param Out
     * dthetadC - derivative of theta wrt the chauchy green tensor
     * \param In
     * matelastic - an elastic material
     * \param In
     * defgrad - the deformation gradient
     * \param In
     * glstrain - the green lagrange strains
     * \param In
     * params - a parameter list as handed in from the element
     * \param In
     * gp - Gauss point
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 06/2015
     */
    void EvaluateGrowth(double* theta, LINALG::Matrix<6, 1>* dthetadC,
        const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* glstrain,
        Teuchos::ParameterList& params, int gp, int eleGID);

    //@}

    //! @name Query  methods
    //@{

    /// Return quick accessible material parameter data
    MAT::PAR::Growth* Parameter() const { return paramsVolumetric_; }

    //@}

   protected:
    /// trace of elastic mandel stress - for postprocessing
    Teuchos::RCP<std::vector<double>> tr_mandel_e_;
    /// elastic fiber stretch - for postprocessing
    Teuchos::RCP<std::vector<double>> lambda_fib_e_;
    /// growth trigger for aniso stress / strain growth laws with constant prescribed trigger
    /// (element-wise, in .dat file!)
    double growthtrig_const_;

   private:
    /// my material parameters
    MAT::PAR::Growth* paramsVolumetric_;

    /// a reference direction (i.e., a fiber direction or a surface normal etc...), for anisotropic
    /// growth laws
    LINALG::Matrix<3, 1> refdir_;

    /// some current direction, used for anisotropic scalar-dependent growth in current radial
    /// direction...
    std::vector<LINALG::Matrix<3, 1>> curdir_;

    /// normal direction, used for anisotropic scalar-dependent growth in current radial
    /// direction...
    std::vector<LINALG::Matrix<3, 1>> curdir_for_update_;

    /// history of growth matrix, needed for anisotropic scalar-dependent growth in current radial
    /// direction...
    std::vector<LINALG::Matrix<3, 3>> F_g_hist_;

  };  // class GrowthVolumetric

}  // namespace MAT

#endif /*GROWTH_H_*/
