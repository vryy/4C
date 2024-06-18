/*----------------------------------------------------------------------*/
/*! \file
\brief This file contains routines for integration point based isotropic and anisotropic volumetric
growth laws.

\level 2

 *----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_GROWTH_HPP
#define FOUR_C_MAT_GROWTH_HPP


#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

/*------------------------------------------------------------------------*/
/* forward declarations */

namespace Mat
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
    class Growth : public Core::Mat::PAR::Parameter
    {
     public:
      /// standard constructor
      Growth(const Core::Mat::PAR::Parameter::Data& matdata);


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
      Teuchos::RCP<Mat::GrowthLaw> growthlaw_;
      //@}

      /// create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;

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
       evaluate(...) and EvaluateNonLinMass(...).

       \author kehl
       \date 6/2015
  */
  class Growth : public So3Material
  {
   public:
    /// construct empty material object
    Growth();

    /// construct the material object given material parameters
    explicit Growth(Mat::PAR::Growth* params);

    //! @name Packing and Unpacking
    //@{

    /*!
      \brief Return unique ParObject id

      \sa Core::Communication::ParObject
    */
    int UniqueParObjectId() const override = 0;

    /*!
      \brief Pack this class so it can be communicated

      \sa Core::Communication::ParObject
      \param data (in/out): char vector to store class information
    */
    void Pack(Core::Communication::PackBuffer& data) const override;

    /*!
      \brief Unpack data from a char vector into this class

      \sa Core::Communication::ParObject
      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    void Unpack(const std::vector<char>& data) override;

    //@}

    /// material type
    Core::Materials::MaterialType MaterialType() const override
    {
      return Core::Materials::m_growth_volumetric;
    }

    /// check if element kinematics and material kinematics are compatible
    void ValidKinematics(Inpar::STR::KinemType kinem) override
    {
      if (kinem != Inpar::STR::KinemType::nonlinearTotLag)
        FOUR_C_THROW("element and material kinematics are not compatible");
    }

    /// return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override = 0;

    /// Setup
    void Setup(int numgp, Input::LineDefinition* linedef) override;

    /// Update
    void Update() override;

    /// Reset time step
    void reset_step() override;

    /// Store history/internal variables
    void StoreHistory(int timestep) override;

    /// Set history/internal variables
    void SetHistory(int timestep) override;


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
    void evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,
        const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        Core::LinAlg::Matrix<6, 1>* stress, Core::LinAlg::Matrix<6, 6>* cmat, const int gp,
        const int eleGID) override = 0;

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
    void EvaluateNonLinMass(const Core::LinAlg::Matrix<3, 3>* defgrd,
        const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        Core::LinAlg::Matrix<6, 1>* linmass_disp, Core::LinAlg::Matrix<6, 1>* linmass_vel,
        const int gp, const int eleGID) override = 0;

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
    void EvaluateElastic(const Core::LinAlg::Matrix<3, 3>* defgrd,
        const Core::LinAlg::Matrix<6, 1>* glstrain, Core::LinAlg::Matrix<6, 1>* stress,
        Core::LinAlg::Matrix<6, 6>* cmat, Teuchos::ParameterList& params, int gp, int eleGID);
    //@}


    //! @name Query  methods
    //@{

    /// Return density
    double Density() const override
    {
      FOUR_C_THROW("growth material needs gauss point data for density!");
      return -1.0;
    }

    // growth law-specific!
    double Density(int gp) const override;

    /// Return whether material has a varying material density
    bool VaryingDensity() const override;

    /// Return quick accessible material parameter data
    Mat::PAR::Growth* Parameter() const override { return params_; }

    /// access to the elastic material
    Teuchos::RCP<Mat::So3Material> Matelastic() const { return matelastic_; }

    // read access to thetaold_
    Teuchos::RCP<const std::vector<double>> ThetaOld() const { return thetaold_; }

    // read access to theta_
    Teuchos::RCP<const std::vector<double>> Theta() const { return theta_; }

    // Return thetaold at Gauss-point
    double ThetaOldAtGp(int gp) const { return ThetaOld()->at(gp); }

    // read access to isinit_
    bool is_init() const { return isinit_; }

    //@}

   protected:
    /// growth stretch
    Teuchos::RCP<std::vector<double>> theta_;

    /// indicates if material is initialized
    bool isinit_;

   private:
    /// my material parameters
    Mat::PAR::Growth* params_;

    /// elastic material
    Teuchos::RCP<Mat::So3Material> matelastic_;

    /// growth stretch at the time step before
    Teuchos::RCP<std::vector<double>> thetaold_;

    /// map to store and recall history data
    std::map<int, std::vector<double>> histdata_;
  };  // class Growth


  /*! \class GrowthVolumetricType
       \brief ParObjectType instance

       \sa Core::Communication::ParObjectType

    \author kehl
    \date 6/2015
  */
  class GrowthVolumetricType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "GrowthVolumetricType"; }

    static GrowthVolumetricType& Instance() { return instance_; };

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

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
    explicit GrowthVolumetric(Mat::PAR::Growth* params);

    //! @name Packing and Unpacking
    //@{

    /*!
      \brief Return unique ParObject id

      \sa Core::Communication::ParObject
    */
    int UniqueParObjectId() const override
    {
      return GrowthVolumetricType::Instance().UniqueParObjectId();
    }

    /*!
      \brief Pack this class so it can be communicated

      \sa Core::Communication::ParObject
      \param data (in/out): char vector to store class information
    */
    void Pack(Core::Communication::PackBuffer& data) const override;

    /*!
      \brief Unpack data from a char vector into this class

      \sa Core::Communication::ParObject
      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    void Unpack(const std::vector<char>& data) override;

    //@}

    /// return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new GrowthVolumetric(*this));
    }

    /// Setup
    void Setup(int numgp, Input::LineDefinition* linedef) override;

    /// Update
    void Update() override;

    /// Return names of visualization data
    void VisNames(std::map<std::string, int>& names) override;

    /// Return visualization data
    bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID) override;

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
    void evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,
        const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        Core::LinAlg::Matrix<6, 1>* stress, Core::LinAlg::Matrix<6, 6>* cmat, const int gp,
        const int eleGID) override;

    /// calculate stresses and elastic material tangent (both in Voigt notation)
    void GetSAndCmatdach(const double theta, const Core::LinAlg::Matrix<3, 3>* defgrd,
        Core::LinAlg::Matrix<6, 1>* stress, Core::LinAlg::Matrix<6, 6>* cmatdach,
        Teuchos::ParameterList& params, int gp, int eleGID);

    /// Function which reads in the given fiber value due to the FIBER1 nomenclature
    void ReadFiber(Input::LineDefinition* linedef, std::string specifier,
        Core::LinAlg::Matrix<3, 1>& fiber_vector);

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
    void EvaluateNonLinMass(const Core::LinAlg::Matrix<3, 3>* defgrd,
        const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        Core::LinAlg::Matrix<6, 1>* linmass_disp, Core::LinAlg::Matrix<6, 1>* linmass_vel,
        const int gp, const int eleGID) override;

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
    void EvaluateGrowth(double* theta, Core::LinAlg::Matrix<6, 1>* dthetadC,
        const Core::LinAlg::Matrix<3, 3>* defgrd, const Core::LinAlg::Matrix<6, 1>* glstrain,
        Teuchos::ParameterList& params, int gp, int eleGID);

    //@}

    //! @name Query  methods
    //@{

    /// Return quick accessible material parameter data
    Mat::PAR::Growth* Parameter() const override { return params_volumetric_; }

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
    Mat::PAR::Growth* params_volumetric_;

    /// a reference direction (i.e., a fiber direction or a surface normal etc...), for anisotropic
    /// growth laws
    Core::LinAlg::Matrix<3, 1> refdir_;

    /// some current direction, used for anisotropic scalar-dependent growth in current radial
    /// direction...
    std::vector<Core::LinAlg::Matrix<3, 1>> curdir_;

    /// normal direction, used for anisotropic scalar-dependent growth in current radial
    /// direction...
    std::vector<Core::LinAlg::Matrix<3, 1>> curdir_for_update_;

    /// history of growth matrix, needed for anisotropic scalar-dependent growth in current radial
    /// direction...
    std::vector<Core::LinAlg::Matrix<3, 3>> f_g_hist_;

  };  // class GrowthVolumetric

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
