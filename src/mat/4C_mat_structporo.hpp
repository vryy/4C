/*----------------------------------------------------------------------*/
/*! \file
 \brief wrapper for structure material of porous media


\level 2
 *-----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_STRUCTPORO_HPP
#define FOUR_C_MAT_STRUCTPORO_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  class StructPoro;

  namespace PAR
  {
    class PoroLaw;

    class StructPoro : public Core::Mat::PAR::Parameter
    {
      friend class Mat::StructPoro;

     public:
      //! standard constructor
      StructPoro(const Core::Mat::PAR::Parameter::Data& matdata);

      //! create material instance of matching type with my parameters
      Teuchos::RCP<Core::Mat::Material> create_material() override;

      //! @name material parameters
      //!@{

      //! material ID of sub-material
      int matid_;

      //! poro law ID
      int poro_law_ID_;

      //! initial porosity
      double init_porosity_;

      //!@}

      //! implementation of porosity law
      PoroLaw* poro_law_;
    };

  }  // namespace PAR

  class StructPoroType : public Core::Communication::ParObjectType
  {
   public:
    std::string Name() const override { return "StructPoroType"; }

    static StructPoroType& Instance() { return instance_; }

    Core::Communication::ParObject* Create(const std::vector<char>& data) override;

   private:
    static StructPoroType instance_;
  };

  /*----------------------------------------------------------------------*/
  //! Wrapper for StructPoro material
  //!
  //! This object exists (several times) at every element

  /*!
    The idea is to use any material formulation within the poro framework.
    Therefore, a poro material wraps the 'real' material and holds it as
    a private member. For most evaluation routines it will just call this material.
    In addition it provides poro specific functions, as giving the constitutive law
    for the porosity.

    Main methods of this material are the compute_porosity(...) methods, providing
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
   */
  class StructPoro : public So3Material
  {
   public:
    //! construct empty material object
    StructPoro();

    //! construct the material object given material parameters
    explicit StructPoro(Mat::PAR::StructPoro* params);

    //! @name Packing and Unpacking

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of parobject.H (this file) and should return it in this method.
     */
    int UniqueParObjectId() const override
    {
      return StructPoroType::Instance().UniqueParObjectId();
    }

    /*!
     \brief Pack this class so it can be communicated

     Resizes the vector data and stores all information of a class in it.
     The first information to be stored in data has to be the
     unique parobject id delivered by UniqueParObjectId() which will then
     identify the exact class on the receiving processor.

     \param data (in/out): char vector to store class information
     */
    void Pack(Core::Communication::PackBuffer& data) const override;

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
    void Unpack(const std::vector<char>& data) override;

    //!@}

    //! material type
    Core::Materials::MaterialType MaterialType() const override
    {
      return Core::Materials::m_structporo;
    }

    //! poro law type
    virtual Core::Materials::MaterialType PoroLawType() const;

    //! return inverse bulk modulus (=compressibility)
    double InvBulkModulus() const;

    //! check if element kinematics and material kinematics are compatible
    void ValidKinematics(Inpar::STR::KinemType kinem) override { mat_->ValidKinematics(kinem); }

    //! return material
    Teuchos::RCP<Core::Mat::Material> GetMaterial() const { return mat_; }

    //! return material ID
    int MatID() const { return params_->matid_; }

    //! return porosity average (for post processing only!)
    double PorosityAv() const;

    //! return initial porosity
    double InitPorosity() const { return params_->init_porosity_; }

    //! return time derivative of reference porosity (only nonzero with reaction)
    virtual double ref_porosity_time_deriv() const { return 0.0; }

    //! compute current porosity and save it
    virtual void compute_porosity(Teuchos::ParameterList& params,  //!< (i) element parameter list
        double press,                                              //!< (i) pressure at gauss point
        double J,          //!< (i) determinant of jacobian at gauss point
        int gp,            //!< (i) number of current gauss point
        double& porosity,  //!< (o) porosity at gauss point
        double* dphi_dp,   //!< (o) first derivative of porosity w.r.t. pressure at gauss point
        double* dphi_dJ,   //!< (o) first derivative of porosity w.r.t. jacobian at gauss point
        double*
            dphi_dJdp,  //!< (o) derivative of porosity w.r.t. pressure and jacobian at gauss point
        double* dphi_dJJ,  //!< (o) second derivative of porosity w.r.t. jacobian at gauss point
        double* dphi_dpp,  //!< (o) second derivative of porosity w.r.t. pressure at gauss point
        bool save = true);

    //! compute current porosity and save it
    void compute_porosity(Teuchos::ParameterList& params,  //!< (i) element parameter list
        double press,                                      //!< (i) pressure at gauss point
        double J,          //!< (i) determinant of jacobian at gauss point
        int gp,            //!< (i) number of current gauss point
        double& porosity,  //!< (o) porosity at gauss point
        bool save = true);

    //! compute current surface porosity and save it
    void ComputeSurfPorosity(Teuchos::ParameterList& params,  //!< (i) element parameter list
        double press,                                         //!< (i) pressure at gauss point
        double J,           //!< (i) determinant of jacobian at gauss point
        const int surfnum,  //!< (i) number of surface
        int gp,             //!< (i) number of current gauss point
        double& porosity,   //!< (o) porosity at gauss point
        double* dphi_dp,    //!< (o) first derivative of porosity w.r.t. pressure at gauss point
        double* dphi_dJ,    //!< (o) first derivative of porosity w.r.t. jacobian at gauss point
        double*
            dphi_dJdp,  //!< (o) derivative of porosity w.r.t. pressure and jacobian at gauss point
        double* dphi_dJJ,  //!< (o) second derivative of porosity w.r.t. jacobian at gauss point
        double* dphi_dpp,  //!< (o) second derivative of porosity w.r.t. pressure at gauss point
        bool save = true);

    //! compute current surface porosity and save it
    void ComputeSurfPorosity(Teuchos::ParameterList& params,  //!< (i) element parameter list
        double press,                                         //!< (i) pressure at gauss point
        double J,           //!< (i) determinant of jacobian at gauss point
        const int surfnum,  //!< (i) number of surface
        int gp,             //!< (i) number of current gauss point
        double& porosity,   //!< (o) porosity at gauss point
        bool save = true);

    //! return copy of this material object
    Teuchos::RCP<Core::Mat::Material> Clone() const override
    {
      return Teuchos::rcp(new StructPoro(*this));
    }

    //! Initialize internal variables
    virtual void poro_setup(int numgp,  //!< number of Gauss points
        Input::LineDefinition* linedef);

    /*!
     * @brief Calculate coupling part of homogenized 2 Piola-Kirchhoff stress (3D)
     *
     * @param[in] defgrad       deformation gradient
     * @param[in] press         pressure at gauss point
     * @param[out] couplstress  coupling stress at gauss point
     */
    void CouplStress(const Core::LinAlg::Matrix<3, 3>& defgrd, const double& press,
        Core::LinAlg::Matrix<6, 1>& couplstress) const;

    /*!
     * @brief Calculate coupling part of homogenized 2 Piola-Kirchhoff stress (2D)
     *
     * @param[in] defgrad       deformation gradient
     * @param[in] press         pressure at gauss point
     * @param[out] couplstress  coupling stress at gauss point
     */
    void CouplStress(const Core::LinAlg::Matrix<2, 2>& defgrd, const double& press,
        Core::LinAlg::Matrix<4, 1>& couplstress) const;

    //! evaluate constitutive relation for porosity and compute derivatives
    virtual void constitutive_derivatives(Teuchos::ParameterList& params,  //!< (i) parameter list
        double press,        //!< (i) fluid pressure at gauss point
        double J,            //!< (i) Jacobian determinant at gauss point
        double porosity,     //!< (i) porosity at gauss point
        double* dW_dp,       //!< (o) derivative of potential w.r.t. pressure
        double* dW_dphi,     //!< (o) derivative of potential w.r.t. porosity
        double* dW_dJ,       //!< (o) derivative of potential w.r.t. jacobian
        double* dW_dphiref,  //!< (o) derivative of potential w.r.t. reference porosity
        double* W            //!< (o) inner potential
    );

    //! evaluate constitutive relation for porosity and compute derivatives using reference porosity
    void constitutive_derivatives(Teuchos::ParameterList& params,  //!< (i) parameter list
        double press,        //!< (i) fluid pressure at gauss point
        double J,            //!< (i) Jacobian determinant at gauss point
        double porosity,     //!< (i) porosity at gauss point
        double refporosity,  //!< (i) porosity at gauss point
        double* dW_dp,       //!< (o) derivative of potential w.r.t. pressure
        double* dW_dphi,     //!< (o) derivative of potential w.r.t. porosity
        double* dW_dJ,       //!< (o) derivative of potential w.r.t. jacobian
        double* dW_dphiref,  //!< (o) derivative of potential w.r.t. reference porosity
        double* W            //!< (o) inner potential
    );

    //! Return quick accessible material parameter data
    Core::Mat::PAR::Parameter* Parameter() const override { return params_; }

    //! @name Evaluation methods

    void Evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,
        const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        Core::LinAlg::Matrix<6, 1>* stress, Core::LinAlg::Matrix<6, 6>* cmat, int gp,
        int EleID) override
    {
      mat_->Evaluate(defgrd, glstrain, params, stress, cmat, gp, EleID);
    }

    void StrainEnergy(const Core::LinAlg::Matrix<6, 1>& glstrain, double& psi, const int gp,
        const int EleID) override
    {
      mat_->StrainEnergy(glstrain, psi, gp, EleID);
    }

    void evaluate_cauchy_n_dir_and_derivatives(const Core::LinAlg::Matrix<3, 3>& defgrd,
        const Core::LinAlg::Matrix<3, 1>& n, const Core::LinAlg::Matrix<3, 1>& dir,
        double& cauchy_n_dir, Core::LinAlg::Matrix<3, 1>* d_cauchyndir_dn,
        Core::LinAlg::Matrix<3, 1>* d_cauchyndir_ddir, Core::LinAlg::Matrix<9, 1>* d_cauchyndir_dF,
        Core::LinAlg::Matrix<9, 9>* d2_cauchyndir_dF2,
        Core::LinAlg::Matrix<9, 3>* d2_cauchyndir_dF_dn,
        Core::LinAlg::Matrix<9, 3>* d2_cauchyndir_dF_ddir, int gp, int eleGID,
        const double* concentration, const double* temp, double* d_cauchyndir_dT,
        Core::LinAlg::Matrix<9, 1>* d2_cauchyndir_dF_dT) override
    {
      mat_->evaluate_cauchy_n_dir_and_derivatives(defgrd, n, dir, cauchy_n_dir, d_cauchyndir_dn,
          d_cauchyndir_ddir, d_cauchyndir_dF, d2_cauchyndir_dF2, d2_cauchyndir_dF_dn,
          d2_cauchyndir_dF_ddir, gp, eleGID, concentration, temp, d_cauchyndir_dT,
          d2_cauchyndir_dF_dT);
    }

    //!@}

    //! Return material density (if provided by the specific material)
    double Density() const override;
    virtual double DensitySolidPhase() const;

    //! @name Handling of Gauss point data. Here, the poro material just calls the underlying
    //! material

    void Setup(int numgp, Input::LineDefinition* linedef) override
    {
      // setup the underlying material
      // Note: poro material itself is setup when calling poro_setup()
      mat_->Setup(numgp, linedef);
    }

    void Update() override { mat_->Update(); }

    void reset_step() override { mat_->reset_step(); }

    //!@}

    //! @name Visualization methods

    void VisNames(std::map<std::string, int>& names) override;

    bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID) override;

    //!@}

   protected:
    //! compute current porosity and save it
    void compute_porosity(
        const double& refporosity,  //!< (i) initial/reference porosity at gauss point
        const double& press,        //!< (i) pressure at gauss point
        const double& J,            //!< (i) determinant of jacobian at gauss point
        const int& gp,              //!< (i) number of current gauss point
        double& porosity,           //!< (o) porosity at gauss point
        double* dphi_dp,  //!< (o) first derivative of porosity w.r.t. pressure at gauss point
        double* dphi_dJ,  //!< (o) first derivative of porosity w.r.t. jacobian at gauss point
        double*
            dphi_dJdp,  //!< (o) derivative of porosity w.r.t. pressure and jacobian at gauss point
        double* dphi_dJJ,      //!< (o) second derivative of porosity w.r.t. jacobian at gauss point
        double* dphi_dpp,      //!< (o) second derivative of porosity w.r.t. pressure at gauss point
        double* dphi_dphiref,  //!< (o) derivative of porosity w.r.t. reference porosity (only
                               //!< nonzero with reaction)
        bool save = true);

    //! my material parameters
    Mat::PAR::StructPoro* params_;

    //! actual material
    Teuchos::RCP<Mat::So3Material> mat_;

    //! porosity at gauss points
    Teuchos::RCP<std::vector<double>> porosity_;

    //! porosity at gauss points of surface element
    Teuchos::RCP<std::map<int, std::vector<double>>> surf_porosity_;

    //! flag indicating initialization of attributes
    bool is_initialized_;
  };

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
