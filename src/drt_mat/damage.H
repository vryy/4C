/*----------------------------------------------------------------------*/
/*! \file
\brief Contains the functions to establish local material law /
       stress-strain law for isotropic material following nonlinear isotropic
       von Mises plasticity and a linear elastic material law
       (St.Venant Kirchhoff).

       isotropic hardening
       - describing the nonlinear (piecewise) hardening curve via \f$ sigma_y
          \ and \f$strainbar_p_ref \ from input file

       geometric linear, for small strains including isotropic ductile damage
       - simplified Lemaitre model only considers isotropic hardening

       ductile isotropic damage
       - elasticity-damage coupling, cf. de Souza Neto Compuatational plasticity Chapt. 12
       - Lemaitre's elastoplastic damage theory (simplified version without
         kinematic hardening)
       - isotropic hardening internal variable \f$ R\
       - damage variable \f$ D\

       ductile isotropic damage according to the book of de Souza Neto et al
       "Computational methods for plasticity", chapter 12

       example input line:
       MAT 1 MAT_Struct_Damage YOUNG 206.9 NUE 0.29 DENS 0.0 SAMPLENUM 2
         SIGMA_Y 0.45 0.65 EPSBAR_P 0.0 1.0 DAMDEN 0.0035 DAMEXP 1.0
         DAMTHRESHOLD 1.0e-06 KINHARD 17 KINHARD_REC 21 TOL 1.0e-6

\level 2

\maintainer Jan Schnabel
*/
/*----------------------------------------------------------------------*
 | definitions                                               dano 09/13 |
 *----------------------------------------------------------------------*/
#ifndef DAMAGE_H_
#define DAMAGE_H_

/*----------------------------------------------------------------------*
 | headers                                                   dano 09/13 |
 *----------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"


namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    //! material parameters for small strain elasto-plastic material with
    //! isotropic ductile damage
    class Damage : public Parameter
    {
     public:
      //! standard constructor
      Damage(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! destructor
      virtual ~Damage() { ; }

      //! @name material parameters
      //@{

      //! Young's modulus
      const double youngs_;
      //! Possion's ratio
      const double poissonratio_;
      //! mass density
      const double density_;
      //! isotropic hardening curve
      const std::vector<double> sigma_y_;
      //! accumulated plastic strain
      const std::vector<double> strainbar_p_ref_;
      //! damage evolution law denominator r
      const double damden_;
      //! damage evolution law exponent s
      const double damexp_;
      //! damage threshold
      const double epsbarD_;
      //! stress-like kinematic hardening modulus
      const double kinhard_;
      //! scalar-valued kinematic hardening modulus
      const double kinhard_rec_;
      // saturation hardening at reference temperature
      const double sathardening_;
      //! hardening exponent
      const double hardexpo_;
      //! tolerance for local Newton iteration
      const double abstol_;

      //@}

      //! create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class Damage

  }  // namespace PAR


  class DamageType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "DamageType"; }

    static DamageType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static DamageType instance_;

  };  // DamageType


  /*----------------------------------------------------------------------*/
  //! wrapper for elasto-plastic material including damage
  class Damage : public So3Material
  {
   public:
    //! construct empty material object
    Damage();

    //! construct the material object given material parameters
    explicit Damage(MAT::PAR::Damage* params);

    //! @name Packing and Unpacking

    //!  \brief Return unique ParObject id
    //!
    //!  every class implementing ParObject needs a unique id defined at the
    //!  top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const { return DamageType::Instance().UniqueParObjectId(); }

    //!  \brief Pack this class so it can be communicated
    //!
    //!  Resizes the vector data and stores all information of a class in it.
    //!  The first information to be stored in data has to be the
    //!  unique parobject id delivered by UniqueParObjectId() which will then
    //!  identify the exact class on the receiving processor.
    //!
    //!  \param data (in/out): char vector to store class information
    virtual void Pack(
        DRT::PackBuffer& data  //!<  data (i/o): char vector to store class information
        ) const;

    //!  \brief Unpack data from a char vector into this class
    //!
    //!  The vector data contains all information to rebuild the
    //!  exact copy of an instance of a class on a different processor.
    //!  The first entry in data has to be an integer which is the unique
    //!  parobject id defined at the top of this file and delivered by
    //!  UniqueParObjectId().
    //!
    //!  \param data (in) : vector storing all data to be unpacked into this
    //!  instance.
    virtual void Unpack(const std::vector<char>&
            data  //!< (i) : vector storing all data to be unpacked into this instance.
    );

    //@}

    //! material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_elpldamage; }

    //! return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new Damage(*this)); }

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem)
    {
      if (!(kinem == INPAR::STR::kinem_linear))
        dserror("element and material kinematics are not compatible");
    }

    //! initialise internal stress variables
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    //! update internal stress variables
    virtual void Update();

    //! evaluate material
    virtual void Evaluate(const LINALG::Matrix<3, 3>* defgrd,  //!< deformation gradient
        const LINALG::Matrix<NUM_STRESS_3D, 1>* linstrain,     //!< linear total strains
        Teuchos::ParameterList& params,                        //!< parameter list for communication
        LINALG::Matrix<NUM_STRESS_3D, 1>* stress,              //!< 2nd PK-stress
        LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>* cmat,    //!< material stiffness matrix
        int gp,                                                ///< Gauss point
        int eleGID);

    //! evaluate material using Lemaitre material model
    virtual void EvaluateFullLemaitre(const LINALG::Matrix<3, 3>* defgrd,  //!< deformation gradient
        const LINALG::Matrix<NUM_STRESS_3D, 1>* linstrain,                 //!< linear total strains
        Teuchos::ParameterList& params,                      //!< parameter list for communication
        LINALG::Matrix<NUM_STRESS_3D, 1>* stress,            //!< 2nd PK-stress
        LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>* cmat,  //!< material stiffness matrix
        const int eleGID);                                   //!< element GID

    //! evaluate material using the simplified Lemaitre model
    virtual void EvaluateSimplifiedLemaitre(
        const LINALG::Matrix<3, 3>* defgrd,                  //!< deformation gradient
        const LINALG::Matrix<NUM_STRESS_3D, 1>* linstrain,   //!< linear total strains
        Teuchos::ParameterList& params,                      //!< parameter list for communication
        LINALG::Matrix<NUM_STRESS_3D, 1>* stress,            //!< 2nd PK-stress
        LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>* cmat,  //!< material stiffness matrix
        const int eleGID);                                   //!< element GID

    //! computes stress
    void Stress(const double p,                             //!< volumetric stress tensor
        const LINALG::Matrix<NUM_STRESS_3D, 1>& devstress,  //!< deviatoric stress tensor
        LINALG::Matrix<NUM_STRESS_3D, 1>& stress            //!< 2nd PK-stress
    );

    //! calculate relative/over stress
    void RelStress(const LINALG::Matrix<NUM_STRESS_3D, 1>& devstrain,  //!< deviatoric stress tensor
        const LINALG::Matrix<NUM_STRESS_3D, 1>& beta,                  //!< back stress tensor
        LINALG::Matrix<NUM_STRESS_3D, 1>& eta                          //!< relative stress
    );

    //! computes isotropic elasticity tensor in matrix notion for 3d
    void SetupCmat(LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat  //!< elastic material tangent
    );

    //! computes isotropic elastoplastic tensor in matrix notion for 3d
    void SetupCmatElastoPlastic(
        LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,  //!< elasto-plastic material tangent
        int eleID,
        double Dgamma,         //!< plastic multiplier
        double G,              //!< shear modulus
        double kappa,          //!< bulk modulus
        double p_tilde,        // undamaged pressure
        double q_tilde,        // undamaged trial von Mises equivalent stress
        double energyrelrate,  // damage energy release rate
        double Ytan,           // derivative of engergy release rate Ytan w.r.t. Dgamma
        double sigma_y, double Hiso, LINALG::Matrix<NUM_STRESS_3D, 1> Nbar,
        int gp,             //!< current Gauss point
        bool damevolution,  //!< flag indicating if damage threshold is passed, i.e. if damage
                            //!< evolves
        double heaviside    //!< heaviside function: decide if loading/unloading
    );

    //! computes isotropic elastoplastic damage tensor in matrix nation for 3d
    void SetupCmatElastoPlasticFullLemaitre(LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>&
                                                cmat,  //!< elasto-plastic tangent modulus (out)
        LINALG::Matrix<NUM_STRESS_3D, 1> N_tilde,      //!< flow vector
        LINALG::Matrix<NUM_STRESS_3D, 1>& stress,      //!< total stress
        double heaviside,                              //!< heaviside-function
        double Dgamma,                                 //!< plastic multiplier
        double s_N,  //!< s_N = 2 * G * (Dgamma / omega_{n+1})^2 * s / r * (Y_{n+1} / r)^(s-1) * (1
                     //!< + nu) / E * s^tilde_{n+1} : N^tilde_{n+1}
        double
            g,  //!< g = (2 * G * (Dgamma / omega_{n+1}) + a * Dgamma / (1 + b * Dgamma)) / q_tilde
        double h_alg,       //!< h_alg = lengthy expression defined in Newton Raphson section
        double G,           //!< shear modulus
        double dkappa_dR,   //!< derivative of the hardening curve
        double kappa,       //!< bulk modulus
        double Hkin,        //!< kinematic hardening variable 1 (kinematic hardening modulus)
        double Hkin_rec,    //!< kinematic hardening variable 2 (saturation effect)
        double Nbetaold,    //!< Nbetaold = (N_tilde : beta_n)
        int gp,             // current Gauss point
        double qbar_tilde,  // effective relative eta^{~}
        double y,           // (-Y / r)^s
        LINALG::Matrix<NUM_STRESS_3D, 1> dy_dsigma_tilde,  // dy/dsigma_tilde
        LINALG::Matrix<NUM_STRESS_3D, 1> b_NbetaoldN  // beta_n - 2/3 . (N_tilde : beta_n) . N_tilde
    );

    //! calculates the derivative of GetSigmaYAtStrainbarnp() w.r.t. astrain_{n+1}
    //! and returns the isotropic hardening modulus
    double GetIsoHardAtStrainbarnp(
        const double strainbar_p  //!< current accumulated strain
                                  //!< or if damage!=0: isotropic hardening internal variable
    );

    //! calculates current yield stress from (sigma_y-astrain)-samples which
    //! describe a piecewise constant function
    double GetSigmaYAtStrainbarnp(
        const double strainbar_p  //!< current accumulated strain, in case of dependent hardening
                                  //!< or if damage!=0: isotropic hardening internal variable
    );

    //! return density
    virtual double Density() const { return params_->density_; }

    //! return accumulated strain at Gauss points
    //! use the old vector (last_) for postprocessing
    //! Output is called after(!!) Update, so that newest values are included in
    //! the old history vectors last_, while the current history vectors curr_
    //! are reset
    double AccumulatedStrain(int gp  //!< current Gauss point
        ) const
    {
      return strainbarpllast_->at(gp);
    }

    //! return damaged isotropic hardening variable R at Gauss points
    //! in undamaged case: corresponding to the accumulated strain strainbar_p
    //! use the old vector (last_) for postprocessing
    //! Output is called after(!!) Update, so that newest values are included in
    //! the old history vectors last_, while the current history vectors curr_
    //! are reset
    double IsotropicHardeningVariable(int gp  //!< current Gauss point
        ) const
    {
      return isohardvarlast_->at(gp);
    }

    //! return scalar-valued isotropic damage variable at Gauss points
    //! use the old vector (last_) for postprocessing
    //! Output is called after(!!) Update, so that newest values are included in
    //! the old history vectors last_, while the current history vectors curr_
    //! are reset
    double IsotropicDamage(int gp  //!< current Gauss point
        ) const
    {
      return damagelast_->at(gp);
    }

    //! check if history variables are already initialised
    bool Initialized() const { return (isinit_ and (strainplcurr_ != Teuchos::null)); }

    //! return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    //! return names of visualization data
    virtual void VisNames(std::map<std::string, int>& names);

    //! return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

   private:
    //! my material parameters
    MAT::PAR::Damage* params_;

    //! plastic history vector
    //! old plastic strain at t_n
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>
        strainpllast_;  //!< \f${\varepsilon}^p_{n}\f$
    //! current plastic strain at t_n+1
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>
        strainplcurr_;  //!< \f${\varepsilon}^p_{n+1}\f$
    //! old back stress at t_n
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>
        backstresslast_;  //!< \f${\beta}_{n}\f$
    //! current back stress at t_n+1
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>
        backstresscurr_;  //!< \f${\beta}_{n+1}\f$
    //! old accumulated plastic strain at t_n
    Teuchos::RCP<std::vector<double>> strainbarpllast_;  //!< \f${\varepsilon}^p_{n}\f$
    //! current accumulated plastic strain at t_n+1
    Teuchos::RCP<std::vector<double>> strainbarplcurr_;  //!< \f${\varepsilon}^p_{n+1}\f$
    //! old damaged accumulated plastic strain corresponding to isotropic hardening at t_n
    Teuchos::RCP<std::vector<double>> isohardvarlast_;  //!< \f${R}_{n+1}\f$
    //! current damaged accumulated plastic strain corresponding to isotropic hardening at t_n+1
    Teuchos::RCP<std::vector<double>> isohardvarcurr_;  //!< \f${R}_{n}\f$
    //! old damage internal state variable D at t_n
    Teuchos::RCP<std::vector<double>> damagelast_;  //!< \f${D}_{n}\f$
    //! current damage internal state variable D at t_n+1
    Teuchos::RCP<std::vector<double>> damagecurr_;  //!< \f${D}_{n+1}\f$

    //! indicator if #Initialize routine has been called
    bool isinit_;
    //! indicator if material has started to be plastic
    bool plastic_step_;

  };  // class Damage : public Material
}  // namespace MAT


/*----------------------------------------------------------------------*/
#endif  // DAMAGE_H_
