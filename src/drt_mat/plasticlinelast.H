/*----------------------------------------------------------------------*/
/*!
\brief Contains the functions to establish local material law /
       stress-strain law for isotropic material for a 3D hex element
       following ideal von Mises plasticity and a linear elastic material law
       (St.Venant Kirchhoff).

       small strains

       plastic materail with kinematic hardening
        - independent yield stress level of degree of plastification
        - constant uniaxial yield stress sigma_y = const.

       extend to linear isotropic hardening
        - yield stress no longer constant, depends on level of accumulated
          plastic strain, i.e. \f$ \sigma_y \,=\, \sigma_y(\bar{\epsilon}_p)\f$

       example input line:
       MAT 1 MAT_Struct_PlasticLinElast YOUNG 206.9 NUE 0.29 DENS 0.0
         YIELD 0.45 ISOHARD 0.0 KINHARD 0.0 TOL 1.0e-6

\level 2
<pre>
\maintainer Matthias Mayr
</pre>
*/
/*----------------------------------------------------------------------*
 | definitions                                               dano 04/11 |
 *----------------------------------------------------------------------*/
#ifndef PLASTICLINELAST_H_
#define PLASTICLINELAST_H_


/*----------------------------------------------------------------------*
 | headers                                                   dano 04/11 |
 *----------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"


namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    //! material parameters for linear elasto-plastic material
    class PlasticLinElast : public Parameter
    {
     public:
      //! standard constructor
      PlasticLinElast(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! destructor
      virtual ~PlasticLinElast() { ; }

      //! @name material parameters
      //@{

      //! Young's modulus
      const double youngs_;
      //! Possion's ratio
      const double poissonratio_;
      //! mass density
      const double density_;
      //! initial yield stress (constant)
      const double yield_;
      //! linear isotropic hardening modulus
      const double isohard_;
      //! linear kinematic hardening modulus
      const double kinhard_;
      //! tolerance for local Newton iteration
      const double abstol_;

      //@}

      //! create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class PlasticLinElast

  }  // namespace PAR


  class PlasticLinElastType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "PlasticLinElastType"; }

    static PlasticLinElastType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static PlasticLinElastType instance_;

  };  // PlasticLinElastType


  /*----------------------------------------------------------------------*/
  //! wrapper for linear elasto-plastic material
  class PlasticLinElast : public So3Material
  {
   public:
    //! construct empty material object
    PlasticLinElast();

    //! construct the material object given material parameters
    explicit PlasticLinElast(MAT::PAR::PlasticLinElast* params);

    //! @name Packing and Unpacking

    //!  \brief Return unique ParObject id
    //!
    //!  every class implementing ParObject needs a unique id defined at the
    //!  top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const
    {
      return PlasticLinElastType::Instance().UniqueParObjectId();
    }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_pllinelast; }

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem)
    {
      if (!(kinem == INPAR::STR::kinem_linear))
        dserror("element and material kinematics are not compatible");
    }

    //! return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new PlasticLinElast(*this));
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
        const int eleGID);

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
        double Dgamma,                                       //!< plastic multiplier
        double G,                                            //!< shear modulus
        double q,                                            //!< effective stress
        const LINALG::Matrix<NUM_STRESS_3D, 1> Nbar,         // unit flow vector
        double heaviside,  //!< heaviside function: decide if loading/unloading
        double Hiso,       //!< isotropic hardening modulus
        double Hkin        //!< kinematic hardening modulus
    );

    //! computes isotropic elastoplastic tensor in matrix notion for 3d
    void SetupCmatElastoPlastic2(
        LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,  //!< elasto-plastic material tangent
        double Dgamma,                                       //!< plastic multiplier
        double q,                                            //!< effective stress
        LINALG::Matrix<NUM_STRESS_3D, 1> unitflow            //!< unit flow vector
    );

    void SetupContinuumCmatElastoPlastic(
        LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat,  //!< elasto-plastic material tangent
        double Dgamma,                                       //!< plastic multiplier
        double q,                                            //!< effective stress
        LINALG::Matrix<NUM_STRESS_3D, 1> unitflow            //!< unit flow vector
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

    //! check if history variables are already initialised
    bool Initialized() const { return (isinit_ and (strainplcurr_ != Teuchos::null)); }

    //! return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    //! return names of visualization data
    virtual void VisNames(std::map<std::string, int>& names);

    //! return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

    //! finite difference check for debugging purposes
    void FDCheck(LINALG::Matrix<NUM_STRESS_3D, 1>& stress,  //!< updated stress sigma_n+1
        LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>&
            cmatFD,  //!< material tangent calculated with FD of stresses
        LINALG::Matrix<NUM_STRESS_3D, 1>& beta,          //!< updated back stresses
        double p,                                        //!< volumetric stress
        const LINALG::Matrix<NUM_STRESS_3D, 1>& strain,  //!< elastic strain vector
        double Dgamma,                                   //!< plastic multiplier
        double G,                                        //!< shear modulus
        double qbar,                          //!< elastic trial von Mises effective stress
        double kappa,                         //!< bulk modulus
        LINALG::Matrix<NUM_STRESS_3D, 1>& N,  // flow vector
        double heaviside                      //!< Heaviside function
    );

   private:
    //! my material parameters
    MAT::PAR::PlasticLinElast* params_;

    //! plastic history vector
    //! old plastic strain at t_n
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>
        strainpllast_;  //!< \f${\varepsilon}^p_{n}\f$
    //! current plastic strain at t_n+1
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>
        strainplcurr_;  //!< \f${\varepsilon}^p_{n+1}\f$
    //! old accumulated plastic strain at t_n
    Teuchos::RCP<std::vector<double>> strainbarpllast_;  //!< \f${\varepsilon}^p_{n}\f$
    //! current accumulated plastic strain at t_n+1
    Teuchos::RCP<std::vector<double>> strainbarplcurr_;  //!< \f${\varepsilon}^p_{n+1}\f$
    //! old back stress at t_n
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>
        backstresslast_;  //!< \f${\beta}_{n}\f$
    //! current back stress at t_n+1
    Teuchos::RCP<std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>>
        backstresscurr_;  //!< \f${\beta}_{n+1}\f$

    //! indicator if #Initialize routine has been called
    bool isinit_;
    //! indicator if material has started to be plastic
    bool plastic_step_;

  };  // class PlasticLinElast : public Material
}  // namespace MAT


/*----------------------------------------------------------------------*/
#endif  // PLASTICLINELAST_H_
