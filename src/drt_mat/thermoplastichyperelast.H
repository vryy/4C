/*----------------------------------------------------------------------*/
/*!
\file thermoplastichyperelast.H
\brief Contains the functions to establish local material law /
       stress-strain law for isotropic material following finite strain
       von-Mises plasticity with nonlinear isotropic hardening and general
       hyperelasticity (for the time being: NeoHooke).

       implementation is based on
       Simo and Miehe: "Associative coupled thermoplasticity at finite strains:
       Formulation, numerical analysis and implementation", in Computer Methods
       in Applied Mechanics and Engineering, 98:41–104, 1992.

       geometrically nonlinear, for finite strains, rate-independent, thermo-plastic

       example input line:
       [mm,ms,kg,K]
       MAT 1 MAT_Struct_ThrPlasticHyperElast YOUNG 206.9 NUE 0.29 DENS 7.8e-6
       CTE 1e-5 INITTEMP 293 YIELD 0.45 ISOHARD 0.12924 SATHARDENING 0.715
       HARDEXPO 16.93 YIELDSOFT 0.002 HARDSOFT 0.002 DISSFACT 0.9

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*
 | definitions                                               dano 03/13 |
 *----------------------------------------------------------------------*/
#ifndef THERMO_PLAST_HYPERELAST_H_
#define THERMO_PLAST_HYPERELAST_H_

/*----------------------------------------------------------------------*
 | headers                                                   dano 03/13 |
 *----------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    //! material parameters for neo-Hooke
    class ThermoPlasticHyperElast : public Parameter
    {
    public:

      //! standard constructor
      ThermoPlasticHyperElast(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! destructor
      virtual ~ThermoPlasticHyperElast() { ; }

      //! @name material parameters
      //@{

      //! Young's modulus
      const double youngs_;
      //! Possion's ratio
      const double poissonratio_;
      //! mass density
      const double density_;
      //! coefficient of thermal expansion
      const double cte_;
      //! initial, reference temperature
      const double inittemp_;
      //! initial yield stress (constant) at reference temperature
      const double yield_;
      //! linear isotropic hardening modulus at reference temperature
      const double isohard_;
      // saturation hardening at reference temperature
      const double sathardening_;
      //! hardening exponent
      const double hardexpo_;
      //! yield stress softening
      const double yieldsoft_;
      //! hardening softening
      const double hardsoft_;

      //@}

      //! create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class ThermoPlasticHyperElast

  }  // namespace PAR


  class ThermoPlasticHyperElastType : public DRT::ParObjectType
  {
  public:

    std::string Name() const { return "ThermoPlasticHyperElastType"; }

    static ThermoPlasticHyperElastType & Instance() { return instance_; };

    virtual DRT::ParObject* Create( const std::vector<char> & data );

  private:

    static ThermoPlasticHyperElastType instance_;

  };  // class ThermoPlasticHyperElastType

  /*----------------------------------------------------------------------*/
  //! wrapper for finite strain elasto-plastic material
  class ThermoPlasticHyperElast : public So3Material
  {
  public:

    //! construct empty material object
    ThermoPlasticHyperElast();

    //! construct the material object given material parameters
    explicit ThermoPlasticHyperElast(MAT::PAR::ThermoPlasticHyperElast* params);

    //! @name Packing and Unpacking

    /*!
    \brief Return unique ParObject id

    every class implementing ParObject needs a unique id defined at the
    top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
    {
      return ThermoPlasticHyperElastType::Instance().UniqueParObjectId();
    }

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

    //! @name Access methods

    //! material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_thermoplhyperelast; }

    //! return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new ThermoPlasticHyperElast( *this ) ); }

    //! Young's modulus
    double Youngs() const { return params_->youngs_; }

    //! Poisson's ratio
    double PoissonRatio() const { return params_->poissonratio_; }

    //! density
    virtual double Density() const { return params_->density_; }

    //! shear modulus
    double ShearMod() const { return 0.5*params_->youngs_/(1.0+params_->poissonratio_); }

    //! yield stress
    virtual double Yield() const { return params_->yield_; }

    //! isotropic hardening modulus
    virtual double IsoHard() const { return params_->isohard_; }

    //! flow stress softening
    virtual double FlowStressSoft() const { return params_->yieldsoft_; }

    //! hardening softening
    virtual double HardSoft() const { return params_->hardsoft_; }

    //! saturation hardening
    virtual double SatHardening() const { return params_->sathardening_; }

    //! coefficient of thermal expansion
    virtual double CTE() const { return params_->cte_; }

    //! coefficient of thermal expansion
    virtual double HardExpo() const { return params_->hardexpo_; }

    //! initial, reference temperature
    virtual double InitTemp() const { return params_->inittemp_; }


    //! return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    //! return accumulated strain at Gauss points
    //! use the old vector (last_) for postprocessing
    //! Output is called after(!!) Update, so that newest values are included in
    //! the old history vectors last_, while the current history vectors curr_
    //! are reset
    double AccumulatedStrain(int gp)
    const { return ( accplstrainlast_->at(gp) ); }

    //! mechanical dissipation
    double MechDiss(int gp)
    const { return ( mechdisscurr_->at(gp) ); }

    //! linearisation of the mechanical dissipation w.r.t. displacements d_{n+1}
    LINALG::Matrix<6,1> MechDiss_kTd(int gp)
    const { return ( mechdiss_kTd_->at(gp) ); }

    // TODO wo ist in cmat die Temperatur enthalten??
    //! linearisation of material tangent w.r.t. temperatures T_{n+1}
    LINALG::Matrix<6,1> CMat_kdT(int gp)
    const { return ( Cmat_kdT_->at(gp) ); }

    // Couple Dissipation K
    double MechDiss_kTT(int gp)
    const { return ( mechdiss_kTT_curr_->at(gp) ); }

    //! check if history variables are already initialized
    bool Initialized() const
    {
      return ( isinit_ and (defgrdcurr_!=Teuchos::null) );
    }

    //! return names of visualization data
    virtual void VisNames(std::map<std::string,int>& names);

    //! return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp);

    //@}

    //! @name Evaluation methods

    //! initialise internal stress variables
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    //! update internal stress variables
    virtual void Update();

    //! evaluate material law
    virtual void Evaluate(
      const LINALG::Matrix<3,3>* defgrd,  //!< input deformation gradient for multiplicative sp
      const LINALG::Matrix<6,1>* glstrain,  //!< input Green-Lagrange strain (redundant with defo
                                            //   but used for neo-hooke evaluation; maybe remove
      Teuchos::ParameterList& params,  //!< input paramter list (e.g. Young's, ...)
      LINALG::Matrix<6,1>* stress,  //!< output (mandatory) second Piola-Kirchhoff stress
      LINALG::Matrix<6,6>* cmat  //!< output (mandatory) material stiffness matrix
      );

    //! TODO
    void SetupCmatElastoPlastic(
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat,  // elasto-plastic tangent modulus (out)
      double Dgamma,  // plastic multiplier
      double Hiso_temp,  // shear modulus
      double sigma_y0infty_temp,  // elastic trial von Mises effective stress
      double sigma_y0_temp,
      double mubar,
      double hardexpo,
      double alphanew,
      double q_trial,
      const LINALG::Matrix<3,3>& defgrd,
      LINALG::Matrix<3,3> invdefgrdcurr,
      LINALG::Matrix<3,3> n,
      double J,
      double kappa,
      LINALG::Matrix<3,3> tmp,
      int gp
      );

    //! main 3D material call to determine stress and constitutive tensor ctemp
    //  originally method of fourieriso with const!!!
    void Evaluate(
      const LINALG::Matrix<1,1>& Ntemp,  //!< temperature of element
      LINALG::Matrix<NUM_STRESS_3D,1>& ctemp,  //!< temperature dependent material tangent
      LINALG::Matrix<NUM_STRESS_3D,1>& stresstemp,  //!< temperature dependent stress term
      const LINALG::Matrix<3,3>& defgrd,
      Teuchos::ParameterList& params
      );
      
    //! computes temperature dependent isotropic elasticity tensor in matrix
    //! notion for 3d
    void SetupCthermo(
      LINALG::Matrix<6,1>& ctemp,  //!< temperature dependent material tangent
      const LINALG::Matrix<3,3>& defgrd,
      Teuchos::ParameterList& params
      );
    

    //! calculates stress-temperature modulus
    double STModulus();

    //@}

  private:

    //! my material parameters
    MAT::PAR::ThermoPlasticHyperElast* params_;

    //! plastic history variables
    //! old (i.e. at t_n)  deformation gradient at each Gauss-point
    Teuchos::RCP<std::vector<LINALG::Matrix<3,3> > > defgrdlast_;
    //! current (i.e. at t_n+1) deformation gradient at each Gauss-point
    Teuchos::RCP<std::vector<LINALG::Matrix<3,3> > > defgrdcurr_;

    //! old (i.e. at t_n) elastic, isochoric right Cauchy-Green tensor
    Teuchos::RCP<std::vector<LINALG::Matrix<3,3> > > bebarlast_;
    //! current (i.e. at t_n+1) elastic, isochoric right Cauchy-Green tensor
    Teuchos::RCP<std::vector<LINALG::Matrix<3,3> > > bebarcurr_;
    
    //! old (i.e. at t_n) accumulated plastic strain
    Teuchos::RCP<std::vector<double> > accplstrainlast_;
    //! current (i.e. at t_n+1) accumulated plastic strain
    Teuchos::RCP<std::vector<double> > accplstraincurr_;

    //! old (i.e. at t_n) mechanical dissipation
    Teuchos::RCP<std::vector<double> > mechdisslast_;
    //! current (i.e. at t_n+1) mechanical dissipation
    Teuchos::RCP<std::vector<double> > mechdisscurr_;

    //! current (i.e. at t_n+1) linearised mechanical dissipation
    Teuchos::RCP<std::vector<double> > mechdiss_kTT_curr_;

    //! current (i.e. at t_n+1) coupling term: mechanical dissipation structure linearization
    Teuchos::RCP<std::vector<LINALG::Matrix<6,1> > > mechdiss_kTd_;
    //! current (i.e. at t_n+1) coupling term: mechanical dissipation temperature linearization
    Teuchos::RCP<std::vector<LINALG::Matrix<6,1> > > Cmat_kdT_;

    //! indicator if #Initialize routine has been called
    bool isinit_;

  };  // class ThermoPlasticHyperElast

}  // namespace MAT

/*----------------------------------------------------------------------*/
#endif  // THERMO_PLAST_HYPERELAST_H_
