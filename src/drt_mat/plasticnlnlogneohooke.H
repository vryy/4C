/*----------------------------------------------------------------------*/
/*!
\brief Contains the functions to establish local material law /
       stress-strain law for isotropic material following large strain
       von-Mises plasticity with linear isotropic hardening
       and logarithmic hyperelastic material (i.e. linear relation
       between Kirchhoff-stress and logarithmic strain; also known as Hencky
       material model).
       The principal stress-based implementation follows
       Bonet and Wood: "Nonlinear continuum mechanics for finite element analysis.",
       Cambridge University Press, Cambridge, 2008

       geometrically nonlinear, finite strains, visco-plastic

       example input line:
       MAT 1 MAT_Struct_PlasticNlnLogNeoHooke YOUNG 206.9 NUE 0.29 DENS 0.0
         YIELD 0.45 ISOHARD 0.12924 SATHARDENING 0.715 HARDEXPO 16.93 VISC 1.0 RATE_DEPENDENCY 0.1

\level 2

\maintainer Matthias Mayr

*/
/*----------------------------------------------------------------------*
 | definitions                                                          |
 *----------------------------------------------------------------------*/
#ifndef PLASTICNLNLOGNEOHOOKE_H_
#define PLASTICNLNLOGNEOHOOKE_H_

/*----------------------------------------------------------------------*
 | headers                                                              |
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
    class PlasticNlnLogNeoHooke : public Parameter
    {
     public:
      //! standard constructor
      PlasticNlnLogNeoHooke(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! destructor
      virtual ~PlasticNlnLogNeoHooke() { ; }

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
      //! saturation yield stress
      const double infyield_;
      //! nonlinear hardening exponent
      const double hardexp_;
      //! viscosity
      const double visc_;
      //! rate dependency
      const double rate_dependency_;

      //@}

      //! create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class PlasticNlnLogNeoHooke

  }  // namespace PAR


  class PlasticNlnLogNeoHookeType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "PlasticNlnLogNeoHookeType"; }

    static PlasticNlnLogNeoHookeType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static PlasticNlnLogNeoHookeType instance_;

  };  // class PlasticNlnLogNeoHookeType

  /*----------------------------------------------------------------------*/
  //! wrapper for finite strain elasto-plastic material
  class PlasticNlnLogNeoHooke : public So3Material
  {
   public:
    //! construct empty material object
    PlasticNlnLogNeoHooke();

    //! construct the material object given material parameters
    explicit PlasticNlnLogNeoHooke(MAT::PAR::PlasticNlnLogNeoHooke* params);

    //! @name Packing and Unpacking

    /*!
    \brief Return unique ParObject id

    every class implementing ParObject needs a unique id defined at the
    top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const
    {
      return PlasticNlnLogNeoHookeType::Instance().UniqueParObjectId();
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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_plnlnlogneohooke; }

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem)
    {
      if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    //! return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new PlasticNlnLogNeoHooke(*this));
    }

    //! density
    virtual double Density() const { return params_->density_; }

    //! return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    //! return accumulated strain at Gauss points
    //! use the old vector (last_) for postprocessing
    //! Output is called after(!!) Update, so that newest values are included in
    //! the old history vectors last_, while the current history vectors curr_
    //! are reset
    double AccumulatedStrain(int gp) const { return (accplstrainlast_->at(gp)); }

    //! check if history variables are already initialized
    bool Initialized() const { return (isinit_ and (accplstraincurr_ != Teuchos::null)); }

    //! return names of visualization data
    virtual void VisNames(std::map<std::string, int>& names);

    //! return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

    /// Return whether the material requires the deformation gradient for its evaluation
    virtual bool NeedsDefgrd() { return true; };

    //@}

    //! @name Evaluation methods

    //! initialise internal stress variables
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    //! update internal stress variables
    virtual void Update();

    //! evaluate material law
    virtual void Evaluate(
        const LINALG::Matrix<3, 3>* defgrd,    //!< input deformation gradient for multiplicative sp
        const LINALG::Matrix<6, 1>* glstrain,  //!< input Green-Lagrange strain (redundant with defo
                                               //   but used for neo-hooke evaluation; maybe remove
        Teuchos::ParameterList& params,        //!< input paramter list (e.g. Young's, ...)
        LINALG::Matrix<6, 1>* stress,          //!< output (mandatory) second Piola-Kirchhoff stress
        LINALG::Matrix<6, 6>* cmat,            //!< output (mandatory) material stiffness matrix
        const int eleGID);

    //@}

   private:
    //! my material parameters
    MAT::PAR::PlasticNlnLogNeoHooke* params_;

    //! inverse right cauchy green of plastic strain
    Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> invplrcglast_;
    //! inverse right cauchy green of plastic strain
    Teuchos::RCP<std::vector<LINALG::Matrix<3, 3>>> invplrcgcurr_;

    //! old (i.e. at t_n) accumulated palstic strain
    Teuchos::RCP<std::vector<double>> accplstrainlast_;
    //! current (i.e. at t_n+1) accumulated palstic strain
    Teuchos::RCP<std::vector<double>> accplstraincurr_;

    //! indicator if #Initialize routine has been called
    bool isinit_;

  };  // class PlasticNlnLogNeoHooke

}  // namespace MAT


/*----------------------------------------------------------------------*/
#endif  // PLASTICNLNLOGNEOHOOKE_H_
