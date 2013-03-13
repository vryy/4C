/*----------------------------------------------------------------------*/
/*!
\file damage.H
\brief TODO Contains the functions to establish local material law /
       stress-strain law for isotropic material for a 3D hex element
       following ideal von Mises plasticity and a linear elastic material law
       (St.Venant Kirchhoff).

       TODO plastic materail with kinematic hardening
        - independent yield stress level of degree of plastification
        - constant uniaxial yield stress sigma_y = const.

       TODO geometric linear, for small strains

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef DAMAGE_H_
#define DAMAGE_H_

#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    //! material parameters for lasto-plastic material
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
      //! initial yield stress (constant)
      const double yield_;
      //! linear kinematic hardening modulus
      const double kinhard_;
      // tolerance for local Newton iteration
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

  static DamageType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static DamageType instance_;

};  // DamageType


  /*----------------------------------------------------------------------*/
  //! wrapper for elasto-plastic material including damage
  class Damage : public Material
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
    virtual void Pack(DRT::PackBuffer& data) const;

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
    virtual void Unpack(const std::vector<char>& data);

    //@}

    //! material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_elpldamage; }

    //! return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new Damage( *this ) ); }

    //! initialise internal stress variables
    void Setup(
      const int numgp //!< number of Gauss points
      );

    //! update internal stress variables
    void Update();

    //! reset internal stress variables
    void Reset();

    //! evaluate material
    void Evaluate(
      const LINALG::Matrix<NUM_STRESS_3D,1>& linstrain,
      LINALG::Matrix<NUM_STRESS_3D,1>& plstrain, //!< return the plastic strains to the element
//      LINALG::Matrix<1,1>& accplstrain, //!< return accumulated plastic strains to the element
      const int gp,  //!< current Gaus point
      Teuchos::ParameterList& params,  //!< parameter list for communication
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat,  //!< material stiffness matrix
      LINALG::Matrix<NUM_STRESS_3D,1>& stress  //!< 2nd PK-stress
      );

    //! computes stress
    void Stress(
      const double p,  //!< volumetric stress tensor
      const LINALG::Matrix<NUM_STRESS_3D,1>& devstress,  //!< deviatoric stress tensor
      LINALG::Matrix<NUM_STRESS_3D,1>& stress //!< 2nd PK-stress
      );

    //! calculate relative/over strees
    void RelStress(
      const LINALG::Matrix<NUM_STRESS_3D,1>& devstrain,  //!< deviatoric stress tensor
      const LINALG::Matrix<NUM_STRESS_3D,1>& beta,  //!< back stress tensor
      LINALG::Matrix<NUM_STRESS_3D,1>& eta //!< relative stress
      );

    //! computes isotropic elasticity tensor in matrix notion for 3d
    void SetupCmat(
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat  //!< elastic material tangent
      );

    //! computes isotropic elastoplastic tensor in matrix notion for 3d
    void SetupCmatElastoPlastic(
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat,  //!< elasto-plastic material tangent
      double Dgamma,  //!< plastic multiplier
      double G,  //!< shear modulus
      double q,  //!< effective stress
      const LINALG::Matrix<NUM_STRESS_3D,1> flowvector,  //!< flow vector
      const LINALG::Matrix<NUM_STRESS_3D,1> devstress,  //!< relative stress
      double heaviside  //!< heaviside function: decide if loading/unloading
      );

    //! Return density
    virtual double Density() const { return params_->density_; }
    
    //! Return accumulated strain at Gauss points
    //! use the old vector (last_) for postprocessing
    //! Output is called after(!!) Update, so that newest values are included in
    //! the old history vectors last_, while the current history vectors curr_
    //! are reset
    LINALG::Matrix<1,1> AccumulatedStrain(
      int gp
      )
    const { return ( strainbarpllast_->at(gp) ); }

    //! Check if history variables are already initialized
    bool Initialized() const
    {
      return ( isinit_ and (strainplcurr_!=Teuchos::null) );
    }

    //! Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    //! Flag plastic step was called
    bool plastic_step;

  private:

    //! my material parameters
    MAT::PAR::Damage* params_;

    //! plastic history vector
    // old plastic strain at t_n
    Teuchos::RCP<std::vector< LINALG::Matrix<NUM_STRESS_3D,1> > > strainpllast_;  //!< \f${\varepsilon}^p_{n}\f$
    // current plastic strain at t_n+1
    Teuchos::RCP<std::vector< LINALG::Matrix<NUM_STRESS_3D,1> > > strainplcurr_;  //!< \f${\varepsilon}^p_{n+1}\f$
    // old accumulated plastic strain at t_n
    Teuchos::RCP<std::vector< LINALG::Matrix<1,1> > > strainbarpllast_;  //!< \f${\varepsilon}^p_{n}\f$
    // current accumulated plastic strain at t_n+1
    Teuchos::RCP<std::vector< LINALG::Matrix<1,1> > > strainbarplcurr_;  //!< \f${\varepsilon}^p_{n+1}\f$
    // old back stress at t_n
    Teuchos::RCP<std::vector< LINALG::Matrix<NUM_STRESS_3D,1> > > backstresslast_;  //!< \f${\beta}_{n}\f$
    // current back stress at t_n+1
    Teuchos::RCP<std::vector< LINALG::Matrix<NUM_STRESS_3D,1> > > backstresscurr_;  //!< \f${\beta}_{n+1}\f$


    //! indicator if #Initialize routine has been called
    bool isinit_;

  };  // class Damage : public Material
}  // namespace MAT


/*----------------------------------------------------------------------*/
#endif  // DAMAGE_H_


