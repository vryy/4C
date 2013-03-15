/*----------------------------------------------------------------------*/
/*!
\file plastichyperelast.H
\brief Contains the functions to establish local material law /
       stress-strain law for isotropic material following large strain
       von-Mises plasticity with combined isotropic / kinematic
       hardening and general hyperelasticity (for the time being: NeoHooke).

       geometrically nonlinear, for finite strains

       example input line:
       MAT 1 MAT_Struct_PlasticHyperElast YOUNG 206.9 NUE 0.29 DENS 0.0
         YIELD 0.45 ISOHARD 0.0 KINHARD 0.0


<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef PLASTICHYPERELAST_H_
#define PLASTICHYPERELAST_H_

#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for neo-Hooke
    class PlasticHyperElast : public Parameter
    {
    public:

      /// standard constructor
      PlasticHyperElast(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~PlasticHyperElast() { ; }

      /// @name material parameters
      //@{

      /// Young's modulus
      const double youngs_;
      /// Possion's ratio
      const double poissonratio_;
      /// mass density
      const double density_;
      //! initial yield stress (constant)
      const double yield_;
      //! linear isotropic hardening modulus
      const double isohard_;
      //! linear kinematic hardening modulus
      const double kinhard_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class PlasticHyperElast

  }  // namespace PAR


  class PlasticHyperElastType : public DRT::ParObjectType
  {
  public:

    std::string Name() const { return "PlasticHyperElastType"; }

    static PlasticHyperElastType & Instance() { return instance_; };

    virtual DRT::ParObject* Create( const std::vector<char> & data );

  private:

    static PlasticHyperElastType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for finite strain elasto-plastic material
  class PlasticHyperElast : public So3Material
  {
  public:

    /// construct empty material object
    PlasticHyperElast();

    /// construct the material object given material parameters
    explicit PlasticHyperElast(MAT::PAR::PlasticHyperElast* params);

    //! @name Packing and Unpacking

    /*!
    \brief Return unique ParObject id

    every class implementing ParObject needs a unique id defined at the
    top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return PlasticHyperElastType::Instance().UniqueParObjectId(); }

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

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_plhyperelast; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new PlasticHyperElast( *this ) ); }

    /// Young's modulus
    double Youngs() const { return params_->youngs_; }

    /// Poisson's ratio
    double PoissonRatio() const { return params_->poissonratio_; }

    /// density
    virtual double Density() const { return params_->density_; }

    /// shear modulus
    double ShearMod() const { return 0.5*params_->youngs_/(1.0+params_->poissonratio_); }

    /// yield stress
    virtual double Yield() const { return params_->yield_; }

    /// isotropic hardening modulus
    virtual double IsoHard() const { return params_->isohard_; }

    /// kinematic hardening modulus
    virtual double KinHard() const { return params_->kinhard_; }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    //@}

    //! @name Evaluation methods

    /// evaluate material law
    virtual void Evaluate(const LINALG::Matrix<3,3>* defgrd,
                  const LINALG::Matrix<6,1>* glstrain,
                  Teuchos::ParameterList& params,
                  LINALG::Matrix<6,1>* stress,
                  LINALG::Matrix<6,6>* cmat);

    /// evaluate strain energy function
    virtual void StrainEnergy(const LINALG::Matrix<6,1>& glstrain,
                      double& psi);

    //@}

  private:

    /// my material parameters
    MAT::PAR::PlasticHyperElast* params_;

  };
}

/*----------------------------------------------------------------------*/
#endif  // PLASTICHYPERELAST_H_
