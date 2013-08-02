/*----------------------------------------------------------------------*/
/*!
\file plasticneohooke.H
\brief Containing all plastic material parameters needed for
       the semi-smooth Newton method for plasticity.
       This is NOT a stand-alone material with proper evaluate
       routines. A hyperelastic material is needed to desribe
       the elastic part.


<pre>
Maintainer: Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef PLASTICSEMISMOOTH_H_
#define PLASTICSEMISMOOTH_H_


#include "matpar_parameter.H"
//#include "material.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Teuchos_ParameterList.hpp>


namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for de St. Venant--Kirchhoff
    class PlasticSemiSmooth : public Parameter
    {
    public:

      /// standard constructor
    	PlasticSemiSmooth(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~PlasticSemiSmooth() { ; }

      /// @name material parameters
      //@{

  	/// initial yield stress
  	const double yield_;
  	/// linear isotropic hardening
  	const double isohard_;
  	/// exponent for nonlinear iso hardening
  	const double expisohard_;
  	/// inf yield stress for nonlinear iso hardening
  	const double infyield_;
  	/// linear kinematic hardening
  	const double kinhard_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class StVenantKirchhoff

  }  // namespace PAR

class PlasticSemiSmoothType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "PlasticSemiSmoothType"; }

  static PlasticSemiSmoothType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static PlasticSemiSmoothType instance_;
};

  /*----------------------------------------------------------------------*/
  /// Wrapper for St.-Venant-Kirchhoff material
  class PlasticSemiSmooth : public So3Material
  {
  public:

    /// construct empty material object
	  PlasticSemiSmooth();

    /// construct the material object given material parameters
    explicit PlasticSemiSmooth(MAT::PAR::PlasticSemiSmooth* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return PlasticSemiSmoothType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_plsemismooth; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new PlasticSemiSmooth( *this ) ); }


    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    // get linear isotropic hardening coefficient
    double InitYield() {return params_->yield_;}

    // get linear isotropic hardening coefficient
    double IsoHard() {return params_->isohard_;}

    // get nonlinear isotropic hardening exponent
    double ExpIsoHard() {return params_->expisohard_;}

    // get linear isotropic hardening coefficient
    double InfYield() {return params_->infyield_;}

    // get linear kinematic hardening coefficient
    double KinHard() {return params_->kinhard_;}
    //@}

    //! @name Evaluation methods

    /// evaluate material law
    virtual void Evaluate(const LINALG::Matrix<3,3>* defgrd,
                          const LINALG::Matrix<6,1>* glstrain,
                          Teuchos::ParameterList& params,
                          LINALG::Matrix<6,1>* stress,
                          LINALG::Matrix<6,6>* cmat)
    {dserror("no evaluate routine for this material. It should come along with a hyperelastic material");return;}

    //@}

  private:

    /// my material parameters
    MAT::PAR::PlasticSemiSmooth* params_;

  };
}

#endif
