
/*----------------------------------------------------------------------*/
/*!
\file Itskov.H
\brief
This file contains the routines required for incompressible Itskov material law
with a penalty-function acocrding to Balzani et al.

Itskov, M.; Ehret, A. & Mavrilas, D.
A polyconvex anisotropic strain energy function for soft collagenous tissues
Biomechanics and Modeling in Mechanobiology, 2006, 5, 17-26

Ehret, A. & Itskov, M.
A polyconvex hyperelastic model for fiber-reinforced materials in application to soft tissues
Journal of Materials Science, 2007, 42, 8853-8863


The input line should read:
MAT 1 MAT_ITSKOV  MU_GS 0.0 MU_FIBERS 2.026255352 ALPHA 20.0 BETA 20.0 EPSILON 100.0 GAMMA 10.0 C 0.0 DENS 1.0E-6

C gives the desired deviation from full incompressibility. If the Value for C
is 0.0, the Penaltyparameters EPSILON and GAMMA are used as given in the input
file. Otherwise the Penaltyparameter GAMMA is used as given and EPSILON is
computed to fullfill the desired rate of incompressibility.
<pre>
Maintainer: Robert Metzke
            metzke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef ITSKOV_H_
#define ITSKOV_H_


#include "material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_container.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{
  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /*!
	\brief material parameters for Itskov material law

	 MAT 1 MAT_ITSKOV  MU_GS 0.0 MU_FIBERS 2.026255352 ALPHA 20.0 BETA 20.0 EPSILON 100.0 GAMMA 10.0 C 0.0 DENS 1.0E-6

	*/
    class Itskov : public Parameter
    {
    public:

      /// standard constructor
      Itskov(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~Itskov() { ; }

      /// @name material parameters
      //@{

      /// alpha parameter fibers
      const double alpha_;
	  /// beta parameter fibers
      const double beta_;
	  /// mu fibers
      const double mu_fibers_;
	  /// mu ground substance
      const double mu_GS_;
	  /// penalty parameter epsilon
      const double epsilon_;
	  /// penalty parameter gamma
      const double gamma_;
      /// I3 deviation parameter
	  const double comp_;
      /// mass density
      const double density_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class CHARMM

  }  // namespace PAR

class ItskovType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "ItskovType"; }

  static ItskovType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static ItskovType instance_;
};

  /*----------------------------------------------------------------------*/
  /*!
  \brief polyconvex Itskov material law including adaptive penalty approach

  Here the polyconvex framework on strain and tension of Itskov is implemented using
   the same material formulation as for the Holzapfel lung model. Meaning that there
   exists some ground substane in a Neo-Hookean formulation and a exponential function
   to capture the nature of collagen. Both formulations are isotropic. Anisotropic
   formulations exist as external code, but have been not included. Feel free to pop by
   for a coffee to get them.

  Adaptive Penaly Approach: In order to avoid problems in the high strain regime (and
   therefore the collagen dominant strain regime) a new method was developed to adaptiveley
   change the penalty parameters according to a maximum deviation from the ideal compressible
   case.

  The input line should read like that:

  MAT 1 MAT_ITSKOV  MU_GS 0.0 MU_FIBERS 2.026255352 ALPHA 20.0 BETA 20.0 EPSILON 100.0 GAMMA 10.0 C 0.0 DENS 1.0E-6

  References:

  Itskov, M.; Ehret, A. & Mavrilas, D.
  A polyconvex anisotropic strain energy function for soft collagenous tissues
  Biomechanics and Modeling in Mechanobiology, 2006, 5, 17-26

  Ehret, A. & Itskov, M.
  A polyconvex hyperelastic model for fiber-reinforced materials in application to soft tissues
  Journal of Materials Science, 2007, 42, 8853-8863

  \author Robert Metzke
  \date 03/09
  */
  class Itskov : public Material
  {
  public:

    /// construct empty material object
    Itskov();

    /// construct the material object given material parameters
    explicit Itskov(MAT::PAR::Itskov* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ItskovType::Instance().UniqueParObjectId(); }

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

    /// material type
	virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_itskov; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new Itskov( *this ) ); }

	/// Computation of stresses and constitutive tensor
    void Evaluate(const LINALG::Matrix<6,1>& glstrain,
    			const int gp, const int ele_ID,
    			DRT::Container& data_,
    			const double time,
				LINALG::Matrix<6,6>& cmat,
				LINALG::Matrix<6,1>& stress);

    virtual double Density() const { return params_->density_; }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

  private:

    /// Itskov material parameters
    MAT::PAR::Itskov* params_;

    void calc_factors (	LINALG::Matrix<3,3>&C,
				LINALG::Matrix<3,3>&Cinv,
				LINALG::Matrix<3,3>&L,
				double alpha,
				double beta,
				double mu,
				double *factors);

    void calc_Itskov (	LINALG::Matrix<3,3>&C,
						LINALG::Matrix<3,3>&Cinv,
						LINALG::Matrix<3,3>&L,
						LINALG::Matrix<3,3>&PK2_,
						LINALG::Matrix<6,6>&cmat_,
						double *factors);

  };
}

#endif /*ITSKOV_H_*/
