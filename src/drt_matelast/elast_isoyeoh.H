/*----------------------------------------------------------------------*/
/*!
\file elast_isoyeoh.H
\brief
This file contains the routines required to calculate the isochoric contribution 
of a Yeoh-type material

<pre>
Maintainer: Sophie Rausch & Thomas Kloeppel
            {rausch,kloeppel}@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
#ifndef ELAST_ISOYEOH_H_
#define ELAST_ISOYEOH_H_
#ifdef CCADISCRET

#include "../drt_inpar/inpar_material.H"
#include "../drt_mat/matpar_parameter.H"
#include "elast_summand.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{
namespace ELAST
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for isochoric contribution of a neo-Hooke material
    ///
    /// <h3>Input line</h3>
    /// MAT 1 ELAST_IsoYeoh C1 100 C2 0 C3 200 
    class IsoYeoh : public MAT::PAR::Parameter
    {
    public:

      /// standard constructor
      IsoYeoh(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~IsoYeoh() { ; }

      /// @name material parameters
      //@{

      /// Shear modulus
      const double c1_;
      const double c2_;
      const double c3_;

      //@}

    };  // class IsoYeoh

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Isochoric neo-Hooke material according to [1]
  ///
  /// This is a compressible, hyperelastic, isotropic material 
  /// of the most simple kind.
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] J Bonet and RD Wood, "Nonlinear continuum mechanics for finite 
  ///          element analysis", Cambridge, 1997.
  /// <li> [2] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
  /// </ul>
  ///
  /// \author bborn
  /// \date 04/09
  class IsoYeoh : public Summand
  {

  public:

    /// empty constructor
    IsoYeoh();

    /// constructor with given material parameters
    IsoYeoh(MAT::ELAST::PAR::IsoYeoh* params);

    /// @name Packing and Unpacking
    //@{

    /// Return unique ParObject id
    ///
    ///  every class implementing ParObject needs a unique id defined at the
    ///  top of drt_parobject.H (this file) and should return it in this method.
//    virtual int UniqueParObjectId() const { return ParObject_IsoYeoh; }

    /// Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by UniqueParObjectId() which will then
    /// identify the exact class on the receiving processor.
//    virtual void Pack(
//      std::vector<char>& data  ///< (in/out): char vector to store class information
//      ) const;

    /// \brief Unpack data from a char vector into this class
    ///
    /// The vector data contains all information to rebuild the
    /// exact copy of an instance of a class on a different processor.
    /// The first entry in data has to be an integer which is the unique
    /// parobject id defined at the top of this file and delivered by
    /// UniqueParObjectId().
    ///
//    virtual void Unpack(
//      const std::vector<char>& data  ///< vector storing all data to be unpacked into this
//      );

    //@}

    /// @name Access material constants
    //@{

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::mes_isoyeoh; }

    //@}

    /// Add principal coeffiencts.
    /// Nothing to do for this (decoupled) material!
    void AddCoefficientsPrincipal(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<8,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& prinv  ///< principal invariants of right Cauchy-Green tensor
    )
    {
      havecoefficients = havecoefficients or false;

      return;
    };
    
    /// Add modified coeffiencts.
    void AddCoefficientsModified(
        bool& havecoefficients,  ///< non-zero coefficients of this type are provided
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<5,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& modinv  ///< principal invariants of right Cauchy-Green tensor
    );

  private:

    /// my material parameters
    MAT::ELAST::PAR::IsoYeoh* params_;

  };

}  // namespace ELAST
}  // namespace MAT

#endif  // CCADISCRET
#endif  // LOGNEOHOOKE_H_
