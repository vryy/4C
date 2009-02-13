/*!----------------------------------------------------------------------
\file lung_ogden.H
\brief
This file contains the routines required for isotropic nearly
incompressible lung tissue

The input line should read

MAT 1 MAT_Struct_LungOgden C 1000. K1 1500. K2 8.5 KAPPA 1000 BETA -2.0 DENS 0.001

<pre>
Maintainer: Lena Wiechert & Sophie Rausch
            {wiechert,rausch}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*----------------------------------------------------------------------*/
#ifndef LUNDOGDEN_H_
#define LUNGOGDEN_H_
#ifdef CCADISCRET

#include "material.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for lung Ogden
    class LungOgden : public Parameter
    {
    public:

      /// standard constructor
      LungOgden(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~LungOgden() { ; }

      /// @name material parameters (these need to be writable for inverse analysis)
      //@{

      double kappa_;
      double beta_;
      double c_;
      double k1_;
      double k2_;
      double dens_;

      //@}

    };  // class LungOgden

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  class LungOgden : public Material
  {

  public:

    // empty constructor
    LungOgden();

    // constructor with given material parameters
    LungOgden(MAT::PAR::LungOgden* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_LungOgden; }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    virtual void Pack(vector<char>& data) const;

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
    virtual void Unpack(const vector<char>& data);

    //@}

    /// @name Access methods
    //@{
    double Density() const { return params_->dens_; }
    double C() const { return params_->c_; }
    double K1() const { return params_->k1_; }
    double K2() const { return params_->k2_; }
    double Kappa() const { return params_->kappa_; }
    double Beta() const { return params_->beta_; }

    //@}

    /// @name Write methods
    //@{
    void SetC(const double c) { params_->c_ = c; }
    void SetK1(const double k1) { params_->k1_ = k1; }
    void SetK2(const double k2) { params_->k2_ = k2; }
    //@}

    // material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_lung_ogden; }

    // THE material routine
    void Evaluate(LINALG::Matrix<6,1>* glstrain,
		  LINALG::Matrix<6,6>* cmat,
		  LINALG::Matrix<6,1>* stress);

  private:

    /// my material parameters
    MAT::PAR::LungOgden* params_;
  };
}

#endif  // CCADISCRET
#endif  // LUNGOGDEN_H_
