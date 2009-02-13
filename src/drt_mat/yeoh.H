/*!----------------------------------------------------------------------
\file yeoh.H
\brief

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*----------------------------------------------------------------------*/
#ifndef YEOH_H_
#define YEOH_H_
#ifdef CCADISCRET

#include "material.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{
  
  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters
    class Yeoh : public Parameter
    {
    public:

      /// standard constructor
      Yeoh(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~Yeoh() { ; }

      /// @name material parameters
      //@{
      const double c1_;
      const double c2_;
      const double c3_;
      const double kap_;
      const double density_;
      //@}

    };  // class Yeoh

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Wrapper for Yeoh material
  class Yeoh : public Material
  {
  public:

    /// construct empty material object
    Yeoh();

    /// construct the material object given material parameters
    explicit Yeoh(MAT::PAR::Yeoh* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_Yeoh; }

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

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_yeoh; }

    /*!
    \brief Yeoh type nearly incompressible, hyperelastic 3D material law.

    The underlying strain-energy function in decoupled form is given by:
    
    W = c_1 (I1_bar - 3)  +  c_2 (I1_bar - 3)^2  +  c_3 (I1_bar-3)^3 + kappa * (J - 1)^2, 
    where I1_bar is the first invariant of the modified Cauchy-Green tensor C_bar = J^{-2/3} C.
    
    
    For references see Holzapfel p. 248f, 255, 263

    */
    void Evaluate(
            const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
            LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> * cmat,
            LINALG::Matrix<NUM_STRESS_3D,1> * stress);

    /// return density
    double Density() const { return params_->density_; }

  private:

    /// my material parameters
    MAT::PAR::Yeoh* params_;

  };
}

#endif /*CCADISCRET*/
#endif /*YEOH_H_*/
