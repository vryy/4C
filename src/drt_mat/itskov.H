
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

#ifdef CCADISCRET

#include "material.H"
#include "../drt_lib/drt_container.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{
  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for neo-Hooke
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

      /// mass density
      const double alpha_;
      const double beta_;
      const double mu_fibers_;
      const double mu_GS_;
      const double epsilon_;
      const double gamma_;
      const double comp_;
      const double density_;

      //@}

    };  // class CHARMM

  }  // namespace PAR
  
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
    virtual int UniqueParObjectId() const { return ParObject_Itskov; }

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
        virtual INPAR::MAT::MaterialType MaterialType() const { return
INPAR::MAT::m_itskov; }
    
    void Evaluate(const LINALG::Matrix<6,1>& glstrain,
    			const int gp, const int ele_ID, 
    			DRT::Container& data_,
    			const double time,
                  LINALG::Matrix<6,6>& cmat,
                  LINALG::Matrix<6,1>& stress);
                  
    /*void Evaluate(const LINALG::Matrix<6,1> *glstrain,
            		const int gp, const int ele_ID, const double time,
                  LINALG::Matrix<6,6> *cmat,
                  LINALG::Matrix<6,1> *stress);*/

    inline double Density() const { return params_->density_; }


  private:

    /// my material parameters
    MAT::PAR::Itskov* params_;

    void calc_factors (	LINALG::Matrix<3,3>&C, 
				LINALG::Matrix<3,3>&Cinv, 
				LINALG::Matrix<3,3>&L, 
				DOUBLE alpha, 
				DOUBLE beta, 
				DOUBLE mu, 
				DOUBLE *factors);

   void calc_Itskov (	LINALG::Matrix<3,3>&C, 
						LINALG::Matrix<3,3>&Cinv, 
						LINALG::Matrix<3,3>&L, 
						LINALG::Matrix<3,3>&PK2_, 
						LINALG::Matrix<6,6>&cmat_, 
						DOUBLE *factors);
  
  };
}

#endif /*CCADISCRET*/
#endif /*ITSKOV_H_*/
