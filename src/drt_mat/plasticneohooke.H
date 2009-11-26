/*----------------------------------------------------------------------*/
/*!
\file plasticneohooke.H
\brief Contains the functions to establish local material law /
       stress-strain law for isotropic material for a 3D hex element
       following rate-independent von Mises plasticity with combined
       (nonlinear) isotropic and kinematic hardening and a nearly
       incompressible Neo Hookean material law.
       
       Refer also to the master thesis of A.A. van der Stelt (2009).

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef PLASTICNEOHOOKE_H_
#define PLASTICNEOHOOKE_H_

#ifdef CCADISCRET

#include "matpar_parameter.H"
#include "material.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

using namespace Teuchos;  

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for plastic neo-Hooke
    class PlasticNeoHooke : public Parameter
    {
    public:

      /// standard constructor
      PlasticNeoHooke(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~PlasticNeoHooke() { ; }

      /// @name material parameters
      //@{

      /// Young's modulus
      const double youngs_;
      /// Possion's ratio
      const double poissonratio_;
      /// mass density
      const double density_;
      /// isotropic hardening
      const double isohard_;
      /// yield stress
      const double yield_;
      /// inf yield stress for nonlinear iso hardening
      const double infyield_;
      /// exponent for nonlinear iso hardening
      const double exp_;
      /// kinematic hardening
      const double kinhard_;
      
      //@}

    };  // class PlasticNeoHooke

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Wrapper for Plastic Neo-Hookean Material
  class PlasticNeoHooke : public Material 
  {
  public:

    /// construct empty material object
    PlasticNeoHooke(); 

    /// construct the material object given material parameters
    explicit PlasticNeoHooke(MAT::PAR::PlasticNeoHooke* params); 

    //! @name Packing and Unpacking

    /*!
    \brief Return unique ParObject id

    every class implementing ParObject needs a unique id defined at the
    top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_PlasticNeoHooke; } 

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

    \param data (in) : vector storing all data to be unpacked into this instance.
    */
    virtual void Unpack(const vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_plneohooke; }
    
    /// Initialize internal stress variables
    void Setup
    (
      const int numgp ///< number of Gauss points
    );
    
    /// Update internal stress variables
    void Update();
    
    /// Reset internal stress variables
    void Reset();
    
    /// Evaluate material
    void Evaluate
    (
      const LINALG::Matrix<3,3>* defgrd,                  ///< deformation gradient 
      const int gp,                                       ///< current Gaus point
      Teuchos::ParameterList& params,                     ///< parameter list for communication 
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,  ///< material stiffness matrix
      LINALG::Matrix<NUM_STRESS_3D,1>* stress             ///< 2nd PK-stress
    );
    
    /// Return density
    double Density() const { return params_->density_; }

    /// Check if history variables are already initialized
    bool Initialized() const { return (isinit_ && (histplasticrcgcurr_!=Teuchos::null)); }

  private:

    /// my material parameters
    MAT::PAR::PlasticNeoHooke* params_;

    /// plastic history deformation
    Teuchos::RCP<vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histplasticrcgcurr_; /// current plastic right cauchy-green tensor
    Teuchos::RCP<vector<LINALG::Matrix<NUM_STRESS_3D,1> > > histplasticrcglast_; /// last converged state of plastic right cauchy-green tensor
    Teuchos::RCP<vector<LINALG::Matrix<1,1> > > histeplasticscurr_;              /// current equivalent plastic strain
    Teuchos::RCP<vector<LINALG::Matrix<1,1> > > histeplasticslast_;              /// last converged state of equivalent plastic strain
    Teuchos::RCP<vector<LINALG::Matrix<3,1> > > histplasticbstresscurr_;         /*added*/
    Teuchos::RCP<vector<LINALG::Matrix<3,1> > > histplasticbstresslast_;         /*added*/
    Teuchos::RCP<vector<LINALG::Matrix<3,1> > > lamdacurr_;                      /*added*/
    Teuchos::RCP<vector<LINALG::Matrix<3,1> > > lamdalast_;                      /*added*/
    Teuchos::RCP<vector<LINALG::Matrix<3,3> > > nsetcurr_;                       /*added*/
    Teuchos::RCP<vector<LINALG::Matrix<3,3> > > nsetlast_;                       /*added*/

    /// indicator if #Initialize routine has been called
    bool isinit_;
    
  };
}

#endif /*CCADISCRET*/
#endif /*PLASTICNEOHOOKE_H_*/
