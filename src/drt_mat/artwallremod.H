/*!----------------------------------------------------------------------
\file artwallremod.H
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*----------------------------------------------------------------------*/
#ifndef ARTWALLREMOD_H_
#define ARTWALLREMOD_H_
#ifdef CCADISCRET

#include "material.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_so3/so_hex8.H"
#include "../drt_so3/so_sh8.H"
#include "../drt_io/io_gmsh.H"
extern struct _FILES  allfiles;


using namespace Teuchos;

namespace MAT
{

  /// Wrapper for Mooney-Rivlin material
  class ArtWallRemod : public Material
  {
  public:

    /// construct empty material object
    ArtWallRemod();

    /// construct the material object given material parameters
    explicit ArtWallRemod(MATERIAL* matdata);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_ArtWallRemod; }

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
    virtual MATERIAL_TYP MaterialType() const { return m_artwallremod; }
    
    void Initialize(const int numgp, const int eleid) ;
    
    void Evaluate(const Epetra_SerialDenseVector* glstrain,
                  const int gp,
                  Teuchos::ParameterList& params,
                  Epetra_SerialDenseMatrix* cmat,
                  Epetra_SerialDenseVector* stress, 
                  int Id);

    double Density() {return matdata_->m.artwallremod->density;}
    
    bool Initialized() const
      { return isinit_&&(gamma_!=Teuchos::null);}
    
    RCP<vector<vector<double> > > Getlambdas() const {return lambda_;}
    

    void Update(const Epetra_SerialDenseMatrix& Phi,
        const Epetra_SerialDenseVector& lambda,
        const vector<double> rem_toggle,
        const int gp, const double r0, const double decay);

    
  private:

    /// my material parameters
    MATERIAL* matdata_;
    
    /// angle of alignment at gp
    Teuchos::RCP<vector<double> > gamma_;
    /// eigenvectors per gp
    Teuchos::RCP<vector<Epetra_SerialDenseMatrix> > phi_;
    /// eigenvalues per gp
    Teuchos::RCP<vector<vector<double> > > lambda_;
    /// current stress
    Teuchos::RCP<vector<Epetra_SerialDenseMatrix> > stresses_;
    
    // remember time
    Teuchos::RCP<vector<double> > mytime_;
    
    bool isinit_;
    
  };

}

#endif /*CCADISCRET*/
#endif /*ARTWALLREMOD_H_*/
