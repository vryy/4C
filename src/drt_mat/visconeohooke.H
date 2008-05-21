/*!----------------------------------------------------------------------
\file visconeohooke.H
\brief

<pre>
Maintainer: Moritz Frenzel & Thomas Kloeppel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*----------------------------------------------------------------------*/
#ifndef VISCONEOHOOKE_H_
#define VISCONEOHOOKE_H_
#ifdef CCADISCRET

#include "material.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Teuchos_RefCountPtr.hpp>

using namespace Teuchos;

namespace MAT
{

  /// Wrapper for Mooney-Rivlin material
  class ViscoNeoHooke : public Material
  {
  public:

    /// construct empty material object
    ViscoNeoHooke();

    /// construct the material object given material parameters
    explicit ViscoNeoHooke(MATERIAL* matdata);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_ViscoNeoHooke; }

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
    virtual MATERIAL_TYP MaterialType() const { return m_visconeohooke; }
    
    void Initialize(const int numgp) ;
    
    void Update();
    
    void Evaluate(const Epetra_SerialDenseVector* glstrain,
                  const int gp,
                  const double dt,
                  Epetra_SerialDenseMatrix* cmat,
                  Epetra_SerialDenseVector* stress);

    double Density();
    
    bool Initialized() { return (histstresscurr_!=Teuchos::null);}
    
  private:

    /// my material parameters
    MATERIAL* matdata_;
    
    /// visco history stresses
    Teuchos::RCP<vector<Epetra_SerialDenseVector> > histstresscurr_;
    Teuchos::RCP<vector<Epetra_SerialDenseVector> > histstresslast_;
    Teuchos::RCP<vector<Epetra_SerialDenseVector> > artstresscurr_;
    Teuchos::RCP<vector<Epetra_SerialDenseVector> > artstresslast_;
    
    bool isinit_;
    
  };
}

#endif /*CCADISCRET*/
#endif /*MOONEYRIVLIN_H_*/
