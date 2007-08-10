#ifndef STVENANTKIRCHHOFF_H_
#define STVENANTKIRCHHOFF_H_

#ifdef CCADISCRET

#include "material.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{

  /// Wrapper for St.-Venant-Kirchhoff material
  class StVenantKirchhoff : public Material
  {
  public:

    /// construct empty material object
    StVenantKirchhoff();

    /// construct the material object given material parameters
    explicit StVenantKirchhoff(MATERIAL* matdata);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_StVenantKirchhoff; }

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
    virtual MATERIAL_TYP MaterialType() const { return m_stvenant; }
    
    //evaluates stresses for 3d
    void Evaluate(const Epetra_SerialDenseVector* glstrain,
                  Epetra_SerialDenseMatrix* cmat,
                  Epetra_SerialDenseVector* stress);
                  
    // computes isotropic elasticity tensor in matrix notion for 3d
    void SetupCmat(Epetra_SerialDenseMatrix* cmat);
    
    // computes isotropic elasticity tensor in matrix notion for 2d
    // plane strain, rotational symmetry
    void SetupCmat2d(Epetra_SerialDenseMatrix* cmat);

    double Density();                  


  private:

    /// my material parameters
    MATERIAL* matdata_;
  
  };
}

#endif /*CCADISCRET*/
#endif /*STVENANTKIRCHHOFF_H_*/
