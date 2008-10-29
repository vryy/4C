/*----------------------------------------------------------------------*/
/*!
\file biocell.H
\brief 
<pre>
Maintainer: Robert Metzke
            metzke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef BIOCELL_H_
#define BIOCELL_H_

#ifdef CCADISCRET

#include "material.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{

  /// Wrapper for Biological Cell Material
  class BioCell : public Material
  {
  public:

    /// construct empty material object
    BioCell();

    /// construct the material object given material parameters
    explicit BioCell(MATERIAL* matdata);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_BioCell; }

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
    virtual MATERIAL_TYP MaterialType() const { return m_biocell; }
    
    void Evaluate(const LINALG::FixedSizeSerialDenseMatrix<6,1>* glstrain,
            LINALG::FixedSizeSerialDenseMatrix<6,6>* cmat,
            LINALG::FixedSizeSerialDenseMatrix<6,1>* stress);
    

    inline double Density() { return matdata_->m.biocell->density; }


  private:

    /// my material parameters
    MATERIAL* matdata_;
    
    void InverseTensor(
                  const Epetra_SerialDenseMatrix& M,
                  Epetra_SerialDenseMatrix& Minv,
		  const double I3);
    void TensorProduct (
                  const LINALG::SerialDenseMatrix& A,
                  const LINALG::SerialDenseMatrix& B,
                  LINALG::SerialDenseMatrix& AB);
  
  };
}

#endif /*CCADISCRET*/
#endif /*BIOCELL_H_*/
