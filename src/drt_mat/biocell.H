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

#include "matpar_parameter.H"
#include "material.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>


namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters
    class BioCell : public Parameter
    {
    public:

      /// standard constructor
      BioCell(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~BioCell() { ; }

      /// @name material parameters
      //@{

      const double density_;

      //@}

    };  // class BioCell

  }  // namespace PAR

  /*----------------------------------------------------------------------*/
  /// Wrapper for Biological Cell Material
  class BioCell : public Material
  {
  public:

    /// construct empty material object
    BioCell();

    /// construct the material object given material parameters
    explicit BioCell(MAT::PAR::BioCell* params);

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_biocell; }
    
    void Evaluate(const LINALG::Matrix<6,1>* glstrain,
            LINALG::Matrix<6,6>* cmat,
            LINALG::Matrix<6,1>* stress);
    

    inline double Density() const { return params_->density_; }


  private:

    /// my material parameters
    MAT::PAR::BioCell* params_;
    
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
