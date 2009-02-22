/*!----------------------------------------------------------------------
\file charmm.H
\brief

<pre>
Maintainer: Robert Metzke
            metzke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
*----------------------------------------------------------------------*/
#ifndef CHARMM_H_
#define CHARMM_H_
#ifdef CCADISCRET

#include "material.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"


namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for neo-Hooke
    class CHARMM : public Parameter
    {
    public:

      /// standard constructor
      CHARMM(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~CHARMM() { ; }

      /// @name material parameters
      //@{

      /// mass density
      const double density_;

      //@}

    };  // class CHARMM

  }  // namespace PAR

  class CHARMM : public Material
  {
  public:

    /// construct empty material object
    CHARMM();

    /// construct the material object given material parameters
    explicit CHARMM(MAT::PAR::CHARMM* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_CHARMM; }

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
INPAR::MAT::m_charmm; }

    /// Initialize internal stress variables
    void Setup
    (
      DRT::Container& data_  ///< Container for History
    );

    void Evaluate(const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
                  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
                  LINALG::Matrix<NUM_STRESS_3D,1>* stress,
                  const int ele_ID,
                  const int gp,
                  DRT::Container& data_,
                  const double time,
                  const LINALG::SerialDenseMatrix& xrefe,
                  const LINALG::SerialDenseMatrix& xcurr);

    void EvalStrain( const bool& origin,
                  const LINALG::SerialDenseMatrix& xrefe,
                  const LINALG::SerialDenseMatrix& xcurr,
                  LINALG::Matrix<3,3>& C,
                  LINALG::SerialDenseVector& lambda );
    
    Epetra_SerialDenseVector charmmfileapi(
                  const double STARTD,
                  const double ENDD,
                  LINALG::SerialDenseVector direction,
                  LINALG::SerialDenseVector& charmm_result);
    
    Epetra_SerialDenseMatrix inverseTensor(
                  const Epetra_SerialDenseMatrix M,
		  const double I3);


    /// density
    inline double Density() const { return params_->density_; }


  private:

    /// my material parameters
    MAT::PAR::CHARMM* params_;
  
  };
}

#endif // CCADISCRET
#endif // CHARMM_H_

