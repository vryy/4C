/*!----------------------------------------------------------------------
\file viscoanisotropic.H
\brief

<pre>
Maintainer: Moritz Frenzel & Thomas Kloeppel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*----------------------------------------------------------------------*/
#ifndef VISCOANISOTROPIC_H_
#define VISCOANISOTROPIC_H_
#ifdef CCADISCRET

#include "material.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>
#include "../drt_lib/linalg_fixedsizematrix.H"

using namespace Teuchos;

namespace MAT
{

  /// Wrapper for Visco-NeoHooke material
  class ViscoAnisotropic : public Material
  {
  public:

    /// construct empty material object
    ViscoAnisotropic();

    /// construct the material object given material parameters
    explicit ViscoAnisotropic(MATERIAL* matdata);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_ViscoAnisotropic; }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.
      This material contains history variables, which are packed for restart purposes.

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
      History data is unpacked in restart.

      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    virtual void Unpack(const vector<char>& data);

    //@}

    /// material type
    virtual MATERIAL_TYP MaterialType() const { return m_viscoanisotropic; }

    /// Setup and Initialize internal stress variables
    void Setup
    (
      const int numgp  ///< number of Gauss points
    );

    /// Update internal stress variables
    void Update();

    /// Reset internal stress variables
    void Reset(){
      // do nothing,
      // because #histstresscurr_ and #artstresscurr_ are recomputed anyway at every iteration
      // based upon #histstresslast_ and #artstresslast_ untouched within time step
      return;
    };

    /// Evaluate material
    void Evaluate
    (
      const Epetra_SerialDenseVector* glstrain,  ///<green lagrange strain
      const int gp,   ///< number of Gauss points
      Teuchos::ParameterList& params,  ///< parameter list for communication
      Epetra_SerialDenseMatrix* cmat,  ///< material stiffness matrix
      Epetra_SerialDenseVector* stress  ///< 2nd PK-stress
    );

    /// Return density
    double Density()
    {
      return matdata_->m.viscoanisotropic->density;
    };

    /// Check if history variables are already initialized
    bool Initialized() const
      { return isinit_&&(histstresscurr_!=Teuchos::null);}

  private:

    /// my material parameters
    MATERIAL* matdata_;

    // visco history stresses for every gausspoint and every stress type
    Teuchos::RCP<vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> > > histstresscurr_; ///< current stress
    Teuchos::RCP<vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> > > histstresslast_; ///< stress of last converged state
    Teuchos::RCP<vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> > > artstresscurr_;  ///< current artificial stress
    Teuchos::RCP<vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> > > artstresslast_;  ///< artificial stress in last converged state

    bool isinit_;  ///< indicates if material is initialized

  };
}

#endif
#endif
