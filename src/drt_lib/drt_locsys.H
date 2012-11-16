/*----------------------------------------------------------------------*/
/*!
\file drt_locsys.H

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef DRT_LOCSYS_H
#define DRT_LOCSYS_H

#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_discret.H"

namespace LINALG {
  class SparseMatrix;
}

namespace DRT {
namespace UTILS {

/*!
\brief Class controlling local coordinate systems on points, lines and surfaces
       and supplying all necessary transformation methods for parallel vectors
       and matrices.

\author popp (popp@lnm.mw.tum.de)
*/
class LocsysManager
{
public:

  //! @name Enums and Friends

  //@}


  /*!
   * \brief Standard Constructor
   *
   * \param discret (in): A discretization containing locsys boundary conditions
   * \param transformleftonly (in): Only a tranformation from the left is going
   *                                to be applied on sysmatrix, if true
   *
   */
  explicit LocsysManager(DRT::Discretization& discret, const bool transformleftonly = false);

  /*!
   * \brief Set-up stuff common to all constructors
   */
  void Setup();

  /*!
   *\brief Destructor
   *
   */
  virtual ~LocsysManager() {};

  /*!
   *\brief Print this Manager
   *
   */
  virtual void Print(ostream& os) const;

  /*!
   * \brief Get Epetra communicator
   *
   */
  inline const Epetra_Comm& Comm() const { return discret_.Comm();};

  //! @name Access methods

  /*!
   * \brief Get discretization
   *
   */
  inline DRT::Discretization& Discret() const { return discret_;};

  /*!
   * \brief Get problem dimension
   *
   */
  inline const int& Dim() {return dim_;};

  /*!
   * \brief Get local system conditions
   *
   */
  inline vector<DRT::Condition*> Conditions() const { return locsysconds_;};

  /*!
   * \brief Get a specific local system condition
   *
   */
  inline DRT::Condition* Conditions(int k) const
  {
    if (k>=numlocsys_)
    {
      dserror("Invalid vector index");
      return NULL;
    }
    else
    {
      return locsysconds_[k];
    }
  };

  /*!
   * \brief Get number of local system conditions
   *
   */
  inline int NumLocsys() const { return numlocsys_;};

  /*!
   * \brief Get types of local system conditions
   *
   */
  inline vector<DRT::Condition::ConditionType> TypeLocsys() const { return typelocsys_;};

  /*!
   * \brief Get type of a specific local system condition
   *
   */
  inline DRT::Condition::ConditionType TypeLocsys(int k) const
  { if (k>=numlocsys_) dserror("Invalid vector index");
    return typelocsys_[k];};

  /*!
   * \brief Retrieve the global transformation matrix
   *
   */
  Teuchos::RCP<const LINALG::SparseMatrix> Trafo() const { return trafo_; }

  //@}

  //! @name Evaluation methods

  /*!
   * \brief Apply forward transformation of linear system of equations
   *
   * This method transform the matrix #sysmat from global co-ordinate
   * into local co-ordinate systems, i.e.
   *   \f[ \tilde{K} = Q \cdot K \cdot Q^T \f]
   * in which \f$K\f$ is the globally  and \f$\tilde{K}\f$ the locally oriented matrix,
   * respectively. The transformation matrix #trafo_ is denoted \f$Q : D \mapsto \tilde{D}\f$.
   * The similar thing is done for the right-hand-side vector:
   *   \f[ \tilde{R} = Q \cdot R \f]
   *
   * In case of #transformleftonly_==true the matrix is only transformed from
   * the left. This yields
   *   \f[ \hat{K} = Q \cdot K \f]
   */
  void RotateGlobalToLocal(
    RCP<LINALG::SparseMatrix> sysmat,  ///< systemmatrix, will be transformed
    RCP<Epetra_Vector> rhs  ///< right-hand-side vector, will be transformed
    );

  /*!
   * \brief Apply forward transformation of a single matrix
   *
   */
  void RotateGlobalToLocal(RCP<LINALG::SparseMatrix> sysmat);

  /*!
   * \brief Apply forward transformation of a single vector
   *
   */
  void RotateGlobalToLocal(RCP<Epetra_Vector> vec);

  /*!
   * \brief Apply backward transformation of result and linear system of equations
   *
   */
  void RotateLocalToGlobal(RCP<Epetra_Vector> result, RCP<LINALG::SparseMatrix> sysmat, RCP<Epetra_Vector> rhs);

  /*!
   * \brief Apply backward transformation of a single vector
   *
   */
  void RotateLocalToGlobal(RCP<Epetra_Vector> vec);

  /*!
	 * \brief Apply backward transformation of a matrix
	 *
	 */
  void RotateLocalToGlobal(RCP<LINALG::SparseMatrix> sysmat);

  //@}

protected:

  // don't want = operator and cctor
  LocsysManager operator = (const LocsysManager& old);
  LocsysManager(const LocsysManager& old);

  /// @name Methods for Local coordinate systems varying in space

  /*!

  \brief Construct rotation matrix around given axis

  \param x     (i) first coordinate of rotation axis
  \param y     (i) second coordinate of rotation axis
  \param z     (i) third coordinate of rotation axis
  \param angle (i) rotation angle

  */
  void SetAxisRotation(const double x,const double y,const double z,const double angle)
  {
    // norm axis vector

    const double norm=sqrt(x*x+y*y+z*z);

    if(norm<1e-9)
    {
      dserror("rotation axis specified incorrectly (length near zero)\n");
    }

    LINALG::Matrix<3,1> n;
    n(0)=x/norm;
    n(1)=y/norm;
    n(2)=z/norm;

    // set up rotation
    const double sa  =sin(angle);
    const double ca  =cos(angle);
    const double omca=1-ca;

    /*

 +--                                            |                                             |                                             --+
 |                                              |                                             |                                               |
 |             cos(a)    + n_x*n_x*(1.0-cos(a)) | n_x*n_y*(1.0-cos(a)) +     n_z*sin(a)       | n_x*n_z*(1.0-cos(a)) +     n_y*sin(a)         |
 |                                              |                                             |                                               |
 |                                              |                                             |                                               |
 |  n_y*n_x*(1.0-cos(a)) +     n_z*sin(a)       |            cos(a)    + n_y*n_y*(1.0-cos(a)) | n_y*n_z*(1.0-cos(a)) -     n_x*sin(a)         |
 |                                              |                                             |                                               |
 |                                              |                                             |                                               |
 |  n_z*n_x*(1.0-cos(a)) -     n_y*sin(a)       | n_z*n_y*(1.0-cos(a)) -     n_x*sin(a)       |            cos(a)    + n_z*n_z*(1.0-cos(a))   |
 |                                              |                                             |                                               |
 +--                                            |                                             |                                             --+

    */

    R_(0,0)=            ca + n(0)*n(0)*omca;
    R_(1,0)=n(1)*n(0)*omca +      n(2)*sa;
    R_(2,0)=n(2)*n(0)*omca -      n(1)*sa;

    R_(0,1)=n(0)*n(1)*omca -      n(2)*sa;
    R_(1,1)=            ca + n(1)*n(1)*omca;
    R_(2,1)=n(2)*n(1)*omca +      n(0)*sa;

    R_(0,2)=n(0)*n(2)*omca +      n(1)*sa;
    R_(1,2)=n(1)*n(2)*omca -      n(0)*sa;
    R_(2,2)=            ca + n(2)*n(2)*omca;

    return;
  }
  /*
  \brief Rotate a given vector using rotation matrix R_

  \param vec (i/o) vector to rotate
  */
  void LocalRotation(LINALG::Matrix<3,1> & vec)
  {
    LINALG::Matrix<3,1> old(vec);

    vec.Multiply(R_,old);

    return;
  }
  // @}

  /// @name Private Attributes
  //@{

  /// current discretization
  DRT::Discretization&                     discret_;

  /// transform system matrix only by left-multiplying
  bool                                     transformleftonly_;

  /// type of locsys
  enum Type
  {
    def,
    functionevaluation,
    originradialsliding
  };
  std::vector<DRT::UTILS::LocsysManager::Type>  type_;

  /// id of locsys condition
  std::vector<int>                         id_;

  /// problem dimension
  int                                      dim_;

  /// local system conditions
  std::vector<DRT::Condition*>             locsysconds_;

  /// number of local systems
  int                                      numlocsys_;

  /// Rotation matrix for coordinate systems varying in space
  LINALG::Matrix<3,3>                      R_;


  /// types of local systems conditions
  std::vector<DRT::Condition::ConditionType> typelocsys_;

  /// local system normals / tangents
  Epetra_SerialDenseMatrix                 normals_;
  Epetra_SerialDenseMatrix                 tangents_;
  Epetra_SerialDenseMatrix								 origins_;
  Epetra_SerialDenseMatrix                 thirddir_;

  /// assignment of local systems to nodes
  RCP<Epetra_Vector>                       locsystoggle_;

  /// maps containing the DOFs affected by locsys
  RCP<Epetra_Map>                          locsysdofmap_;

  /// store OriginRadialSliding locsys IDs
  std::vector<int>                         radslideids_;

  /// Transformation matrix which maps globally oriented components
  /// into locally oriented components (dubbed 'forward' transformation)
  ///
  /// Even if this matrix bears the general name transformation matrix here,
  /// one should be aware of that it is treated like a rotational
  /// transformation matrix. This means: The inverse (or 'backward') mapping
  /// is simply implemented as its transpose rather than a general inverse.
  RCP<LINALG::SparseMatrix>                trafo_;

  /// Transformation 'sub'-matrix with non-identity entries
  ///
  /// This is actually not a sub-matrix, but a global matrix
  /// with nil entries at (a lot) places
  RCP<LINALG::SparseMatrix>                subtrafo_;

  //@}

}; // class LocsysManager
} // namespace UTILS
} // namespace DRT

// << operator
ostream& operator << (ostream& os, const DRT::UTILS::LocsysManager& manager);


#endif  // #ifndef DRT_LOCSYS_H
