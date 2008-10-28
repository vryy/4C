/*!----------------------------------------------------------------------
\file drt_locsys.H

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef DRT_LOCSYS_H
#define DRT_LOCSYS_H

#include "Teuchos_RefCountPtr.hpp"
#include "../drt_lib/drt_discret.H"

using namespace std;
using namespace Teuchos;

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
  \brief Standard Constructor

  \param discret (in): A discretization containing locsys boundary conditions
  
  */
  explicit LocsysManager(DRT::Discretization& discret);

  /*!
  \brief Destructor

  */
  virtual ~LocsysManager() {};

  /*!
  \brief Print this Manager

  */
  virtual void Print(ostream& os) const;

  /*!
  \brief Get Epetra communicator

  */
  inline const Epetra_Comm& Comm() const { return discret_.Comm();};
    
  //! @name Access methods
  
  /*!
  \brief Get discretization

  */
  inline DRT::Discretization& Discret() const { return discret_;};
  
  /*!
  \brief Get problem dimension
    
  */
  inline const int& Dim() {return dim_;};
    
  /*!
  \brief Get local system conditions

  */
  inline vector<DRT::Condition*> Conditions() const { return locsysconds_;};
  
  /*!
  \brief Get a specific local system condition

  */
  inline DRT::Condition* Conditions(int k) const
  { if (k>=numlocsys_) dserror("Invalid vector index");
    else return locsysconds_[k];};
                                                   
  /*!
  \brief Get number of local system conditions

  */
  inline int NumLocsys() const { return numlocsys_;};
  
  /*!
  \brief Get types of local system conditions

  */
  inline vector<DRT::Condition::ConditionType> TypeLocsys() const { return typelocsys_;};
  
  /*!
  \brief Get type of a specific local system condition

  */
  inline DRT::Condition::ConditionType TypeLocsys(int k) const
  { if (k>=numlocsys_) dserror("Invalid vector index");
    return typelocsys_[k];};
    
  //@}
  
  //! @name Evaluation methods
  
  /*!
  \brief Apply forward transformation of linear system of equations

  */
  void RotateGlobalToLocal(RCP<LINALG::SparseMatrix> sysmat, RCP<Epetra_Vector> rhs);
  
  /*!
  \brief Apply forward transformation of a single vector

  */
  void RotateGlobalToLocal(RCP<Epetra_Vector> vec);
    
  /*!
  \brief Apply backward transformation of result and linear system of equations

  */
  void RotateLocalToGlobal(RCP<Epetra_Vector> result, RCP<LINALG::SparseMatrix> sysmat, RCP<Epetra_Vector> rhs);
  
  /*!
  \brief Apply backward transformation of a single vector

  */
  void RotateLocalToGlobal(RCP<Epetra_Vector> vec);
    
  //@}
  
  
protected:

  // don't want = operator and cctor
  LocsysManager operator = (const LocsysManager& old);
  LocsysManager(const LocsysManager& old);

  //! name@ Private Attributes
  
  // current discretization
  DRT::Discretization&                     discret_;
  
  // problem dimension
  int                                      dim_;
  
  // local system conditions
  vector<DRT::Condition*>                  locsysconds_;
  
  // number of local systems
  int                                      numlocsys_;
  
  // types of local systems
  vector<DRT::Condition::ConditionType>    typelocsys_;
  
  // local system normals / tangents
  Epetra_SerialDenseMatrix                 normals_;
  Epetra_SerialDenseMatrix                 tangents_;
  Epetra_SerialDenseMatrix                 thirddir_;
  
  // assignment of local systems to nodes
  RCP<Epetra_Vector>                       locsystoggle_;
  
  // transformation matrix
  RCP<LINALG::SparseMatrix>                trafo_;
  
  //@}
  
}; // class LocsysManager
} // namespace UTILS
} // namespace DRT

// << operator
ostream& operator << (ostream& os, const DRT::UTILS::LocsysManager& manager);


#endif  // #ifndef DRT_LOCSYS_H
#endif  // #ifdef CCADISCRET
