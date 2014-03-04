/*!----------------------------------------------------------------------
\file beam3tospherecontact.H

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>
*----------------------------------------------------------------------*/
#ifndef BEAMTOSPHEREC_H
#define BEAMTOSPHEREC_H

#include "beam3contact_defines.H"
#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "../drt_lib/drt_condition.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_sparsematrix.H"

namespace CONTACT
{

/*!
 \brief contact element for contact between a 3D beam element and a rigid sphere

 \author meier (meier@lnm.mw.tum.de)
 */


class Beam3tospherecontact
{
public:

  //! @name Friends
  
  //no friend classes defined
  
  //@}

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param pdiscret (in): the problem discretization
  \param cdiscret (in): the beam contact discretization
  \param dofoffset (in): offset of dof between pdiscret and cdiscret
  \param element1 (in): first element of contact pair
  \param element2 (in): second element of contact pair
  \param ele1pos (in): nodal coordinates of first element
  \param ele2pos (in): nodal coordinates of second element
  
  */
  Beam3tospherecontact(const DRT::Discretization& pdiscret,
               const DRT::Discretization& cdiscret,
               const int& dofoffset,
               DRT::Element* element1,
               DRT::Element* element2);

  /*!
  \brief Copy Constructor

  Makes a deep copy of this contact element pair
  
  */
  Beam3tospherecontact(const Beam3tospherecontact& old);

  /*!
  \brief Destructor
  
  */
  virtual ~Beam3tospherecontact() {};
  
  //@}
  
  //! @name Access methods
  
  /*!
  \brief Get problem discretization

  */
  inline const DRT::Discretization& ProblemDiscret() const { return pdiscret_;};

  /*!
  \brief Get beam contact discretization

  */
  inline const DRT::Discretization& ContactDiscret() const { return cdiscret_;};

  /*!
  \brief Get offset of dofs between cdiscret and pdiscret

  */
  inline const int& DofOffset() const { return dofoffset_;};

  /*!
  \brief Get first element
  
  */
  inline const DRT::Element* Element1() { return element1_;};
  
  /*!
  \brief Get first element
    
  */
  inline const DRT::Element* Element2() { return element2_;};

  /*!
  \brief Get gap of this contact pair
    
  */
  const double GetGap() { return gap_;};
  
  /*!
  \brief Get current augmented lagrange multiplier 'lmuzawa_'
  
  */
  const double Getlmuzawa() { return lmuzawa_; };
  
  /*! 
  \brief Get flag ndicating whether contact is active (true) or inactive (false)
  
  */
  bool GetContactFlag() { return contactflag_; };
    
  /*!
  \brief Get coordinates of contact point of element1 and element2
  
  */
  const Epetra_SerialDenseVector GetX1() { return x1_; };
  const Epetra_SerialDenseVector GetX2() { return x2_; };

  //@}
  
  //! @name Public evaluation methods
  
  /*!
  \brief Evaluate this contact element pair
  */
  bool Evaluate(LINALG::SparseMatrix& stiffmatrix, Epetra_Vector& fint, double& pp);
  
  /*!
  \brief Reset Uzawa-based Lagrange multiplier
  */
  void Resetlmuzawa();
  
  /*!
    \brief Update Uzawa-based Lagrange multiplier
  */
  void Updatelmuzawa(const double& currentpp);
  
  /* 
  \brief Update nodal coordinates of both elements
  */
  void UpdateElePos(Epetra_SerialDenseMatrix& newele1pos,
                    Epetra_SerialDenseMatrix& newele2pos);
  
  //@}

private:
  
  //! @name member variables

  //! reference to problem discretization
  const DRT::Discretization& pdiscret_;

  //! reference to beam contact discretization
  const DRT::Discretization& cdiscret_;
  
  //! dof offset between pdiscret and cdiscret
  const int& dofoffset_;

  //! first element of contact pair
  DRT::Element* element1_;
  
  //! second element of contact pair
  DRT::Element* element2_;
  
  //! current node coordinates of the two elements
  // TODO: Set to const if finite differences are not needed anymore
  Epetra_SerialDenseMatrix ele1pos_;
  Epetra_SerialDenseMatrix ele2pos_;

  //! variable to identify first call of a pair (for initializing)
  bool firstcall_;

  //! Uzawa-based Lagrange multiplier
  double lmuzawa_;

  //! gap function
  double gap_;
  
  //! flag indicating contact (active/inactive)
  bool contactflag_;

  //! coordinates of contact point on center lines of beams
  Epetra_SerialDenseVector x1_;
  Epetra_SerialDenseVector x2_;
  
  //! parameter value of contact point
  double xicontact_;

  //! normal vector of current time step
  Epetra_SerialDenseVector normal_;

  //@}
 
  //! @name Private evaluation methods
  
  /*!
  \brief Find contact point via closest point projection
  */
  void ClosestPointProjection();
  
  /*!
  \brief Utility method for CPP (evaluate nonlinear function f)
  */
  void EvaluateNewtonF(double& f,
                       const std::vector<double>& x1, const std::vector<double>& x2,
                       const std::vector<double>& dx1, const double& norm);
  
  /*!
  \brief Utility method for CPP (evaluate Jacobian of nonlinear function f)
  */
  void EvaluateNewtonGradF(double& df,
                           const std::vector<double>& x1, const std::vector<double>& x2,
                           const std::vector<double>& dx1,const std::vector<double>& ddx1,
                           const double& norm);
  
  /*!
  \brief Evaluate and assemble contact forces
  */
  void EvaluateFcContact(const double& pp, const double& gap,
                         const std::vector<double>& normal, Epetra_Vector& fint,
                         Epetra_SerialDenseVector funct1, const int numnode1);
  
  /*!
  \brief Evaluate and assemble contact stiffness
  */
  void EvaluateStiffcContact(const double& pp, const double& gap,
                             const std::vector<double>& normal, const double& norm,
                             LINALG::SparseMatrix& stiffmatrix,
                             const std::vector<double>& x1, const std::vector<double>& x2,
                             const std::vector<double>& dx1,const std::vector<double>& ddx1,
                             const Epetra_SerialDenseVector& funct1,
                             const Epetra_SerialDenseMatrix& deriv1,
                             const Epetra_SerialDenseMatrix& secondderiv1,
                             const int numnode1, double& XiContact);
  
  /*!
  \brief Compute normal vector in contact point
  */
  void ComputeNormal(std::vector<double>& normal, double& gap, double& norm, 
                     const std::vector<double>& x1, const std::vector<double>& x2);

  /*!
  \brief Evaluate gap function
  */
  void ComputeGap(double& gap, const double& norm);
 
  /*!
  \brief Compute radius of cross section based on moment of inertia
  */
  void ComputeEleRadius(double& radius, const double& moi);

  /*!
  \brief Compute coordinates and their derivatives from the discretization
  */
  void ComputeCoordsAndDerivs(std::vector<double>& x1, std::vector<double>& x2,
                              std::vector<double>& dx1,std::vector<double>& ddx1,
                              const Epetra_SerialDenseVector& funct1,
                              const Epetra_SerialDenseMatrix& deriv1,
                              const Epetra_SerialDenseMatrix& secondderiv1,
                              const int& numnode1);

  /*!
  \brief Get shape functions and their derivatives at eta
  */
  void GetShapeFunctions(Epetra_SerialDenseVector& funct1,
                         Epetra_SerialDenseMatrix& deriv1,
                         Epetra_SerialDenseMatrix& secondderiv1,
                         const double& eta);
  
  /*!
  \brief Compute linearizations of contact point
  */
  void ComputeLinXiAndLinEta(std::vector<double>& delta_xi,
                             const std::vector<double>& x1, const std::vector<double>& x2,
                             const std::vector<double>& dx1,
                             const std::vector<double>& ddx1,
                             const Epetra_SerialDenseVector& funct1,
                             const Epetra_SerialDenseMatrix& deriv1,
                             const std::vector<double>& normal, const double& norm,
                             const int numnode1, const double& XiContact);
  
  /*!
  \brief Compute linearization of gap
  */
  void ComputeLinGap(std::vector<double>& delta_gap,
                     std::vector<double>& delta_xi,
                     const std::vector<double>& x1, const std::vector<double>& x2,
                     const std::vector<double>& dx1,
                     const Epetra_SerialDenseVector& funct1,
                     const double& normdist,
                     const int& numnode1,
                     const std::vector<double>& normal, const double& norm,
                     const double& gap,Epetra_SerialDenseMatrix& delta_x1_minus_x2);
  
  /*!
  \brief Compute linearization of normal
  */
  void ComputeLinNormal(Epetra_SerialDenseMatrix& delta_n, 
                        const std::vector<double>& x1, const std::vector<double>& x2,
                        const double& norm,
                        const int& numnode1,
                        const Epetra_SerialDenseMatrix& delta_x1_minus_x2,
                        const std::vector<double>& normal, 
                        const double& XiContact);
  
  /*!
  \brief Compute normal contact disctance
  */                
  void ComputeDistance(std::vector<double>& distance, double& normdist,
                       const std::vector<double>& normal, const double& norm);
  
  /*!
  \brief Check if contact is active and set flag accordingly
  */
  void CheckContactStatus(double& pp);

  /*!
  \brief Get global dofs of a node

  Internally this method first extracts the dofs of the given node
  in the beam contact discretization (which has its own dofs) and
  then transfers these dofs to their actual GIDs in the underlying
  problem discretization by applying the pre-computed dofoffset_.

  */
  std::vector<int> GetGlobalDofs(DRT::Node* node);

  //@}

}; // class Beam3tospherecontact
} // namespace CONTACT

#endif  // #ifndef BEAMTOSPHEREC_H
