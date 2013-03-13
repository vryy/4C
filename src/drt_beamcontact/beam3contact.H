/*!----------------------------------------------------------------------
\file beam3contact.H

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>
*----------------------------------------------------------------------*/
#ifndef BEAM3C_H
#define BEAM3C_H

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

using namespace Teuchos;

namespace CONTACT
{

/*!
 \brief contact element for contact between two 3D beam elements

 Refer also to the Semesterarbeit of Matthias Mayr, 2010

 \author popp (popp@lnm.mw.tum.de), cyron (cyron@lnm.mw.tum.de)
 */


class B3CNeighbor
{
public:
  DRT::Element* left_neighbor_;      // pointer to the Element on the left side (eta=-1)
  DRT::Element* right_neighbor_;      // pointer to the Element on the right side (eta=1)
  int connecting_node_left_;  // local node-ID of the connecting node of the left neighbor
  int connecting_node_right_;  // local node-ID of the connecting node of the right neighbor


  B3CNeighbor(); // standard constructor


   virtual ~B3CNeighbor() {}; // destructor




}; // class B3CNeighbor



class Beam3contact 
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
  Beam3contact(const DRT::Discretization& pdiscret,
               const DRT::Discretization& cdiscret,
               const int& dofoffset,
               DRT::Element* element1,
               DRT::Element* element2,
               const Epetra_SerialDenseMatrix ele1pos,
               const Epetra_SerialDenseMatrix ele2pos);

  /*!
  \brief Copy Constructor

  Makes a deep copy of this contact element pair
  
  */
  Beam3contact(const Beam3contact& old);

  /*!
  \brief Destructor
  
  */
  virtual ~Beam3contact() {};
  
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
  \brief Closest Point Projection before Evaluation can start.

  */
  bool PreEvaluation(int beams_smoothing, std::map<int,LINALG::Matrix<3,1> >& currentpositions);


  /*!
  \brief Evaluate this contact element pair

  */
  bool Evaluate(LINALG::SparseMatrix& stiffmatrix,
                Epetra_Vector& fint, double& pp, bool ngf,
                std::map<std::pair<int,int>, RCP<Beam3contact> >& contactpairmap,
                int beams_smoothing);
  
  /*!
  \brief Reset Uzawa-based Lagrange multiplier
    
  */
  void Resetlmuzawa();
  

  /*!
    \brief Change the sign of the normal vector

    */
    void InvertNormal();



  /*!
    \brief Get local ID of the right boundary node

  */

  void GetBoundaryNode(int& nright, int nnode);



  /*!
    \brief Get length of an element

  */

  void GetEleLength(Epetra_SerialDenseMatrix& elepos, int& nright, double& length);



  /*!
    \brief Shift current normal vector to old normal vector

  */

  void ShiftNormal();


  /*!
    \Check if there is a difference between new and old gap definition e.g. if the beams centerlines have allready passed or not

  */

  bool GetNewGapStatus();



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
  
  //! current node tangents of both elements
  Epetra_SerialDenseMatrix ele1tangent_;
  Epetra_SerialDenseMatrix ele2tangent_;


  //! variable to check if old or modified gap function
  bool ngf_;

  //! variable to identify first call of a pair (for initializing)
  bool firstcall_;

  //! Uzawa-based Lagrange multiplier
  double lmuzawa_;
  
  //! gap function
 double oldgap_;

  //! gap function
  double gap_;
  
  //! flag indicating contact (active/inactive)
  bool contactflag_;
  
  //! flag indicating if elements are colinear or not
  bool elementscolinear_;

  //! coordinates of contact point on center lines of beams
  Epetra_SerialDenseVector x1_;
  Epetra_SerialDenseVector x2_;
  

  //! parameter values of contact point
  Epetra_SerialDenseVector xicontact_;


  //! neighbor elements of element 1

  RCP<B3CNeighbor> neighbor1_;

  //! neighbor elements of element 2

  RCP<B3CNeighbor> neighbor2_;

  //! normal vector of current time step
  Epetra_SerialDenseVector normal_;

  //! normal vector of last time step
  Epetra_SerialDenseVector normal_old_;











  //@}
 
  //! @name Private evaluation methods
  
  /*!
  \brief Find contact point via closest point projection
  
  */
  void ClosestPointProjection(int beams_smoothing);

  /*!
   \brief determine the interpolated nodal tangent vectors of both elements for smoothing

   */

  void CalculateNodalTangents(std::map<int,LINALG::Matrix<3,1> >& currentpositions);

  
  /*!
  \brief Utility method for CPP (evaluate nonlinear function f)
         
  */
  void EvaluateNewtonF(std::vector<double>& f,
                       const std::vector<double>& x1, const std::vector<double>& x2,
                       const std::vector<double>& dx1, const std::vector<double>& dx2,
                       const std::vector<double>& t1, const std::vector<double>& t2,
                       const double& norm, int beams_smoothing);
  
  /*!
  \brief Utility method for CPP (evaluate Jacobian of nonlinear function f)
  
  */
  void EvaluateNewtonGradF(LINALG::Matrix<2,2>& df,
                           const std::vector<double>& x1, const std::vector<double>& x2,
                           const std::vector<double>& dx1, const std::vector<double>& dx2,
                           const std::vector<double>& ddx1, const std::vector<double>& ddx2,
                           const std::vector<double>& t1, const std::vector<double>& t2,
                           const std::vector<double>& dt1, const std::vector<double>& dt2,
                           const double& norm, int beams_smoothing);
  
  /*!
  \brief Evaluate and assemble contact forces
  
  */
  void EvaluateFcContact(const double& pp, const double& gap,
                         const std::vector<double>& normal, Epetra_Vector& fint,
                         Epetra_SerialDenseVector funct1, Epetra_SerialDenseVector funct2, 
                         const int numnode1, const int numnode2);
  
  /*!
  \brief Evaluate and assemble contact stiffness
  
  */
  void EvaluateStiffcContact(const double& pp, const double& gap,
                             const std::vector<double>& normal, const double& norm,
                             LINALG::SparseMatrix& stiffmatrix,
                             const std::vector<double>& x1, const std::vector<double>& x2,
                             const std::vector<double>& dx1, const std::vector<double>& dx2,
                             const std::vector<double>& ddx1, const std::vector<double>& ddx2,
                             const Epetra_SerialDenseVector& funct1,
                             const Epetra_SerialDenseVector& funct2,  
                             const Epetra_SerialDenseMatrix& deriv1,
                             const Epetra_SerialDenseMatrix& deriv2,  
                             const Epetra_SerialDenseMatrix& secondderiv1,
                             const Epetra_SerialDenseMatrix& secondderiv2, 
                             const int numnode1, const int numnode2,
                             const std::vector<double>& XiContact,
                             int beams_smoothing);
  

   /*!
    \brief determine the neighbour elements of an element

    */
  void DetermineNeigbours(DRT::Element* element1,
                          DRT::Element* element2);



  /*!
  \brief Compute normal vector in contact point
          
  */
  void ComputeNormal(std::vector<double>& normal, double& gap, double& norm, 
                     const std::vector<double>& x1, const std::vector<double>& x2);
  


  /*!
  \brief Compute skalar product between old normal vector and new normal vector

  */
  double Computeskalar();

  /*!
  \brief Compute sign function of an arbitrary number

  */
  int Signum(double val);

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
                              std::vector<double>& dx1, std::vector<double>& dx2,
                              std::vector<double>& ddx1, std::vector<double>& ddx2,
                              const Epetra_SerialDenseVector& funct1,
                              const Epetra_SerialDenseVector& funct2,
                              const Epetra_SerialDenseMatrix& deriv1,
                              const Epetra_SerialDenseMatrix& deriv2,
                              const Epetra_SerialDenseMatrix& secondderiv1,
                              const Epetra_SerialDenseMatrix& secondderiv2,
                              const int& numnode1, const int& numnode2);
  

  /*!
  \brief Compute contact point tangents and their derivatives

  */
  void ComputeTangentsAndDerivs(std::vector<double>& t1, std::vector<double>& t2,
                                std::vector<double>& dt1, std::vector<double>& dt2,
                                const Epetra_SerialDenseVector& funct1,
                                const Epetra_SerialDenseVector& funct2,
                                const Epetra_SerialDenseMatrix& deriv1,
                                const Epetra_SerialDenseMatrix& deriv2,
                                const int& numnode1, const int& numnode2);



  /*!
  \brief Get derivatives of shape functions at the nodes

  */

  void GetNodalDerivatives(Epetra_SerialDenseMatrix& deriv1,
                           int node,
                           const int nnode,
                           double length,
                           const DRT::Element::DiscretizationType distype);


  /*!
  \brief Get shape functions and their derivatives at eta
          
  */
  void GetShapeFunctions(Epetra_SerialDenseVector& funct1,
                         Epetra_SerialDenseVector& funct2,
                         Epetra_SerialDenseMatrix& deriv1,
                         Epetra_SerialDenseMatrix& deriv2,
                         Epetra_SerialDenseMatrix& secondderiv1,
                         Epetra_SerialDenseMatrix& secondderiv2,
                         const std::vector<double>& eta);
  
  /*!
  \brief Compute linearizations of contact point
          
  */
  void ComputeLinXiAndLinEta(std::vector<double>& delta_xi, std::vector<double>& delta_eta,
                             const std::vector<double>& x1, const std::vector<double>& x2,
                             const std::vector<double>& dx1, const std::vector<double>& dx2,
                             const std::vector<double>& ddx1, const std::vector<double>& ddx2,
                             const Epetra_SerialDenseVector& funct1,
                             const Epetra_SerialDenseVector& funct2, 
                             const Epetra_SerialDenseMatrix& deriv1,
                             const Epetra_SerialDenseMatrix& deriv2, 
                             const std::vector<double>& normal, const double& norm,
                             const int numnode1, const int numnode2,
                             const std::vector<double>& XiContact,
                             int beams_smoothing);
  
  /*!
  \brief Compute linearization of gap
          
  */
  void ComputeLinGap(std::vector<double>& delta_gap,
                     std::vector<double>& delta_xi, std::vector<double>& delta_eta,
                     const std::vector<double>& x1, const std::vector<double>& x2,
                     const std::vector<double>& dx1, const std::vector<double>& dx2,
                     const Epetra_SerialDenseVector& funct1,
                     const Epetra_SerialDenseVector& funct2, 
                     const double& normdist,
                     const int& numnode1, const int& numnode2,
                     const std::vector<double>& normal, const double& norm,
                     const double& gap,Epetra_SerialDenseMatrix& delta_x1_minus_x2,
                     int beams_smoothing);
  
  /*!
  \brief Compute linearization of normal
          
  */
  void ComputeLinNormal(Epetra_SerialDenseMatrix& delta_n, 
                        const std::vector<double>& x1, const std::vector<double>& x2,
                        const double& norm,
                        const int& numnode1, const int& numnode2,
                        const Epetra_SerialDenseMatrix& delta_x1_minus_x2,
                        const std::vector<double>& normal, 
                        const std::vector<double>& XiContact,
                        int beams_smoothing);
  
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
   
  //! @name Finite difference checks
  
  /*!
  \brief FD check of closest point projection (delta_xi, delta_eta)
          
  */                
  void FDCheckCPP(const int& numnode1, const int& numnode2, 
                  const std::vector<double>& delta_xi,
                  const std::vector<double>& delta_eta,
                  const std::vector<double>& XiContact,
                  int beams_smoothing);
  



  /*!
  \brief FD check of local Newton Iteration for CCP

  */
  void FDCheckNewtonCPP(const int& numnode1, const int& numnode2,
         const std::vector<double>& eta,int beams_smoothing);



  /*!
  \brief FD check of of gap function (delta_gap)
          
  */
  void FDCheckLinGap(const int& numnode1, const int& numnode2, 
                     const std::vector<double>& normal,
                     const double& norm,
                     const std::vector<double>& delta_gap,
                     const double& gap,
                     int beams_smoothing);
  
  /*!
  \brief FD check of contact normal vector (delta_n)
          
  */
  void FDCheckLinNormal(const int& numnode1, const int& numnode2, 
                        const Epetra_SerialDenseMatrix& delta_n,
                        const std::vector<double>& normal,
                        int beams_smoothing);
  
  /*!
  \brief FD check of contact stiffness matrices (stiffc1, stiffc2)
          
  */
  void FDCheckStiffc(const int& numnode1, const int& numnode2, 
                     const Epetra_SerialDenseMatrix& stiffc1,
                     const Epetra_SerialDenseMatrix& stiffc2, 
                     const double& pp, const std::vector<double>& normal,
                     const double& gap,
                     const Epetra_SerialDenseVector& funct1,
                     const Epetra_SerialDenseVector& funct2,
                     int beams_smoothing);

  //@}

}; // class Beam3contact
} // namespace CONTACT

#endif  // #ifndef BEAM3C_H
