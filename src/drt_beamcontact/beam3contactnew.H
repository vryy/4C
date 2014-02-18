/*!----------------------------------------------------------------------
\file beam3contactnew.H

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>
*----------------------------------------------------------------------*/
#ifndef BEAM3CNEW_H
#define BEAM3CNEW_H

#include "beam3contact_utils.H"
#include "beam3contact_tangentsmoothing.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_sparsematrix.H"
#include "beam3contact_defines.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "Sacado.hpp"

#ifdef NEWBEAMCONTACT

namespace CONTACT
{

/*!
 \brief contact element for contact between two 3D beam elements

 Refer also to the Semesterarbeit of Matthias Mayr, 2010

 \author meier (meier@lnm.mw.tum.de)
 */

class Beam3contactinterface
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
  Beam3contactinterface() {}

  /*!
  \brief Destructor
  */
  virtual ~Beam3contactinterface() {}

  //@}

  //! @name Access methods
  /*!
  \brief Get problem discretization
  */
  virtual const DRT::Discretization& ProblemDiscret() const =0;

  /*!
  \brief Get beam contact discretization
  */
  virtual const DRT::Discretization& ContactDiscret() const =0;

  /*!
  \brief Get offset of dofs between cdiscret and pdiscret
  */
  virtual const int& DofOffset() const =0;

  /*!
  \brief Get first element
  */
  virtual const DRT::Element* Element1()=0;
  //inline const DRT::Element* Element1() { return element1_;};

  /*!
  \brief Get first element
  */
  virtual const DRT::Element* Element2()=0;

  /*!
  \brief Get gap of this contact pair
  */
  virtual const double GetGap()=0;

  /*!
  \brief Get current augmented lagrange multiplier 'lmuzawa_'
  */
  virtual const double Getlmuzawa()=0;

  /*!
  \brief Get flag ndicating whether contact is active (true) or inactive (false)
  */
  virtual bool GetContactFlag()=0;

  /*!
  \brief Get coordinates of contact point of element1 and element2
  */
  virtual const Epetra_SerialDenseVector GetX1()=0;

  virtual const Epetra_SerialDenseVector GetX2()=0;

  /*!
    \Check, if there is a difference between the result of the new and old gap definition, i.e. if the beams centerlines have already crossed or not.
  */
  virtual bool GetNewGapStatus()=0;

  /*!
  \brief Get flag indicating whether the nodal values of one element had been shifted due to r1=r2
  */
  virtual bool GetShiftStatus()=0;
  //@}


  //! @name Public evaluation methods
  /*!
  \brief Evaluate this contact element pair
  */
  virtual bool Evaluate(LINALG::SparseMatrix& stiffmatrix,
                        Epetra_Vector& fint,
                        const double& pp)=0;

  //! return appropriate internal implementation class (acts as a simple factory)
  Teuchos::RCP<Beam3contactinterface> Impl(const int numnodes,
                                     const int numnodalvalues,
                                     const DRT::Discretization& pdiscret,
                                     const DRT::Discretization& cdiscret,
                                     const int& dofoffset,
                                     DRT::Element* element1,
                                     DRT::Element* element2,
                                     Teuchos::ParameterList beamcontactparams);

  /*!
  \brief Reset Uzawa-based Lagrange multiplier
  */
  virtual void Resetlmuzawa()=0;


  /*!
  \brief Change the sign of the normal vector: This has to be done at the end of a time step when the remainig penetration is
  larger that the sum of the beam radii (R1+R2). Otherwise, the beams could cross in the next time step when the new gap function
  definition (ngf_=true) for slender beams is applied!
  */
  virtual void InvertNormal()=0;

  /*!
    \brief Shift current normal vector to old normal vector at the end of a time step. This is necessary when the new gap function
    definition (ngf_=true) for slender beams is applied!
  */
  virtual void ShiftNormal()=0;

  /*!
    \brief Update Uzawa-based Lagrange multiplier
    */
  virtual void Updatelmuzawa(const double& currentpp)=0;

  /*
  \brief Update nodal coordinates of both elements at the beginning of a new time step!
  */
  virtual void UpdateElePos(Epetra_SerialDenseMatrix& newele1pos,
                            Epetra_SerialDenseMatrix& newele2pos)=0;

  /*
  \brief Update interpolated nodal tangents for tangent smoothing
  */
  virtual void UpdateEleSmoothTangents(std::map<int,LINALG::Matrix<3,1> >& currentpositions)=0;

}; // class Beam3contactinterface




template<const int numnodes , const int numnodalvalues>
class Beam3contactnew: public Beam3contactinterface
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
  Beam3contactnew(  const DRT::Discretization& pdiscret,
                    const DRT::Discretization& cdiscret,
                    const int& dofoffset,
                    DRT::Element* element1,
                    DRT::Element* element2,
                    Teuchos::ParameterList beamcontactparams);

  /*!
  \brief Copy Constructor
  Makes a deep copy of this contact element pair
  */
  Beam3contactnew(const Beam3contactnew& old);

  /*!
  \brief Destructor
  */
  virtual ~Beam3contactnew() {};
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
  const double GetGap() { return BEAMCONTACT::CastToDouble(gap_); };

  /*!
  \brief Get current augmented lagrange multiplier 'lmuzawa_'
  */
  const double Getlmuzawa() { return BEAMCONTACT::CastToDouble(lmuzawa_); };

  /*!
  \brief Get flag indicating whether contact is active (true) or inactive (false)
  */
  bool GetContactFlag() { return contactflag_; };

  /*!
  \brief Get coordinates of contact point of element1 and element2
  */
  const Epetra_SerialDenseVector GetX1()
  {
    Epetra_SerialDenseVector r1;
    r1.Resize(3);
    for (int i=0;i<3;i++)
      r1(i)=BEAMCONTACT::CastToDouble(r1_(i));

    return r1;
  };

  const Epetra_SerialDenseVector GetX2()
  {
    Epetra_SerialDenseVector r2;
    r2.Resize(3);
    for (int i=0;i<3;i++)
      r2(i)=BEAMCONTACT::CastToDouble(r2_(i));

    return r2;
  };

  /*!
  \brief Get flag indicating whether the nodal values of one element had been shifted due to r1=r2
  */
  bool GetShiftStatus() {return shiftnodalvalues_;};

  /*!
    \Check, if there is a difference between the result of the new and old gap definition, i.e. if the beams centerlines have already crossed or not.
  */
  bool GetNewGapStatus();
  //@}


  //! @name Public evaluation methods
  /*!
  \brief Evaluate this contact element pair
  */
  bool Evaluate(LINALG::SparseMatrix& stiffmatrix,
                Epetra_Vector& fint,
                const double& pp);

  /*!
  \brief Reset Uzawa-based Lagrange multiplier
  */
  void Resetlmuzawa();

  /*!
  \brief Change the sign of the normal vector: This has to be done at the end of a time step when the remainig penetration is
  larger that the sum of the beam radii (R1+R2). Otherwise, the beams could cross in the next time step when the new gap function
  definition (ngf_=true) for slender beams is applied!
  */
  void InvertNormal();

  /*!
    \brief Shift current normal vector to old normal vector at the end of a time step. This is necessary when the new gap function
    definition (ngf_=true) for slender beams is applied!
  */
  void ShiftNormal();

  /*!
    \brief Update Uzawa-based Lagrange multiplier
    */
  void Updatelmuzawa(const double& currentpp);

  /*
  \brief Update nodal coordinates of both elements at the beginning of a new time step!
  */
  void UpdateElePos(Epetra_SerialDenseMatrix& newele1pos,
                    Epetra_SerialDenseMatrix& newele2pos);

  /*
  \brief Update interpolated nodal tangents for tangent smoothing
  */
  void UpdateEleSmoothTangents(std::map<int,LINALG::Matrix<3,1> >& currentpositions);


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

  //const int numnode1_ = NumNodes1(element1_);

  //! current node coordinates of the two elements
  LINALG::TMatrix<TYPE,3*numnodes*numnodalvalues,1> ele1pos_;
  LINALG::TMatrix<TYPE,3*numnodes*numnodalvalues,1> ele2pos_;

  //! variable to check if old or modified gap function
  bool ngf_;

  //! variable to check which smoothing type should be applied
  int smoothing_;

  //! sgn(normal*normal_old)
  double sgn_;

  //! variable to identify first call of a pair (for initializing)
  bool firstcall_;

  //! Uzawa-based Lagrange multiplier
  TYPE lmuzawa_;

  //! gap function according to original (ngf_==false) or modified (ngf_==true) definition
  TYPE gap_;

  //! gap function according to original definition
  TYPE gap_original_;

  //! flag indicating contact (active/inactive)
  bool contactflag_;

  //! flag indicating if elements are collinear or not
  bool elementscolinear_;

  //! flag indicating if elements share the same contact point, i.e. r1_=r2_ --> evaluation not possible
  bool elementscrossing_;

  //! flag indicating if the element nodal positions have been shifted in order to avoid r1_=r2_
  bool shiftnodalvalues_;

  //! coordinates of contact point on center lines of beams
  LINALG::TMatrix<TYPE, 3, 1> r1_;
  LINALG::TMatrix<TYPE, 3, 1> r2_;

  //! parameter values of contact point
  TYPE xi1_;
  TYPE xi2_;

  //! normal vector of current time step
  LINALG::TMatrix<TYPE, 3, 1>  normal_;

  //! normal vector of last time step
  LINALG::TMatrix<TYPE, 3, 1> normal_old_;

  //! neighbor elements of element 1
  Teuchos::RCP<B3CNeighbor> neighbors1_;

  //! neighbor elements of element 2
  Teuchos::RCP<B3CNeighbor> neighbors2_;

  //! averaged nodal tangents, necessary for smoothed tangent fields of C^0 Reissner beams
  LINALG::Matrix<3*numnodes,1> nodaltangentssmooth1_;
  LINALG::Matrix<3*numnodes,1> nodaltangentssmooth2_;

  //@}

  //! @name Private evaluation methods

  /*!
  \brief Evaluate contact forces
  */
  void EvaluateFcContact( const double& pp,
                          Epetra_Vector* fint,
                          const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                          const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2,
                          LINALG::TMatrix<TYPE, 3*numnodes*numnodalvalues, 1>* fc1_FAD = NULL,
                          LINALG::TMatrix<TYPE, 3*numnodes*numnodalvalues, 1>* fc2_FAD = NULL);

  /*!
  \brief Evaluate contact stiffness
  */
  void EvaluateStiffcContact(const double& pp,
                             const TYPE& norm_delta_r,
                             const LINALG::TMatrix<TYPE, 3, 1>& delta_r,
                             LINALG::SparseMatrix& stiffmatrix,
                             const LINALG::TMatrix<TYPE, 3, 1>& r1,
                             const LINALG::TMatrix<TYPE, 3, 1>& r2,
                             const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                             const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                             const LINALG::TMatrix<TYPE, 3, 1>& r1_xixi,
                             const LINALG::TMatrix<TYPE, 3, 1>& r2_xixi,
                             const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                             const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2,
                             const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1_xi,
                             const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xi,
                             const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1_xixi,
                             const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xixi);

  /*!
  \brief Compute linearizations of contact point parameter coordinates xi and eta
  */
  void ComputeLinXiAndLinEta( LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>&  delta_xi,
                              LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_eta,
                              const LINALG::TMatrix<TYPE, 3, 1>& delta_r,
                              const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                              const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                              const LINALG::TMatrix<TYPE, 3, 1>& r1_xixi,
                              const LINALG::TMatrix<TYPE, 3, 1>& r2_xixi,
                              const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                              const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2,
                              const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1_xi,
                              const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xi);

  /*!
  \brief Compute linearization of gap g
  */
  void ComputeLinGap(LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_gap,
                     const LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>&  delta_xi,
                     const LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_eta,
                     const LINALG::TMatrix<TYPE, 3, 1>& delta_r,
                     const TYPE& norm_delta_r,
                     const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                     const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                     const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                     const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2);

  /*!
  \brief Compute linearization of normal vector n
  */
  void ComputeLinNormal(LINALG::TMatrix<TYPE, 3, 2*3*numnodes*numnodalvalues>& delta_normal,
                        const LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>&  delta_xi,
                        const LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_eta,
                        const TYPE& norm_delta_r,
                        const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                        const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                        const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                        const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2);

  /*!
  \brief Find contact point via closest point projection
  */
  void ClosestPointProjection();

  /*!
  \brief Calculate shape function values for given parameter values
  */
  void GetShapeFunctions( LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                          LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2,
                          LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1_xi,
                          LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xi,
                          LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1_xixi,
                          LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xixi,
                          const TYPE& eta1,
                          const TYPE& eta2);

  /*!
  \brief Assemble the shape functions into corresponding matrices
  */
  void AssembleShapefunctions(const LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues>& N_i,
                              const LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues>& N_i_xi,
                              const LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues>& N_i_xixi,
                              LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N,
                              LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N_xi,
                              LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N_xixi);

  /*!
  \brief Compute coordinates and their derivatives from the discretization
  */
  void ComputeCoordsAndDerivs(LINALG::TMatrix<TYPE,3,1>& r1,
                              LINALG::TMatrix<TYPE,3,1>& r2,
                              LINALG::TMatrix<TYPE,3,1>& r1_xi,
                              LINALG::TMatrix<TYPE,3,1>& r2_xi,
                              LINALG::TMatrix<TYPE,3,1>& r1_xixi,
                              LINALG::TMatrix<TYPE,3,1>& r2_xixi,
                              const LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N1,
                              const LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N2,
                              const LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N1_xi,
                              const LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N2_xi,
                              const LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N1_xixi,
                              const LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N2_xixi);

  /*!
  \brief Utility method for CPP (evaluate nonlinear function f)
  */
  void EvaluateOrthogonalityCondition(LINALG::TMatrix<TYPE,2,1>& f,
                                      const LINALG::TMatrix<TYPE,3,1>& delta_r,
                                      const double norm_delta_r,
                                      const LINALG::TMatrix<TYPE,3,1>& r1_xi,
                                      const LINALG::TMatrix<TYPE,3,1>& r2_xi,
                                      const LINALG::TMatrix<TYPE,3,1>& t1,
                                      const LINALG::TMatrix<TYPE,3,1>& t2);

  /*!
  \brief Utility method for CPP (evaluate Jacobian of nonlinear function f)
  */
  void EvaluateLinOrthogonalityCondition( LINALG::TMatrix<TYPE,2,2>& df,
                                          LINALG::TMatrix<TYPE,2,2>& dfinv,
                                          const LINALG::TMatrix<TYPE,3,1>& delta_r,
                                          const double norm_delta_r,
                                          const LINALG::TMatrix<TYPE,3,1>& r1_xi,
                                          const LINALG::TMatrix<TYPE,3,1>& r2_xi,
                                          const LINALG::TMatrix<TYPE,3,1>& r1_xixi,
                                          const LINALG::TMatrix<TYPE,3,1>& r2_xixi,
                                          const LINALG::TMatrix<TYPE,3,1>& t1,
                                          const LINALG::TMatrix<TYPE,3,1>& t2,
                                          const LINALG::TMatrix<TYPE,3,1>& t1_xi,
                                          const LINALG::TMatrix<TYPE,3,1>& t2_xi);

  /*!
  \brief Compute normal vector and gap function at contact point
  */
  void ComputeNormal(LINALG::TMatrix<TYPE,3,1>& delta_r,
                     TYPE& norm_delta_r);

  /*!
  \brief Check, if we have contact or not (gap < 0 ???)
  */
  void CheckContactStatus(const double& pp);

  /*!
    \brief Get global dofs of a node

    Internally this method first extracts the dofs of the given node
    in the beam contact discretization (which has its own dofs) and
    then transfers these dofs to their actual GIDs in the underlying
    problem discretization by applying the pre-computed dofoffset_.
    */
    std::vector<int> GetGlobalDofs(const DRT::Node* node);

    /*!
      \brief These method shifts the nodal positions applied within the beam contact framework py a small pre-defined
      amount in order to enable contact evaluation in the case of two identical contact points, i.e r1=r2
    */
    void ShiftNodalPositions();

    void FADCheckLinXiAndLinEta(const LINALG::TMatrix<TYPE, 3, 1>& delta_r,
                                const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                                const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                                const LINALG::TMatrix<TYPE, 3, 1>& r1_xixi,
                                const LINALG::TMatrix<TYPE, 3, 1>& r2_xixi,
                                const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                                const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2,
                                const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1_xi,
                                const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xi);

    void FADCheckLinOrthogonalityCondition(const LINALG::TMatrix<TYPE, 3, 1>& delta_r,
                                           const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                                           const LINALG::TMatrix<TYPE, 3, 1>& r2_xi);

  //@}

}; // class Beam3contactnew
} // namespace CONTACT

#endif //#ifdef NEWBEAMCONTACT

#endif  // #ifndef BEAM3CNEW_H
