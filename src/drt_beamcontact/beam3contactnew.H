/*----------------------------------------------------------------------*/
/*!
\file beam3contactnew.H

\brief contact element for contact between two 3D beam elements

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/

#ifndef BEAM3CNEW_H
#define BEAM3CNEW_H

#include "../drt_beaminteraction/beam3contact_utils.H"
#include "beam3contactinterface.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "Sacado.hpp"

#include "../drt_beaminteraction/beam3contact_defines.H"
#include "../drt_beaminteraction/beam3contact_tangentsmoothing.H"
#include "../headers/FAD_utils.H"

#include "../drt_lib/standardtypes_cpp.H"

namespace CONTACT
{

/*!
 \brief contact element for contact between two 3D beam elements

 Refer also to the Semesterarbeit of Matthias Mayr, 2010

 \author meier (meier@lnm.mw.tum.de)
 */

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
                    const std::map<int,int>& dofoffsetmap,
                    DRT::Element* element1,
                    DRT::Element* element2,
                    Teuchos::ParameterList& beamcontactparams);

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
  inline const std::map<int,int>& DofOffset() const { return dofoffsetmap_;};

  /*!
  \brief Get first element
  */
  inline const DRT::Element* Element1() { return element1_;};

  /*!
  \brief Get first element
  */
  inline const DRT::Element* Element2() { return element2_;};

  /*!
  \brief Get number of contact points on this element pair
  */
  int GetNumCps(){return 1;};
  int GetNumGps(){return 0;};
  int GetNumEps(){return 0;};

  /*!
  \brief Get vector of type declarations (0=closest point contact, 1=gauss point contact, 2= end point contact) of all contact pairs
  */
  std::vector<int> GetContactType()
  {

    std::vector<int > types(1, 0);

    return types;
  };

  /*!
  \brief Get gap of this contact pair
  */
  std::vector<double> GetGap()
  {
    std::vector<double> gaps(1,FADUTILS::CastToDouble(gap_));

    return gaps;
  };

  /*!
  \brief Get contact force of this contact pair
  */
  std::vector<double> GetContactForce()
  {
    std::vector<double> forces(1,FADUTILS::CastToDouble(fp_));

    return forces;
  };

  /*!
  \brief Get contact angle of this contact pair
  */
  std::vector<double> GetContactAngle()
  {
    double angle = 0.0;
    double cosangle = FADUTILS::CastToDouble(tangentproduct_);
    if(cosangle<1.0)
      angle=acos(cosangle); //returns an angle \in [0;pi/2] since scalarproduct \in [0;1.0]
    else
      angle=0;
    std::vector<double> angles(1,angle);

    return angles;
  };

  /*!
  \brief Get closest point of this contact pair
  */
  std::vector<std::pair<double,double> > GetClosestPoint()
  {
    std::pair<double,double> cp(0.0,0.0);
    cp.first=FADUTILS::CastToDouble(xi1_);
    cp.second=FADUTILS::CastToDouble(xi2_);
    std::vector<std::pair<double,double> > cps(1,cp);

    return cps;
  };

  /*!
  \brief Return number of individual contact segments on element pair
  */
  std::pair<int,int> GetNumSegments()
  {
    return std::make_pair(1,1);
  };

  /*!
  \brief Return ids of active segments
  */
  std::vector<std::pair<int,int> > GetSegmentIds()
  {
    std::vector<std::pair<int,int> > ids(1,std::make_pair(0,0));

    return ids;
  };

  /*!
  \brief Get flag indicating whether contact is active (true) or inactive (false)
  */
  bool GetContactFlag() { return contactflag_; };

  /*!
  \brief Get coordinates of contact point of element1
  */
  const std::vector<LINALG::Matrix<3,1> > GetX1()
  {
    std::vector<LINALG::Matrix<3,1> > r1(1,LINALG::Matrix<3,1>(true));

    for(int j=0;j<3;j++)
      r1[0](j)=FADUTILS::CastToDouble(r1_(j));

    return r1;
  };

  /*!
  \brief Get coordinates of contact point of element2
  */
  const std::vector<LINALG::Matrix<3,1> > GetX2()
  {
    std::vector<LINALG::Matrix<3,1> > r2(1,LINALG::Matrix<3,1>(true));

    for(int j=0;j<3;j++)
      r2[0](j)=FADUTILS::CastToDouble(r2_(j));

    return r2;
  };

  /*!
  \brief Get normal vector
  */
  const Epetra_SerialDenseVector GetNormal()
  {
    Epetra_SerialDenseVector normal(3);

    if(GetNewGapStatus()==true)
    {
      for(int i=0;i<3;i++)
        normal(i)=-FADUTILS::CastToDouble(normal_(i));
    }
    else
    {
      for(int i=0;i<3;i++)
        normal(i)=FADUTILS::CastToDouble(normal_(i));
    }

    return normal;
  }

  /*!
  \brief Get flag indicating whether the nodal values of one element had been shifted due to r1=r2
  */
  bool GetShiftStatus() {return shiftnodalvalues_;};

  /*!
    \Check, if there is a difference between the result of the new and old gap definition, i.e. if the beams centerlines have already crossed or not.
  */
  bool GetNewGapStatus();
  //@}

  /*!
    \Get energy of penalty contact.
  */
  double GetEnergy()
  {
    TYPE energy=0.5*pp_*gap_*gap_;
    return FADUTILS::CastToDouble(energy);
  };

  /*!
    Not needed for beam3contactnew variables
  */
  double GetUnscaledPerpEnergy()
  {
    return 0.0;
  };

  /*!
    Not needed for beam3contactnew variables
  */
  double GetUnscaledParallelEnergy()
  {
    return 0.0;
  };

  /*!
  \brief Get normal vector of last time step
  */
  LINALG::TMatrix<TYPE, 3, 1>* GetNormalOld()
  {
    if(firsttimestep_)
    {
      dserror("Vector normal_old_requested but not available in the first time step the pair has been found: Choose larger search radius!!!");
    }

    //Only deliver normal_old_ when a valid closest point projection for last time step existed!
    if (FADUTILS::Norm(xi1_old_) < 1.0+XIETATOL and FADUTILS::Norm(xi2_old_) < 1.0+XIETATOL and !oldcppunconverged_)
    {
      return &normal_old_;
    }
    else
      return NULL;
  };

  /*!
  \brief Check, if it is the first time step the element is in contact
  */
  bool FirstTimeStep()
  {
    return firsttimestep_;
  };
  //@}

  /** \brief print this beam contact element pair to screen
   *
   *  \author grill
   *  \date 05/16 */
  void Print() const {};

  //! @name Public evaluation methods
  /*!
  \brief Evaluate this contact element pair
  */
  bool Evaluate(LINALG::SparseMatrix& stiffmatrix,
                Epetra_Vector& fint,
                const double& pp,
                std::map<std::pair<int,int>, Teuchos::RCP<Beam3contactinterface > >& contactpairmap,
                Teuchos::ParameterList& timeintparams,
                bool fdcheck=false);

  /*!
  \brief Change the sign of the normal vector: This has to be done at the end of a time step when the remainig penetration is
  larger that the sum of the beam radii (R1+R2). Otherwise, the beams could cross in the next time step when the new gap function
  definition (ngf_=true) for slender beams is applied!
  */
  void InvertNormal();

  /*!
    \brief Update of class variables at the end of a time step
  */
  void UpdateClassVariablesStep();

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
  const std::map<int,int>& dofoffsetmap_;

  //! first element of contact pair
  DRT::Element* element1_;

  //! second element of contact pair
  DRT::Element* element2_;

  //! beam contact parameter list
  Teuchos::ParameterList& bcparams_;

  //! current node coordinates of the two elements
  LINALG::TMatrix<TYPE,3*numnodes*numnodalvalues,1> ele1pos_;
  LINALG::TMatrix<TYPE,3*numnodes*numnodalvalues,1> ele2pos_;

  LINALG::TMatrix<double,3*numnodes*numnodalvalues,1> ele1pos_old_;
  LINALG::TMatrix<double,3*numnodes*numnodalvalues,1> ele2pos_old_;

  LINALG::TMatrix<double,3*numnodes*numnodalvalues,1> ele1pos_lastiter_;
  LINALG::TMatrix<double,3*numnodes*numnodalvalues,1> ele2pos_lastiter_;

  //! variable to check which smoothing type should be applied
  int smoothing_;

  //! sgn(normal*normal_old)
  double sgn_;

  //! variable to identify first call of a pair within considered time step (for initializing)
  bool firstcallofstep_;

  //! variable to identify the first time step in which the pair was found
  bool firsttimestep_;

  //! gap function according to original (ngf_==false) or modified (ngf_==true) definition
  TYPE gap_;

  //! gap function according to original definition
  TYPE gap_original_;

  //! flag indicating contact (active/inactive)
  bool contactflag_;

  //! flag indicating (active/inactive) damping forces
  bool dampingcontactflag_;

  //! flag indicating if contact was active in the last time step
  bool oldcontactflag_;

  //! flag indicating whether an element has at least once been active during the current time step
  bool currentlyincontact_;

  //! flag indicating if elements are collinear or not
  bool elementscolinear_;

  //! flag indicating if elements share the same contact point, i.e. r1_=r2_ --> evaluation not possible
  bool elementscrossing_;

  //! flag indicating if the element nodal positions have been shifted in order to avoid r1_=r2_
  bool shiftnodalvalues_;

  //! coordinates of contact point on center lines of beams
  LINALG::TMatrix<TYPE, 3, 1> r1_;
  LINALG::TMatrix<TYPE, 3, 1> r2_;
  LINALG::TMatrix<TYPE, 3, 1> r1_old_;
  LINALG::TMatrix<TYPE, 3, 1> r2_old_;

  //! derivative of contact point coordinates
  LINALG::TMatrix<TYPE, 3, 1> r1_xi_;
  LINALG::TMatrix<TYPE, 3, 1> r2_xi_;
  LINALG::TMatrix<TYPE, 3, 1> r1_xi_old_;
  LINALG::TMatrix<TYPE, 3, 1> r2_xi_old_;

  //! parameter values of contact point
  TYPE xi1_;
  TYPE xi2_;

  //! parameter values of contact point of last time step
  double xi1_old_;
  double xi2_old_;

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

  //! current penalty parameter
  double pp_;

  //! scalar value of current contact force, standard case: fp_ = -pp_ * gap_
  TYPE fp_;

  //! derivative of scalar value of current contact force, standard case: fp_ = -pp_
  TYPE dfp_;

  //! scalar value of current damping force and its derivative, standard case: fdamp_ = -d_ * fd_ e.g. = -d_ * gap_t =
  TYPE fd_;
  TYPE dfd_;
  TYPE d_;
  TYPE dd_;

  //!current Newton iteration
  int iter_;

  //!current time step
  int numstep_;

  //!time step size
  double dt_;

  //!flag indicating when one beam slides over the end point of the second beam
  bool beamendcontactopened_;

  //!flag indicating that a pair is sorted out due to almost parallelism
  bool beamsalmostparallel_;

  //!flag indicating that closest point projection did not converge
  bool cppunconverged_;

  //!flag indicating that closest point projection of last time step did not converge
  bool oldcppunconverged_;

  //! initial element lengths
  double ele1length_;
  double ele2length_;

  //! Has the pair asked for the vector normal_old_ of an neighbor element during the current iteration step?
  bool neighbornormalrequired_;

  //! Scalar product of normalized tangents at contact point (needed in order to exclude almost parallel elements)
  TYPE tangentproduct_;

  //! Cross-section radius of beam 1
  double radius1_;

  //! Cross-section radius of beam 2
  double radius2_;

  //! Additive extrusion value of bounding box applied in contact search
  double searchboxinc_;

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
  \brief Evaluate algorithmic forces that make convergence more robust / faster
  */
  void EvaluateAlgorithmicForce(const double& pp,
                                Epetra_Vector* fint,
                                const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                                const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2,
                                LINALG::TMatrix<TYPE, 3*numnodes*numnodalvalues, 1>* fc1_FAD = NULL,
                                LINALG::TMatrix<TYPE, 3*numnodes*numnodalvalues, 1>* fc2_FAD = NULL);

  /*!
  \brief Evaluate algorithmic stiffnesses that make convergence more robust / faster
  */
  void EvaluateAlgorithmicStiff(const double& pp,
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
  \brief Compute linearization of time derivative g_t of gap g
  */
  void ComputeLinGapt(LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_gap_t,
                      const LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>&  delta_xi,
                      const LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_eta,
                      const LINALG::TMatrix<TYPE, 3, 2*3*numnodes*numnodalvalues>& delta_n,
                      const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                      const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                      const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                      const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2,
                      const LINALG::TMatrix<TYPE, 3, 1>& r1_old,
                      const LINALG::TMatrix<TYPE, 3, 1>& r2_old,
                      const LINALG::TMatrix<TYPE, 3, 1>& r1_xi_old,
                      const LINALG::TMatrix<TYPE, 3, 1>& r2_xi_old);

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
  \brief Calculate scalar contact force and linearization according to chosen penalty law
  */
  void CalcPenaltyLaw();

  /*!
  \brief Calculate scalar damping contribution to contact force and linearization
  */
  void CalcDampingLaw();

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
  \brief Compute coordinates of contact points of last time step from the discretization
  */
  void ComputeOldCoordsAndDerivs(LINALG::TMatrix<TYPE,3,1>& r1_old,
                              LINALG::TMatrix<TYPE,3,1>& r2_old,
                              LINALG::TMatrix<TYPE,3,1>& r1_xi_old,
                              LINALG::TMatrix<TYPE,3,1>& r2_xi_old,
                              const LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N1,
                              const LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N2,
                              const LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N1_xi,
                              const LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N2_xi);

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
                     TYPE& norm_delta_r,
                     std::map<std::pair<int,int>, Teuchos::RCP<Beam3contactinterface > >& contactpairmap);

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
      \brief These method shifts the nodal positions applied within the beam contact framework by a small pre-defined
      amount in order to enable contact evaluation in the case of two identical contact points, i.e r1=r2
    */
    void ShiftNodalPositions();

    /*!
      \brief Get the vector normal_old_ from the neighbor element
    */
    void GetNeighborNormalOld(std::map<std::pair<int,int>, Teuchos::RCP<Beam3contactinterface > >& contactpairmap);

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

    /*!
      \brief Check, if the current contact point is at one of the two ends of a physical beam
    */
    void CheckBoundaryContact();

    /*!
      \brief Get jacobi factor of beam element
    */
    double GetJacobi(DRT::Element* element1);


    /*!
      \brief Update of class variables at the end of a Newton step
    */
    void UpdateClassVariablesIter();

    /*!
      \brief Set class variables at the beginning of a Newton step
    */
    void SetClassVariables(const double& pp, Teuchos::ParameterList timeintparams);

  //@}

}; // class Beam3contactnew
} // namespace CONTACT

#endif  // #ifndef BEAM3CNEW_H
