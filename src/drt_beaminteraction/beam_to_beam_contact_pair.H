/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_to_beam_contact_pair.H

\brief One beam contact pair (two beam elements) consisting of several contact segments

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef BEAMTOBEAMCONTACTPAIR_H
#define BEAMTOBEAMCONTACTPAIR_H

#include "beam_contact_pair.H"

// Todo get rid of header inclusions
#include "beam3contact_utils.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "Sacado.hpp"
#include "../headers/FAD_utils.H"

#include "../drt_inpar/inpar_beamcontact.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_beam3/beam3_base.H"

#include "../drt_lib/standardtypes_cpp.H"

#include <Epetra_SerialDenseVector.h>

#include "beam3contact_defines.H"
#include "beam3contact_tangentsmoothing.H"

//forward declaration ...
namespace LINALG {
  class SerialDenseVector;
  class SerialDenseMatrix;
}
namespace BEAMINTERACTION
{

template <unsigned int numnodes, unsigned int numnodalvalues>
class BeamToBeamContactVariables;

/*!
 \brief class for contact between two 3D beam elements
 */
template<unsigned int numnodes, unsigned int numnodalvalues>
class BeamToBeamContactPair: public BeamContactPair
{
public:

  //! @name Friends
  //no friend classes defined
  //@}

  //! @name Constructors and destructors and related methods
  /*!
  \brief Standard Constructor
  */
  BeamToBeamContactPair();

  /*!
  \brief Destructor
  */
  virtual ~BeamToBeamContactPair() {};

  //! Setup
  virtual void Setup();
  //@}

  //! @name Derived methods from base class
  /*!
  \brief things that need to be done in a separate loop before the actual evaluation loop
         over all contact pairs
  */
  void PreEvaluate() override;

  /*!
  \brief Evaluate this contact element pair, return value indicates whether pair is active,
         i.e. non-zero values for force and stiffmat are returned
  */
  virtual bool Evaluate(
      LINALG::SerialDenseVector* forcevec1,
      LINALG::SerialDenseVector* forcevec2,
      LINALG::SerialDenseMatrix* stiffmat11,
      LINALG::SerialDenseMatrix* stiffmat12,
      LINALG::SerialDenseMatrix* stiffmat21,
      LINALG::SerialDenseMatrix* stiffmat22);

  /*
  \brief Update state of translational nodal DoFs (absolute positions and tangents) of both elements
  */
  virtual void ResetState(
      const std::vector<double>& centerline_dofvec_ele1,
      const std::vector<double>& centerline_dofvec_ele2);

  /** \brief print information about this beam contact element pair to screen
   *
   *  \author grill
   *  \date 05/16 */
  virtual void Print(std::ostream& out) const;


  /** \brief print this beam contact element pair to screen
   *
   *  \author grill
   *  \date 12/16 */
  virtual void PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const;
  //@}

  //! @name Access methods

  /*!
  \brief Get flag indicating whether contact is active (true) or inactive (false)
  */
  virtual inline bool GetContactFlag() const
  {
    // The element pair is assumed to be active when we have at least one active contact point
    return (cpvariables_.size()+gpvariables_.size()+epvariables_.size());
  }

  /*!
  \brief Get number of active contact point pairs on this element pair
  */
  virtual unsigned int GetNumAllActiveContactPointPairs() const
  {
    return (unsigned int) (cpvariables_.size() + gpvariables_.size() + epvariables_.size());
  }

  /*!
  \brief Get coordinates of all active contact points on element1 and element2
  */
  virtual void GetAllActiveContactPointCoordsElement1(
      std::vector<LINALG::TMatrix<double,3,1> >& coords ) const;

  virtual void GetAllActiveContactPointCoordsElement2(
      std::vector<LINALG::TMatrix<double,3,1> >& coords ) const;

  /*!
  \brief Get all (scalar) contact forces of this contact pair
  */
  virtual void GetAllActiveContactForces( std::vector<double> & forces ) const;

  /*!
  \brief Get all (scalar) gap values of this contact pair
  */
  virtual void GetAllActiveContactGaps( std::vector<double> & gaps ) const;

  /*!
  \brief Get energy of penalty contact.
  */
  double GetEnergy() const override;

  //@}

private:

  //! @name member variables

  //! current node coordinates of the two elements
  LINALG::TMatrix<TYPE,3*numnodes*numnodalvalues,1> ele1pos_;
  LINALG::TMatrix<TYPE,3*numnodes*numnodalvalues,1> ele2pos_;

  //! neighbor elements of element 1
  Teuchos::RCP<CONTACT::B3CNeighbor> neighbors1_;

  //! neighbor elements of element 2
  Teuchos::RCP<CONTACT::B3CNeighbor> neighbors2_;

  //! cross section radius of first beam
  double R1_;

  //! cross section radius of second beam
  double R2_;

  //! Maximal gap at which a contact can become active
  double maxactivegap_;

  //! Maximal distance between a real segment on beam element 1 and its straight approximation
  double maxsegdist1_;

  //! Maximal distance between a real segment on beam element 2 and its straight approximation
  double maxsegdist2_;

  //! Number of segments on element1
  int numseg1_;

  //! Number of segments on element2
  int numseg2_;

  //! bound for search of large angle contact segment pairs
  double deltalargeangle_;

  //! bound for search of small angle contact segment pairs
  double deltasmallangle_;

  //! Indicates if the left / right node of the slave element 1 coincides with the endpoint of the physical beam (true) or not (false)
  std::pair<bool,bool> boundarynode1_;

  //! Indicates if the left / right node of the master element 2 coincides with the endpoint of the physical beam (true) or not (false)
  std::pair<bool,bool> boundarynode2_;

  //! Variables stored at the closest points of the large-angle-contact algorithm
  std::vector<Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues> > > cpvariables_;

  //! Variables stored at the Gauss points of the small-angle-contact algorithm
  std::vector<Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues> > > gpvariables_;

  //! Variables stored at the end points of the endpoint-contact algorithm
  std::vector<Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues> > > epvariables_;

  //@}

  //! @name Private evaluation methods

  /*!
  \brief Get active large angle pairs
  */
  void GetActiveLargeAnglePairs( std::vector<LINALG::TMatrix<double,3,1> >& endpoints1,
                                 std::vector<LINALG::TMatrix<double,3,1> >& endpoints2,
                                 std::map<std::pair<int,int>,LINALG::TMatrix<double,3,1> >& closelargeanglesegments,
                                 const double& pp);

  /*!
  \brief Evaluate active large angle pairs
  */
  void EvaluateActiveLargeAnglePairs(
      LINALG::SerialDenseVector* forcevec1,
      LINALG::SerialDenseVector* forcevec2,
      LINALG::SerialDenseMatrix* stiffmat11,
      LINALG::SerialDenseMatrix* stiffmat12,
      LINALG::SerialDenseMatrix* stiffmat21,
      LINALG::SerialDenseMatrix* stiffmat22);

  /*!
  \brief Get active small angle pairs
  */
  void GetActiveSmallAnglePairs(std::map<std::pair<int,int>,LINALG::TMatrix<double,3,1> >& closesmallanglesegments,
                                std::pair<int,int>* iminmax = NULL,
                                std::pair<bool,bool>* leftrightsolutionwithinsegment = NULL,
                                std::pair<double,double>* eta1_leftrightboundary = NULL);

  /*!
  \brief Evaluate active small angle pairs
  */
  void EvaluateActiveSmallAnglePairs(
      LINALG::SerialDenseVector* forcevec1,
      LINALG::SerialDenseVector* forcevec2,
      LINALG::SerialDenseMatrix* stiffmat11,
      LINALG::SerialDenseMatrix* stiffmat12,
      LINALG::SerialDenseMatrix* stiffmat21,
      LINALG::SerialDenseMatrix* stiffmat22,
      std::pair<int,int>* iminmax = NULL,
      std::pair<bool,bool>* leftrightsolutionwithinsegment = NULL,
      std::pair<double,double>* eta1_leftrightboundary = NULL);

  /*!
  \brief Get active endpoint pairs
  */
  void GetActiveEndPointPairs(std::vector<std::pair<int,int> >& closeendpointsegments, const double pp);

  /*!
  \brief Evaluate active endpoint pairs
  */
  void EvaluateActiveEndPointPairs(
      LINALG::SerialDenseVector* forcevec1,
      LINALG::SerialDenseVector* forcevec2,
      LINALG::SerialDenseMatrix* stiffmat11,
      LINALG::SerialDenseMatrix* stiffmat12,
      LINALG::SerialDenseMatrix* stiffmat21,
      LINALG::SerialDenseMatrix* stiffmat22);

  /*!
  \brief Find segments close to each other
  */
  void GetCloseSegments( const std::vector<LINALG::TMatrix<double,3,1> >& endpoints1,
                         const std::vector<LINALG::TMatrix<double,3,1> >& endpoints2,
                         std::map<std::pair<int,int>,LINALG::TMatrix<double,3,1> >& closesmallanglesegments,
                         std::map<std::pair<int,int>,LINALG::TMatrix<double,3,1> >& closelargeanglesegments,
                         std::vector<std::pair<int,int> >& closeendpointsegments,
                         double maxactivedist);

  /*!
  \brief Find contact point via closest point projection
  */
  bool ClosestPointProjection(double& eta_left1,
                              double& eta_left2,
                              double& l1,
                              double& l2,
                              LINALG::TMatrix<double,3,1>& segmentdata,
                              std::pair<TYPE,TYPE>& solutionpoints,
                              int segid1, int segid2);

  /*!
  \brief Find closest point eta2_master on a line for a given slave point eta1_slave
  */
  bool PointToLineProjection(double& eta1_slave,
                             double& eta_left2,
                             double& l2,
                             double& eta2_master,
                             double& gap,
                             double& alpha,
                             bool& pairactive,
                             bool smallanglepair,
                             bool invertpairs=false,
                             bool orthogonalprojection=false);

  /*!
  \brief Determine minimal distance and contact angle for unconverged segment pair
  */
  void CheckUnconvergedSegmentPair( double& eta_left1,
                                    double& eta_left2,
                                    double& l1,
                                    double& l2,
                                    double& eta1_min,
                                    double& eta2_min,
                                    double& g_min,
                                    double& alpha_g_min,
                                    bool& pointtolinesolfound);

  /*!
  \brief Subdivide elements into segments for CPP
  */
  double CreateSegments(const DRT::Element* ele,
                        std::vector<LINALG::TMatrix<double,3,1> >& endpoints_final,
                        int& numsegment,
                        int i);

  /*!
  \brief Get maximal gap at which a contact can become active
  */
  double GetMaxActiveDist();

  /*!
  \brief Check, if segments are fine enough
  */
  bool CheckSegment(LINALG::TMatrix<double,3,1>& r1,
                    LINALG::TMatrix<double,3,1>& t1,
                    LINALG::TMatrix<double,3,1>& r2,
                    LINALG::TMatrix<double,3,1>& t2,
                    LINALG::TMatrix<double,3,1>& rm,
                    double& segdist);

  /*!
  \brief Calculate scalar contact force
  */
  void CalcPenaltyLaw(Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues> > variables);

  /*!
  \brief Calculate angle-dependent penalty scale factor for large-angle-contact
  */
  void CalcPerpPenaltyScaleFac( Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues> > cpvariables,
                            LINALG::TMatrix<TYPE,3,1>& r1_xi,
                            LINALG::TMatrix<TYPE,3,1>& r2_xi,
                            const double shiftangle1,
                            const double shiftangle2);

  /*!
   * Todo redundant ?!
  \brief Calculate angle-dependent penalty scale factor for small-angle-contact
  */
  void CalcParPenaltyScaleFac(Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues> > gpvariables,
                              LINALG::TMatrix<TYPE,3,1>& r1_xi,
                              LINALG::TMatrix<TYPE,3,1>& r2_xi,
                              const double shiftangle1,
                              const double shiftangle2);

  /*!
   \brief Compute contact forces
   */
  void EvaluateFcContact(
      LINALG::SerialDenseVector& forcevec1,
      LINALG::SerialDenseVector& forcevec2,
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
      Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues> > variables,
      const double& intfac,
      bool cpp,
      bool gp,
      bool fixedendpointxi,
      bool fixedendpointeta,
      LINALG::TMatrix<TYPE, 3*numnodes*numnodalvalues, 1>* fc1_FAD=NULL,
      LINALG::TMatrix<TYPE, 3*numnodes*numnodalvalues, 1>* fc2_FAD=NULL);

  /*!
  \brief Evaluate contact stiffness
  */
  void EvaluateStiffcContact(
      LINALG::SerialDenseMatrix& stiffmat11,
      LINALG::SerialDenseMatrix& stiffmat12,
      LINALG::SerialDenseMatrix& stiffmat21,
      LINALG::SerialDenseMatrix& stiffmat22,
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
      const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xixi,
      Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues> > variables,
      const double& intfac,
      bool cpp,
      bool gp,
      bool fixedendpointxi,
      bool fixedendpointeta);

#ifdef ENDPOINTSEGMENTATION
  /*!
  \brief FAD-based Evaluation of contact stiffness in case of ENDPOINTSEGMENTATION
  */
  void EvaluateStiffcContactIntSeg(LINALG::SparseMatrix& stiffmatrix,
                                   const LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_xi_bound,
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
                                   Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues> > cpvariables,
                                   const double& intfac,
                                   const double& d_xi_ele_d_xi_bound,
                                   TYPE signed_jacobi_interval);
#endif

  /*!
  \brief Linearizations of contact point
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
  \brief Lin. of contact point coordinate eta with fixed xi
  */
  void ComputeLinEtaFixXi(LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_eta,
                          const LINALG::TMatrix<TYPE, 3, 1>& delta_r,
                          const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                          const LINALG::TMatrix<TYPE, 3, 1>& r2_xixi,
                          const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                          const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2,
                          const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xi);

  /*!
  \brief Lin. of contact point coordinate xi with fixed eta
  */
  void ComputeLinXiFixEta(LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_eta,
                          const LINALG::TMatrix<TYPE, 3, 1>& delta_r,
                          const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                          const LINALG::TMatrix<TYPE, 3, 1>& r1_xixi,
                          const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2,
                          const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                          const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1_xi);

  /*!
  \brief Compute linearization of integration interval bounds (necessary in case of ENDPOINTSEGMENTATION)
  */
  void ComputeLinXiBound(LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_xi_bound,
                         TYPE& eta1_bound,
                         TYPE eta2);

  /*!
  \brief Compute linearization of gap
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
  \brief Compute linearization of cosine of contact angle
  */
  void ComputeLinCosContactAngle(LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_coscontactangle,
                                 LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_xi,
                                 LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_eta,
                                 const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                                 const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                                 const LINALG::TMatrix<TYPE, 3, 1>& r1_xixi,
                                 const LINALG::TMatrix<TYPE, 3, 1>& r2_xixi,
                                 const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1_xi,
                                 const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xi);

  /*!
  \brief Compute linearization of normal vector
  */
  void ComputeLinNormal(LINALG::TMatrix<TYPE, 3, 2*3*numnodes*numnodalvalues>& delta_normal,
                        const LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>&  delta_xi,
                        const LINALG::TMatrix<TYPE, 2*3*numnodes*numnodalvalues, 1>& delta_eta,
                        const LINALG::TMatrix<TYPE, 3, 1>& delta_r,
                        const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                        const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                        const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                        const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2);

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
  \brief Calculate one specified shape function value / derivative for given parameter value and element
  */
  void GetShapeFunctions( LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N,
                          const TYPE& eta,
                          int deriv,
                          const DRT::Element* ele) const;

  /*!
  \brief compute coordinate at given curve point
  */
  LINALG::TMatrix<TYPE,3,1> r(const TYPE& eta, const DRT::Element* ele) const;

  /*!
  \brief compute derivative at given curve point
  */
  LINALG::TMatrix<TYPE,3,1> r_xi(const TYPE& eta, const DRT::Element* ele);

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
                                          const LINALG::TMatrix<TYPE,3,1>& t2_xi,
                                          bool& elementscolinear);

  /*!
  \brief Evaluate orthogonality cond. of point to line projeciton
  */
  void EvaluatePTLOrthogonalityCondition( TYPE& f,
                                          const LINALG::TMatrix<TYPE,3,1>& delta_r,
                                          const double norm_delta_r,
                                          const LINALG::TMatrix<TYPE,3,1>& r1_xi,
                                          const LINALG::TMatrix<TYPE,3,1>& r2_xi,
                                          bool orthogonalprojection);

  /*!
  \brief Evaluate Jacobian df of PTLOrthogonalityCondition
  */
  bool EvaluateLinPTLOrthogonalityCondition(TYPE& df,
                                            const LINALG::TMatrix<TYPE,3,1>& delta_r,
                                            const double norm_delta_r,
                                            const LINALG::TMatrix<TYPE,3,1>& r1_xi,
                                            const LINALG::TMatrix<TYPE,3,1>& r2_xi,
                                            const LINALG::TMatrix<TYPE,3,1>& r2_xixi,
                                            bool orthogonalprojection);

  /*!
  \brief Compute normal vector and gap function at contact point
  */
  void ComputeNormal( LINALG::TMatrix<TYPE, 3, 1>& r1,
                      LINALG::TMatrix<TYPE, 3, 1>& r2,
                      LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                      LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                      Teuchos::RCP<BeamToBeamContactVariables<numnodes, numnodalvalues> > variables,
                      int contacttype);

  /*!
  \brief Check, if we have contact or not (e.g. gap < gmax [e.g. gmax=0]?)
  */
  bool CheckContactStatus(const double& gap);

  /*!
    \brief Get jacobi factor of beam element
  */
  double GetJacobi(const DRT::Element* element1);

  /** \brief get Jacobi factor of beam element at xi \in [-1;1]
   *
   *  \author grill
   *  \date 06/16 */
  inline double GetJacobiAtXi(const DRT::Element* element1, const double& xi)
  {
    const DRT::ELEMENTS::Beam3Base* ele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(element1);

    if (ele==NULL)
      dserror("Dynamic cast to Beam3Base failed");

    return ele->GetJacobiFacAtXi(xi);
  }

  /*!
    \brief clear class variables at the beginning of a Newton step
  */
  void ClearClassVariables();

  /*!
    \brief Linearization-check of coordinates xi and eta via FAD
  */
  void FADCheckLinXiAndLinEta(const LINALG::TMatrix<TYPE, 3, 1>& delta_r,
                              const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                              const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                              const LINALG::TMatrix<TYPE, 3, 1>& r1_xixi,
                              const LINALG::TMatrix<TYPE, 3, 1>& r2_xixi,
                              const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1,
                              const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2,
                              const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N1_xi,
                              const LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues>& N2_xi);

  /*!
    \brief Linearization-check for local Newton in CPP via FAD
  */
  void FADCheckLinOrthogonalityCondition( const LINALG::TMatrix<TYPE, 3, 1>& delta_r,
                                          const double& norm_delta_r,
                                          const LINALG::TMatrix<TYPE, 3, 1>& r1_xi,
                                          const LINALG::TMatrix<TYPE, 3, 1>& r2_xi,
                                          const LINALG::TMatrix<TYPE, 3, 1>& t1,
                                          const LINALG::TMatrix<TYPE, 3, 1>& t2);

//  /*!
//    \brief FD-Check of stiffness matrix
//  */
//  void FDCheck( LINALG::SparseMatrix& stiffmatrix,
//                Epetra_Vector& fint,
//                const double& pp,
//                std::map<std::pair<int,int>, Teuchos::RCP<BeamContactPair> >& contactpairmap,
//                Teuchos::ParameterList& timeintparams,
//                bool fdcheck);

  //@}

};
}

#endif
