/*----------------------------------------------------------------------------*/
/*!
\file beam3contact.H

\brief One beam contact pair (two beam elements) consisting of several contact segments

\level 2

\maintainer Christoph Meier
*/
/*----------------------------------------------------------------------------*/

#ifndef BEAM3C_H
#define BEAM3C_H

#include "beam3contact_utils.H"
#include "beam3contactvariables.H"
#include "beam3contact_tangentsmoothing.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_sparsematrix.H"
#include "beam3contact_defines.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "Sacado.hpp"
#include "../headers/FAD_utils.H"
#include "beam3contactinterface.H"
#include "../drt_inpar/inpar_beamcontact.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_beam3/beam3_base.H"

#include "../drt_lib/standardtypes_cpp.H"

namespace CONTACT
{

/*!
 \brief contact element for contact between two 3D beam elements

 Refer also to the Semesterarbeit of Matthias Mayr, 2010

 \author meier (meier@lnm.mw.tum.de)
 */

template<const int numnodes , const int numnodalvalues>
class Beam3contact: public Beam3contactinterface
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
  Beam3contact(  const DRT::Discretization& pdiscret,
                    const DRT::Discretization& cdiscret,
                    const std::map<int,int>& dofoffsetmap,
                    DRT::Element* element1,
                    DRT::Element* element2,
                    Teuchos::ParameterList& beamcontactparams);

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
  \brief Get number of standard large angle/small angle/endpoint contact points on this element pair
  */
  int GetNumCps()
  {
    return cpvariables_.size();
  };

  int GetNumGps()
  {
    return gpvariables_.size();
  };

  int GetNumEps()
  {
    return epvariables_.size();
  };

  /*!
  \brief Get vector of type declarations (0=closest point contact, 1=gauss point contact, 2= end point contact) of all contact pairs
  */
  std::vector<int> GetContactType()
  {
    int size1=cpvariables_.size();
    int size2=gpvariables_.size();
    int size3=epvariables_.size();
    std::vector<int > types(size1+size2+size3, 0);

    for (int i=0;i<size1;i++)
    {
      types[i] = 0;
    }

    for (int i=size1;i<size2+size1;i++)
    {
      types[i] = 1;
    }

    for (int i=size1+size2;i<size1+size2+size3;i++)
    {
      types[i] = 2;
    }

    return types;
  };

  /*!
  \brief Get vector of all gaps of this contact pair
  */
  std::vector<double> GetGap()
  {
    int size1=cpvariables_.size();
    int size2=gpvariables_.size();
    int size3=epvariables_.size();
    std::vector<double > gaps(size1+size2+size3, 0.0);

    for (int i=0;i<size1;i++)
    {
      gaps[i] = FADUTILS::CastToDouble(cpvariables_[i]->GetGap());
    }

    for (int i=size1;i<size2+size1;i++)
    {
      gaps[i] = FADUTILS::CastToDouble(gpvariables_[i-size1]->GetGap());
    }

    for (int i=size1+size2;i<size1+size2+size3;i++)
    {
      gaps[i] = FADUTILS::CastToDouble(epvariables_[i-size1-size2]->GetGap());
    }

    return gaps;
  };

  /*!
  \brief Get vector of all contact forces of this contact pair
  */
  std::vector<double> GetContactForce()
  {
    int size1=cpvariables_.size();
    int size2=gpvariables_.size();
    int size3=epvariables_.size();
    std::vector<double > f(size1+size2+size3, 0.0);

    for (int i=0;i<size1;i++)
    {
      f[i] = FADUTILS::CastToDouble(cpvariables_[i]->Getfp()*cpvariables_[i]->GetPPfac());
    }

    for (int i=size1;i<size2+size1;i++)
    {
      f[i] = FADUTILS::CastToDouble(gpvariables_[i-size1]->Getfp()*gpvariables_[i-size1]->GetPPfac());
    }

    for (int i=size1+size2;i<size1+size2+size3;i++)
    {
      f[i] = FADUTILS::CastToDouble(epvariables_[i-size1-size2]->Getfp()*epvariables_[i-size1-size2]->GetPPfac());
    }

    return f;
  };

  /*!
  \brief Get vector of all contact angles of this contact pair
  */
  std::vector<double> GetContactAngle()
  {
    int size1=cpvariables_.size();
    int size2=gpvariables_.size();
    int size3=epvariables_.size();
    std::vector<double > angles(size1+size2+size3, 0.0);

    for (int i=0;i<size1;i++)
    {
      angles[i] = FADUTILS::CastToDouble(cpvariables_[i]->GetAngle());
    }

    for (int i=size1;i<size2+size1;i++)
    {
      angles[i] = FADUTILS::CastToDouble(gpvariables_[i-size1]->GetAngle());
    }

    for (int i=size1+size2;i<size1+size2+size3;i++)
    {
      angles[i] = FADUTILS::CastToDouble(epvariables_[i-size1-size2]->GetAngle());
    }

    return angles;
  };

  /*!
  \brief Get vector of all closest points of this contact pair
  */
  std::vector<std::pair<double,double> > GetClosestPoint()
  {
    int size1=cpvariables_.size();
    int size2=gpvariables_.size();
    int size3=epvariables_.size();
    std::vector<std::pair<double,double> > cps(size1+size2+size3,std::make_pair(0.0,0.0));

    for (int i=0;i<size1;i++)
    {
      double xi=FADUTILS::CastToDouble(cpvariables_[i]->GetCP().first);
      double eta=FADUTILS::CastToDouble(cpvariables_[i]->GetCP().second);
      cps[i]=std::make_pair(xi,eta);
    }

    for (int i=size1;i<size2+size1;i++)
    {
      double xi=FADUTILS::CastToDouble(gpvariables_[i-size1]->GetCP().first);
      double eta=FADUTILS::CastToDouble(gpvariables_[i-size1]->GetCP().second);
      cps[i]=std::make_pair(xi,eta);
    }

    for (int i=size1+size2;i<size1+size2+size3;i++)
    {
      double xi=FADUTILS::CastToDouble(epvariables_[i-size1-size2]->GetCP().first);
      double eta=FADUTILS::CastToDouble(epvariables_[i-size1-size2]->GetCP().second);
      cps[i]=std::make_pair(xi,eta);
    }

    return cps;
  };

  /*!
  \brief Return number of individual contact segments on element pair
  */
  std::pair<int,int> GetNumSegments()
  {
    return std::make_pair(numseg1_,numseg2_);
  };

  /*!
  \brief Return ids of active segments
  */
  std::vector<std::pair<int,int> > GetSegmentIds()
  {
    int size1=cpvariables_.size();
    int size2=gpvariables_.size();
    int size3=epvariables_.size();
    std::vector<std::pair<int,int> > ids(size1+size2+size3,std::make_pair(1,1));

    for (int i=0;i<size1;i++)
    {
      ids[i]=cpvariables_[i]->GetSegIds();
    }

    for (int i=size1;i<size2+size1;i++)
    {
      ids[i]=gpvariables_[i-size1]->GetSegIds();
    }

    for (int i=size1+size2;i<size1+size2+size3;i++)
    {
      ids[i]=epvariables_[i-size1-size2]->GetSegIds();
    }

    return ids;
  };

  /*!
  \brief Get current augmented lagrange multiplier 'lmuzawa_'
  */
  double Getlmuzawa()
  {
    return 0.0;
  };

  /*!
  \brief Get flag indicating whether contact is active (true) or inactive (false)
  */
  bool GetContactFlag()
  {
    //The element pair is assumed to be active when we have at least one active contact point
    return (cpvariables_.size()+gpvariables_.size()+epvariables_.size());
  };

  /*!
  \brief Get coordinates of contact point of element1
  */
  const std::vector<LINALG::Matrix<3,1> > GetX1()
  {
    int size1=cpvariables_.size();
    int size2=gpvariables_.size();
    int size3=epvariables_.size();
    std::vector<LINALG::Matrix<3,1> > r1(size1+size2+size3, LINALG::Matrix<3,1>(true));

    for (int i=0;i<size1;i++)
    {
      TYPE eta1 = cpvariables_[i]->GetCP().first;
      for(int j=0;j<3;j++)
        r1[i](j) = FADUTILS::CastToDouble(r(eta1,element1_)(j));
    }

    for (int i=size1;i<size2+size1;i++)
    {
      TYPE eta1 = gpvariables_[i-size1]->GetCP().first;
      for(int j=0;j<3;j++)
        r1[i](j) = FADUTILS::CastToDouble(r(eta1,element1_)(j));
    }

    for (int i=size1+size2;i<size1+size2+size3;i++)
    {
      TYPE eta1 = epvariables_[i-size1-size2]->GetCP().first;
      for(int j=0;j<3;j++)
        r1[i](j) = FADUTILS::CastToDouble(r(eta1,element1_)(j));
    }

    return r1;
  };

  /*!
  \brief Get coordinates of contact point of element2
  */
  const std::vector<LINALG::Matrix<3,1> > GetX2()
  {
    int size1=cpvariables_.size();
    int size2=gpvariables_.size();
    int size3=epvariables_.size();
    std::vector<LINALG::Matrix<3,1> > r2(size1+size2+size3, LINALG::Matrix<3,1>(true));

    for (int i=0;i<size1;i++)
    {
      TYPE eta2 = cpvariables_[i]->GetCP().second;
      for(int j=0;j<3;j++)
        r2[i](j) = FADUTILS::CastToDouble(r(eta2,element2_)(j));
    }

    for (int i=size1;i<size2+size1;i++)
    {
      TYPE eta2 = gpvariables_[i-size1]->GetCP().second;
      for(int j=0;j<3;j++)
        r2[i](j) = FADUTILS::CastToDouble(r(eta2,element2_)(j));
    }

    for (int i=size1+size2;i<size1+size2+size3;i++)
    {
      TYPE eta2 = epvariables_[i-size1-size2]->GetCP().second;
      for(int j=0;j<3;j++)
        r2[i](j) = FADUTILS::CastToDouble(r(eta2,element2_)(j));
    }

    return r2;
  };


  //TODO:
  /*!
  \brief Get normal vector
  */
  const Epetra_SerialDenseVector GetNormal()
  {
    Epetra_SerialDenseVector normal(3);

    for(int i=0;i<3;i++)
      normal(i)=0.0;

    return normal;
  }

  /*!
  \brief Get flag indicating whether the nodal values of one element had been shifted due to r1=r2
         Since this is only possible for beam3contactnew elements but not for beam3contact elements
         we always return false within this class.
  */
  bool GetShiftStatus() {return false;};

  /*!
    \Check, if there is a difference between the result of the new and old gap definition, i.e. if the beams centerlines have already crossed or not.
    Since this is only possible for beam3contactnew elements but not for beam3contact elements
    we always return false within this class.
  */
  bool GetNewGapStatus()
  {
    return false;
  };
  //@}

  /*!
    \Get energy of penalty contact.
  */
  double GetEnergy()
  {
    if(     DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(bcparams_,"BEAMS_PENALTYLAW") != INPAR::BEAMCONTACT::pl_lp
        and DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(bcparams_,"BEAMS_PENALTYLAW") != INPAR::BEAMCONTACT::pl_qp
        and DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(bcparams_,"BEAMS_PENALTYLAW") != INPAR::BEAMCONTACT::pl_lpqp)
      dserror("Contact Energy calculation not implemented for the chosen penalty law!");


    double energy=0.0;

    for (int i=0;i<(int)cpvariables_.size();i++)
    {
      double ppfac = FADUTILS::CastToDouble(cpvariables_[i]->GetPPfac());
      double e = -cpvariables_[i]->GetIntegratedEnergy();
      energy += ppfac*e;
    }

    for (int i=0;i<(int)gpvariables_.size();i++)
    {
      double ppfac = FADUTILS::CastToDouble(gpvariables_[i]->GetPPfac());
      double e = -gpvariables_[i]->GetIntegratedEnergy();
      energy += ppfac*e;
    }

    for (int i=0;i<(int)epvariables_.size();i++)
    {
      double ppfac = FADUTILS::CastToDouble(epvariables_[i]->GetPPfac());
      double e = -epvariables_[i]->GetIntegratedEnergy();
      energy += ppfac*e;
    }

    return energy;
  };

  /*!
    \Get energy of perp penalty contact without transition factor contribution.
  */
  double GetUnscaledPerpEnergy()
  {
    if(     DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(bcparams_,"BEAMS_PENALTYLAW") != INPAR::BEAMCONTACT::pl_lp
        and DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(bcparams_,"BEAMS_PENALTYLAW") != INPAR::BEAMCONTACT::pl_qp
        and DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(bcparams_,"BEAMS_PENALTYLAW") != INPAR::BEAMCONTACT::pl_lpqp)
      dserror("Contact Energy calculation not implemented for the chosen penalty law!");


    double energy=0.0;

    for (int i=0;i<(int)cpvariables_.size();i++)
    {
      double e = -cpvariables_[i]->GetIntegratedEnergy();
      energy += e;
    }

    return energy;
  };

  /*!
    \Get energy of parallel penalty contact without transition factor contribution.
  */
  double GetUnscaledParallelEnergy()
  {
    if(     DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(bcparams_,"BEAMS_PENALTYLAW") != INPAR::BEAMCONTACT::pl_lp
        and DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(bcparams_,"BEAMS_PENALTYLAW") != INPAR::BEAMCONTACT::pl_qp
        and DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(bcparams_,"BEAMS_PENALTYLAW") != INPAR::BEAMCONTACT::pl_lpqp)
      dserror("Contact Energy calculation not implemented for the chosen penalty law!");

    double energy=0.0;

    for (int i=0;i<(int)gpvariables_.size();i++)
    {
      double e = -gpvariables_[i]->GetIntegratedEnergy();
      energy += e;
    }

    return energy;
  };

  //TODO
  /*!
    \We don't need this method for beam3contact elements!
  */
  LINALG::TMatrix<TYPE, 3, 1>* GetNormalOld()
  {
      return NULL;
  };

  //TODO
  /*!
    \We don't need this method for beam3contact elements!
  */
  bool FirstTimeStep()
  {
    return false;
  };
  //@}

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
  \brief Reset Uzawa-based Lagrange multiplier
  */
  void Resetlmuzawa();

  /*!
  \brief Change the sign of the normal vector: This has to be done at the end of a time step when the remainig penetration is
  larger that the sum of the beam radii (R1+R2). Otherwise, the beams could cross in the next time step when the new gap function
  definition (ngf_=true) for slender beams is applied!
  */
  void InvertNormal(){dserror("Function not implemented!");};

  //TODO
  /*!
    \brief We don't need this method for beam3contact elements!
  */
  void UpdateClassVariablesStep(){};

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

  /** \brief print information about this beam contact element pair to screen
   *
   *  \author grill
   *  \date 05/16 */
  void Print() const;

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

  //! neighbor elements of element 1
  Teuchos::RCP<B3CNeighbor> neighbors1_;

  //! neighbor elements of element 2
  Teuchos::RCP<B3CNeighbor> neighbors2_;

  //! averaged nodal tangents, necessary for smoothed tangent fields of C^0 Reissner beams
  LINALG::Matrix<3*numnodes,1> nodaltangentssmooth1_;
  LINALG::Matrix<3*numnodes,1> nodaltangentssmooth2_;

  //!current Newton iteration
  int iter_;

  //!current time step
  int numstep_;

  //! cross section radius of first beam
  const double R1_;

  //! cross section radius of second beam
  const double R2_;

  //! Maximal gap at which a contact can become active
  const double maxactivegap_;

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
  std::vector<Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > > cpvariables_;

  //! Variables stored at the Gauss points of the small-angle-contact algorithm
  std::vector<Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > > gpvariables_;

  //! Variables stored at the end points of the endpoint-contact algorithm
  std::vector<Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > > epvariables_;

  //@}

  //! @name Private evaluation methods

  /*!
  \brief Get active large angle pairs
  */
  void GetActiveLargeAnglePairs( std::vector<LINALG::TMatrix<double,3,1> >& endpoints1,
                                 std::vector<LINALG::TMatrix<double,3,1> >& endpoints2,
                                 std::map<std::pair<int,int>,LINALG::TMatrix<double,3,1> >& closelargeanglesegments,
                                 const double pp);

  /*!
  \brief Evaluate active large angle pairs
  */
  void EvaluateActiveLargeAnglePairs(LINALG::SparseMatrix& stiffmatrix,
                                     Epetra_Vector& fint);

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
  void EvaluateActiveSmallAnglePairs(LINALG::SparseMatrix& stiffmatrix,
                                     Epetra_Vector& fint,
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
  void EvaluateActiveEndPointPairs(LINALG::SparseMatrix& stiffmatrix,
                                     Epetra_Vector& fint);

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
  double CreateSegments(DRT::Element* ele,
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
  void CalcPenaltyLaw(Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > variables);

  /*!
  \brief Calculate angle-dependent penalty scale factor for small-angle-contact
  */
  void CalcPerpPenaltyScaleFac( Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > cpvariables,
                            LINALG::TMatrix<TYPE,3,1>& r1_xi,
                            LINALG::TMatrix<TYPE,3,1>& r2_xi,
                            const double shiftangle1,
                            const double shiftangle2);

  /*!
  \brief Calculate angle-dependent penalty scale factor for large-angle-contact
  */
  void CalcParPenaltyScaleFac(Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > gpvariables,
                              LINALG::TMatrix<TYPE,3,1>& r1_xi,
                              LINALG::TMatrix<TYPE,3,1>& r2_xi,
                              const double shiftangle1,
                              const double shiftangle2);

  /*!
   \brief Compute contact forces
   */
  void EvaluateFcContact( Epetra_Vector* fint,
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
                          Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > variables,
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
  void EvaluateStiffcContact(LINALG::SparseMatrix& stiffmatrix,
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
                             Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > variables,
                             const double& intfac,
                             bool cpp,
                             bool gp,
                             bool fixedendpointxi,
                             bool fixedendpointeta);

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
                                   Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > cpvariables,
                                   const double& intfac,
                                   const double& d_xi_ele_d_xi_bound,
                                   TYPE signed_jacobi_interval);

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
                          DRT::Element* ele);

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
  \brief Assemble shape functions for one given matrix
  */
  void AssembleShapefunctions(const LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues>& N_i,
                              LINALG::TMatrix<TYPE,3,3*numnodes*numnodalvalues>& N);

  /*!
  \brief compute coordinate at given curve point
  */
  LINALG::TMatrix<TYPE,3,1> r(const TYPE& eta, DRT::Element* ele);

  /*!
  \brief compute derivative at given curve point
  */
  LINALG::TMatrix<TYPE,3,1> r_xi(const TYPE& eta, DRT::Element* ele);

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
                      Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues> > variables,
                      int contacttype);

  /*!
  \brief Check, if we have contact or not (e.g. gap < gmax [e.g. gmax=0]?)
  */
  bool CheckContactStatus(const double& gap);

  /*!
  \brief Check, if we have contact or not (e.g. gap < gdmax?)
  */
  bool CheckDampingStatus(const double& gap);

  /*!
  \brief Get global dofs of a node

  Internally this method first extracts the dofs of the given node
  in the beam contact discretization (which has its own dofs) and
  then transfers these dofs to their actual GIDs in the underlying
  problem discretization by applying the pre-computed dofoffset_.
  */
  std::vector<int> GetGlobalDofs(const DRT::Node* node);

  /*!
    \brief Get jacobi factor of beam element
  */
  double GetJacobi(DRT::Element* element1);

  /** \brief get Jacobi factor of beam element at xi \in [-1;1]
   *
   *  \author grill
   *  \date 06/16 */
  inline double GetJacobiAtXi(DRT::Element* element1, const double& xi)
  {
    const DRT::ELEMENTS::Beam3Base* ele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(element1);

    if (ele==NULL)
      dserror("Dynamic cast to Beam3Base failed");

    return ele->GetJacobiFacAtXi(xi);
  }

  /*!
    \brief Set class variables at the beginning of a Newton step
  */
  void SetClassVariables(Teuchos::ParameterList& timeintparams);

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

  /*!
    \brief FD-Check of stiffness matrix
  */
  void FDCheck( LINALG::SparseMatrix& stiffmatrix,
                Epetra_Vector& fint,
                const double& pp,
                std::map<std::pair<int,int>, Teuchos::RCP<Beam3contactinterface > >& contactpairmap,
                Teuchos::ParameterList& timeintparams,
                bool fdcheck);

  //@}

}; // class Beam3contact
} // namespace CONTACT

#endif  // #ifndef BEAM3C_H
