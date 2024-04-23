/*----------------------------------------------------------------------------*/
/*! \file

\brief One beam contact pair (two beam elements) consisting of several contact segments

\level 2

*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_BEAMCONTACT_BEAM3CONTACT_HPP
#define FOUR_C_BEAMCONTACT_BEAM3CONTACT_HPP

#include "baci_config.hpp"

#include "baci_beam3_base.hpp"
#include "baci_beamcontact_beam3contactinterface.hpp"
#include "baci_beamcontact_beam3contactvariables.hpp"
#include "baci_beaminteraction_beam_to_beam_contact_defines.hpp"
#include "baci_beaminteraction_beam_to_beam_contact_tangentsmoothing.hpp"
#include "baci_beaminteraction_beam_to_beam_contact_utils.hpp"
#include "baci_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "baci_inpar_beamcontact.hpp"
#include "baci_inpar_contact.hpp"
#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_linalg_sparsematrix.hpp"
#include "baci_utils_fad.hpp"

#include <Sacado.hpp>

FOUR_C_NAMESPACE_OPEN


namespace CONTACT
{
  /*!
   \brief contact element for contact between two 3D beam elements

   Refer also to the Semesterarbeit of Matthias Mayr, 2010

   */

  template <const int numnodes, const int numnodalvalues>
  class Beam3contact : public Beam3contactinterface
  {
   public:
    //! @name Friends
    // no friend classes defined
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
    Beam3contact(const DRT::Discretization& pdiscret, const DRT::Discretization& cdiscret,
        const std::map<int, int>& dofoffsetmap, DRT::Element* element1, DRT::Element* element2,
        Teuchos::ParameterList& beamcontactparams);


    //@}


    //! @name Access methods
    /*!
    \brief Get problem discretization
    */
    inline const DRT::Discretization& ProblemDiscret() const override { return pdiscret_; };

    /*!
    \brief Get beam contact discretization
    */
    inline const DRT::Discretization& ContactDiscret() const override { return cdiscret_; };

    /*!
    \brief Get offset of dofs between cdiscret and pdiscret
    */
    inline const std::map<int, int>& DofOffset() const override { return dofoffsetmap_; };

    /*!
    \brief Get first element
    */
    inline const DRT::Element* Element1() override { return element1_; };

    /*!
    \brief Get first element
    */
    inline const DRT::Element* Element2() override { return element2_; };

    /*!
    \brief Get number of standard large angle/small angle/endpoint contact points on this element
    pair
    */
    int GetNumCps() override { return cpvariables_.size(); };

    int GetNumGps() override { return gpvariables_.size(); };

    int GetNumEps() override { return epvariables_.size(); };

    /*!
    \brief Get vector of type declarations (0=closest point contact, 1=gauss point contact, 2= end
    point contact) of all contact pairs
    */
    std::vector<int> GetContactType() override
    {
      int size1 = cpvariables_.size();
      int size2 = gpvariables_.size();
      int size3 = epvariables_.size();
      std::vector<int> types(size1 + size2 + size3, 0);

      for (int i = 0; i < size1; i++)
      {
        types[i] = 0;
      }

      for (int i = size1; i < size2 + size1; i++)
      {
        types[i] = 1;
      }

      for (int i = size1 + size2; i < size1 + size2 + size3; i++)
      {
        types[i] = 2;
      }

      return types;
    };

    /*!
    \brief Get vector of all gaps of this contact pair
    */
    std::vector<double> GetGap() override
    {
      int size1 = cpvariables_.size();
      int size2 = gpvariables_.size();
      int size3 = epvariables_.size();
      std::vector<double> gaps(size1 + size2 + size3, 0.0);

      for (int i = 0; i < size1; i++)
      {
        gaps[i] = CORE::FADUTILS::CastToDouble(cpvariables_[i]->GetGap());
      }

      for (int i = size1; i < size2 + size1; i++)
      {
        gaps[i] = CORE::FADUTILS::CastToDouble(gpvariables_[i - size1]->GetGap());
      }

      for (int i = size1 + size2; i < size1 + size2 + size3; i++)
      {
        gaps[i] = CORE::FADUTILS::CastToDouble(epvariables_[i - size1 - size2]->GetGap());
      }

      return gaps;
    };

    /*!
    \brief Get vector of all contact forces of this contact pair
    */
    std::vector<double> GetContactForce() override
    {
      int size1 = cpvariables_.size();
      int size2 = gpvariables_.size();
      int size3 = epvariables_.size();
      std::vector<double> f(size1 + size2 + size3, 0.0);

      for (int i = 0; i < size1; i++)
      {
        f[i] = CORE::FADUTILS::CastToDouble(cpvariables_[i]->Getfp() * cpvariables_[i]->GetPPfac());
      }

      for (int i = size1; i < size2 + size1; i++)
      {
        f[i] = CORE::FADUTILS::CastToDouble(
            gpvariables_[i - size1]->Getfp() * gpvariables_[i - size1]->GetPPfac());
      }

      for (int i = size1 + size2; i < size1 + size2 + size3; i++)
      {
        f[i] = CORE::FADUTILS::CastToDouble(
            epvariables_[i - size1 - size2]->Getfp() * epvariables_[i - size1 - size2]->GetPPfac());
      }

      return f;
    };

    /*!
    \brief Get vector of all contact angles of this contact pair
    */
    std::vector<double> GetContactAngle() override
    {
      int size1 = cpvariables_.size();
      int size2 = gpvariables_.size();
      int size3 = epvariables_.size();
      std::vector<double> angles(size1 + size2 + size3, 0.0);

      for (int i = 0; i < size1; i++)
      {
        angles[i] = CORE::FADUTILS::CastToDouble(cpvariables_[i]->GetAngle());
      }

      for (int i = size1; i < size2 + size1; i++)
      {
        angles[i] = CORE::FADUTILS::CastToDouble(gpvariables_[i - size1]->GetAngle());
      }

      for (int i = size1 + size2; i < size1 + size2 + size3; i++)
      {
        angles[i] = CORE::FADUTILS::CastToDouble(epvariables_[i - size1 - size2]->GetAngle());
      }

      return angles;
    };

    /*!
    \brief Get vector of all closest points of this contact pair
    */
    std::vector<std::pair<double, double>> GetClosestPoint() override
    {
      int size1 = cpvariables_.size();
      int size2 = gpvariables_.size();
      int size3 = epvariables_.size();
      std::vector<std::pair<double, double>> cps(size1 + size2 + size3, std::make_pair(0.0, 0.0));

      for (int i = 0; i < size1; i++)
      {
        double xi = CORE::FADUTILS::CastToDouble(cpvariables_[i]->GetCP().first);
        double eta = CORE::FADUTILS::CastToDouble(cpvariables_[i]->GetCP().second);
        cps[i] = std::make_pair(xi, eta);
      }

      for (int i = size1; i < size2 + size1; i++)
      {
        double xi = CORE::FADUTILS::CastToDouble(gpvariables_[i - size1]->GetCP().first);
        double eta = CORE::FADUTILS::CastToDouble(gpvariables_[i - size1]->GetCP().second);
        cps[i] = std::make_pair(xi, eta);
      }

      for (int i = size1 + size2; i < size1 + size2 + size3; i++)
      {
        double xi = CORE::FADUTILS::CastToDouble(epvariables_[i - size1 - size2]->GetCP().first);
        double eta = CORE::FADUTILS::CastToDouble(epvariables_[i - size1 - size2]->GetCP().second);
        cps[i] = std::make_pair(xi, eta);
      }

      return cps;
    };

    /*!
    \brief Return number of individual contact segments on element pair
    */
    std::pair<int, int> GetNumSegments() override { return std::make_pair(numseg1_, numseg2_); };

    /*!
    \brief Return ids of active segments
    */
    std::vector<std::pair<int, int>> GetSegmentIds() override
    {
      int size1 = cpvariables_.size();
      int size2 = gpvariables_.size();
      int size3 = epvariables_.size();
      std::vector<std::pair<int, int>> ids(size1 + size2 + size3, std::make_pair(1, 1));

      for (int i = 0; i < size1; i++)
      {
        ids[i] = cpvariables_[i]->GetSegIds();
      }

      for (int i = size1; i < size2 + size1; i++)
      {
        ids[i] = gpvariables_[i - size1]->GetSegIds();
      }

      for (int i = size1 + size2; i < size1 + size2 + size3; i++)
      {
        ids[i] = epvariables_[i - size1 - size2]->GetSegIds();
      }

      return ids;
    };

    /*!
    \brief Get flag indicating whether contact is active (true) or inactive (false)
    */
    bool GetContactFlag() override
    {
      // The element pair is assumed to be active when we have at least one active contact point
      return (cpvariables_.size() + gpvariables_.size() + epvariables_.size());
    };

    /*!
    \brief Get coordinates of contact point of element1
    */
    std::vector<CORE::LINALG::Matrix<3, 1>> GetX1() override
    {
      int size1 = cpvariables_.size();
      int size2 = gpvariables_.size();
      int size3 = epvariables_.size();
      std::vector<CORE::LINALG::Matrix<3, 1>> r1(
          size1 + size2 + size3, CORE::LINALG::Matrix<3, 1>(true));

      for (int i = 0; i < size1; i++)
      {
        TYPE eta1 = cpvariables_[i]->GetCP().first;
        for (int j = 0; j < 3; j++) r1[i](j) = CORE::FADUTILS::CastToDouble(r(eta1, element1_)(j));
      }

      for (int i = size1; i < size2 + size1; i++)
      {
        TYPE eta1 = gpvariables_[i - size1]->GetCP().first;
        for (int j = 0; j < 3; j++) r1[i](j) = CORE::FADUTILS::CastToDouble(r(eta1, element1_)(j));
      }

      for (int i = size1 + size2; i < size1 + size2 + size3; i++)
      {
        TYPE eta1 = epvariables_[i - size1 - size2]->GetCP().first;
        for (int j = 0; j < 3; j++) r1[i](j) = CORE::FADUTILS::CastToDouble(r(eta1, element1_)(j));
      }

      return r1;
    };

    /*!
    \brief Get coordinates of contact point of element2
    */
    std::vector<CORE::LINALG::Matrix<3, 1>> GetX2() override
    {
      int size1 = cpvariables_.size();
      int size2 = gpvariables_.size();
      int size3 = epvariables_.size();
      std::vector<CORE::LINALG::Matrix<3, 1>> r2(
          size1 + size2 + size3, CORE::LINALG::Matrix<3, 1>(true));

      for (int i = 0; i < size1; i++)
      {
        TYPE eta2 = cpvariables_[i]->GetCP().second;
        for (int j = 0; j < 3; j++) r2[i](j) = CORE::FADUTILS::CastToDouble(r(eta2, element2_)(j));
      }

      for (int i = size1; i < size2 + size1; i++)
      {
        TYPE eta2 = gpvariables_[i - size1]->GetCP().second;
        for (int j = 0; j < 3; j++) r2[i](j) = CORE::FADUTILS::CastToDouble(r(eta2, element2_)(j));
      }

      for (int i = size1 + size2; i < size1 + size2 + size3; i++)
      {
        TYPE eta2 = epvariables_[i - size1 - size2]->GetCP().second;
        for (int j = 0; j < 3; j++) r2[i](j) = CORE::FADUTILS::CastToDouble(r(eta2, element2_)(j));
      }

      return r2;
    };


    // TODO:
    /*!
    \brief Get normal vector
    */
    CORE::LINALG::SerialDenseVector GetNormal() override
    {
      CORE::LINALG::SerialDenseVector normal(3);

      for (int i = 0; i < 3; i++) normal(i) = 0.0;

      return normal;
    }

    /*!
    \brief Get flag indicating whether the nodal values of one element had been shifted due to r1=r2
           Since this is only possible for beam3contactnew elements but not for beam3contact
    elements we always return false within this class.
    */
    bool GetShiftStatus() override { return false; };

    /*!
      \Check, if there is a difference between the result of the new and old gap definition, i.e. if
      the beams centerlines have already crossed or not. Since this is only possible for
      beam3contactnew elements but not for beam3contact elements we always return false within this
      class.
    */
    bool GetNewGapStatus() override { return false; };
    //@}

    /*!
      \Get energy of penalty contact.
    */
    double GetEnergy() override
    {
      if (CORE::UTILS::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(
              bcparams_, "BEAMS_PENALTYLAW") != INPAR::BEAMCONTACT::pl_lp and
          CORE::UTILS::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(
              bcparams_, "BEAMS_PENALTYLAW") != INPAR::BEAMCONTACT::pl_qp and
          CORE::UTILS::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(
              bcparams_, "BEAMS_PENALTYLAW") != INPAR::BEAMCONTACT::pl_lpqp)
        FOUR_C_THROW("Contact Energy calculation not implemented for the chosen penalty law!");


      double energy = 0.0;

      for (int i = 0; i < (int)cpvariables_.size(); i++)
      {
        double ppfac = CORE::FADUTILS::CastToDouble(cpvariables_[i]->GetPPfac());
        double e = -cpvariables_[i]->GetIntegratedEnergy();
        energy += ppfac * e;
      }

      for (int i = 0; i < (int)gpvariables_.size(); i++)
      {
        double ppfac = CORE::FADUTILS::CastToDouble(gpvariables_[i]->GetPPfac());
        double e = -gpvariables_[i]->GetIntegratedEnergy();
        energy += ppfac * e;
      }

      for (int i = 0; i < (int)epvariables_.size(); i++)
      {
        double ppfac = CORE::FADUTILS::CastToDouble(epvariables_[i]->GetPPfac());
        double e = -epvariables_[i]->GetIntegratedEnergy();
        energy += ppfac * e;
      }

      return energy;
    };

    /*!
      \Get energy of perp penalty contact without transition factor contribution.
    */
    double GetUnscaledPerpEnergy() override
    {
      if (CORE::UTILS::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(
              bcparams_, "BEAMS_PENALTYLAW") != INPAR::BEAMCONTACT::pl_lp and
          CORE::UTILS::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(
              bcparams_, "BEAMS_PENALTYLAW") != INPAR::BEAMCONTACT::pl_qp and
          CORE::UTILS::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(
              bcparams_, "BEAMS_PENALTYLAW") != INPAR::BEAMCONTACT::pl_lpqp)
        FOUR_C_THROW("Contact Energy calculation not implemented for the chosen penalty law!");


      double energy = 0.0;

      for (int i = 0; i < (int)cpvariables_.size(); i++)
      {
        double e = -cpvariables_[i]->GetIntegratedEnergy();
        energy += e;
      }

      return energy;
    };

    /*!
      \Get energy of parallel penalty contact without transition factor contribution.
    */
    double GetUnscaledParallelEnergy() override
    {
      if (CORE::UTILS::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(
              bcparams_, "BEAMS_PENALTYLAW") != INPAR::BEAMCONTACT::pl_lp and
          CORE::UTILS::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(
              bcparams_, "BEAMS_PENALTYLAW") != INPAR::BEAMCONTACT::pl_qp and
          CORE::UTILS::IntegralValue<INPAR::BEAMCONTACT::PenaltyLaw>(
              bcparams_, "BEAMS_PENALTYLAW") != INPAR::BEAMCONTACT::pl_lpqp)
        FOUR_C_THROW("Contact Energy calculation not implemented for the chosen penalty law!");

      double energy = 0.0;

      for (int i = 0; i < (int)gpvariables_.size(); i++)
      {
        double e = -gpvariables_[i]->GetIntegratedEnergy();
        energy += e;
      }

      return energy;
    };

    // TODO
    /*!
      \We don't need this method for beam3contact elements!
    */
    CORE::LINALG::Matrix<3, 1, TYPE>* GetNormalOld() override { return nullptr; };

    // TODO
    /*!
      \We don't need this method for beam3contact elements!
    */
    bool FirstTimeStep() override { return false; };
    //@}

    //! @name Public evaluation methods
    /*!
    \brief Evaluate this contact element pair
    */
    bool Evaluate(CORE::LINALG::SparseMatrix& stiffmatrix, Epetra_Vector& fint, const double& pp,
        std::map<std::pair<int, int>, Teuchos::RCP<Beam3contactinterface>>& contactpairmap,
        Teuchos::ParameterList& timeintparams, bool fdcheck = false) override;

    /*!
    \brief Change the sign of the normal vector: This has to be done at the end of a time step when
    the remainig penetration is larger that the sum of the beam radii (R1+R2). Otherwise, the beams
    could cross in the next time step when the new gap function definition (ngf_=true) for slender
    beams is applied!
    */
    void InvertNormal() override { FOUR_C_THROW("Function not implemented!"); };

    // TODO
    /*!
      \brief We don't need this method for beam3contact elements!
    */
    void UpdateClassVariablesStep() override{};

    /*
    \brief Update nodal coordinates of both elements at the beginning of a new time step!
    */
    void UpdateElePos(CORE::LINALG::SerialDenseMatrix& newele1pos,
        CORE::LINALG::SerialDenseMatrix& newele2pos) override;

    /*
    \brief Update interpolated nodal tangents for tangent smoothing
    */
    void UpdateEleSmoothTangents(
        std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions) override;

    /** \brief print information about this beam contact element pair to screen
     *
     *  \author grill
     *  \date 05/16 */
    void Print() const override;

    //@}
   private:
    //! @name member variables

    //! reference to problem discretization
    const DRT::Discretization& pdiscret_;

    //! reference to beam contact discretization
    const DRT::Discretization& cdiscret_;

    //! dof offset between pdiscret and cdiscret
    const std::map<int, int>& dofoffsetmap_;

    //! first element of contact pair
    DRT::Element* element1_;

    //! second element of contact pair
    DRT::Element* element2_;

    //! beam contact parameter list
    Teuchos::ParameterList& bcparams_;

    //! current node coordinates of the two elements
    CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPE> ele1pos_;
    CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPE> ele2pos_;

    //! neighbor elements of element 1
    Teuchos::RCP<BEAMINTERACTION::B3CNeighbor> neighbors1_;

    //! neighbor elements of element 2
    Teuchos::RCP<BEAMINTERACTION::B3CNeighbor> neighbors2_;

    //! averaged nodal tangents, necessary for smoothed tangent fields of C^0 Reissner beams
    CORE::LINALG::Matrix<3 * numnodes, 1> nodaltangentssmooth1_;
    CORE::LINALG::Matrix<3 * numnodes, 1> nodaltangentssmooth2_;

    //! current Newton iteration
    int iter_;

    //! current time step
    int numstep_;

    //! cross section radius of first beam
    const double r1_;

    //! cross section radius of second beam
    const double r2_;

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

    //! Indicates if the left / right node of the slave element 1 coincides with the endpoint of the
    //! physical beam (true) or not (false)
    std::pair<bool, bool> boundarynode1_;

    //! Indicates if the left / right node of the master element 2 coincides with the endpoint of
    //! the physical beam (true) or not (false)
    std::pair<bool, bool> boundarynode2_;

    //! Variables stored at the closest points of the large-angle-contact algorithm
    std::vector<Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues>>> cpvariables_;

    //! Variables stored at the Gauss points of the small-angle-contact algorithm
    std::vector<Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues>>> gpvariables_;

    //! Variables stored at the end points of the endpoint-contact algorithm
    std::vector<Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues>>> epvariables_;

    //@}

    //! @name Private evaluation methods

    /*!
    \brief Get active large angle pairs
    */
    void GetActiveLargeAnglePairs(std::vector<CORE::LINALG::Matrix<3, 1, double>>& endpoints1,
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& endpoints2,
        std::map<std::pair<int, int>, CORE::LINALG::Matrix<3, 1, double>>& closelargeanglesegments,
        const double pp);

    /*!
    \brief Evaluate active large angle pairs
    */
    void EvaluateActiveLargeAnglePairs(
        CORE::LINALG::SparseMatrix& stiffmatrix, Epetra_Vector& fint);

    /*!
    \brief Get active small angle pairs
    */
    void GetActiveSmallAnglePairs(
        std::map<std::pair<int, int>, CORE::LINALG::Matrix<3, 1, double>>& closesmallanglesegments,
        std::pair<int, int>* iminmax = nullptr,
        std::pair<bool, bool>* leftrightsolutionwithinsegment = nullptr,
        std::pair<double, double>* eta1_leftrightboundary = nullptr);

    /*!
    \brief Evaluate active small angle pairs
    */
    void EvaluateActiveSmallAnglePairs(CORE::LINALG::SparseMatrix& stiffmatrix, Epetra_Vector& fint,
        std::pair<int, int>* iminmax = nullptr,
        std::pair<bool, bool>* leftrightsolutionwithinsegment = nullptr,
        std::pair<double, double>* eta1_leftrightboundary = nullptr);


    /*!
    \brief Get active endpoint pairs
    */
    void GetActiveEndPointPairs(
        std::vector<std::pair<int, int>>& closeendpointsegments, const double pp);

    /*!
    \brief Evaluate active endpoint pairs
    */
    void EvaluateActiveEndPointPairs(CORE::LINALG::SparseMatrix& stiffmatrix, Epetra_Vector& fint);

    /*!
    \brief Find segments close to each other
    */
    void GetCloseSegments(const std::vector<CORE::LINALG::Matrix<3, 1, double>>& endpoints1,
        const std::vector<CORE::LINALG::Matrix<3, 1, double>>& endpoints2,
        std::map<std::pair<int, int>, CORE::LINALG::Matrix<3, 1, double>>& closesmallanglesegments,
        std::map<std::pair<int, int>, CORE::LINALG::Matrix<3, 1, double>>& closelargeanglesegments,
        std::vector<std::pair<int, int>>& closeendpointsegments, double maxactivedist);

    /*!
    \brief Find contact point via closest point projection
    */
    bool ClosestPointProjection(double& eta_left1, double& eta_left2, double& l1, double& l2,
        CORE::LINALG::Matrix<3, 1, double>& segmentdata, std::pair<TYPE, TYPE>& solutionpoints,
        int segid1, int segid2);

    /*!
    \brief Find closest point eta2_master on a line for a given slave point eta1_slave
    */
    bool PointToLineProjection(double& eta1_slave, double& eta_left2, double& l2,
        double& eta2_master, double& gap, double& alpha, bool& pairactive, bool smallanglepair,
        bool invertpairs = false, bool orthogonalprojection = false);

    /*!
    \brief Determine minimal distance and contact angle for unconverged segment pair
    */
    void CheckUnconvergedSegmentPair(double& eta_left1, double& eta_left2, double& l1, double& l2,
        double& eta1_min, double& eta2_min, double& g_min, double& alpha_g_min,
        bool& pointtolinesolfound);

    /*!
    \brief Subdivide elements into segments for CPP
    */
    double CreateSegments(DRT::Element* ele,
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& endpoints_final, int& numsegment, int i);

    /*!
    \brief Get maximal gap at which a contact can become active
    */
    double GetMaxActiveDist();

    /*!
    \brief Check, if segments are fine enough
    */
    bool CheckSegment(CORE::LINALG::Matrix<3, 1, double>& r1,
        CORE::LINALG::Matrix<3, 1, double>& t1, CORE::LINALG::Matrix<3, 1, double>& r2,
        CORE::LINALG::Matrix<3, 1, double>& t2, CORE::LINALG::Matrix<3, 1, double>& rm,
        double& segdist);

    /*!
    \brief Calculate scalar contact force
    */
    void CalcPenaltyLaw(Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues>> variables);

    /*!
    \brief Calculate angle-dependent penalty scale factor for small-angle-contact
    */
    void CalcPerpPenaltyScaleFac(
        Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues>> cpvariables,
        CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi, CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const double shiftangle1, const double shiftangle2);

    /*!
    \brief Calculate angle-dependent penalty scale factor for large-angle-contact
    */
    void CalcParPenaltyScaleFac(
        Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues>> gpvariables,
        CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi, CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const double shiftangle1, const double shiftangle2);

    /*!
     \brief Compute contact forces
     */
    void EvaluateFcContact(Epetra_Vector* fint, const CORE::LINALG::Matrix<3, 1, TYPE>& r1,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2, const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi,
        Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues>> variables,
        const double& intfac, bool cpp, bool gp, bool fixedendpointxi, bool fixedendpointeta,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPE>* fc1_FAD = nullptr,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPE>* fc2_FAD = nullptr);

    /*!
    \brief Evaluate contact stiffness
    */
    void EvaluateStiffcContact(CORE::LINALG::SparseMatrix& stiffmatrix,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1, const CORE::LINALG::Matrix<3, 1, TYPE>& r2,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xixi,
        Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues>> variables,
        const double& intfac, bool cpp, bool gp, bool fixedendpointxi, bool fixedendpointeta);

    /*!
    \brief FAD-based Evaluation of contact stiffness in case of ENDPOINTSEGMENTATION
    */
    void EvaluateStiffcContactIntSeg(CORE::LINALG::SparseMatrix& stiffmatrix,
        const CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi_bound,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1, const CORE::LINALG::Matrix<3, 1, TYPE>& r2,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi,
        Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues>> cpvariables,
        const double& intfac, const double& d_xi_ele_d_xi_bound, TYPE signed_jacobi_interval);

    /*!
    \brief Linearizations of contact point
    */
    void ComputeLinXiAndLinEta(
        CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi,
        CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_eta,
        const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi);

    /*!
    \brief Lin. of contact point coordinate eta with fixed xi
    */
    void ComputeLinEtaFixXi(
        CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_eta,
        const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi);

    /*!
    \brief Lin. of contact point coordinate xi with fixed eta
    */
    void ComputeLinXiFixEta(
        CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi);

    /*!
    \brief Compute linearization of integration interval bounds (necessary in case of
    ENDPOINTSEGMENTATION)
    */
    void ComputeLinXiBound(
        CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi_bound,
        TYPE& eta1_bound, TYPE eta2);

    /*!
    \brief Compute linearization of gap
    */
    void ComputeLinGap(CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_gap,
        const CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi,
        const CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_eta,
        const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r, const TYPE& norm_delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2);

    /*!
    \brief Compute linearization of cosine of contact angle
    */
    void ComputeLinCosContactAngle(
        CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_coscontactangle,
        CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi,
        CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_eta,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi);

    /*!
    \brief Compute linearization of normal vector
    */
    void ComputeLinNormal(
        CORE::LINALG::Matrix<3, 2 * 3 * numnodes * numnodalvalues, TYPE>& delta_normal,
        const CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_xi,
        const CORE::LINALG::Matrix<2 * 3 * numnodes * numnodalvalues, 1, TYPE>& delta_eta,
        const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2);

    /*!
    \brief Calculate shape function values for given parameter values
    */
    void GetShapeFunctions(CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xixi,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xixi, const TYPE& eta1,
        const TYPE& eta2);

    /*!
    \brief Calculate one specified shape function value / derivative for given parameter value and
    element
    */
    void GetShapeFunctions(CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N,
        const TYPE& eta, int deriv, DRT::Element* ele);

    /*!
    \brief Assemble the shape functions into corresponding matrices
    */
    void AssembleShapefunctions(const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N_i,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N_i_xi,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N_i_xixi,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N_xi,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N_xixi);

    /*!
    \brief Assemble shape functions for one given matrix
    */
    void AssembleShapefunctions(const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N_i,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N);

    /*!
    \brief compute coordinate at given curve point
    */
    CORE::LINALG::Matrix<3, 1, TYPE> r(const TYPE& eta, DRT::Element* ele);

    /*!
    \brief compute derivative at given curve point
    */
    CORE::LINALG::Matrix<3, 1, TYPE> r_xi(const TYPE& eta, DRT::Element* ele);

    /*!
    \brief Compute coordinates and their derivatives from the discretization
    */
    void ComputeCoordsAndDerivs(CORE::LINALG::Matrix<3, 1, TYPE>& r1,
        CORE::LINALG::Matrix<3, 1, TYPE>& r2, CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi, CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xixi);

    /*!
    \brief Compute coordinates of contact points of last time step from the discretization
    */
    void ComputeOldCoordsAndDerivs(CORE::LINALG::Matrix<3, 1, TYPE>& r1_old,
        CORE::LINALG::Matrix<3, 1, TYPE>& r2_old, CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi_old,
        CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi_old,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi);

    /*!
    \brief Utility method for CPP (evaluate nonlinear function f)
    */
    void EvaluateOrthogonalityCondition(CORE::LINALG::Matrix<2, 1, TYPE>& f,
        const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r, const double norm_delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi, const CORE::LINALG::Matrix<3, 1, TYPE>& t1,
        const CORE::LINALG::Matrix<3, 1, TYPE>& t2);

    /*!
    \brief Utility method for CPP (evaluate Jacobian of nonlinear function f)
    */
    void EvaluateLinOrthogonalityCondition(CORE::LINALG::Matrix<2, 2, TYPE>& df,
        CORE::LINALG::Matrix<2, 2, TYPE>& dfinv, const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r,
        const double norm_delta_r, const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi, const CORE::LINALG::Matrix<3, 1, TYPE>& t1,
        const CORE::LINALG::Matrix<3, 1, TYPE>& t2, const CORE::LINALG::Matrix<3, 1, TYPE>& t1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& t2_xi, bool& elementscolinear);

    /*!
    \brief Evaluate orthogonality cond. of point to line projeciton
    */
    void EvaluatePTLOrthogonalityCondition(TYPE& f, const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r,
        const double norm_delta_r, const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi, bool orthogonalprojection);

    /*!
    \brief Evaluate Jacobian df of PTLOrthogonalityCondition
    */
    bool EvaluateLinPTLOrthogonalityCondition(TYPE& df,
        const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r, const double norm_delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi, bool orthogonalprojection);

    /*!
    \brief Compute normal vector and gap function at contact point
    */
    void ComputeNormal(CORE::LINALG::Matrix<3, 1, TYPE>& r1, CORE::LINALG::Matrix<3, 1, TYPE>& r2,
        CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi, CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        Teuchos::RCP<Beam3contactvariables<numnodes, numnodalvalues>> variables, int contacttype);

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

      if (ele == nullptr) FOUR_C_THROW("Dynamic cast to Beam3Base failed");

      return ele->GetJacobiFacAtXi(xi);
    }

    /*!
      \brief Set class variables at the beginning of a Newton step
    */
    void SetClassVariables(Teuchos::ParameterList& timeintparams);

    /*!
      \brief Linearization-check of coordinates xi and eta via FAD
    */
    void FADCheckLinXiAndLinEta(const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xixi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xixi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N1_xi,
        const CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues, TYPE>& N2_xi);

    /*!
      \brief Linearization-check for local Newton in CPP via FAD
    */
    void FADCheckLinOrthogonalityCondition(const CORE::LINALG::Matrix<3, 1, TYPE>& delta_r,
        const double& norm_delta_r, const CORE::LINALG::Matrix<3, 1, TYPE>& r1_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& r2_xi, const CORE::LINALG::Matrix<3, 1, TYPE>& t1,
        const CORE::LINALG::Matrix<3, 1, TYPE>& t2);

    /*!
      \brief FD-Check of stiffness matrix
    */
    void FDCheck(CORE::LINALG::SparseMatrix& stiffmatrix, Epetra_Vector& fint, const double& pp,
        std::map<std::pair<int, int>, Teuchos::RCP<Beam3contactinterface>>& contactpairmap,
        Teuchos::ParameterList& timeintparams, bool fdcheck);

    //@}

  };  // class Beam3contact
}  // namespace CONTACT

FOUR_C_NAMESPACE_CLOSE

#endif
