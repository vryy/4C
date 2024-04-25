/*----------------------------------------------------------------------*/
/*! \file

\brief interface class for templated classes beam3contact and beam3contactnew

\level 2

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_BEAMCONTACT_BEAM3CONTACTINTERFACE_HPP
#define FOUR_C_BEAMCONTACT_BEAM3CONTACTINTERFACE_HPP

#include "4C_config.hpp"

#include "4C_beaminteraction_beam_to_beam_contact_defines.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_sparsematrix.hpp"

FOUR_C_NAMESPACE_OPEN


namespace CONTACT
{
  /*!
   \brief interface class for templated classes beam3contact and beam3contactnew

   */

  class Beam3contactinterface
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
    Beam3contactinterface() {}

    /*!
    \brief Destructor
    */
    virtual ~Beam3contactinterface() = default;
    //@}

    //! @name Access methods
    /*!
    \brief Get problem discretization
    */
    virtual const DRT::Discretization& ProblemDiscret() const = 0;

    /*!
    \brief Get beam contact discretization
    */
    virtual const DRT::Discretization& ContactDiscret() const = 0;

    /*!
    \brief Get offset of dofs between cdiscret and pdiscret
    */
    virtual const std::map<int, int>& DofOffset() const = 0;

    /*!
    \brief Get first element
    */
    virtual const DRT::Element* Element1() = 0;
    // inline const DRT::Element* Element1() { return element1_;};

    /*!
    \brief Get first element
    */
    virtual const DRT::Element* Element2() = 0;

    /*!
    \brief Get number of contact points on this element pair
    */
    virtual int GetNumCps() = 0;
    virtual int GetNumGps() = 0;
    virtual int GetNumEps() = 0;

    /*!
    \brief Get vector of type declarations (0=closest point contact, 1=gauss point contact, 2= end
    point contact) of all contact pairs
    */
    virtual std::vector<int> GetContactType() = 0;

    /*!
    \brief Get vector of all gaps of this contact pair
    */
    virtual std::vector<double> GetGap() = 0;

    /*!
    \brief Get vector of all contact forces of this contact pair
    */
    virtual std::vector<double> GetContactForce() = 0;

    /*!
    \brief Get vector of all contact angles of this contact pair
    */
    virtual std::vector<double> GetContactAngle() = 0;

    /*!
    \brief Get vector of all closest points of this contact pair
    */
    virtual std::vector<std::pair<double, double>> GetClosestPoint() = 0;

    /*!
    \brief Return number of individual contact segments on element pair
    */
    virtual std::pair<int, int> GetNumSegments() = 0;

    /*!
    \brief Return ids of active segments
    */
    virtual std::vector<std::pair<int, int>> GetSegmentIds() = 0;

    /*!
    \brief Get flag ndicating whether contact is active (true) or inactive (false)
    */
    virtual bool GetContactFlag() = 0;

    /*!
    \brief Get coordinates of contact point of element1 and element2
    */
    virtual std::vector<CORE::LINALG::Matrix<3, 1>> GetX1() = 0;

    virtual std::vector<CORE::LINALG::Matrix<3, 1>> GetX2() = 0;

    virtual CORE::LINALG::SerialDenseVector GetNormal() = 0;

    virtual CORE::LINALG::Matrix<3, 1, TYPE>* GetNormalOld() = 0;

    /*!
      \Check, if there is a difference between the result of the new and old gap definition, i.e. if
      the beams centerlines have already crossed or not.
    */
    virtual bool GetNewGapStatus() = 0;

    /*!
      \Get energy of penalty contact.
    */
    virtual double GetEnergy() = 0;

    /*!
      \Get energy of perp penalty contact without transition factor contribution.
    */
    virtual double GetUnscaledPerpEnergy() = 0;

    /*!
      \Get energy of parallel penalty contact without transition factor contribution.
    */
    virtual double GetUnscaledParallelEnergy() = 0;

    virtual bool FirstTimeStep() = 0;

    /*!
    \brief Get flag indicating whether the nodal values of one element had been shifted due to r1=r2
    */
    virtual bool GetShiftStatus() = 0;
    //@}

    /** \brief print this beam contact element pair to screen
     *
     *  \author grill
     *  \date 05/16 */
    virtual void Print() const = 0;


    //! @name Public evaluation methods
    /*!
    \brief Evaluate this contact element pair
    */
    virtual bool Evaluate(CORE::LINALG::SparseMatrix& stiffmatrix, Epetra_Vector& fint,
        const double& pp,
        std::map<std::pair<int, int>, Teuchos::RCP<Beam3contactinterface>>& contactpairmap,
        Teuchos::ParameterList& timeintparams, bool fdcheck = false) = 0;

    //! return appropriate internal implementation class (acts as a simple factory)
    static Teuchos::RCP<Beam3contactinterface> Impl(const int numnodes, const int numnodalvalues,
        const DRT::Discretization& pdiscret, const DRT::Discretization& cdiscret,
        const std::map<int, int>& dofoffsetmap, DRT::Element* element1, DRT::Element* element2,
        Teuchos::ParameterList& beamcontactparams);

    /*!
    \brief Change the sign of the normal vector: This has to be done at the end of a time step when
    the remainig penetration is larger that the sum of the beam radii (R1+R2). Otherwise, the beams
    could cross in the next time step when the new gap function definition (ngf_=true) for slender
    beams is applied!
    */
    virtual void InvertNormal() = 0;

    /*!
      \brief Update of class variables at the end of a time step
    */
    virtual void UpdateClassVariablesStep() = 0;

    /*
    \brief Update nodal coordinates of both elements at the beginning of a new time step!
    */
    virtual void UpdateElePos(CORE::LINALG::SerialDenseMatrix& newele1pos,
        CORE::LINALG::SerialDenseMatrix& newele2pos) = 0;

    /*
    \brief Update interpolated nodal tangents for tangent smoothing
    */
    virtual void UpdateEleSmoothTangents(
        std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions) = 0;

  };  // class Beam3contactinterface
}  // namespace CONTACT

FOUR_C_NAMESPACE_CLOSE

#endif
