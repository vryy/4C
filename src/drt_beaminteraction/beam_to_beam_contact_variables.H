/*----------------------------------------------------------------------------*/
/*!

\brief one beam contact segment living on an element pair

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/

#ifndef BEAMTOBEAMCONTACTVARIABLES_H
#define BEAMTOBEAMCONTACTVARIABLES_H


#include "beam3contact_utils.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../headers/FAD_utils.H"
#include "Sacado.hpp"

#include "../drt_beaminteraction/beam_contact_pair.H"

#include "../drt_lib/standardtypes_cpp.H"
#include "beam3contact_defines.H"
#include "beam3contact_tangentsmoothing.H"

namespace BEAMINTERACTION
{
  /*!
   \brief Todo
   */

  template <unsigned int numnodes, unsigned int numnodalvalues>
  class BeamToBeamContactVariables
  {
   public:
    //! @name Friends
    // no friend classes defined
    //@}

    //! @name Constructors and destructors and related methods
    /*!
    \brief Standard Constructor
    */
    BeamToBeamContactVariables(std::pair<TYPE, TYPE>& closestpoint, std::pair<int, int>& segids,
        std::pair<int, int>& intids, const double& pp, TYPE jacobi);

    /*!
    \brief Destructor
    */
    virtual ~BeamToBeamContactVariables(){};
    //@}


    //! @name Access methods

    /*!
    \brief Set closest point
    */

    void SetCP(std::pair<TYPE, TYPE> cp) { closestpoint_ = cp; };

    /*!
    \brief Get closest point
    */
    std::pair<TYPE, TYPE> GetCP() const { return closestpoint_; };

    /*!
    \brief Get Segment Ids
    */
    std::pair<int, int> GetSegIds() const { return segids_; };

    /*!
    \brief Get Integration Ids
    */
    std::pair<int, int> GetIntIds() const { return intids_; };

    /*!
    \brief Get jacobi factor
    */
    TYPE GetJacobi() const { return jacobi_; };

    /*!
    \brief Set gap
    */
    void SetGap(TYPE gap) { gap_ = gap; };

    /*!
    \brief Get gap
    */
    TYPE GetGap() const { return gap_; };

    /*!
    \brief Set gap
    */
    void SetNormal(LINALG::Matrix<3, 1, TYPE> normal) { normal_ = normal; };

    /*!
    \brief Get gap
    */
    LINALG::Matrix<3, 1, TYPE> GetNormal() const { return normal_; };

    /*!
    \brief Get penalty parameter
    */
    double GetPP() const { return pp_; };

    /*!
    \brief Set penalty force
    */
    void Setfp(TYPE fp) { fp_ = fp; };

    /*!
    \brief Get penalty force
    */
    TYPE Getfp() const { return fp_; };

    /*!
    \brief Set derivative of penalty force
    */
    void Setdfp(TYPE dfp) { dfp_ = dfp; };

    /*!
    \brief Get pre-factor for penalty parameter
    */
    TYPE GetPPfac() const { return ppfac_; };

    /*!
    \brief Set pre-factor for penalty parameter
    */
    void SetPPfac(TYPE ppfac) { ppfac_ = ppfac; };

    /*!
    \brief Get linearization of pre-factor for penalty parameter
    */
    TYPE GetDPPfac() const { return dppfac_; };

    /*!
    \brief Set linearization of pre-factor for penalty parameter
    */
    void SetDPPfac(TYPE dppfac) { dppfac_ = dppfac; };

    /*!
    \brief Get derivative of penalty force
    */
    TYPE Getdfp() const { return dfp_; };

    /*!
    \brief Set penalty energy
    */
    void SetEnergy(TYPE e) { energy_ = e; };

    /*!
    \brief Get penalty energy
    */
    TYPE GetEnergy() const { return energy_; };

    /*!
    \brief Set length integrated penalty energy
    */
    void SetIntegratedEnergy(double inte) { integratedenergy_ = inte; };

    /*!
    \brief Get length integrated penalty energy
    */
    double GetIntegratedEnergy() const { return integratedenergy_; };

    /*!
    \brief Set contact angle
    */
    void SetAngle(double angle) { angle_ = angle; };

    /*!
    \brief Get contact angle
    */
    double GetAngle() const { return angle_; };

    //@}


    //! @name Print methods

    /** \brief print information to screen
     *
     *  \author grill
     *  \date 12/16 */
    inline void Print(std::ostream& out) const
    {
      out << "\nInstance of BeamToBeamContactVariables (SegmentIds " << segids_.first << " & "
          << segids_.second << "):";
      out << "\ngap= " << FADUTILS::CastToDouble(gap_);
      out << "\nangle= " << angle_ / M_PI * 180.0 << " degree";
      out << "\nclosest point coords: " << FADUTILS::CastToDouble(closestpoint_.first) << " "
          << FADUTILS::CastToDouble(closestpoint_.second);

      out << "\n";
    }

    /** \brief print information to screen
     *
     *  \author grill
     *  \date 12/16 */
    inline void PrintSummaryInOneLine(std::ostream& out) const
    {
      out << std::setw(9) << std::left << std::setprecision(2) << closestpoint_.first
          << std::setw(9) << std::left << std::setprecision(2) << closestpoint_.second
          << std::setw(9) << std::left << std::setprecision(3) << angle_ / M_PI * 180.0
          << std::setw(12) << std::left << std::scientific << gap_ << std::setw(12) << std::left
          << std::scientific << fp_ << std::setprecision(6)
          << std::resetiosflags(std::ios::scientific) << std::right;
    }

    //@}
   private:
    // closest point coordinates
    std::pair<TYPE, TYPE> closestpoint_;

    // element local Ids of considered segments
    std::pair<int, int> segids_;

    // stores numgp and number of integration interval (only necessary for small-angle contact)
    std::pair<int, int> intids_;

    // jacobi factor for integration (only necessary for line contact)
    TYPE jacobi_;

    // gap function
    TYPE gap_;

    // normal vector
    LINALG::Matrix<3, 1, TYPE> normal_;

    // penalty parameter
    double pp_;

    // penalty parameter
    TYPE ppfac_;

    // linearization of penalty parameter
    TYPE dppfac_;

    // penalty force
    TYPE fp_;

    // derivative of penalty force with respect to gap: d(fp_)/d(gap_)
    TYPE dfp_;

    // penalty energy of collocation point / Gauss point
    TYPE energy_;

    // length integrated penalty energy of collocation point (in this case identical to energy_) /
    // Gauss point
    double integratedenergy_;

    // contact angle
    double angle_;
    //@}
  };
}  // namespace BEAMINTERACTION

#endif
