/*----------------------------------------------------------------------------*/
/*!
\file beam3contactvariables.H

\brief One beam contact segment living on an element pair

\level 2

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/
#ifndef BEAM3CV_H
#define BEAM3CV_H

#include "beam3contact_utils.H"
#include "beam3contact_tangentsmoothing.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_sparsematrix.H"
#include "beam3contact_defines.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "Sacado.hpp"
#include "beam3contactinterface.H"

#include "../drt_lib/standardtypes_cpp.H"

namespace CONTACT
{

/*!
 \brief contact element for contact between two 3D beam elements

 Refer also to the Semesterarbeit of Matthias Mayr, 2010

 \author meier (meier@lnm.mw.tum.de)
 */

template<const int numnodes , const int numnodalvalues>
class Beam3contactvariables
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
  Beam3contactvariables(std::pair<TYPE,TYPE>& closestpoint,
                        std::pair<int,int>& segids,
                        std::pair<int,int>& intids,
                        const double& pp,
                        TYPE jacobi);

  /*!
  \brief Destructor
  */
  virtual ~Beam3contactvariables() {};
  //@}


  //! @name Access methods

  /*!
  \brief Set closest point
  */

  void SetCP(std::pair<TYPE,TYPE> cp) { closestpoint_=cp; };

  /*!
  \brief Get closest point
  */
  std::pair<TYPE,TYPE> GetCP() { return closestpoint_; };

  /*!
  \brief Get Segment Ids
  */
  std::pair<int,int> GetSegIds() { return segids_; };

  /*!
  \brief Get Integration Ids
  */
  std::pair<int,int> GetIntIds() { return intids_; };

  /*!
  \brief Get jacobi factor
  */
  TYPE GetJacobi() { return jacobi_; };

  /*!
  \brief Set gap
  */
  void SetGap(TYPE gap) { gap_=gap; };

  /*!
  \brief Get gap
  */
  TYPE GetGap() { return gap_; };

  /*!
  \brief Set gap
  */
  void SetNormal(LINALG::TMatrix<TYPE,3,1> normal) { normal_=normal; };

  /*!
  \brief Get gap
  */
  LINALG::TMatrix<TYPE,3,1> GetNormal() { return normal_; };

  /*!
  \brief Get penalty parameter
  */
  double GetPP() { return pp_; };

  /*!
  \brief Set penalty force
  */
  void Setfp(TYPE fp) { fp_=fp; };

  /*!
  \brief Get penalty force
  */
  TYPE Getfp() { return fp_; };

  /*!
  \brief Set derivative of penalty force
  */
  void Setdfp(TYPE dfp) { dfp_=dfp; };

  /*!
  \brief Get pre-factor for penalty parameter
  */
  TYPE GetPPfac() { return ppfac_; };

  /*!
  \brief Set pre-factor for penalty parameter
  */
  void SetPPfac(TYPE ppfac) { ppfac_=ppfac; };

  /*!
  \brief Get linearization of pre-factor for penalty parameter
  */
  TYPE GetDPPfac() { return dppfac_; };

  /*!
  \brief Set linearization of pre-factor for penalty parameter
  */
  void SetDPPfac(TYPE dppfac) { dppfac_=dppfac; };

  /*!
  \brief Get derivative of penalty force
  */
  TYPE Getdfp() { return dfp_; };

  /*!
  \brief Set penalty energy
  */
  void SetEnergy(TYPE e) { energy_=e; };

  /*!
  \brief Get penalty energy
  */
  TYPE GetEnergy() { return energy_; };

  /*!
  \brief Set length integrated penalty energy
  */
  void SetIntegratedEnergy(double inte) { integratedenergy_=inte; };

  /*!
  \brief Get length integrated penalty energy
  */
  double GetIntegratedEnergy() { return integratedenergy_; };

  /*!
  \brief Set contact angle
  */
  void SetAngle(double angle) { angle_=angle; };

  /*!
  \brief Get contact angle
  */
  double GetAngle() { return angle_; };

  //@}


  //@}
private:

  //closest point coordinates
  std::pair<TYPE,TYPE> closestpoint_;

  //element local Ids of considered segments
  std::pair<int,int> segids_;

  //stores numgp and number of integration interval (only necessary for small-angle contact)
  std::pair<int,int> intids_;

  //jacobi factor for integration (only necessary for line contact)
  TYPE jacobi_;

  //gap function
  TYPE gap_;

  //normal vector
  LINALG::TMatrix<TYPE,3,1> normal_;

  //penalty parameter
  double pp_;

  //penalty parameter
  TYPE ppfac_;

  //linearization of penalty parameter
  TYPE dppfac_;

  //penalty force
  TYPE fp_;

  //derivative of penalty force with respect to gap: d(fp_)/d(gap_)
  TYPE dfp_;

  //penalty energy of collocation point / Gauss point
  TYPE energy_;

  //length integrated penalty energy of collocation point (in this case identical to energy_) / Gauss point
  double integratedenergy_;

  //contact angle
  double angle_;
  //@}

}; // class Beam3contactvariables
} // namespace CONTACT

#endif  // #ifndef BEAM3CV_H
