/*----------------------------------------------------------------------------*/
/*!
\file beam3_base.H

\brief base class for all beam elements

\level 2

\maintainer Christoph Meier
*/
/*----------------------------------------------------------------------------*/

#ifndef BEAM3_BASE_H_
#define BEAM3_BASE_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"

#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_lib/drt_element.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_integration.H"

/* number of linear subsegments used for visualization of beam centerline
 * in vtu format (see post_drt_vtu_writer.cpp) */
#define BEAMSVTUVISUALSUBSEGMENTS 10

// forward declaration ...
namespace STR {
namespace ELEMENTS {
  class ParamsInterface;
} // namespace ELEMENTS
} // namespace STR

namespace STATMECH {
  class ParamsInterface;
}

namespace DRT {
namespace ELEMENTS {

//! base class for all beam elements
class Beam3Base : public DRT::Element
{

public:

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param etype (in): Type of element
  \param owner (in): owner processor of the element
  */
  Beam3Base(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element
  */
  Beam3Base(const Beam3Base& old);


  //! Destructor
  virtual ~Beam3Base(){};

  /** \brief set the parameter interface ptr for the solid elements
   *
   *  \param p (in): Parameter list coming from the time integrator.
   *
   *  \author hiermeier
   *  \date 04/16 */
  virtual void SetParamsInterfacePtr(const Teuchos::ParameterList& p);

  virtual void SetStatMechParamsInterfacePtr();

  /** \brief returns true if the parameter interface is defined and initialized, otherwise false
   *
   *  \author hiermeier
   *  \date 04/16 */
  virtual inline bool IsParamsInterface() const
  { return (not interface_ptr_.is_null()); }

  /** \brief get access to the parameter interface pointer
   *
   *  \author hiermeier
   *  \date 04/16 */
  virtual Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> ParamsInterfacePtr();
  virtual Teuchos::RCP<STATMECH::ParamsInterface> StatMechParamsInterfacePtr() const;

  //!computes the number of different random numbers required in each time step for generation of stochastic forces
  virtual int HowManyRandomNumbersINeed() {return 0;};

  /** \brief get access to the element reference length
   *        (i.e. arc-length in stress-free configuration)
   *
   *  \author grill
   *  \date 05/16 */
  virtual double RefLength() const = 0;

  /** \brief get access to the first principal area moment of inertia
   *
   *  \author grill
   *  \date 08/16 */
  virtual const double& Iyy() const = 0;

  /** \brief get number of nodes used for centerline interpolation
   *
   *  \author grill
   *  \date 05/16 */
  virtual int NumCenterlineNodes() const = 0;

  /** \brief return flag whether Hermite polynomials are applied for centerline interpolation
   *
   *  \author grill
   *  \date 05/16 */
  virtual bool HermiteCenterlineInterpolation() const = 0;

  /** \brief return GIDs of all additive DoFs for a given node
   *
   *  \author grill
   *  \date 07/16 */
  std::vector<int> GetAdditiveDofGIDs(const DRT::Discretization& discret, const DRT::Node& node) const;

  /** \brief return GIDs of all non-additive, i.e. rotation pseudo vector DoFs for a given node
   *
   *  \author grill
   *  \date 07/16 */
  std::vector<int> GetRotVecDofGIDs(const DRT::Discretization& discret, const DRT::Node& node) const;

  /** \brief add indices of those DOFs of a given node that are positions
   *
   *  \author grill
   *  \date 07/16 */
  virtual void PositionDofIndices(std::vector<int>& posdofs, const DRT::Node& node) const = 0;

  /** \brief add indices of those DOFs of a given node that are tangents (in the case of Hermite interpolation)
   *
   *  \author grill
   *  \date 07/16 */
  virtual void TangentDofIndices(std::vector<int>& tangdofs, const DRT::Node& node) const = 0;

  /** \brief add indices of those DOFs of a given node that are rotation DOFs (non-additive rotation vectors)
   *
   *  \author grill
   *  \date 07/16 */
  virtual void RotationVecDofIndices(std::vector<int>& rotvecdofs, const DRT::Node& node) const = 0;

  /** \brief add indices of those DOFs of a given node that are 1D rotation DOFs
   *         (planar rotations are additive, e.g. in case of relative twist DOF of beam3k with rotvec=false)
   *
   *  \author grill
   *  \date 07/16 */
  virtual void Rotation1DDofIndices(std::vector<int>& twistdofs, const DRT::Node& node) const = 0;

  /** \brief add indices of those DOFs of a given node that represent norm of tangent vector
   *         (additive, e.g. in case of beam3k with rotvec=true)
   *
   *  \author grill
   *  \date 07/16 */
  virtual void TangentLengthDofIndices(std::vector<int>& tangnormdofs, const DRT::Node& node) const = 0;

  /** \brief get Jacobi factor ds/dxi(xi) at xi \in [-1;1]
   *
   *  \author grill
   *  \date 06/16 */
  virtual double GetJacobiFacAtXi(const double& xi) const = 0;

  /** \brief Get number of degrees of freedom of a single node
   *
   *  \author eichinger
   *  \date 08/16 */
  virtual int NumDofPerNode(const DRT::Node& node) const
  {
    dserror("not implemented");
    return -1;
  }

  /** \brief get centerline position at xi \in [-1,1] (element parameter space) in stress-free reference configuration
   *
   *  \author grill
   *  \date 06/16 */
  void GetRefPosAtXi(LINALG::Matrix<3,1>& refpos, const double& xi) const;

  /** \brief get unit tangent vector in reference configuration at i-th node of beam element (element-internal numbering)
   *
   *  \author grill
   *  \date 06/16 */
  virtual void GetRefTangentAtNode(LINALG::Matrix<3,1>& Tref_i, const int& i) const
  {
    // ToDo make pure virtual and add/generalize implementations in beam eles
    dserror("not implemented");
  }

  /** \brief get centerline position at xi \in [-1,1] (element parameter space)
   *
   *  \author grill
   *  \date 06/16 */
  virtual void GetPosAtXi(LINALG::Matrix<3,1>& pos,
                          const double& xi,
                          const std::vector<double>& disp_totlag) const
  {
    // ToDo make pure virtual and add/generalize implementations in beam eles
    dserror("not implemented");
  }

  /** \brief get triad at xi \in [-1,1] (element parameter space)
   *
   *  \author grill
   *  \date 07/16 */
  virtual void GetTriadAtXi(LINALG::Matrix<3,3>&       triad,
                            const double&              xi,
                            const std::vector<double>& psi_totlag) const
  {
    // ToDo make pure virtual and add/generalize implementations in beam eles
    // so far only implemented for beam3r
    dserror("not implemented");
  }

protected:
  /** \brief get access to the interface
   *
   *  \author hiermeier
   *  \date 04/16 */
  inline STR::ELEMENTS::ParamsInterface& ParamsInterface()
  { if (not IsParamsInterface())
      dserror("The interface ptr is not set!");
    return *interface_ptr_;
  }

  inline STATMECH::ParamsInterface& StatMechParamsInterface() const
  {
    return *sm_interface_ptr_;
  }

  //! shifts nodes so that proper evaluation is possible even in case of periodic boundary conditions
  virtual void UnShiftNodePosition(std::vector<double>& disp,
                                   unsigned int nnode); //!<element disp vector

private:
  /** \brief interface ptr
   *
   *  data exchange between the element and the time integrator. */
  Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_;

 Teuchos::RCP<STATMECH::ParamsInterface> sm_interface_ptr_;

  /*!
  \brief Default Constructor must not be called
  */
  Beam3Base();

};  // class Beam3Base

} // namespace ELEMENTS
} // namespace DRT

#endif /* BEAM3_BASE_H_ */
