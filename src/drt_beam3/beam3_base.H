/*----------------------------------------------------------------------------*/
/*!
\file beam3_base.H

\brief base class for all beam elements

\level 3

\maintainer Christoph Meier
*/
/*----------------------------------------------------------------------------*/

#ifndef BEAM3_BASE_H_
#define BEAM3_BASE_H_

#include "Teuchos_RCP.hpp"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"

#include "../drt_lib/drt_elementtype.H"

#include "../drt_inpar/inpar_structure.H"
#include "../drt_lib/drt_element.H"

// forward declaration ...
namespace STR {
namespace ELEMENTS {
  class ParamsInterface;
} // namespace ELEMENTS
} // namespace STR

namespace DRT {
namespace ELEMENTS {

//! A wrapper for structural elements
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

  /** \brief get access to the element reference length
   *
   *  \author grill
   *  \date 05/16 */
  virtual double GetRefLength() const = 0;      // TODO maybe skip "Get" and rename to RefLength()

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

protected:
  /** \brief get access to the interface
   *
   *  \author hiermeier
   *  \date 04/16 */
  inline STR::ELEMENTS::ParamsInterface& ParamsInterface()
  { if (not IsParamsInterface())
      dserror("The interface ptr is not set!");
    return *interface_ptr_; }

private:
  /** \brief interface ptr
   *
   *  data exchange between the element and the time integrator. */
  Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_;

  /*!
  \brief Default Constructor must not be called
  */
  Beam3Base();

};  // class Beam3Base

} // namespace ELEMENTS
} // namespace DRT

#endif /* BEAM3_BASE_H_ */
