/*!----------------------------------------------------------------------
\file drt_potential_container.H

\brief A container class which stores all necessary data for a potential element

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef POTENTIALCONTAINER_H_
#define POTENTIALCONTAINER_H_

#include <string>
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_SerialDenseMatrix.h"

#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_element.H"


namespace POTENTIAL
{
/*!
  \brief A set of degrees of freedom for potential discretization
 */
class PotentialElementContainer : public DRT::ParObject
{
public:

  /*!
    \brief Standard Constructor

   */
  PotentialElementContainer ();


  /*!
        \brief Standard Constructor

   */
  PotentialElementContainer (
      const int                               gid,
      const int                               owner,
      const DRT::Element::DiscretizationType  distype,
      const int                               body_label,
      const LINALG::SerialDenseMatrix&        xyz_e,
      const LINALG::SerialDenseMatrix&        XYZ_e,
      const std::vector<int>&      	          lm);


  /*!
        \brief Copy Constructor

   */
  PotentialElementContainer(const PotentialElementContainer& old);
  
  
  //! Assignment operator
  // explicit PotentialElementContainer& PotentialElementContainer operator = (const PotentialElementContainer& old);


  /*!
    \brief Destructor

   */
  virtual ~PotentialElementContainer () { return; }
  
  //! unique par object
  virtual int UniqueParObjectId() const {return ParObject_Pec;}
  
  //! packs the potential element container in a vector<char>
  virtual void Pack(vector<char>& data) const;
  
  //! unpacks the potential element container from a vector<char>
  virtual void Unpack(const vector<char>& data);

  //! unpacks the potential element container from a vector<char>
  void Unpack(const vector<char>& data, int& position);

  //! returns global Id of an element
  int Id() {return gid_;}
  
  //! returns owner of an element
  int Owner() {return owner_;}
  
  //! returns discretization type of an element
  DRT::Element::DiscretizationType Shape() {return distype_;}
  
  //! returns body label of an element
  int Body_label(){return body_label_;}
  
  //! returns number of dofs of an element
  int NumDofs() {return (int) lm_.size();}
  
  //! returns number of dofs of an element
  int NumNode() {return xyz_e_.M()*xyz_e_.N();}
  
  //! returns spatial configuration
  LINALG::SerialDenseMatrix GetSpatialConfiguration(){return xyz_e_;}

  //! returns reference configuration
  LINALG::SerialDenseMatrix GetReferenceConfiguration(){return XYZ_e_;}
  
  //! returns global degrees of freedom
  std::vector<int> GetLm(){return lm_;}


private:

  //! global Id of an element
  int                                   gid_;

  //! owner of an element
  int                                   owner_;

  //! shape of an element
  DRT::Element::DiscretizationType      distype_;

  //! label of an element
  int                                   body_label_;

  //! nodal values in spatial configuration of an element
  LINALG::SerialDenseMatrix             xyz_e_;

  //! nodal values in reference configuration of an element
  LINALG::SerialDenseMatrix              XYZ_e_;

  //! global dof ids
  std::vector<int>       		             lm_;

  // assignment operator
  //virtual PotentialElementContainer operator = (const PotentialElementContainer& old);


  }; // class PotentialElementContainer
} //namespace

#endif  /*POTENTIALCONTAINER_H_*/
#endif  // #ifdef CCADISCRET











