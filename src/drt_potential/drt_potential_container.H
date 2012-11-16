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
#ifndef POTENTIALCONTAINER_H_
#define POTENTIALCONTAINER_H_

#include <string>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_SerialDenseMatrix.h"

#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_parobjectfactory.H"


namespace POTENTIAL
{

class PecType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "PecType"; }

  static PecType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static PecType instance_;
};


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
      const DRT::Element::DiscretizationType  distype,
      const int                               body_label,
      const double                            beta,
      const Epetra_SerialDenseMatrix&         xyz_e,
      const Epetra_SerialDenseMatrix&         XYZ_e,
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
  virtual ~PotentialElementContainer ();

  //! unique par object
  virtual int UniqueParObjectId() const {return PecType::Instance().UniqueParObjectId(); }

  //! packs the potential element container in a vector<char>
  virtual void Pack(DRT::PackBuffer& data) const;

  //! unpacks the potential element container from a vector<char>
  virtual void Unpack(const std::vector<char>& data);

  //! unpacks the potential element container from a vector<char>
  void Unpack(const std::vector<char>& data, std::vector<char>::size_type& position);

  //! print
  void Print();

  //! returns global Id of an element
  int Id() {return gid_;}

  //! returns discretization type of an element
  DRT::Element::DiscretizationType Shape() {return distype_;}

  //! returns body label of an element
  int Body_label(){return body_label_;}

  //! returns atom density of an element
  double Beta(){return beta_;}

  //! returns number of dofs of an element
  int NumDofs() {return (int) lm_.size();}

  //! returns number of dofs of an element
  int NumNode() {return xyz_e_.M();}

  //! returns spatial configuration
  Epetra_SerialDenseMatrix GetSpatialConfiguration(){return xyz_e_;}

  //! returns reference configuration
  Epetra_SerialDenseMatrix GetReferenceConfiguration(){return XYZ_e_;}

  //! returns global degrees of freedom
  std::vector<int> GetLm(){return lm_;}


private:

  //! global Id of an element
  int                                         gid_;

  //! shape of an element
  DRT::Element::DiscretizationType            distype_;

  //! label of an element
  int                                         body_label_;

  //! atom density
  double                                      beta_;

  //! nodal values in spatial configuration of an element
  Epetra_SerialDenseMatrix                    xyz_e_;

  //! nodal values in reference configuration of an element
  Epetra_SerialDenseMatrix                    XYZ_e_;

  //! global dof ids
  std::vector<int>                            lm_;

  // assignment operator
  //virtual PotentialElementContainer operator = (const PotentialElementContainer& old);


  }; // class PotentialElementContainer
} //namespace

#endif  /*POTENTIALCONTAINER_H_*/











