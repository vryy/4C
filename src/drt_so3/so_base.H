/*----------------------------------------------------------------------*/
/*!
 \file so_base.H

\brief a common base class for all solid elements

\level 2

\maintainer Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271

 *----------------------------------------------------------------------*/


#ifndef SO_BASE_H_
#define SO_BASE_H_

#include "../drt_inpar/inpar_structure.H"
#include "../drt_lib/drt_element.H"

// forward declaration ...
namespace STR {
namespace ELEMENTS {
  class ParamsInterface;
} // namespace ELEMENTS
} // namespace STR
namespace MAT {
  class So3Material;
} // namespace MAT
namespace DRT {
namespace ELEMENTS {

//! A wrapper for structural elements
class So_base : public DRT::Element
{

public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  \param id    (in): A globally unique element id
  \param owner (in): owner processor of the element
  */
  So_base(int id, int owner);

  /*!
  \brief Copy Constructor

  Makes a deep copy of a Element

  */
  So_base(const So_base& old);

  //! Destructor
  virtual ~So_base(){};

  //@}

  /*!
  \brief Pack this class so it can be communicated

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Pack(DRT::PackBuffer& data) const;

  /*!
  \brief Unpack data from a char vector into this class

  \ref Pack and \ref Unpack are used to communicate this element

  */
  virtual void Unpack(const std::vector<char>& data);

  // get the kinematic type from the element
  INPAR::STR::KinemType KinematicType() const
  {
    return kintype_;
  }

  // get the kinematic type from the element
  void SetKinematicType(INPAR::STR::KinemType kintype)
  {
    kintype_=kintype;
  }

  /*!
  \brief Does this element use EAS?

  ToDo: This function can be declared as pure virtual and each concrete derived
        class has to implement this function. This can be done during the up-coming
        cleaning procedure.                                      hiermeier 09/15
  */
  virtual bool HaveEAS() const { return false; };

  /*!
  \brief Return the material of this element

  Note: The input parameter nummat is not the material number from input file
        as in SetMaterial(int matnum), but the number of the material within
        the vector of materials the element holds

  \param nummat (in): number of requested material
  */
  virtual Teuchos::RCP<MAT::So3Material> SolidMaterial(int nummat=0) const;

  /// evaluate Cauchy stress at given point in parameter space
  virtual void GetCauchyAtXi(
      const LINALG::Matrix<3,1>& xi,
      const std::vector<double>& disp,
      const LINALG::Matrix<3,1>& n,
      const LINALG::Matrix<3,1>& t,
      double& sigma_nt,
      Epetra_SerialDenseMatrix* dsntdd,
      Epetra_SerialDenseMatrix* d2sntdd2,
      Epetra_SerialDenseMatrix* d2sntDdDn,
      Epetra_SerialDenseMatrix* d2sntDdDt,
      Epetra_SerialDenseMatrix* d2sntDdDpxi,
      LINALG::Matrix<3,1>* dsntdn,
      LINALG::Matrix<3,1>* dsntdt,
      LINALG::Matrix<3,1>* dsntdpxi,
      const std::vector<double>* temp=NULL,
      Epetra_SerialDenseMatrix* dsntdT=NULL,
      Epetra_SerialDenseMatrix* d2sntDdDT=NULL)
  { dserror("not implemented for chosen solid element"); }

  /// evaluate Cauchy stress at given point in parameter space
  virtual void GetCauchyAtXi(
      const LINALG::Matrix<2,1>& xi,
      const std::vector<double>& disp,
      const LINALG::Matrix<2,1>& n,
      const LINALG::Matrix<2,1>& t,
      double& sigma_nt,
      Epetra_SerialDenseMatrix* dsntdd,
      Epetra_SerialDenseMatrix* d2sntdd2,
      Epetra_SerialDenseMatrix* d2sntDdDn,
      Epetra_SerialDenseMatrix* d2sntDdDt,
      Epetra_SerialDenseMatrix* d2sntDdDpxi,
      LINALG::Matrix<2,1>* dsntdn,
      LINALG::Matrix<2,1>* dsntdt,
      LINALG::Matrix<2,1>* dsntdpxi,
      const std::vector<double>* temp=NULL,
      Epetra_SerialDenseMatrix* dsntdT=NULL,
      Epetra_SerialDenseMatrix* d2sntDdDT=NULL
      )
  { dserror("not implemented for chosen solid element"); }

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

protected:
  /** \brief get access to the interface
   *
   *  \author hiermeier
   *  \date 04/16 */
  inline DRT::ELEMENTS::ParamsInterface& ParamsInterface()
  { if (not IsParamsInterface())
      dserror("The interface ptr is not set!");
    return *interface_ptr_; }

  /** \brief get access to the structure interface
   *
   *  \author vuong
   *  \date 11/16 */
  STR::ELEMENTS::ParamsInterface& StrParamsInterface();

protected:
  //! kinematic type
  INPAR::STR::KinemType kintype_;

private:
  /** \brief interface ptr
   *
   *  data exchange between the element and the time integrator. */
  Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> interface_ptr_;

private:
  /*!
  \brief Default Constructor must not be called

  */
  So_base();

};  // class So_base

}
}

#endif /* SO_BASE_H_ */
