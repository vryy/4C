/*!----------------------------------------------------------------------
\file beam3tosolidmeshtying.H

\brief meshtying element for meshtying between a 3D beam end a 2D surface (belonging to a 3D solid) element

\level 3

\maintainer Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238

*----------------------------------------------------------------------*/
#ifndef BEAM3TSMT_H
#define BEAM3TSMT_H

#include "beam3contact_defines.H"
#include "beam3contact_utils.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_sparsematrix.H"

namespace CONTACT
{

/*!
 \brief meshtying element for meshtying between a 3D beam and several 2D surface (belonging to 3D solid) elements

 \author popp (popp@lnm.mw.tum.de)
 */

class Beam3tosolidmeshtyinginterface
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
  \param element1 (in): first element (beam) of meshtying group
  \param element2 (in): list of second elements (solid surface) of meshtying group
  */
  Beam3tosolidmeshtyinginterface() {}

  /*!
  \brief Destructor
  */
  virtual ~Beam3tosolidmeshtyinginterface() {}

  //@}

  //! @name Access methods
  /*!
  \brief Get problem discretization
  */
  virtual const DRT::Discretization& ProblemDiscret() const =0;

  /*!
  \brief Get beam contact discretization
  */
  virtual const DRT::Discretization& ContactDiscret() const =0;

  /*!
  \brief Get offset of dofs between cdiscret and pdiscret
  */
  virtual const std::map<int,int>& DofOffset() const =0;

  /*!
  \brief Get first element (beam)
  */
  virtual const DRT::Element* Element1()=0;

  /*!
  \brief Get list of second elements (solid surface)
  */
  virtual const std::vector<DRT::Element*> Element2()=0;

  /*!
  \brief Get flag ndicating whether contact is active (true) or inactive (false)
  */
  virtual bool GetContactFlag()=0;


  //! @name Public evaluation methods

  /*!
  \brief Evaluate this contact element pair
  */
  virtual bool Evaluate(LINALG::SparseMatrix& stiffmatrix,
                        Epetra_Vector& fint,
                        const double& pp)=0;

  /*
  \brief Update nodal coordinates (and tangents) of all elements
  */
  virtual void UpdateElePos(Epetra_SerialDenseMatrix& newele1pos,
                            std::vector<Epetra_SerialDenseMatrix>& newele2pos)=0;

  //! return appropriate internal implementation class (acts as a simple factory)
  static Teuchos::RCP<Beam3tosolidmeshtyinginterface> Impl(const int numnodessol,
                                                           const int numnodes,
                                                           const int numnodalvalues,
                                                           const DRT::Discretization& pdiscret,
                                                           const DRT::Discretization& cdiscret,
                                                           const std::map<int,int>& dofoffsetmap,
                                                           DRT::Element* element1,
                                                           std::vector<DRT::Element*> element2,
                                                           Teuchos::ParameterList beamcontactparams);

  //@}

}; // class Beam3tosolidmeshtyinginterface




template<const int numnodessol, const int numnodes, const int numnodalvalues>
class Beam3tosolidmeshtying: public Beam3tosolidmeshtyinginterface
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
  \param element1 (in): first element (beam) of meshtying group
  \param element2 (in): list of second elements (solid surface) of meshtying group
  */
  Beam3tosolidmeshtying(const DRT::Discretization& pdiscret,
                        const DRT::Discretization& cdiscret,
                        const std::map<int,int>& dofoffsetmap,
                        DRT::Element* element1,
                        std::vector<DRT::Element*> element2,
                        Teuchos::ParameterList beamcontactparams);

  /*!
  \brief Copy Constructor
  Makes a deep copy of this contact element pair
  */
  Beam3tosolidmeshtying(const Beam3tosolidmeshtying& old);

  /*!
  \brief Destructor
  */
  virtual ~Beam3tosolidmeshtying() {};
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
  \brief Get first element (beam)
  */
  inline const DRT::Element* Element1() { return element1_;};

  /*!
  \brief Get list of second elements (solid surface)
  */
  inline const std::vector<DRT::Element*> Element2() { return element2_;};

  /*!
  \brief Get flag indicating whether contact is active (true) or inactive (false)
  */
  bool GetContactFlag() { return contactflag_; };

  //@}


  //! @name Public evaluation methods
  /*!
  \brief Evaluate this contact element pair
  */
  bool Evaluate(LINALG::SparseMatrix& stiffmatrix,
                Epetra_Vector& fint,
                const double& pp);

  /*
  \brief Update nodal coordinates (and tangents) of all elements
  */
  void UpdateElePos(Epetra_SerialDenseMatrix& newele1pos,
                    std::vector<Epetra_SerialDenseMatrix>& newele2pos);
  //@}

private:

  //! @name member variables

  //! reference to problem discretization
  const DRT::Discretization& pdiscret_;

  //! reference to beam contact discretization
  const DRT::Discretization& cdiscret_;

  //! dof offset between pdiscret and cdiscret
  const std::map<int,int>& dofoffsetmap_;

  //! first element (beam) of meshtying group
  DRT::Element* element1_;

  //! list of second elements (solid surface) of meshtying group
  std::vector<DRT::Element*> element2_;

  //! current nodal positions (and tangents) of the two elements
  LINALG::TMatrix<TYPEBTS,3*numnodes*numnodalvalues,1> ele1pos_;
  std::vector<LINALG::TMatrix<TYPEBTS,3*numnodessol,1> > ele2pos_;

  //! reference nodal positions (and tangents) of the two elements
  LINALG::TMatrix<TYPEBTS,3*numnodes*numnodalvalues,1> ele1posref_;
  std::vector<LINALG::TMatrix<TYPEBTS,3*numnodessol,1> > ele2posref_;

  //! flag indicating contact (active/inactive)
  bool contactflag_;

  //@}

  //! @name Private evaluation methods

  //@}

}; // class Beam3tosolidmeshtying
} // namespace CONTACT

#endif  // #ifndef BEAM3TSMT_H
