/*----------------------------------------------------------------------*/
/*!
\file beam3tospherepotential.H

\brief One beam-to-sphere potential-based interacting pair

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/

#ifndef BEAM3TOSPHEREPOT_H
#define BEAM3TOSPHEREPOT_H

#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include <Sacado.hpp>

#include "../drt_beaminteraction/beam3contact_defines.H"
#include "../drt_beaminteraction/beam3contact_tangentsmoothing.H"
#include "../drt_lib/standardtypes_cpp.H"

// forward declarations
namespace DRT {
namespace ELEMENTS {
  class Rigidsphere;
  class Beam3Base;
}
}

namespace CONTACT
{

/*!
 \brief class for potential-based interaction of a 3D beam element and a rigid sphere
 */

class Beam3tospherepotentialinterface
{
public:

  //! @name Friends
  //no friend classes defined
  //@}


  Beam3tospherepotentialinterface() {}

  /*!
  \brief Destructor
  */
  virtual ~Beam3tospherepotentialinterface() {}

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
  \brief Get first element
  */
  virtual const DRT::Element* Element1()=0;

  /*!
  \brief Get second element
  */
  virtual const DRT::Element* Element2()=0;


  //! @name Public evaluation methods
  /*!
  \brief Evaluate this contact element pair
  */
  virtual bool Evaluate(LINALG::SparseMatrix& stiffmatrix,
                        Epetra_Vector& fint,
                        const std::vector<DRT::Condition*> chargeconds,
                        const double k,
                        const double m) =0;

  //! return appropriate internal implementation class (acts as a simple factory)
  static Teuchos::RCP<Beam3tospherepotentialinterface> Impl(const int numnodes,
                                     const int numnodalvalues,
                                     const DRT::Discretization& pdiscret,
                                     const DRT::Discretization& cdiscret,
                                     const std::map<int,int>& dofoffsetmap,
                                     DRT::Element* element1,
                                     DRT::Element* element2,
                                     Teuchos::ParameterList beampotparams);

  /*
  \brief Update nodal coordinates of both elements at the beginning of a new time step!
  */
  virtual void UpdateElePos(Epetra_SerialDenseMatrix& newele1pos,
                            Epetra_SerialDenseMatrix& newele2pos)=0;

}; // class Beam3tospherepotentialinterface



template<const int numnodes , const int numnodalvalues>
class Beam3tospherepotential: public CONTACT::Beam3tospherepotentialinterface
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
  Beam3tospherepotential(  const DRT::Discretization& pdiscret,
                    const DRT::Discretization& cdiscret,
                    const std::map<int,int>& dofoffsetmap,
                    DRT::Element* element1,
                    DRT::Element* element2,
                    Teuchos::ParameterList beampotparams);

  /*!
  \brief Copy Constructor
  Makes a deep copy of this contact element pair
  */
  Beam3tospherepotential(const Beam3tospherepotential& old);

  /*!
  \brief Destructor
  */
  virtual ~Beam3tospherepotential() {};
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
  \brief Get first element
  */
  inline const DRT::Element* Element1() { return element1_;};

  /*!
  \brief Get first element
  */
  inline const DRT::Element* Element2() { return element2_;};


  //! @name Public evaluation methods
  /*!
  \brief Evaluate this contact element pair
  */

  /*!
  \brief Evaluate this contact element pair
  */
  bool Evaluate(LINALG::SparseMatrix& stiffmatrix,
                Epetra_Vector& fint,
                const std::vector<DRT::Condition*> chargeconds,
                const double k,
                const double m);

  /*
  \brief Update nodal coordinates of both elements at the beginning of a new time step!
  */
  void UpdateElePos(Epetra_SerialDenseMatrix& newele1pos,
                    Epetra_SerialDenseMatrix& newele2pos);

  //@}
private:

  //! @name member variables

  //! reference to problem discretization
  const DRT::Discretization& pdiscret_;

  //! reference to beam contact discretization
  const DRT::Discretization& cdiscret_;

  //! dof offset between pdiscret and cdiscret
  const std::map<int,int>& dofoffsetmap_;

  //! first element of contact pair
  DRT::ELEMENTS::Beam3Base* element1_;

  //! second element of contact pair
  DRT::ELEMENTS::Rigidsphere* element2_;

  //! beam contact parameter list
  Teuchos::ParameterList bpotparams_;

  //! line and point charge condition
  std::vector<DRT::Condition*> chargeconds_;

  //! current node coordinates of the two elements
  LINALG::TMatrix<TYPE,3*numnodes*numnodalvalues,1> ele1pos_;
  LINALG::TMatrix<TYPE,3,1> ele2pos_;

  //! parameters of the applied (point-point) potential law Phi(r)=k_ * r^(-m_)
  double k_;
  double m_;

  //!current Newton iteration
  int iter_;

  //!current time step
  int numstep_;

  //!time step size
  double dt_;

  //! initial beam element length
  double ele1length_;

  //! Cross-section radius of beam 1
  double radius1_;

  //! Cross-section radius of beam 2
  double radius2_;

  //! resulting forces on element 1 and 2
  LINALG::TMatrix<TYPE, 3*numnodes*numnodalvalues, 1> fpot1_;
  LINALG::TMatrix<TYPE, 3, 1> fpot2_;

  //! stiffness contributions
  LINALG::TMatrix<TYPE, 3*numnodes*numnodalvalues, 3*numnodes*numnodalvalues+3> stiffpot1_;
  LINALG::TMatrix<TYPE, 3, 3*numnodes*numnodalvalues+3> stiffpot2_;

  //@}

  //! @name Private evaluation methods

  /*!
  \brief Evaluate forces and stiffness contribution resulting from potential-based interaction
  */
  void EvaluateFpotandStiffpot_LargeSepApprox();

  /*!
  \brief Assemble forces fpot1/2 into fint
  */
  void AssembleFpot(Epetra_Vector& fint);

  /*!
  \brief Assemble stiffness contributions stiffpot1/2 into stiffmatrix
  */
  void AssembleStiffpot(LINALG::SparseMatrix& stiffmatrix);

  /*!
  \brief Calculate shape function values for given parameter values
  */
  void GetShapeFunctions( std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> >& N1_i,
                          std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> >& N1_i_xi,
                          DRT::UTILS::IntegrationPoints1D& gausspoints);

  /*!
  \brief Compute coordinates of centreline points from the discretization
  */
  void ComputeCoords(LINALG::TMatrix<TYPE,3,1>& r,
                     const LINALG::Matrix<1,numnodes*numnodalvalues>& N_i,
                     const LINALG::TMatrix<TYPE,3*numnodes*numnodalvalues,1> elepos);

  /*!
  \brief Get global dofs of a node

  Internally this method first extracts the dofs of the given node
  in the beam contact discretization (which has its own dofs) and
  then transfers these dofs to their actual GIDs in the underlying
  problem discretization by applying the pre-computed dofoffset_.
  */
  std::vector<int> GetGlobalDofs(const DRT::Node* node);

  //@}

}; // class Beam3tospherepotential
} // namespace CONTACT

#endif  // #ifndef BEAM3TOSPHEREPOT_H
