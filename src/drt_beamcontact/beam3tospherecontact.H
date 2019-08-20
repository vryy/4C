/*----------------------------------------------------------------------------*/
/*! \file

\brief class to handle contact between a 3D beam element and a rigid sphere

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/
#ifndef BEAMTOSPHEREC_H
#define BEAMTOSPHEREC_H

#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"

#include "../drt_beaminteraction/beam3contact_defines.H"
#include "../drt_lib/drt_condition.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../headers/FAD_utils.H"

namespace CONTACT
{
  class Beam3tospherecontactinterface
  {
   public:
    //! @name Friends
    // no friend classes defined
    //@}

    //! @name Constructors and destructors and related methods
    /*!
    \brief Standard Constructor

    */
    Beam3tospherecontactinterface() {}

    /*!
    \brief Destructor
    */
    virtual ~Beam3tospherecontactinterface(){};

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

    /*!
    \brief Get second element

    */
    virtual const DRT::Element* Element2() = 0;

    /*!
    \brief Get gap of this contact pair

    */
    virtual double GetGap() = 0;

    /*!
    \brief Get closest point of this contact pair

    */
    virtual double GetClosestPoint() const = 0;

    /*!
    \brief Get flag indicating whether contact is active (true) or inactive (false)

    */
    virtual bool GetContactFlag() = 0;

    /*!
    \brief Get coordinates of contact point of element1 and element2

    */
    virtual const Epetra_SerialDenseVector GetX1() = 0;
    virtual const Epetra_SerialDenseVector GetX2() = 0;

    //@}

    //! @name Public evaluation methods

    //! return appropriate internal implementation class (acts as a simple factory)
    static Teuchos::RCP<Beam3tospherecontactinterface> Impl(const int numnodes,
        const int numnodalvalues, const DRT::Discretization& pdiscret,
        const DRT::Discretization& cdiscret, const std::map<int, int>& dofoffsetmap,
        DRT::Element* element1, DRT::Element* element2);

    /*!
    \brief Evaluate this contact element pair
    */
    virtual bool Evaluate(LINALG::SparseMatrix& stiffmatrix, Epetra_Vector& fint, double& pp) = 0;

    /*
    \brief Update nodal coordinates of both elements
    */
    virtual void UpdateElePos(
        Epetra_SerialDenseMatrix& newele1pos, Epetra_SerialDenseMatrix& newele2pos) = 0;

    virtual void Print() const = 0;
    //@}
  };  // class Beam3tospherecontactinterface


  template <const int numnodes, const int numnodalvalues>
  class Beam3tospherecontact : public Beam3tospherecontactinterface
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
    Beam3tospherecontact(const DRT::Discretization& pdiscret, const DRT::Discretization& cdiscret,
        const std::map<int, int>& dofoffsetmap, DRT::Element* element1, DRT::Element* element2);

    /*!
    \brief Copy Constructor

    Makes a deep copy of this contact element pair

    */
    Beam3tospherecontact(const Beam3tospherecontact& old);

    /*!
    \brief Destructor

    */
    virtual ~Beam3tospherecontact(){};

    //@}

    //! @name Access methods

    /*!
    \brief Get problem discretization
    */
    inline const DRT::Discretization& ProblemDiscret() const { return pdiscret_; };

    /*!
    \brief Get beam contact discretization
    */
    inline const DRT::Discretization& ContactDiscret() const { return cdiscret_; };

    /*!
    \brief Get offset of dofs between cdiscret and pdiscret
    */
    inline const std::map<int, int>& DofOffset() const { return dofoffsetmap_; };

    /*!
    \brief Get first element
    */
    inline const DRT::Element* Element1() { return element1_; };

    /*!
    \brief Get first element
    */
    inline const DRT::Element* Element2() { return element2_; };

    /*!
    \brief Get gap of this contact pair
    */
    double GetGap() { return FADUTILS::CastToDouble(gap_); };

    /*!
    \brief Get closest point of this contact pair
    */
    double GetClosestPoint() const { return FADUTILS::CastToDouble(xicontact_); };

    /*!
    \brief Get flag indicating whether contact is active (true) or inactive (false)
    */
    bool GetContactFlag()
    {
      return (contactflag_ || nodalcontactflag_[0] || nodalcontactflag_[1]);
    };

    /*!
    \brief Get coordinates of contact point of element1 and element2
    */
    const Epetra_SerialDenseVector GetX1()
    {
      Epetra_SerialDenseVector r1;
      r1.Resize(3);
      for (int i = 0; i < 3; i++) r1(i) = FADUTILS::CastToDouble(x1_(i));

      return r1;
    };

    const Epetra_SerialDenseVector GetX2()
    {
      Epetra_SerialDenseVector r2;
      r2.Resize(3);
      for (int i = 0; i < 3; i++) r2(i) = FADUTILS::CastToDouble(x2_(i));

      return r2;
    };

    void Print() const;
    //@}

    //! @name Public evaluation methods

    /*!
    \brief Evaluate this contact element pair
    */
    bool Evaluate(LINALG::SparseMatrix& stiffmatrix, Epetra_Vector& fint, double& pp);

    /*
    \brief Update nodal coordinates of both elements
    */
    void UpdateElePos(Epetra_SerialDenseMatrix& newele1pos, Epetra_SerialDenseMatrix& newele2pos);

    //@}

   private:
    //! @name member variables

    //! reference to problem discretization
    const DRT::Discretization& pdiscret_;

    //! reference to beam contact discretization
    const DRT::Discretization& cdiscret_;

    //! dof offset between pdiscret and cdiscret
    const std::map<int, int>& dofoffsetmap_;

    //! first element of contact pair
    DRT::Element* element1_;

    //! second element of contact pair
    DRT::Element* element2_;

    //! current node coordinates of the two elements
    LINALG::TMatrix<TYPE, 3 * numnodes * numnodalvalues, 1> ele1pos_;
    LINALG::TMatrix<TYPE, 3, 1> ele2pos_;

    //! gap function
    TYPE gap_;

    //! flag indicating contact (active/inactive)
    bool contactflag_;

    //! flag indicating contact of beam end points (nodes) (active/inactive)
    std::vector<bool> nodalcontactflag_;

    //! resulting nodal contact forces on ele 1/2
    LINALG::TMatrix<TYPE, 3 * numnodes * numnodalvalues, 1> fc1_;
    LINALG::TMatrix<TYPE, 3, 1> fc2_;

    //! coordinates of contact point on center lines of beams
    LINALG::TMatrix<TYPE, 3, 1> x1_;
    LINALG::TMatrix<TYPE, 3, 1> x2_;

    //! parameter value of contact point
    TYPE xicontact_;

    //! normal vector of current time step
    LINALG::TMatrix<TYPE, 3, 1> normal_;

    //@}

    //! @name Private evaluation methods

    /*!
    \brief Find contact point via closest point projection
    */
    void ClosestPointProjection();

    /*!
    \brief Utility method for CPP (evaluate nonlinear function f)
    */
    void EvaluateOrthogonalityCondition(TYPE& f, const LINALG::TMatrix<TYPE, 3, 1>& delta_x,
        const double norm_delta_x, const LINALG::TMatrix<TYPE, 3, 1>& dx1);

    /*!
    \brief Utility method for CPP (evaluate Jacobian of nonlinear function f)
    */
    void EvaluateLinOrthogonalityCondition(TYPE& df, LINALG::TMatrix<TYPE, 3, 1>& delta_x,
        const double norm_delta_x, const LINALG::TMatrix<TYPE, 3, 1>& dx1,
        const LINALG::TMatrix<TYPE, 3, 1>& ddx1);

    /*!
    \brief Evaluate and assemble contact forces
    */
    void EvaluateFcContact(const double& pp, const TYPE& gap,
        const LINALG::TMatrix<TYPE, 3, 1>& normal, Epetra_Vector& fint,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i, const bool contactactive);

    /*!
    \brief Evaluate and assemble contact stiffness
    */
    void EvaluateStiffcContact(const double& pp, const TYPE& gap,
        const LINALG::TMatrix<TYPE, 3, 1>& normal, const TYPE& norm,
        LINALG::SparseMatrix& stiffmatrix, const LINALG::TMatrix<TYPE, 3, 1>& x1,
        const LINALG::TMatrix<TYPE, 3, 1>& x2, const LINALG::TMatrix<TYPE, 3, 1>& dx1,
        const LINALG::TMatrix<TYPE, 3, 1>& ddx1,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xi,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xixi, bool activecontact,
        bool linxi = true);

    /*!
    \brief Compute normal vector in contact point
    */
    void ComputeNormal(LINALG::TMatrix<TYPE, 3, 1>& normal, TYPE& gap, TYPE& norm,
        const LINALG::TMatrix<TYPE, 3, 1>& x1, const LINALG::TMatrix<TYPE, 3, 1>& x2);

    /*!
    \brief Evaluate gap function
    */
    void ComputeGap(TYPE& gap, const TYPE& norm);

    /*!
    \brief Compute coordinates and their derivatives from the discretization
    */
    void ComputeCoordsAndDerivs(LINALG::TMatrix<TYPE, 3, 1>& x1, LINALG::TMatrix<TYPE, 3, 1>& x2,
        LINALG::TMatrix<TYPE, 3, 1>& dx1, LINALG::TMatrix<TYPE, 3, 1>& ddx1,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xi,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xixi);

    /*!
    \brief Get shape functions and their derivatives at eta
    */
    void GetShapeFunctions(LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i,
        LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xi,
        LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xixi, const TYPE& xi);

    /*!
    \brief Compute linearizations of contact point
    */
    void ComputeLinXi(LINALG::TMatrix<TYPE, 3 * numnodes * numnodalvalues + 3, 1>& delta_xi,
        const LINALG::TMatrix<TYPE, 3, 1>& x1, const LINALG::TMatrix<TYPE, 3, 1>& x2,
        const LINALG::TMatrix<TYPE, 3, 1>& dx1, const LINALG::TMatrix<TYPE, 3, 1>& ddx1,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xi);

    /*!
    \brief Compute linearization of gap
    */
    void ComputeLinGap(LINALG::TMatrix<TYPE, 3 * numnodes * numnodalvalues + 3, 1>& delta_gap,
        LINALG::TMatrix<TYPE, 3 * numnodes * numnodalvalues + 3, 1>& delta_xi,
        const LINALG::TMatrix<TYPE, 3, 1>& x1, const LINALG::TMatrix<TYPE, 3, 1>& x2,
        const LINALG::TMatrix<TYPE, 3, 1>& dx1,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i, const TYPE& normdist,
        const LINALG::TMatrix<TYPE, 3, 1>& normal, const TYPE& norm, const TYPE& gap,
        LINALG::TMatrix<TYPE, 3, 3 * numnodes * numnodalvalues + 3>& delta_x1_minus_x2);

    /*!
    \brief Compute linearization of normal
    */

    void ComputeLinNormal(LINALG::TMatrix<TYPE, 3, 3 * numnodes * numnodalvalues + 3>& delta_normal,
        const LINALG::TMatrix<TYPE, 3 * numnodes * numnodalvalues + 3, 1>& delta_xi,
        const LINALG::TMatrix<TYPE, 3, 1>& normal, const TYPE& norm_delta_x,
        const LINALG::TMatrix<TYPE, 3, 1>& x1_xi,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i);

    /*!
    \brief Compute normal contact disctance
    */
    void ComputeDistance(LINALG::TMatrix<TYPE, 3, 1>& distance, TYPE& normdist,
        const LINALG::TMatrix<TYPE, 3, 1>& normal, const TYPE& norm);

    /*!
    \brief Check if contact is active and set flag accordingly
    */
    void CheckContactStatus(double& pp);

    /*!
    \brief Get global dofs of a node

    Internally this method first extracts the dofs of the given node
    in the beam contact discretization (which has its own dofs) and
    then transfers these dofs to their actual GIDs in the underlying
    problem discretization by means of the pre-computed dofoffsetmap_.

    */
    std::vector<int> GetGlobalDofs(DRT::Node* node);

    //@}

  };  // class Beam3tospherecontact
}  // namespace CONTACT

#endif  // #ifndef BEAMTOSPHEREC_H
