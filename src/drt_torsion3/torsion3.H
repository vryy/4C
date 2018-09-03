/*----------------------------------------------------------------------------*/
/*!
\file torsion3.H

\brief three dimensional torsion spring element

\level 2

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/
#ifndef TORSION3_H
#define TORSION3_H


#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"

#include "../drt_lib/drt_elementtype.H"
#include "../drt_fem_general/drt_utils_integration.H"

// forward declaration ...
namespace STR
{
  namespace ELEMENTS
  {
    class ParamsInterface;
  }  // namespace ELEMENTS
}  // namespace STR

namespace DRT
{
  namespace ELEMENTS
  {
    class Torsion3Type : public DRT::ElementType
    {
     public:
      std::string Name() const { return "Torsion3Type"; }

      static Torsion3Type& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

      virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
          const double* x0, int numdf, int dimns);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static Torsion3Type instance_;
    };

    /*!
    \brief three dimensional torsion element

    \author cyron (cyron @lnm.mw.tum.de)
    */
    class Torsion3 : public DRT::Element
    {
     public:
      //! @name Friends


      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id    (in): A globally unique element id
      \param etype (in): Type of element
      \param owner (in): owner processor of the element
      */
      Torsion3(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element
      */
      Torsion3(const Torsion3& old);

      /*!
      \brief Destructor
      */
      virtual ~Torsion3();

      /*!
      \brief Deep copy this instance of Torsion3 and return pointer to the copy

      The Clone() method is used by the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed
    .
      */
      virtual DRT::Element* Clone() const;

      /*!
     \brief Get shape type of element
     */
      virtual DiscretizationType Shape() const;


      /*!
      \brief Return unique ParObject id

      Every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H
      */
      virtual int UniqueParObjectId() const { return Torsion3Type::Instance().UniqueParObjectId(); }

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

      virtual DRT::ElementType& ElementType() const { return Torsion3Type::Instance(); }

      //@}

      /*!
      \brief Return number of lines to this element
      */
      virtual int NumLine() const { return 1; }


      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element
      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Lines();


      /*!
      \brief Get number of degrees of freedom of a single node
      */
      virtual int NumDofPerNode(const DRT::Node& node) const
      {
        /*note: this is not necessarily the number of DOF assigned to this node by the
         *discretization finally, but only the number of DOF requested for this node by this
         *element; the discretization will finally assign the maximal number of DOF to this node
         *requested by any element connected to this node*/
        return 3;
      }


      /*!
      \brief Get number of degrees of freedom per element not including nodal degrees of freedom
      */
      virtual int NumDofPerElement() const { return 0; }

      /*!
      \brief Print this element
      */
      virtual void Print(std::ostream& os) const;


      //@}

      //! @name Construction


      /*!
      \brief Read input for this element

      This class implements a dummy of this method that prints a warning and
      returns false. A derived class would read one line from the input file and
      store all necessary information.

      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef);


      //@}


      //! @name Evaluation methods


      /*!
      \brief Evaluate an element

      An element derived from this class uses the Evaluate method to receive commands
      and parameters from some control routine in params and evaluates element matrices and
      vectors accoring to the command in params.

      \note This class implements a dummy of this method that prints a warning and
            returns false.

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param lm (in)            : location vector of this element
      \param elemat1 (out)      : matrix to be filled by element depending on commands
                                  given in params
      \param elemat2 (out)      : matrix to be filled by element depending on commands
                                  given in params
      \param elevec1 (out)      : vector to be filled by element depending on commands
                                  given in params
      \param elevec2 (out)      : vector to be filled by element depending on commands
                                  given in params
      \param elevec3 (out)      : vector to be filled by element depending on commands
                                  given in params
      \return 0 if successful, negative otherwise
      */
      virtual int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);

      /*!
      \brief Evaluate a Neumann boundary condition

      An element derived from this class uses the EvaluateNeumann method to receive commands
      and parameters from some control routine in params and evaluates a Neumann boundary condition
      given in condition

      \note This class implements a dummy of this method that prints a warning and
            returns false.

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param condition (in)     : The condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : Force vector to be filled by element

      \return 0 if successful, negative otherwise
      */
      virtual int EvaluateNeumann(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1 = NULL);

      //@}

      /*! \brief set the parameter interface ptr for the solid elements
       *
       *  \param p (in): Parameter list coming from the time integrator.
       *
       *  \author hiermeier
       *  \date 04/16 */
      virtual void SetParamsInterfacePtr(const Teuchos::ParameterList& p);

      /*! \brief returns true if the parameter interface is defined and initialized, otherwise false
       *
       *  \author hiermeier
       *  \date 04/16 */
      virtual inline bool IsParamsInterface() const { return (not interface_ptr_.is_null()); }

      /*! \brief get access to the parameter interface pointer
       *
       *  \author hiermeier
       *  \date 04/16 */
      virtual Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> ParamsInterfacePtr();

     protected:
      /** \brief get access to the interface
       *
       *  \author hiermeier
       *  \date 04/16 */
      inline STR::ELEMENTS::ParamsInterface& ParamsInterface()
      {
        if (not IsParamsInterface()) dserror("The interface ptr is not set!");
        return *interface_ptr_;
      }

     private:
      //! possible bending potentials
      enum BendingPotential
      {
        quadratic,
        cosine
      };

      /*! \brief interface ptr
       *
       *  data exchange between the element and the time integrator. */
      Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_;

      // container containing all kind of Torsion3 thisandthat
      DRT::Container data_;

      //! Bending potential
      BendingPotential bendingpotential_;

      //! @name Internal calculation methods

      //! calculation of nonlinear stiffness and mass matrix
      void t3_nlnstiffmass(std::vector<double>& disp, Epetra_SerialDenseMatrix* stiffmatrix,
          Epetra_SerialDenseMatrix* massmatrix, Epetra_SerialDenseVector* force);

      //! calculation of elastic energy
      void t3_energy(Teuchos::ParameterList& params, std::vector<double>& disp,
          Epetra_SerialDenseVector* intenergy);


      //@}


      //! @name Methods for Brownian dynamics simulations

      //! shifts nodes so that proper evaluation is possible even in case of periodic boundary
      //! conditions
      template <int nnode, int ndim>                  // number of nodes, number of dimensions
      void NodeShift(Teuchos::ParameterList& params,  //!< parameter list
          std::vector<double>& disp);                 //!< element disp vector

      //@}

      // don't want = operator
      Torsion3& operator=(const Torsion3& old);


    };  // class Torsion3



    // << operator
    std::ostream& operator<<(std::ostream& os, const DRT::Element& ele);


  }  // namespace ELEMENTS
}  // namespace DRT



#endif  // #ifndef TORSION3_H
