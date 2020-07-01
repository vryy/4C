/*----------------------------------------------------------------------*/
/*! \file
\brief One-dimensional artery element


\level 3

*----------------------------------------------------------------------*/
#ifndef ARTERY_H
#define ARTERY_H

#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"


#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_inpar/inpar_bio.H"
#include "../drt_lib/drt_utils_nullspace.H"

// forward declarations
struct _MATERIAL;



namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    class ArteryType : public DRT::ElementType
    {
     public:
      std::string Name() const { return "ArteryType"; }

      static ArteryType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
      {
        numdf = dwele->NumDofPerNode(*(dwele->Nodes()[0]));
        dimns = numdf;
        nv = numdf;
      }

      virtual void ComputeNullSpace(
          DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
      {
        DRT::UTILS::ComputeFluidDNullSpace(dis, ns, x0, numdf, dimns);
      }

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static ArteryType instance_;
    };

    /*!
    \brief A C++ wrapper for the artery element

    */
    class Artery : public DRT::Element
    {
     public:
      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      */
      Artery(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      Artery(const Artery& old);

      /*!
      \brief Deep copy this instance of Artery and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const;

      /*!
      \brief Get shape type of element
      */
      virtual DiscretizationType Shape() const;

      /*!
      \brief Return number of lines of this element
      */
      virtual int NumLine() const
      {
        if (NumNode() == 2)
          return 1;
        else
        {
          dserror("Could not determine number of lines");
          return -1;
        }
      }

      /*!
      \brief Return number of surfaces of this element (always 1)
      */
      virtual int NumSurface() const { return -1; }

      /*!
      \brief Return number of volumes of this element (always 1)
      */
      virtual int NumVolume() const { return -1; }

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      virtual int UniqueParObjectId() const { return ArteryType::Instance().UniqueParObjectId(); }

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

      /*!
       \brief Get vector of Teuchos::RCPs to the lines of this element
       */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Lines();

      /*!
      \brief Destructor

      */
      virtual ~Artery();

      //@}

      //! @name Acess methods


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual DRT::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can redecide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      virtual int NumDofPerNode(const DRT::Node& node) const
      {
        switch (impltype_)
        {
          case INPAR::ARTDYN::impltype_lin_exp:
          {
            return 2;
            break;
          }
          case INPAR::ARTDYN::impltype_pressure_based:
          {
            return 1;
            break;
          }
          default:
          {
            dserror("Defined problem type %d does not exist!!", impltype_);
            break;
          }
        }

        return 0;
      }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual DRT::Element)

      The element decides how many element degrees of freedom it has.
      It can redecide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      virtual int NumDofPerElement() const { return 0; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      virtual ArteryType& ElementType() const { return ArteryType::Instance(); }

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef);

      /**
       * Set diameter in material
       * @param diam: diameter to be set
       */
      void SetDiamInMaterial(const double* diam);

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      An element derived from this class uses the Evaluate method to receive commands
      and parameters from some control routine in params and evaluates element matrices and
      vectors accoring to the command in params.

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param la (in)            : location data for all dofsets of the discretization
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
          LocationArray& la, Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
          Epetra_SerialDenseVector& elevec3);

      int ScatraEvaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
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

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a line Neumann condition on the artery element

      \return 0 if successful, negative otherwise
      */
      virtual int EvaluateDirichlet(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1);

      /*!
      \brief return implementation type (physics)

      */
      INPAR::ARTDYN::ImplType ImplType() { return impltype_; }


      //@}


      //! @name Other

      DRT::UTILS::GaussRule1D GaussRule() const { return gaussrule_; }


      //@}


     private:
      //! implementation type (physics)
      INPAR::ARTDYN::ImplType impltype_;

      //! Gaussrule
      DRT::UTILS::GaussRule1D gaussrule_;

      //! container containing all kind of Artery thisandthat
      DRT::Container data_;


      // internal calculation methods

      // don't want = operator
      Artery& operator=(const Artery& old);


      /// set number of gauss points to element shape default
      DRT::UTILS::GaussRule1D getOptimalGaussrule(const DiscretizationType& distype);

      /*!
       * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are
       * necessary \return boolean indicating higher order status
       */
      bool isHigherOrderElement(
          const DRT::Element::DiscretizationType distype  ///< discretization type
          ) const;


    };  // class Artery


    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================


  }  // namespace ELEMENTS
}  // namespace DRT



#endif  // #ifndef ARTERY_H
