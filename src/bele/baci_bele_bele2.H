/*----------------------------------------------------------------------*/
/*! \file

\brief 2D boundary elment


\level 2
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_BELE_BELE2_H
#define BACI_BELE_BELE2_H


#include "baci_lib_element.H"
#include "baci_lib_elementtype.H"
#include "baci_lib_node.H"
#include "baci_lib_parobjectfactory.H"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>



namespace DRT
{
  // forward declarations
  class Discretization;


  namespace ELEMENTS
  {
    class Bele2Type : public DRT::ElementType
    {
     public:
      std::string Name() const { return "Bele2Type"; }

      static Bele2Type& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

      virtual Teuchos::SerialDenseMatrix<int, double> ComputeNullSpace(
          DRT::Node& node, const double* x0, const int numdof, const int dimnsp);

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static Bele2Type instance_;
    };

    /*!
     * A 2D boundary element
     *
     * It can be used to have a boundary discretization
     * of surface/boundary elements. They can be of any 1d shape (line2,line3)
     *
     * The number of dof per node is set to 2, so we can define displacement vectors by
     * using FillComplete on the boundary discretization.
     *
     */
    class Bele2 : public DRT::Element
    {
     public:
      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor
      */
      explicit Bele2(int id,  ///< A unique global id
          int owner           ///< proc num that owns this element
      );

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      explicit Bele2(const Bele2& old);

      DRT::Element* Clone() const;
      virtual DiscretizationType Shape() const;
      virtual int NumLine() const { return 1; }
      virtual int NumSurface() const { return -1; }
      virtual int NumVolume() const { return -1; }
      virtual std::vector<Teuchos::RCP<DRT::Element>> Lines();
      virtual std::vector<Teuchos::RCP<DRT::Element>> Surfaces();
      virtual int UniqueParObjectId() const { return Bele2Type::Instance().UniqueParObjectId(); }
      virtual void Pack(DRT::PackBuffer& data) const;
      virtual void Unpack(const std::vector<char>& data);

      /*!
      \brief Destructor

      */
      virtual ~Bele2();

      //@}

      //! @name Access methods

      virtual int NumDofPerNode(const DRT::Node&) const { return 2; }
      virtual int NumDofPerElement() const { return 0; }
      void Print(std::ostream& os) const;

      virtual DRT::ElementType& ElementType() const { return Bele2Type::Instance(); }

      //@}

      //! @name Evaluation

      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
          CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseVector& elevec2, CORE::LINALG::SerialDenseVector& elevec3);

      virtual int EvaluateNeumann(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          CORE::LINALG::SerialDenseVector& elevec1,
          CORE::LINALG::SerialDenseMatrix* elemat1 = NULL);

      /// Read input for this element
      virtual bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef);
      //@}

      //! @name Other
      //! does this element have non-zero displacements or not
      //  bool IsMoving() const { return is_moving_; }

      //@}


     private:
      //! action parameters recognized by bele3
      enum ActionType
      {
        none
      };

      /*!
       * \brief check, whether higher order derivatives for shape functions (dxdx, dxdy, ...) are
       * necessary \return boolean indicating higher order status
       */
      bool isHigherOrderElement(const DRT::Element::DiscretizationType distype) const
      {
        bool hoel = true;
        switch (distype)
        {
          case line3:
            hoel = true;
            break;
          case line2:
            hoel = false;
            break;
          default:
            dserror("distype unknown!");
        }
        return hoel;
      };

      //! don't want = operator
      Bele2& operator=(const Bele2& old);


    };  // class Bele2



  }  // namespace ELEMENTS
}  // namespace DRT


#endif
