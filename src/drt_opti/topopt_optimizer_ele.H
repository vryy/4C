/*----------------------------------------------------------------------*/
/*!

\brief element routines of the topology optimizer

\maintainer Martin Kronbichler

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef TOPOPT_OPTIMIZER_ELE_H_
#define TOPOPT_OPTIMIZER_ELE_H_


#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"


namespace DRT
{
  class ParObject;
  class Element;

  namespace INPUT
  {
    class LineDefinition;
  }

  namespace ELEMENTS
  {
    class TopOptType : public DRT::ElementType
    {
     public:
      std::string Name() const { return "TopOptType"; }

      static TopOptType& Instance();

      virtual DRT::ParObject* Create(const std::vector<char>& data);

      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
      {
        dserror("used?");
      }

      virtual void ComputeNullSpace(
          DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
      {
        dserror("used?");
      }

      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

     private:
      static TopOptType instance_;
    };

    /*!
    \brief The topology optimization element

    This class computes the domain parts of the values required for an optimization step (currently
    for topology optimization for flow problems). These are the objective value, the constraint
    value(s) and potentially their derivatives.

    \author Martin Winklmaier (winklmaier@lnm.mw.tum.de)
    */
    class TopOpt : public DRT::Element
    {
     public:
      //@}
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor
      */
      TopOpt(int id,  ///< A unique global id of this element
          int owner   ///< processor id who owns a certain instance of this class
      );

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      TopOpt(const TopOpt& old);

      /*!
      \brief Deep copy this instance of Transport and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const;

      /// Set element material
      /*!
        Material numbers are read from the input file. The element stores
        a corresponding material object. These material objects can be
        anything from very simple (just a little calculation) to highly
        sophisticated with history data. The material is packed and
        unpacked along with its element.

        \param matnum : material number from input file

        \note reimplementation of this method, due to initialising
              numdofpernode_, since the material is known now.
       */
      virtual void SetMaterial(int matnum);

      /*!
      \brief Get shape type of element
      */
      virtual DiscretizationType Shape() const;

      /*!
      \brief set discretization type of element
      */
      virtual void SetDisType(DRT::Element::DiscretizationType shape)
      {
        distype_ = shape;
        return;
      };

      /*!
      \brief Return number of lines of this element
      */
      virtual int NumLine() const;

      /*!
      \brief Return number of surfaces of this element
      */
      virtual int NumSurface() const;

      /*!
      \brief Return number of volumes of this element
      */
      virtual int NumVolume() const;

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element
      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Lines();

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element
      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Surfaces();

      /*!
      \brief Get vector of Teuchos::RCPs to the volume of this element
      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Volumes();

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of this file.
      */
      virtual int UniqueParObjectId() const { return TopOptType::Instance().UniqueParObjectId(); }

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
      \brief Destructor

      */
      virtual ~TopOpt();

      //@}

      //! @name Access methods


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual DRT::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can re-decide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      virtual int NumDofPerNode(const DRT::Node& node) const { return 1; }

      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual DRT::Element)

      The element decides how many element degrees of freedom it has.
      It can re-decide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      virtual int NumDofPerElement() const { return 0; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      /*!
      \brief Return ElementType
      */
      virtual DRT::ElementType& ElementType() const { return TopOptType::Instance(); }

      //@}

      //! @name Input and Creation

      /*!
      \brief Read input for this element
      */
      bool ReadElement(const std::string& eletype, const std::string& distype,
          DRT::INPUT::LineDefinition* linedef);

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element, i.e. call the implementation to evaluate Transport
      element stiffness, mass, internal forces or evaluate errors, statistics or updates
      etc. directly.

      Following implementations of the element are allowed:
      <pre>
        o Evaluation of system matrix and residual for the One-Step-Theta or BDF2
          method

        o Evaluation of system matrix and residual for the stationary fluid solver
      </pre>

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param discretization (in): A reference to the underlying discretization
      \param lm (in):         location vector of this element
      \param elemat1 (out)  : matrix to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elemat2 (out)  : matrix to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elevec1 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);


      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a surfaces Neumann condition on the shell element

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param condition (in)     : The condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If NULL on input,

      \return 0 if successful, negative otherwise
      */
      virtual int EvaluateNeumann(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1 = NULL);


      //@}

      //@}

     private:
      //! action parameters recognized by fluid3
      enum ActionType
      {
        none,
        set_general_optimization_parameter,
        update_general_optimization_parameter,
        compute_values,
        compute_gradients
      };

      //! the element discretization type (shape)
      DRT::Element::DiscretizationType distype_;


      //! don't want = operator
      TopOpt& operator=(const TopOpt& old);

      /*!
       * \brief convert std::string to action type
       */
      ActionType convertStringToActionType(const std::string& action) const;

    };  // class TopOpt


    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================

    class TopOptBoundaryType : public DRT::ElementType
    {
     public:
      std::string Name() const { return "TopOptBoundaryType"; }

      static TopOptBoundaryType& Instance();

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
      {
        dserror("used?");
      }

      virtual void ComputeNullSpace(
          DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
      {
        dserror("used?");
      }

     private:
      static TopOptBoundaryType instance_;
    };


    /*!
    \brief An element representing a boundary element of a Transport element

    \note This is a pure boundary condition element. Its only
          purpose is to evaluate certain boundary conditions that might be
          adjacent to a parent topology optimization element.

    \author Martin Winklmaier (winklmaier@lnm.mw.tum.de)
    */
    class TopOptBoundary : public DRT::FaceElement
    {
     public:
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id : A unique global id
      \param owner: Processor owning this surface
      \param nnode: Number of nodes attached to this element
      \param nodeids: global ids of nodes attached to this element
      \param nodes: the discretizations map of nodes to build ptrs to nodes from
      \param parent: The parent fluid element of this surface
      \param lsurface: the local surface number of this surface w.r.t. the parent element
      */
      TopOptBoundary(int id, int owner, int nnode, const int* nodeids, DRT::Node** nodes,
          DRT::ELEMENTS::TopOpt* parent, const int lsurface);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      TopOptBoundary(const TopOptBoundary& old);

      /*!
      \brief Deep copy this instance of an element and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const;

      /*!
      \brief Get shape type of element
      */
      virtual DiscretizationType Shape() const;

      /*!
      \brief Return number of lines of boundary element
      */
      virtual int NumLine() const;

      /*!
      \brief Return number of surfaces of boundary element
       */
      virtual int NumSurface() const;

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element

      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Lines();

      /*!
      \brief Get vector of Teuchos::RCPs to the surfaces of this element
      */
      virtual std::vector<Teuchos::RCP<DRT::Element>> Surfaces();

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the drt_parobject.H file.
      */
      virtual int UniqueParObjectId() const
      {
        return TopOptBoundaryType::Instance().UniqueParObjectId();
      }

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
      \brief Destructor

      */
      virtual ~TopOptBoundary();

      //@}

      //! @name Access methods


      /*!
      \brief Get number of degrees of freedom of a certain node
             (implements pure virtual DRT::Element)

      The element decides how many degrees of freedom its nodes must have.
      As this may vary along a simulation, the element can re-decide the
      number of degrees of freedom per node along the way for each of it's nodes
      separately.
      */
      virtual int NumDofPerNode(const DRT::Node& node) const
      {
        return ParentElement()->NumDofPerNode(node);
      }
      /*
        //! Return a pointer to the parent element of this boundary element
        virtual DRT::ELEMENTS::TopOpt* ParentElement()
        {return ParentElement();}

        //! Return a boundary element id of this boundary element
        virtual int BeleNumber()
        {
          return lbeleid_;
        }
        */


      /*!
      \brief Get number of degrees of freedom per element
             (implements pure virtual DRT::Element)

      The element decides how many element degrees of freedom it has.
      It can re-decide along the way of a simulation.

      \note Element degrees of freedom mentioned here are dofs that are visible
            at the level of the total system of equations. Purely internal
            element dofs that are condensed internally should NOT be considered.
      */
      virtual int NumDofPerElement() const { return 0; }

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      /*!
      \brief Return ElementType
      */
      virtual DRT::ElementType& ElementType() const { return TopOptBoundaryType::Instance(); }

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate an element

      Evaluate Transport element stiffness, mass, internal forces etc

      \param params (in/out): ParameterList for communication between control routine
                              and elements
      \param discretization (in): A reference to the underlying discretization
      \param lm (in):         location vector of this element
      \param elemat1 (out)  : matrix to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elemat2 (out)  : matrix to be filled by element. If NULL on input,
                              the controling method does not epxect the element to fill
                              this matrix.
      \param elevec1 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec2 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \param elevec3 (out)  : vector to be filled by element. If NULL on input,
                              the controlling method does not epxect the element
                              to fill this vector
      \return 0 if successful, negative otherwise
      */
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
          Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
          Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3);

      //@}

      //! @name Evaluate methods

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a surface Neumann condition on the Transport element

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): A reference to the underlying discretization
      \param condition (in)     : The condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If NULL on input,

      \return 0 if successful, negative otherwise
      */
      virtual int EvaluateNeumann(Teuchos::ParameterList& params,
          DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1 = NULL);

      //@}

     private:
      //! action parameters recognized by fluid3
      enum ActionType
      {
        none,
        set_general_optimization_parameter,
        update_general_optimization_parameter,
        compute_values,
        compute_gradients
      };

      // don't want = operator
      TopOptBoundary& operator=(const TopOptBoundary& old);

      /*!
       * \brief convert std::string to action type
       */
      ActionType convertStringToActionType(const std::string& action) const;

    };  // class TopOptBoundary

  }  // namespace ELEMENTS
}  // namespace DRT


#endif /* TOPOPT_OPTIMIZER_ELE_ */
