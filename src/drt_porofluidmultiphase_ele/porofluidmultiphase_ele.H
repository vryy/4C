/*----------------------------------------------------------------------*/
/*!
 \file porofluidmultiphase_ele.H

 \brief definition of porofluidmultiphase elements

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_H_
#define SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_H_


#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"

namespace DRT
{
  namespace ELEMENTS
  {
    class PoroFluidMultiPhaseType : public DRT::ElementType
    {
     public:
      std::string Name() const { return "PoroFluidMultiPhaseType"; }

      /// singleton instance method
      static PoroFluidMultiPhaseType& Instance();

      /// create an element from data
      virtual DRT::ParObject* Create(const std::vector<char>& data);

      /// create an element from a dat file specifier
      virtual Teuchos::RCP<DRT::Element> Create(
          const std::string eletype, const std::string eledistype, const int id, const int owner);

      /// create an empty element
      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      /// nodal block information to create a null space description
      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

      /// do the null space computation
      virtual void ComputeNullSpace(DRT::Discretization& dis, std::vector<double>& ns,
          const double* x0, int numdf, int dimns);

      /// setup the dat file input line definitions for this type of element
      virtual void SetupElementDefinition(
          std::map<std::string, std::map<std::string, DRT::INPUT::LineDefinition>>& definitions);

      /// initialize element
      virtual int Initialize(DRT::Discretization& dis);

      /// pre-evaluation
      virtual void PreEvaluate(DRT::Discretization& dis, Teuchos::ParameterList& p,
          Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
          Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
          Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
          Teuchos::RCP<Epetra_Vector> systemvector3);

     private:
      /// the actual instance
      static PoroFluidMultiPhaseType instance_;
    };

    /*!
    \brief The PoroFluidMultiPhase element
    */
    class PoroFluidMultiPhase : public DRT::Element
    {
      //! @name Friends
      friend class PoroFluidMultiPhaseType;

     public:
      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor
      */
      PoroFluidMultiPhase(int id,  ///< A unique global id of this element
          int owner                ///< processor id who owns a certain instance of this class
      );

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      PoroFluidMultiPhase(const PoroFluidMultiPhase& old);

      /*!
      \brief Deep copy this instance of PoroFluidMultiPhase and return pointer to the copy

      The Clone() method is used from the virtual base class Element in cases
      where the type of the derived class is unknown and a copy-ctor is needed

      */
      DRT::Element* Clone() const;

      /*!
      \brief Destructor

      */
      virtual ~PoroFluidMultiPhase();

      //@}

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
      virtual int UniqueParObjectId() const
      {
        return PoroFluidMultiPhaseType::Instance().UniqueParObjectId();
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
        if (numdofpernode_ < 1) dserror("NumDofPerNode is < 1");
        return numdofpernode_;
      }

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
      virtual DRT::ElementType& ElementType() const { return PoroFluidMultiPhaseType::Instance(); }

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
      \brief Evaluate an element (multiple dofset version)

      An element derived from this class uses the Evaluate method to receive commands
      and parameters from some control routine in params and evaluates element matrices and
      vectors accoring to the command in params.

      \note This class implements a dummy of this method that prints a warning and
            returns false.

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
      int Evaluate(Teuchos::ParameterList& params, DRT::Discretization& discretization,
          LocationArray& la, Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
          Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
          Epetra_SerialDenseVector& elevec3);

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

     protected:
      virtual void Initialize();

     private:
      //! the element discretization type (shape)
      DRT::Element::DiscretizationType distype_;

      /*!
       * \brief number of dofs per node (for systems of transport equations)
       * (storage neccessary because we dont know the material in the post filters anymore)
       */
      int numdofpernode_;

      //! don't want = operator
      PoroFluidMultiPhase& operator=(const PoroFluidMultiPhase& old);
    };  // class PoroFluidMultiPhase


    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================

    class PoroFluidMultiPhaseBoundaryType : public DRT::ElementType
    {
     public:
      std::string Name() const { return "PoroFluidMultiPhaseBoundaryType"; }

      static PoroFluidMultiPhaseBoundaryType& Instance();

      virtual Teuchos::RCP<DRT::Element> Create(const int id, const int owner);

      virtual void NodalBlockInformation(
          DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
      {
      }

      virtual void ComputeNullSpace(
          DRT::Discretization& dis, std::vector<double>& ns, const double* x0, int numdf, int dimns)
      {
      }

     private:
      static PoroFluidMultiPhaseBoundaryType instance_;
    };


    /*!
    \brief An element representing a boundary element of a PoroFluidMultiPhase element

    \note This is a pure boundary condition element. Its only
          purpose is to evaluate certain boundary conditions that might be
          adjacent to a parent PoroFluidMultiPhase element.
    */
    class PoroFluidMultiPhaseBoundary : public DRT::FaceElement
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
      \param parent: The parent PoroFluidMultiPhase element of this surface
      \param lsurface: the local surface number of this surface w.r.t. the parent element
      */
      PoroFluidMultiPhaseBoundary(int id, int owner, int nnode, const int* nodeids,
          DRT::Node** nodes, DRT::ELEMENTS::PoroFluidMultiPhase* parent, const int lsurface);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element

      */
      PoroFluidMultiPhaseBoundary(const PoroFluidMultiPhaseBoundary& old);

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
        return PoroFluidMultiPhaseBoundaryType::Instance().UniqueParObjectId();
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
      virtual ~PoroFluidMultiPhaseBoundary();

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

      //! Return a pointer to the parent element of this boundary element
      virtual DRT::ELEMENTS::PoroFluidMultiPhase* ParentElement() const
      {
        DRT::Element* parent = DRT::FaceElement::ParentElement();
        // make sure the static cast below is really valid
        dsassert(dynamic_cast<DRT::ELEMENTS::PoroFluidMultiPhase*>(parent) != NULL,
            "Master element is no PoroFluidMultiPhase element");
        return static_cast<DRT::ELEMENTS::PoroFluidMultiPhase*>(parent);
      }


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
      virtual DRT::ElementType& ElementType() const
      {
        return PoroFluidMultiPhaseBoundaryType::Instance();
      }

      //@}

      //! @name Evaluation

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

      //@}

      //! @name Evaluate methods

      /*!
      \brief Evaluate a Neumann boundary condition

      this method evaluates a surface Neumann condition on the PoroFluidMultiPhase element

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
      // don't want = operator
      PoroFluidMultiPhaseBoundary& operator=(const PoroFluidMultiPhaseBoundary& old);

    };  // class PoroFluidMultiPhaseBoundary


  }  // namespace ELEMENTS
}  // namespace DRT


#endif /* SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_H_ */
