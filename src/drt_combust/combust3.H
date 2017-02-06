/*!----------------------------------------------------------------------*
\file combust3.H

\brief wrapper for the combust3 element

  This file contains the element specific service routines like
  Pack, Unpack, NumDofPerNode etc.

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*----------------------------------------------------------------------*/

#ifndef COMBUST3_H
#define COMBUST3_H


#include "../drt_lib/drt_elementtype.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_inpar/inpar_combust.H"
#include "../drt_xfem/physics.H"


namespace COMBUST
{
  class InterfaceHandleCombust;
}

namespace LINALG
{
  class SerialDenseVector;
}

namespace XFEM
{
  class ElementDofManager;
}

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {
    // forward declarations
//    class Combust3Surface;
//    class Combust3Line;
  class Combust3IntFace;

class Combust3Type : public DRT::ElementType
{
public:

  std::string Name() const { return "Combust3Type"; }

  static Combust3Type & Instance();

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const std::string eletype,
                                             const std::string eledistype,
                                             const int id,
                                             const int owner );

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np );

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns );

  virtual void SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions );

private:

  static Combust3Type instance_;
};

    //! A combust3 element with combustion and XFEM
    class Combust3 : public DRT::Element
    {
      public:

        //! @name Friends
//        friend class Combust3Surface;
//        friend class Combust3Line;
        friend class Combust3IntFace;
        //@}

        class MyState
        {
        public:
          MyState(
              const DRT::Discretization&             discretization,
              const std::vector<int>&                lm,
              const bool                             instationary,
              const bool                             genalpha,
              const bool                             gradphi,
              const bool                             gradphi2,
              const bool                             avm3,
              const DRT::ELEMENTS::Combust3*         ele,
              const Epetra_Vector*                   phinp,
              const Epetra_MultiVector*              gradphinp,
              const Epetra_MultiVector*              gradphi2np,
              const Epetra_Vector*                   curvature
              );

          const bool instationary_;
          const bool genalpha_;
          const bool gradphi_;
          const bool gradphi2_;
          const bool avm3_;
          std::vector<double> velnp_;     ///< velocity/pressure unknowns at n+1
          std::vector<double> veln_;      ///< velocity/pressure unknowns at n
          std::vector<double> velnm_;     ///< velocity/pressure unknowns at n-1
          std::vector<double> velaf_;     ///< velocity/pressure unknowns at alpha+f for generalized alpha scheme
          std::vector<double> accn_;      ///< time derivative of unknowns at n
          std::vector<double> accam_;     ///< time derivative of unknowns at alpha+m for generalized alpha scheme
          std::vector<double> fsvelaf_;   ///< fine-scale velocity at n+1
          std::vector<double> phinp_;     ///< G-function unknowns at n+1
          std::vector<double> gradphinp_; ///< smoothed G-function gradient at n+1
          std::vector<double> gradphi2np_; ///< smoothed G-function second derivative at n+1
          std::vector<double> curv_;      ///< curvature based on G-function gradient at n+1
        };

        class MyStateSurface
        {
        public:
          MyStateSurface(
              const DRT::Discretization&             discretization,
              const std::vector<int>&                lm,
              const bool                             instationary,
              const bool                             genalpha,
              const bool                             gradphi,
              const DRT::ELEMENTS::Combust3*         ele,
              const Epetra_Vector*                   phinp
              );

          const bool instationary_;
          const bool genalpha_;
          std::vector<double> velnp_;     ///< velocity/pressure unknowns at n+1
          std::vector<double> velaf_;     ///< velocity/pressure unknowns at alpha+f for generalized alpha scheme
          std::vector<double> phinp_;     ///< G-function unknowns at n+1
        };

        //! @name Constructors and destructors and related methods

        //! Standard Constructor
        explicit Combust3(
            int id,    ///< A unique global id
            int owner  ///< processor id where the element is located
        );

        /*!
         * \brief Copy Constructor
         *
         *  Makes a deep copy of a Element
         */
        Combust3(const Combust3& old);

        DRT::Element* Clone() const;

        virtual DiscretizationType Shape() const;

        virtual int NumLine() const
        {if (NumNode()==8 || NumNode()==20 || NumNode()==27) return 12;
        else if (NumNode()==4 || NumNode()==10) return 6;
        else {
          dserror("Could not determine number of lines");
          return -1;
        }
        }

        virtual int NumSurface() const
        {if (NumNode()==8 || NumNode()==20 || NumNode()==27) return 6;
        else if (NumNode()==4 || NumNode()==10) return 4;
        else if (NumNode()==6 || NumNode()==15 || NumNode()==5) return 5;
        else {
          dserror("Could not determine number of surfaces");
          return -1;
        }
        }

        virtual int NumVolume() const {return 1;}

        virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

        virtual std::vector<Teuchos::RCP<DRT::Element> > Surfaces();

        virtual std::vector<Teuchos::RCP<DRT::Element> > Volumes();

        /*!
        \brief Get Teuchos::RCP to the internal face adjacent to this element as master element and the parent_slave element
        */
        virtual Teuchos::RCP<DRT::Element> CreateFaceElement(
                                                      DRT::Element* parent_slave,          //!< parent slave fluid3 element
                                                      int nnode,                           //!< number of surface nodes
                                                      const int* nodeids,                  //!< node ids of surface element
                                                      DRT::Node** nodes,                   //!< nodes of surface element
                                                      const int lsurface_master,           //!< local surface number w.r.t master parent element
                                                      const int lsurface_slave,            //!< local surface number w.r.t slave parent element
                                                      const std::vector<int>& localtrafomap//! local trafo map
            );

        virtual int UniqueParObjectId() const { return Combust3Type::Instance().UniqueParObjectId(); }

        virtual void Pack(DRT::PackBuffer& data) const;

        virtual void Unpack(const std::vector<char>& data);

        virtual bool Intersected() const { return intersected_; };

        virtual bool Splited() const { return splited_; };

        virtual bool Bisected() const { return bisected_; };

        virtual COMBUST::InterfaceHandleCombust* InterfaceHandle() const { return ih_; };

        virtual Epetra_Vector* Phinp() const { return epetra_phinp_;}

        //! Destructor
        virtual ~Combust3();

        //@}

        //! @name Access methods

        virtual int NumDofPerNode(const DRT::Node& node) const;

        virtual int NumDofPerElement() const;

        Teuchos::RCP<XFEM::ElementDofManager> GetEleDofManager() const { return eleDofManager_; }

        void Print(std::ostream& os) const;

        virtual DRT::ElementType & ElementType() const
        { return Combust3Type::Instance(); }

        //@}

        //! @name Input and Creation

        bool ReadElement(const std::string& eletype,
                         const std::string& distype,
                         DRT::INPUT::LineDefinition* linedef);

        //@}

        //! @name Evaluation

        virtual int Evaluate(Teuchos::ParameterList& params,
            DRT::Discretization&      discretization,
            std::vector<int>&         lm,
            Epetra_SerialDenseMatrix& elemat1,
            Epetra_SerialDenseMatrix& elemat2,
            Epetra_SerialDenseVector& elevec1,
            Epetra_SerialDenseVector& elevec2,
            Epetra_SerialDenseVector& elevec3);

        virtual int EvaluateNeumann(Teuchos::ParameterList& params,
            DRT::Discretization&      discretization,
            DRT::Condition&           condition,
            std::vector<int>&              lm,
            Epetra_SerialDenseVector& elevec1,
            Epetra_SerialDenseMatrix* elemat1 = NULL);

        //@}

        //! @name Other

      private:

        //! action parameters recognized by combust3
        enum ActionType
        {
          none,
          calc_fluid_systemmat_and_residual,
          calc_fluid_stationary_systemmat_and_residual,
          calc_fluid_beltrami_error,
          calc_nitsche_error,
          calc_void_fraction,
          calc_fluid_box_filter,
          calc_smagorinsky_const,
          store_xfem_info,
          get_density,
          integrate_shapefunction,
          reset,
          set_standard_mode
        };

        //! store information about the unknowns of this element
        //!
        //! this dofman will answer to the outside world, how many dofs this element has
        Teuchos::RCP<XFEM::ElementDofManager> eleDofManager_;

        //! if element is in 'standard mode', it will answer without asking the elementdofmanager
        //! -> it will pretend to have 4 dofs per node and no element dofs
        //! needed for postprocessing in the filters (output) and for some initialization procedures
        bool standard_mode_;

        //! store information about the interface handle
        COMBUST::InterfaceHandleCombust* ih_;

        //! store the current phi-vector
        Epetra_Vector* epetra_phinp_;

        //! store the smoothed phi-gradient
        Epetra_MultiVector* gradphi_;
        //! store the smoothed phi2-gradient
        Epetra_MultiVector* gradphi2_;
        //! store the curvature
        Epetra_Vector* curvature_;

        //! boolean indicating whether element has boundary cells
        bool intersected_;
        //! boolean indicating whether element has more than one domain integration cell
        bool splited_;
        //! boolean indicating whether element is completely enriched
        bool bisected_;

        // internal calculation methods

        //! don't want = operator
        Combust3& operator = (const Combust3& old);

        //! evaluate error for beltrami test problem
        void f3_int_beltrami_err(
            std::vector<double>&       myvelnp,
            std::vector<double>&       myprenp,
            Teuchos::RCP<const MAT::Material> material,
            Teuchos::ParameterList&    params
        );

        /*!
          \brief calculate spatial mean void fraction for channel flow
          (requires wall parallel layers of elements)

                                                       rasthofer 06/11
                                                       DA wichmann
          <pre>

          This method calculates the volume fraction per element
          plane in the homogeneous directions.
          All results are communicated via the parameter list!

          </pre>

         */
        void calc_volume_fraction(
          DRT::Discretization&       discretization,
          Teuchos::ParameterList&    params
          );

        /*!
         * \brief get optimal gaussrule for given discretization type
         * \return enum with optimal gausrule for this element
         */
        DRT::UTILS::GaussRule3D getOptimalGaussrule(
            const DiscretizationType&  distype  ///< switch gausrule based on this discretization type
        );

        ActionType convertStringToActionType(const std::string& action) const;

    }; // class Combust3


class Combust3SurfaceType : public DRT::ElementType
{
public:

  std::string Name() const { return "Combust3SurfaceType"; }

  static Combust3SurfaceType & Instance();

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

private:

  static Combust3SurfaceType instance_;
};

    /*!
     * \brief An element representing a surface of a combust3 element
     *
     * \author rasthofer (rasthofer@lnm.mw.tum.de)
     */
    class Combust3Surface : public DRT::FaceElement
    {
      public:

        //friend class Combust3Line;

        //! @name Constructors and destructors and related methods

        //! Standard Constructor
        explicit Combust3Surface(
            int id,                         ///< A unique global id
            int owner,                      ///< Processor owning this surface
            const int nnode,                ///< Number of nodes attached to this element
            const int* nodeids,             ///< global ids of nodes attached to this element
            DRT::Node** nodes,              ///< the discretizations map of nodes to build ptrs to nodes from
            DRT::ELEMENTS::Combust3* parent, ///< The parent fluid element of this surface
            const int lsurface              ///< the local surface number of this surface w.r.t. the parent element
            );

        /*!
         * \brief Copy Constructor
         *
         * Makes a deep copy of a Element
         */
        explicit Combust3Surface(const Combust3Surface& old);

        DRT::Element* Clone() const;

        virtual DiscretizationType Shape() const;

        virtual int NumLine() const
        {if (NumNode()==4 || NumNode()==8 || NumNode()==9) return 4;
        else if (NumNode()==3 || NumNode()==6) return 3;
        else
        {
          dserror("Could not determine number of surfaces");
          return -1;
        }
        }

        virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

        virtual int UniqueParObjectId() const { return Combust3SurfaceType::Instance().UniqueParObjectId(); }

        virtual void Pack(std::vector<char>& data) const;

        virtual void Unpack(const std::vector<char>& data);

        //! Destructor
        virtual ~Combust3Surface();

        //@}

        //! @name Access methods

        virtual int NumDofPerNode(const DRT::Node& node) const
        {return ParentElement()->NumDofPerNode(node);}

        virtual int NumDofPerElement() const { return 0; }


        /*!
         * \brief Return pointer to the parent element
         */
        virtual DRT::ELEMENTS::Combust3* ParentElement() const
        {
          DRT::Element* parent = this->DRT::FaceElement::ParentElement();
          // make sure the static cast below is really valid
          dsassert(dynamic_cast<DRT::ELEMENTS::Combust3*>(parent) != NULL, "Parent element is no combust element");
          return static_cast<DRT::ELEMENTS::Combust3*>(parent);
        }


        void Print(std::ostream& os) const;

        virtual DRT::ElementType & ElementType() const
        { return Combust3SurfaceType::Instance(); }

        //@}

        //! @name Evaluation

        int Evaluate(Teuchos::ParameterList& params,
            DRT::Discretization&      discretization,
            std::vector<int>&         lm,
            Epetra_SerialDenseMatrix& elemat1,
            Epetra_SerialDenseMatrix& elemat2,
            Epetra_SerialDenseVector& elevec1,
            Epetra_SerialDenseVector& elevec2,
            Epetra_SerialDenseVector& elevec3);

        virtual int EvaluateNeumann(Teuchos::ParameterList& params,
            DRT::Discretization&      discretization,
            DRT::Condition&           condition,
            std::vector<int>&         lm,
            Epetra_SerialDenseVector& elevec1,
            Epetra_SerialDenseMatrix* elemat1 = NULL);

        //@}

        /*!
        \brief Return the location vector of this element

        The method computes degrees of freedom this element adresses.
        Degree of freedom ordering is as follows:<br>
        First all degrees of freedom of adjacent nodes are numbered in
        local nodal order, then the element internal degrees of freedom are
        given if present.<br>
        If a derived element has to use a different ordering scheme,
        it is welcome to overload this method as the assembly routines actually
        don't care as long as matrices and vectors evaluated by the element
        match the ordering, which is implicitly assumed.<br>
        Length of the output vector matches number of degrees of freedom
        exactly.<br>
        This version is intended to fill the LocationArray with the dofs
        the element will assemble into. In the standard case these dofs are
        the dofs of the element itself. For some special conditions (e.g.
        the weak dirichlet boundary condtion) a surface element will assemble
        into the dofs of a volume element.<br>

        \note The degrees of freedom returned are not neccessarily only nodal dofs.
              Depending on the element implementation, output might also include
              element dofs.

        \param dis (in)      : the discretization this element belongs to
        \param la (out)      : location data for all dofsets of the discretization
        \param doDirichlet (in): whether to get the Dirichlet flags
        \param condstring (in): Name of condition to be evaluated
        \param condstring (in):  List of parameters for use at element level
        */
        virtual void LocationVector(const Discretization&   dis,
                                    LocationArray&          la,
                                    bool                    doDirichlet,
                                    const std::string&      condstring,
                                    Teuchos::ParameterList& params
                                    ) const;

      private:

        //! action parameters recognized by FluidSurface
        enum ActionType
        {
          none,
          integrate_Shapefunction,
          calc_flux,
          calc_Neumann_inflow
        };

        //! don't want = operator
        Combust3Surface& operator = (const Combust3Surface& old);

        //!  compute kovariant metric tensor G for fluid surface element
        void ComputeMetricTensorForSurface(
            const int                       numnode,
            const Epetra_SerialDenseMatrix& xyze,
            const Epetra_SerialDenseMatrix& deriv,
            LINALG::Matrix<2,2>&            metrictensor,
            double&                         detmetric
            ) const;

        //! integrate shapefunction over surface element
        virtual void IntegrateShapeFunction(
            Teuchos::ParameterList&          params,          ///< ParameterList for communication between control routine and elements
            const DRT::Discretization&       discretization,  ///< vector to be filled by element
            const std::vector<int>&          lm,              ///< location vector of this element
            Epetra_SerialDenseVector&        elevec1,         ///< vector to be filled by element
            const std::vector<double>&       edispnp          ///< current displacement of fluid mesh
            );

        //! integrate shapefunction over surface element
        void IntegrateSurfaceFlow(
            Teuchos::ParameterList&          params,          ///< ParameterList for communication between control routine and elements
            const DRT::Discretization&       discretization,  ///< vector to be filled by element
            const std::vector<int>&          lm,              ///< location vector of this element
            Epetra_SerialDenseVector&        elevec1,         ///< vector to be filled by element
            const std::vector<double>&       edispnp          ///< current displacement of fluid mesh
            );

    }; // class FluidSurface

class Combust3LineType : public DRT::ElementType
{
public:

  std::string Name() const { return "Combust3LineType"; }

  static Combust3LineType & Instance();

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

private:

  static Combust3LineType instance_;
};

    /*!
     * \brief An element representing a line of a combust3 element
     *
     * \author gerstenberger (gerstenberger@lnm.mw.tum.de)
     */
    class Combust3Line : public DRT::FaceElement
    {
      public:

        //! @name Constructors and destructors and related methods

        /*!
  \brief Standard Constructor

  \param id : A unique global id
  \param owner: Processor owning this line
  \param nnode: Number of nodes attached to this element
  \param nodeids: global ids of nodes attached to this element
  \param nodes: the discretizations map of nodes to build ptrs to nodes from
  \param parent: The parent fluid element of this line
  \param lline: the local line number of this line w.r.t. the parent element
         */

        explicit Combust3Line(
            const int id,
            const int owner,
            const int nnode,
            const int* nodeids,
            DRT::Node** nodes,
            DRT::Element* parent,
            const int lline);

        /*!
         * \brief Copy Constructor
         *
         * Makes a deep copy of a Element
         */
        explicit Combust3Line(const Combust3Line& old);

        DRT::Element* Clone() const;

        virtual int UniqueParObjectId() const { return Combust3LineType::Instance().UniqueParObjectId(); }

        virtual void Pack(std::vector<char>& data) const;

        virtual void Unpack(const std::vector<char>& data);

        //! Destructor
        virtual ~Combust3Line();

        //@}

        //! @name Access methods

        virtual DiscretizationType Shape() const;

        virtual int NumDofPerNode(const DRT::Node& node) const
        {return ParentElement()->NumDofPerNode(node);}

        virtual int NumDofPerElement() const { return 0; }

        void Print(std::ostream& os) const;

        virtual DRT::ElementType & ElementType() const
        { return Combust3LineType::Instance(); }

        //! @name Evaluation
        virtual int Evaluate(
            Teuchos::ParameterList&   params,
            DRT::Discretization&      discretization,
            std::vector<int>&   lm,
            Epetra_SerialDenseMatrix& elemat1,
            Epetra_SerialDenseMatrix& elemat2,
            Epetra_SerialDenseVector& elevec1,
            Epetra_SerialDenseVector& elevec2,
            Epetra_SerialDenseVector& elevec3);

        virtual int EvaluateNeumann(
            Teuchos::ParameterList&    params,
            DRT::Discretization&       discretization,
            DRT::Condition&            condition,
            std::vector<int>&          lm,
            Epetra_SerialDenseVector&  elevec1,
            Epetra_SerialDenseMatrix* elemat1 = NULL
            );

        //@}

        /*!
        \brief Return the location vector of this element

        The method computes degrees of freedom this element adresses.
        Degree of freedom ordering is as follows:<br>
        First all degrees of freedom of adjacent nodes are numbered in
        local nodal order, then the element internal degrees of freedom are
        given if present.<br>
        If a derived element has to use a different ordering scheme,
        it is welcome to overload this method as the assembly routines actually
        don't care as long as matrices and vectors evaluated by the element
        match the ordering, which is implicitly assumed.<br>
        Length of the output vector matches number of degrees of freedom
        exactly.<br>
        This version is intended to fill the LocationArray with the dofs
        the element will assemble into. In the standard case these dofs are
        the dofs of the element itself. For some special conditions (e.g.
        the weak dirichlet boundary condtion) a surface element will assemble
        into the dofs of a volume element.<br>

        \note The degrees of freedom returned are not neccessarily only nodal dofs.
              Depending on the element implementation, output might also include
              element dofs.

        \param dis (in)      : the discretization this element belongs to
        \param la (out)      : location data for all dofsets of the discretization
        \param doDirichlet (in): whether to get the Dirichlet flags
        \param condstring (in): Name of condition to be evaluated
        \param condstring (in):  List of parameters for use at element level
        */
        virtual void LocationVector(
            const Discretization&   dis,
            LocationArray&          la,
            bool                    doDirichlet,
            const std::string&      condstring,
            Teuchos::ParameterList& params
            ) const;

      private:

        //! action parameters recognized by FluidLine
        enum ActionType
        {
          none
        };

        // don't want = operator
        Combust3Line& operator = (const Combust3Line& old);

    }; // class Combust3Line


    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================


    /*!
    \brief An element representing an internal face element between two fluid elements

    \note It's only purpose is to evaluate edge based stabilizations for XFEM.
    */
    class Combust3IntFaceType : public DRT::ElementType
    {
    public:

      std::string Name() const { return "Combust3IntFaceType"; }

      static Combust3IntFaceType & Instance();

      virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

      virtual void NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

      virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

    private:

      static Combust3IntFaceType instance_;
    };


    // class Combust3IntFace

    class Combust3IntFace : public DRT::FaceElement
    {
    public:

      //! @name Constructors and destructors and related methods

      //! number of space dimensions
      /*!
      \brief Standard Constructor

      \param id: A unique global id
      \param owner: Processor owning this surface
      \param nnode: Number of nodes attached to this element
      \param nodeids: global ids of nodes attached to this element
      \param nodes: the discretizations map of nodes to build ptrs to nodes from
      \param master_parent: The master parent fluid element of this surface
      \param slave_parent: The slave parent fluid element of this surface
      \param lsurface_master: the local surface number of this surface w.r.t. the master parent element
      \param lsurface_slave: the local surface number of this surface w.r.t. the slave parent element
      \param localtrafomap: transformation map between the local coordinate systems of the face w.r.t the master parent element's face's coordinate system and the slave element's face's coordinate system
      */
      Combust3IntFace( int id,
                     int owner,
                     int nnode,
                     const int* nodeids,
                     DRT::Node** nodes,
                     DRT::ELEMENTS::Combust3* parent_master,
                     DRT::ELEMENTS::Combust3* parent_slave,
                     const int lsurface_master,
                     const int lsurface_slave,
                     const std::vector<int> localtrafomap
                     );

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element
      */
      Combust3IntFace(const Combust3IntFace& old);

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
      \brief Return number of lines of this element
      */
      virtual int NumLine() const
      {
        return DRT::UTILS::getNumberOfElementLines(Shape());
      }

      /*!
      \brief Return number of surfaces of this element
      */
      virtual int NumSurface() const
      {
        return DRT::UTILS::getNumberOfElementSurfaces(Shape());
      }

      /*!
      \brief Get vector of Teuchos::RCPs to the lines of this element
      */
      virtual std::vector<Teuchos::RCP<DRT::Element> > Lines();

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the drt_parobject.H file.
      */
      virtual std::vector<Teuchos::RCP<DRT::Element> > Surfaces();

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of the drt_parobject.H file.
      */
      virtual int UniqueParObjectId() const { return Combust3IntFaceType::Instance().UniqueParObjectId(); }

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
      virtual ~Combust3IntFace();

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
        return std::max(ParentMasterElement()->NumDofPerNode(node),ParentSlaveElement()->NumDofPerNode(node));
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
      \brief create the location vector for patch of master and slave element

      \note All dofs shared by master and slave element are contained only once. Dofs from interface nodes are also included.
      */
      void PatchLocationVector(
          DRT::Discretization & discretization,       ///< discretization
          std::vector<int>&     master_lm,            ///< local map for gdof ids for master element
          std::vector<int>&     slave_lm,             ///< local map for gdof ids for slave element
          std::map<XFEM::PHYSICS::Field,std::vector<int> >&  lm_masterDofPerFieldToPatch, ///< local map between master nodes and nodes in patch
          std::map<XFEM::PHYSICS::Field,std::vector<int> >&  lm_slaveDofPerFieldToPatch,  ///< local map between slave nodes and nodes in patch
          std::map<XFEM::PHYSICS::Field,std::vector<int> >&  patch_components_lm,     ///< rearranged local map for gdof ids for patch of elements
          std::map<XFEM::PHYSICS::Field,std::vector<int> >&  patch_components_lmowner ///< owner of patch
          );

      /*!
      \brief Print this element
      */
      void Print(std::ostream& os) const;

      virtual DRT::ElementType & ElementType() const
      { return Combust3IntFaceType::Instance(); }

      //@}

      //! @name Evaluation

      /*!
      \brief Evaluate element

      \param params (in/out): ParameterList for communication between control routine
                              and elements
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
      int Evaluate(Teuchos::ParameterList&   params,
                   DRT::Discretization&      discretization,
                   std::vector<int>&         lm,
                   Epetra_SerialDenseMatrix& elemat1,
                   Epetra_SerialDenseMatrix& elemat2,
                   Epetra_SerialDenseVector& elevec1,
                   Epetra_SerialDenseVector& elevec2,
                   Epetra_SerialDenseVector& elevec3);

      int Evaluate(Teuchos::ParameterList&            params,               ///< parameter list
                   DRT::Discretization&               discretization,       ///< discretization
                   std::vector<int>&                  lm_master,            ///< master local map
                   std::vector<int>&                  lm_slave,             ///< slave local map
                   std::map<XFEM::PHYSICS::Field,std::vector<int> > &   lm_masterDofPerFieldToPatch, ///< local map between master nodes and nodes in patch
                   std::map<XFEM::PHYSICS::Field,std::vector<int> > &   lm_slaveDofPerFieldToPatch,  ///< local map between slave nodes and nodes in patch
                   std::map<std::pair<XFEM::PHYSICS::Field,XFEM::PHYSICS::Field>,Epetra_SerialDenseMatrix>&  elemat_blocks,   ///< element matrix blocks
                   std::map<XFEM::PHYSICS::Field,Epetra_SerialDenseVector>&  elevec_blocks);   ///< element vector blocks

      //@}

      //! @name Evaluate methods

      /*!
      \brief Evaluate Neumann boundary condition

      \param params (in/out)    : ParameterList for communication between control routine
                                  and elements
      \param discretization (in): reference to the underlying discretization
      \param condition (in)     : condition to be evaluated
      \param lm (in)            : location vector of this element
      \param elevec1 (out)      : vector to be filled by element. If NULL on input,

      \return 0 if successful, negative otherwise
      */
      virtual int EvaluateNeumann(Teuchos::ParameterList&   params,
                                  DRT::Discretization&      discretization,
                                  DRT::Condition&           condition,
                                  std::vector<int>&         lm,
                                  Epetra_SerialDenseVector& elevec1,
                                  Epetra_SerialDenseMatrix* elemat1 = NULL);

      /*!
      \brief return the master parent fluid element
      */
      DRT::ELEMENTS::Combust3* ParentMasterElement() const
      {
        DRT::Element* parent = this->DRT::FaceElement::ParentMasterElement();
        // make sure the static cast below is really valid
        dsassert(dynamic_cast<DRT::ELEMENTS::Combust3*>(parent) != NULL,
            "Master element is no combust element");
        return static_cast<DRT::ELEMENTS::Combust3*>(parent);
      }

      /*!
      \brief return the slave parent fluid element
      */
      DRT::ELEMENTS::Combust3* ParentSlaveElement() const
      {
        DRT::Element* parent = this->DRT::FaceElement::ParentSlaveElement();
        // make sure the static cast below is really valid
        dsassert(dynamic_cast<DRT::ELEMENTS::Combust3*>(parent) != NULL,
            "Master element is no combust element");
        return static_cast<DRT::ELEMENTS::Combust3*>(parent);
      }

      //@}

    private:

      // don't want = operator
      Combust3IntFace& operator = (const Combust3IntFace& old);

      //! action parameters recognized by combust3
      enum ActionType
      {
        none,
        calc_edge_based_stab_terms
      };

      //! convert string to desired action
      ActionType convertStringToActionType(const std::string& action) const;

    }; // class Combust3IntFace

  } // namespace ELEMENTS
} // namespace DRT

#endif // #ifndef COMBUST3_H
