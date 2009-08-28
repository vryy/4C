/*!
\file xfluid3.H

\brief provides the XFluid3 element and related classes

In addition to that, it contains the interface between element call
and integration point loop (depending on the fluid implementation)
as well as some additional service routines (for the evaluation
of errors, turbulence statistics etc.).

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
 */
#ifdef D_FLUID3
#ifdef CCADISCRET
#ifndef DRT_XFLUID3_H
#define DRT_XFLUID3_H


#include "../drt_lib/drt_elementregister.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_xfem/dof_management_element.H"
#include "../drt_xfem/xfem_enums.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"


namespace DRT
{
  // forward declarations
  class Discretization;

  namespace ELEMENTS
  {

    //! A fluid3 element with XFEM capabilities
    class XFluid3 : public DRT::Element
    {
      public:

        //! this struct stores local vectors extracted from the global unknown arrays
        class MyState
        {
        public:
          MyState(
              const DRT::Discretization&      discret,
              const std::vector<int>&         lm,
              const bool                      instat
              );

          const bool instationary;   ///< whether we need unknowns for transient computation
          std::vector<double> velnp; ///< unknowns at n+1
          std::vector<double> veln;  ///< unknowns at n
          std::vector<double> velnm; ///< unknowns at n-1
          std::vector<double> accn;  ///< time derivative of unknowns at n
        };

        //! @name Constructors and destructors and related methods

        //! Standard Constructor
        explicit XFluid3(
            int id,    ///< A unique global id
            int owner  ///< processor id where the element is located
        );

        /*!
         * \brief Copy Constructor
         *
         *  Makes a deep copy of a Element
         */
        XFluid3(const XFluid3& old);

        DRT::Element* Clone() const;

        virtual DiscretizationType Shape() const;

        virtual int NumLine() const
        {
          return DRT::UTILS::getNumberOfElementLines(Shape());
        }

        virtual int NumSurface() const
        {
          return DRT::UTILS::getNumberOfElementSurfaces(Shape());
        }

        virtual int NumVolume() const {return 1;}

        virtual std::vector<RCP<DRT::Element> > Lines();

        virtual std::vector<RCP<DRT::Element> > Surfaces();

        virtual std::vector<RCP<DRT::Element> > Volumes();

        virtual int UniqueParObjectId() const { return ParObject_XFluid3;}

        virtual void Pack(std::vector<char>& data) const;

        virtual void Unpack(const std::vector<char>& data);

        //! Destructor
        virtual ~XFluid3();

        //@}

        //! @name Access methods

        virtual int NumDofPerNode(const DRT::Node& node) const
        {
          if (output_mode_)
          {
            return 4;
          }
          else
          {
            if (eleDofManager_ != Teuchos::null)
            {
              return eleDofManager_->NumDofPerNode(node.Id());
            }
            else
            {
              dserror("no element dof information available!");
              return 0;
            }
          }
        }

        virtual int NumDofPerElement() const
        {
          if (output_mode_)
          {
            return 0;
          }
          else
          {
            if (eleDofManager_ != Teuchos::null)
            {
              return eleDofManager_->NumElemDof();
            }
            else
            {
              dserror("no element dof information available!");
              return 0;
            }
          }
        }

        void Print(ostream& os) const;

        virtual Teuchos::RCP<DRT::ElementRegister> ElementRegister() const;

        //@}

        //! @name Input and Creation

        //bool ReadElement();

        bool ReadElement(const std::string& eletype,
                         const std::string& distype,
                         DRT::INPUT::LineDefinition* linedef);

        //@}

        //! @name Evaluation

        virtual int Evaluate(ParameterList& params,
            DRT::Discretization&      discretization,
            std::vector<int>&         lm,
            Epetra_SerialDenseMatrix& elemat1,
            Epetra_SerialDenseMatrix& elemat2,
            Epetra_SerialDenseVector& elevec1,
            Epetra_SerialDenseVector& elevec2,
            Epetra_SerialDenseVector& elevec3);

        virtual int EvaluateNeumann(ParameterList& params,
            DRT::Discretization&      discretization,
            DRT::Condition&           condition,
            std::vector<int>&              lm,
            Epetra_SerialDenseVector& elevec1,
            Epetra_SerialDenseMatrix* elemat1 = NULL);

        //@}

      private:

        //! action parameters recognized by fluid3
        enum ActionType
        {
          none,
          calc_fluid_systemmat_and_residual,
          calc_linear_fluid,
          calc_fluid_stationary_systemmat_and_residual,
          calc_fluid_beltrami_error,
          store_xfem_info,
          get_density,
          reset,
          set_output_mode,
          integrate_shape
        };

        //! stabilization actions recognized by fluid3 (genalpha implementation)
        enum StabilisationAction
        {
          stabaction_unspecified,
          subscales_quasistatic,
          subscales_time_dependent,
          inertia_stab_drop,
          inertia_stab_keep,
          pstab_assume_inf_sup_stable,
          pstab_use_pspg,
          convective_stab_none,
          convective_stab_supg,
          viscous_stab_none,
          viscous_stab_gls,
          viscous_stab_gls_only_rhs,
          viscous_stab_usfem,
          viscous_stab_usfem_only_rhs,
          continuity_stab_none,
          continuity_stab_yes,
          continuity_stab_td,
          cross_stress_stab_none,
          cross_stress_stab,
          cross_stress_stab_only_rhs,
          reynolds_stress_stab_none,
          reynolds_stress_stab_only_rhs,
          fssgv_no,
          fssgv_artificial_all,
          fssgv_artificial_small,
          fssgv_Smagorinsky_all,
          fssgv_Smagorinsky_small,
          fssgv_mixed_Smagorinsky_all,
          fssgv_mixed_Smagorinsky_small,
          fssgv_scale_similarity
        };


        //! physical turbulence models recognized by fluid3_genalpha implementation
        enum TurbModelAction
        {
          no_model,
          smagorinsky,
          smagorinsky_with_wall_damping,
          dynamic_smagorinsky
        };

        //! map that converts strings to actions (stabilization)
        static map<string,StabilisationAction>          stabstrtoact_;

        //! store information about the unknowns of this element
        //!
        //! this dofman will answer to the outside world, how many dofs this element has
        Teuchos::RCP<XFEM::ElementDofManager>                         eleDofManager_;

        //! store information about the unknowns of this element - only for intersected elements
        //!
        //! this dofman contains also element dofs, that are not part of the global dofs
        //! it is used for integration only. After integration, condensation on the element level
        //! will remove all dofs that are not in eleDofManager_
        Teuchos::RCP<XFEM::ElementDofManager>                         eleDofManager_uncondensed_;

        //! if no dofmanager is present, we want to answer the number of dof question with 4
        //! for this purpose, one has to set the element output mode -> true
        //! when creating the elementdofmanager, output_mode_ is set to false
        bool                             output_mode_;

        //! store information about the intersection handle
        Teuchos::RCP<XFEM::InterfaceHandleXFSI>                       ih_;

        //! store element stress info and corresponding matrices used during condensation
        class DLMInfo
        {
        public:
          //! Kaa inv of old iteration step
          LINALG::SerialDenseMatrix      oldKaainv_;
          //! Kad of old iteration step
          LINALG::SerialDenseMatrix      oldKad_;
          //! fa inv of old iteration step
          LINALG::SerialDenseVector      oldfa_;

          //! this map stores all element stresses (independent of the xfem label)
          LINALG::SerialDenseVector      stressdofs_;

          //! constructor
          explicit DLMInfo(const int nd, const int na);
        private:
          //! disabled standard constructor
          DLMInfo();
          //! disabled copy constructor
          DLMInfo(const DLMInfo&);
        };

        Teuchos::RCP<DLMInfo>            DLM_info_;

        // internal calculation methods

        //! don't want = operator
        XFluid3& operator = (const XFluid3& old);

        //! evaluate error for beltrami test problem
        void f3_int_beltrami_err(
            std::vector<double>&       myvelnp,
            std::vector<double>&       myprenp,
            Teuchos::RCP<const MAT::Material> material,
            ParameterList&             params
        );

        //! compute stresses from previous newton iteration and store them in dlm_info_
        void UpdateOldDLMAndDLMRHS(
            const DRT::Discretization&      discretization,
            const std::vector<int>&         lm,
            MyState&                        mystate
            ) const;

        //! condense (u,p,sigma) matrix to (u,p) matrix and store this iterations sub-matrices
        void CondenseDLMAndStoreOldIterationStep(
            const Epetra_SerialDenseMatrix& elemat1_uncond,
            const Epetra_SerialDenseVector& elevec1_uncond,
            const Epetra_SerialDenseMatrix& Gds_uncond,
            const Epetra_SerialDenseVector& rhsd_uncond,
            Epetra_SerialDenseMatrix&       elemat1,
            Epetra_SerialDenseVector&       elevec1,
            Epetra_SerialDenseMatrix&       Cuu,
            Epetra_SerialDenseMatrix&       Mud,
            Epetra_SerialDenseMatrix&       Mdu,
            Epetra_SerialDenseMatrix&       Cdd,
            Epetra_SerialDenseVector&       rhsd,
            const bool                      monolithic_FSI
            ) const;

        /*!
         * \brief get optimal gaussrule for given discretization type
         * \return enum with optimal gausrule for this element
         */
        DRT::UTILS::GaussRule3D getOptimalGaussrule(
            const DiscretizationType&  distype  ///< switch gausrule based on this discretization type
        );

        /*!
          Calculate matrix and rhs for stationary problem formulation
          */
        void integrateShapefunction(
                const XFEM::AssemblyType&                       assembly_type,
                const DRT::Element*                             ele,     ///< the element those matrix is calculated
                const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,      ///< connection to the interface handler
                const XFEM::ElementDofManager&                  dofman,  ///< dofmanager of the current element
                Epetra_SerialDenseMatrix&                       estif,   ///< element matrix to calculate
                Epetra_SerialDenseVector&                       eforce   ///< element rhs to calculate
                );

        ActionType convertStringToActionType(const string& action) const;

        //! convert string to stabilization action
        StabilisationAction ConvertStringToStabAction (const string& action) const;

    }; // class XFluid3


    //=======================================================================
    //=======================================================================
    //=======================================================================
    //=======================================================================


    /*!
     * \brief A register for xfluid3 element
     *
     * \author gerstenberger (gerstenberger@lnm.mw.tum.de)
     */
    class XFluid3Register : public DRT::ElementRegister
    {
      public:

        //! @name Constructors and destructors and related methods

        //! Standard Constructor
        explicit XFluid3Register(DRT::Element::ElementType etype);

        //! Copy Constructor
        explicit XFluid3Register(const DRT::ELEMENTS::XFluid3Register& old);

        DRT::ELEMENTS::XFluid3Register* Clone() const;

        virtual int UniqueParObjectId() const { return ParObject_XFluid3Register;}

        virtual void Pack(std::vector<char>& data) const;

        virtual void Unpack(const std::vector<char>& data);

        //! Destructor
        virtual ~XFluid3Register();

        //@}

        //! @name Access methods

        void Print(ostream& os) const;

        //@}

        //! @name Construction

        int Initialize(DRT::Discretization& dis);

        //@}

      private:


    }; // class XFluid3Register



    /*!
     * \brief An element representing a surface of a xfluid3 element
     *
     * \author gerstenberger (gerstenberger@lnm.mw.tum.de)
     */
    class XFluid3Surface : public DRT::Element
    {
      public:

        //! @name Constructors and destructors and related methods

        //! Standard Constructor
        explicit XFluid3Surface(
            int id,                         ///< A unique global id
            int owner,                      ///< Processor owning this surface
            const int nnode,                ///< Number of nodes attached to this element
            const int* nodeids,             ///< global ids of nodes attached to this element
            DRT::Node** nodes,              ///< the discretizations map of nodes to build ptrs to nodes from
            DRT::ELEMENTS::XFluid3* parent, ///< The parent fluid element of this surface
            const int lsurface              ///< the local surface number of this surface w.r.t. the parent element
            );

        /*!
         * \brief Copy Constructor
         *
         * Makes a deep copy of a Element
         */
        explicit XFluid3Surface(const XFluid3Surface& old);

        DRT::Element* Clone() const;

        virtual DiscretizationType Shape() const;

        virtual int NumLine() const
        {
          return DRT::UTILS::getNumberOfElementLines(Shape());
        }

        virtual std::vector<RCP<DRT::Element> > Lines();

        virtual int UniqueParObjectId() const { return ParObject_XFluid3Surface;}

        virtual void Pack(std::vector<char>& data) const;

        virtual void Unpack(const std::vector<char>& data);

        //! Destructor
        virtual ~XFluid3Surface();

        //@}

        //! @name Access methods

        virtual int NumDofPerNode(const DRT::Node& node) const
        {return parent_->NumDofPerNode(node);}

        virtual int NumDofPerElement() const { return 0; }

        void Print(ostream& os) const;

        virtual Teuchos::RCP<DRT::ElementRegister> ElementRegister() const { return null; }

        //@}

        //! @name Evaluation

        int Evaluate(ParameterList& params,
            DRT::Discretization&      discretization,
            std::vector<int>&         lm,
            Epetra_SerialDenseMatrix& elemat1,
            Epetra_SerialDenseMatrix& elemat2,
            Epetra_SerialDenseVector& elevec1,
            Epetra_SerialDenseVector& elevec2,
            Epetra_SerialDenseVector& elevec3);

        virtual int EvaluateNeumann(ParameterList& params,
            DRT::Discretization&      discretization,
            DRT::Condition&           condition,
            std::vector<int>&         lm,
            Epetra_SerialDenseVector& elevec1,
            Epetra_SerialDenseMatrix* elemat1 = NULL);

        //@}

      private:

        //! action parameters recognized by Fluid3Surface
        enum ActionType
        {
          none,
          integrate_Shapefunction,
          calc_flow_rate,
          calc_impuls_rate,
          calc_Neumann_inflow
        };

        //! don't want = operator
        XFluid3Surface& operator = (const XFluid3Surface& old);

        //! The parent element of this surface
        DRT::ELEMENTS::XFluid3* parent_;
        //! The local surface number of this surface w.r.t to the parent_ element
        const int                    lsurface_;

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
            ParameterList&                   params,          ///< ParameterList for communication between control routine and elements
            DRT::Discretization&             discretization,  ///< vector to be filled by element
            const std::vector<int>&          lm,              ///< location vector of this element
            Epetra_SerialDenseVector&        elevec1,         ///< vector to be filled by element
            const std::vector<double>&       edispnp          ///< current displacement of fluid mesh
            );

        //! integrate mass flow over surface element
        void IntegrateSurfaceFlowRate(
            ParameterList&                   params,          ///< ParameterList for communication between control routine and elements
            DRT::Discretization&             discretization,  ///< vector to be filled by element
            const std::vector<int>&          lm,              ///< location vector of this element
            Epetra_SerialDenseVector&        elevec1,         ///< vector to be filled by element
            const std::vector<double>&       edispnp          ///< current displacement of fluid mesh
            );

        //! integrate impuls rate over surface element
        void IntegrateSurfaceImpulsRate(
            ParameterList&                   params,          ///< ParameterList for communication between control routine and elements
            DRT::Discretization&             discretization,  ///< vector to be filled by element
            const std::vector<int>&          lm,              ///< location vector of this element
            Epetra_SerialDenseVector&        elevec1,         ///< vector to be filled by element
            const std::vector<double>&       edispnp          ///< current displacement of fluid mesh
            );

        /*!
        \brief apply outflow boundary condition which is necessary for the conservative
        element formulation (since the convective term was partially integrated)

        \param params (in)        : ParameterList for communication between control routine
        \param discretization (in): A reference to the underlying discretization
        \param lm (in)            : location vector of this element
        \param elemat1 (out)      : matrix to be filled by element.
        \param elevec1 (out)      : vector to be filled by element.

        */
        void SurfaceConservativeOutflowConsistency(
          ParameterList&             params,
          DRT::Discretization&       discretization,
          vector<int>&               lm,
          Epetra_SerialDenseMatrix&  elemat1,
          Epetra_SerialDenseVector&  elevec1);

        /*!
        \brief compute potential Neumann inflow

        \param params (in)        : ParameterList for communication between control routine
        \param discretization (in): A reference to the underlying discretization
        \param lm (in)            : location vector of this element
        \param elemat1 (out)      : matrix to be filled by element. If NULL on input,
                                    the controling method does not epxect the element to fill
                                    this matrix.                          and elements
        \param elevec1 (out)      : vector to be filled by element. If NULL on input,

        */
        void NeumannInflow(
          ParameterList&             params,
          DRT::Discretization&       discretization,
          vector<int>&               lm,
          Epetra_SerialDenseMatrix&  elemat1,
          Epetra_SerialDenseVector&  elevec1);

    }; // class Fluid3Surface


    /*!
     * \brief An element representing a line of a xfluid3 element
     *
     * \author gerstenberger (gerstenberger@lnm.mw.tum.de)
     */
    class XFluid3Line : public DRT::Element
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
        explicit XFluid3Line(
            const int id,
            const int owner,
            const int nnode,
            const int* nodeids,
            DRT::Node** nodes,
            DRT::Element* parent,
            const int lline
            );

        /*!
         * \brief Copy Constructor
         *
         * Makes a deep copy of a Element
         */
        explicit XFluid3Line(const XFluid3Line& old);

        DRT::Element* Clone() const;

        virtual int UniqueParObjectId() const { return ParObject_XFluid3Line;}

        virtual void Pack(std::vector<char>& data) const;

        virtual void Unpack(const std::vector<char>& data);

        //! Destructor
        virtual ~XFluid3Line();

        //@}

        //! @name Access methods

        virtual DiscretizationType Shape() const;

        virtual int NumDofPerNode(const DRT::Node& node) const
        {return parent_->NumDofPerNode(node);}

        virtual int NumDofPerElement() const { return 0; }

        void Print(ostream& os) const;

        virtual Teuchos::RCP<DRT::ElementRegister> ElementRegister() const { return null; }

        //! @name Evaluation
        virtual int Evaluate(
            ParameterList& params,
            DRT::Discretization&      discretization,
            std::vector<int>&   lm,
            Epetra_SerialDenseMatrix& elemat1,
            Epetra_SerialDenseMatrix& elemat2,
            Epetra_SerialDenseVector& elevec1,
            Epetra_SerialDenseVector& elevec2,
            Epetra_SerialDenseVector& elevec3);

        virtual int EvaluateNeumann(
            ParameterList&             params,
            DRT::Discretization&       discretization,
            DRT::Condition&            condition,
            std::vector<int>&          lm,
            Epetra_SerialDenseVector&  elevec1,
            Epetra_SerialDenseMatrix* elemat1 = NULL
            );

        //@}

      private:

        //! action parameters recognized by Fluid3Line
        enum ActionType
        {
          none
        };

        //! shape of the element
        DRT::Element::DiscretizationType distype_;
        DiscretizationType MyShape() const;

        // don't want = operator
        XFluid3Line& operator = (const XFluid3Line& old);

        //! The volume parent element of this line
        DRT::Element* parent_;

        //! The local line number of this line w.r.t to the parent_ element
        int                    lline_;

    }; // class XFluid3Line
  } // namespace ELEMENTS
} // namespace DRT




#endif  // #ifndef XFLUID3_H
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
