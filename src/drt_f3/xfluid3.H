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
        struct MyState
        {
          bool instationary;
          vector<double> velnp;
          vector<double> veln;
          vector<double> velnm;
          vector<double> accn;
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

        virtual int NumDofPerNode(const DRT::Node& node) const {return eleDofManager_->NumDofPerNode(node.Id(), 4);}

        virtual int NumDofPerElement() const {return eleDofManager_->NumElemDof();}
        
        void Print(ostream& os) const;
        
        virtual Teuchos::RCP<DRT::ElementRegister> ElementRegister() const;
        
        //@}
        
        //! @name Input and Creation
        
        bool ReadElement();
        
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
            Epetra_SerialDenseVector& elevec1);
        
        //@}
        
        //! @name Other
        
        //! do we have a moving fluid mesh or not
        bool IsAle() const { return is_ale_; }
        
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
          get_density
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
        
        //! flag for euler/ale net algorithm
        bool                                            is_ale_;
        //! container containing all kind of fluid3 this and that
        DRT::Container                                  data_;
        
      private:
        
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
        
        //! store information about the intersection handle
        Teuchos::RCP<XFEM::InterfaceHandleXFSI>                       ih_;

        //! store element stress info and corresponding matrices used during condensation
        class DLMInfo
        {
        public:
          //! store matrizes and rhs of the old iteration step
          LINALG::SerialDenseMatrix      oldKaainv_;
          LINALG::SerialDenseMatrix      oldKad_;
          LINALG::SerialDenseVector      oldfa_;
          
          //! this map stores all element stresses (independent of the xfem label)
          LINALG::SerialDenseVector      stressdofs_;
          
          //! constructor
          explicit DLMInfo(const int nd, const int na);
        private:
          //! no standard constructor
          DLMInfo();
          //! no copy constructor
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
            const struct _MATERIAL*    material,
            ParameterList&             params
        );
        
        //! compute stresses from previous newton iteration and store them in dlm_info_
        void UpdateOldDLMAndDLMRHS(
            const DRT::Discretization&      discretization,
            const std::vector<int>&         lm
            ) const;
        
        //! condense (u,p,sigma) matrix to (u,p) matrix and store this iterations sub-matrices
        void CondenseDLMAndStoreOldIterationStep(
            const Epetra_SerialDenseMatrix& elemat1_uncond,
            const Epetra_SerialDenseVector& elevec1_uncond,
            Epetra_SerialDenseMatrix& elemat1,
            Epetra_SerialDenseVector& elevec1
            ) const;
        
        /*!
         * \brief get optimal gaussrule for given discretization type
         * \return enum with optimal gausrule for this element
         */
        DRT::UTILS::GaussRule3D getOptimalGaussrule(
            const DiscretizationType&  distype  ///< switch gausrule based on this discretization type
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
        {if (NumNode()==4 || NumNode()==8 || NumNode()==9) return 4;
        else if (NumNode()==3 || NumNode()==6) return 3;
        else
        {
          dserror("Could not determine number of surfaces");
          return -1;
        }
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
            Epetra_SerialDenseVector& elevec1);
        
        //@}
        
      private:
        
        //! action parameters recognized by Fluid3Surface
        enum ActionType
        {
          none,
          integrate_Shapefunction,
          calc_flux
        };
        
        //! don't want = operator
        XFluid3Surface& operator = (const XFluid3Surface& old);
        
        //! The parent element of this surface
        DRT::ELEMENTS::XFluid3* parent_;
        //! The local surface number of this surface w.r.t to the parent_ element
        const int                    lsurface_;

        //!  compute kovariant metric tensor G for fluid surface element
        void f3_metric_tensor_for_surface(
            const Epetra_SerialDenseMatrix& xyze,
            const Epetra_SerialDenseMatrix& deriv,
            Epetra_SerialDenseMatrix&       metrictensor,
            double&                         drs);
        
        //! integrate shapefunction over surface element
        virtual void IntegrateShapeFunction(
            ParameterList&                   params,          ///< ParameterList for communication between control routine and elements
            DRT::Discretization&             discretization,  ///< vector to be filled by element
            const std::vector<int>&          lm,              ///< location vector of this element
            Epetra_SerialDenseVector&        elevec1,         ///< vector to be filled by element
            const std::vector<double>&       edispnp          ///< current displacement of fluid mesh
            );
        
        //! integrate shapefunction over surface element
        void IntegrateSurfaceFlow(
            ParameterList&                   params,          ///< ParameterList for communication between control routine and elements
            DRT::Discretization&             discretization,  ///< vector to be filled by element
            const std::vector<int>&          lm,              ///< location vector of this element
            Epetra_SerialDenseVector&        elevec1,         ///< vector to be filled by element
            const std::vector<double>&       edispnp          ///< current displacement of fluid mesh
            );
        
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
            Epetra_SerialDenseVector&  elevec1
            );
        
        //@}
        
      private:
        
        //! action parameters recognized by Fluid3Line
        enum ActionType
        {
          none
        };
        
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
