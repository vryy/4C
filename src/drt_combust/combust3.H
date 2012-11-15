/*!----------------------------------------------------------------------*
\file combust3.H

\brief wrapper for the combust3 element

	This file contains the element specific service routines like
	Pack, Unpack, NumDofPerNode etc.

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/

#ifndef COMBUST3_H
#define COMBUST3_H


#include "../drt_lib/drt_elementtype.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_inpar/inpar_combust.H"


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

class Combust3Type : public DRT::ElementType
{
public:

  std::string Name() const { return "Combust3Type"; }

  static Combust3Type & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

  virtual Teuchos::RCP<DRT::Element> Create( const string eletype,
                                             const string eledistype,
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
              const DRT::ELEMENTS::Combust3*         ele,
              const Epetra_Vector*                   phinp,
              const Epetra_MultiVector*              gradphinp,
              const Epetra_Vector*                   curvature
              );

          const bool instationary_;
          const bool genalpha_;
          const bool gradphi_;
          std::vector<double> velnp_;     ///< velocity/pressure unknowns at n+1
          std::vector<double> veln_;      ///< velocity/pressure unknowns at n
          std::vector<double> velnm_;     ///< velocity/pressure unknowns at n-1
          std::vector<double> velaf_;     ///< velocity/pressure unknowns at alpha+f for generalized alpha scheme
          std::vector<double> accn_;      ///< time derivative of unknowns at n
          std::vector<double> accam_;     ///< time derivative of unknowns at alpha+m for generalized alpha scheme
          std::vector<double> phinp_;     ///< G-function unknowns at n+1
          std::vector<double> gradphinp_; ///< smoothed G-function gradient at n+1
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

        virtual int UniqueParObjectId() const { return Combust3Type::Instance().UniqueParObjectId(); }

        virtual void Pack(DRT::PackBuffer& data) const;

        virtual void Unpack(const std::vector<char>& data);

        virtual bool Bisected() const { return bisected_; };

        virtual bool Touched() const { return touched_; };

        virtual COMBUST::InterfaceHandleCombust* InterfaceHandle() const { return ih_; };

        virtual Epetra_Vector* Phinp() const { return epetra_phinp_;}

        //! Destructor
        virtual ~Combust3();

        //@}

        //! @name Access methods

        virtual int NumDofPerNode(const DRT::Node& node) const;

        virtual int NumDofPerElement() const;

        Teuchos::RCP<XFEM::ElementDofManager> GetEleDofManager(){ return eleDofManager_; }

        void Print(ostream& os) const;

        virtual DRT::ElementType & ElementType() const
        { return Combust3Type::Instance(); }

        //@}

        //! @name Input and Creation

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
          calc_turbulence_statistics,
          calc_fluid_box_filter,
          calc_smagorinsky_const,
          store_xfem_info,
          get_density,
          integrate_shapefunction,
          reset,
          set_standard_mode
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
        static std::map<string,StabilisationAction>          stabstrtoact_;

        //! store information about the unknowns of this element
        //!
        //! this dofman will answer to the outside world, how many dofs this element has
        Teuchos::RCP<XFEM::ElementDofManager> eleDofManager_;

        //! store information about the unknowns of this element - only for intersected elements
        //!
        //! this dofman contains also element dofs, that are not part of the global dofs
        //! it is used for integration only. After integration, condensation on the element level
        //! will remove all dofs that are not in eleDofManager_
        Teuchos::RCP<XFEM::ElementDofManager> eleDofManager_uncondensed_;

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

        //! store the curvature
        Epetra_Vector* curvature_;

        //! boolean indicating whether element is bisected by the interface
        bool bisected_;

        //! boolean indicating whether element is trisected by the interface
        bool trisected_;

        //! boolean indicating whether element is touched by the interface at a whole 2D face
        bool touched_;

        //! store element stress info and corresponding matrices used during condensation
        class DLMInfo
        {
        public:
          //! Kaa inv of old iteration step
          Teuchos::RCP<LINALG::SerialDenseMatrix>      oldKaainv_;
          //! Kad of old iteration step
          Teuchos::RCP<LINALG::SerialDenseMatrix>      oldKad_;
          //! fa inv of old iteration step
          Teuchos::RCP<LINALG::SerialDenseVector>      oldfa_;
          //! this map stores all element stresses (independent of the xfem label)
          Teuchos::RCP<LINALG::SerialDenseVector>      stressdofs_;

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
        Combust3& operator = (const Combust3& old);

        //! evaluate error for beltrami test problem
        void f3_int_beltrami_err(
            std::vector<double>&       myvelnp,
            std::vector<double>&       myprenp,
            Teuchos::RCP<const MAT::Material> material,
            ParameterList&             params
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
          DRT::Discretization&      discretization,
          const vector<double>&     solution,
          ParameterList&            params
          );

        //! compute stresses from previous newton iteration and store them in dlm_info_
        void UpdateOldDLMAndDLMRHS(
            const DRT::Discretization& discretization,
            const std::vector<int>&    lm,
            MyState&                   mystate
        ) const;

        //! condense (u,p,sigma) matrix to (u,p) matrix and store this iterations sub-matrices
        void CondenseElementStressAndStoreOldIterationStep(
            const Epetra_SerialDenseMatrix& elemat1_uncond,
            const Epetra_SerialDenseVector& elevec1_uncond,
            Epetra_SerialDenseMatrix&       elemat1,
            Epetra_SerialDenseVector&       elevec1
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

    }; // class Combust3


class Combust3SurfaceType : public DRT::ElementType
{
public:

  std::string Name() const { return "Combust3SurfaceType"; }

  static Combust3SurfaceType & Instance() { return instance_; };

  virtual Teuchos::RCP<DRT::Element> Create( const int id, const int owner );

  virtual void NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np ) {}

  virtual void ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns ) {}

private:

  static Combust3SurfaceType instance_;
};

    /*!
     * \brief An element representing a surface of a combust3 element
     *
     * \author gerstenberger (gerstenberger@lnm.mw.tum.de)
     */
    class Combust3Surface : public DRT::Element
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
        {return parent_->NumDofPerNode(node);}

        virtual int NumDofPerElement() const { return 0; }

        void Print(ostream& os) const;

        virtual DRT::ElementType & ElementType() const
        { return Combust3SurfaceType::Instance(); }

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

        //! The parent element of this surface
        DRT::ELEMENTS::Combust3* parent_;
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
            const DRT::Discretization&       discretization,  ///< vector to be filled by element
            const std::vector<int>&          lm,              ///< location vector of this element
            Epetra_SerialDenseVector&        elevec1,         ///< vector to be filled by element
            const std::vector<double>&       edispnp          ///< current displacement of fluid mesh
            );

        //! integrate shapefunction over surface element
        void IntegrateSurfaceFlow(
            ParameterList&                   params,          ///< ParameterList for communication between control routine and elements
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

  static Combust3LineType & Instance() { return instance_; };

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
    class Combust3Line : public DRT::Element
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
        {return parent_->NumDofPerNode(node);}

        virtual int NumDofPerElement() const { return 0; }

        void Print(ostream& os) const;

        virtual DRT::ElementType & ElementType() const
        { return Combust3LineType::Instance(); }

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

        //! action parameters recognized by FluidLine
        enum ActionType
        {
          none
        };

        // don't want = operator
        Combust3Line& operator = (const Combust3Line& old);

        //! The volume parent element of this line
        DRT::Element* parent_;

        //! The local line number of this line w.r.t to the parent_ element
        int                    lline_;

    }; // class Combust3Line

  } // namespace ELEMENTS
} // namespace DRT

#endif // #ifndef COMBUST3_H
