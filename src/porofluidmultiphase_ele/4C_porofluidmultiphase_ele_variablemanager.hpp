/*----------------------------------------------------------------------*/
/*! \file
 \brief variable manager class for poro multiphase fluid element

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROFLUIDMULTIPHASE_ELE_VARIABLEMANAGER_HPP
#define FOUR_C_POROFLUIDMULTIPHASE_ELE_VARIABLEMANAGER_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_element.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_porofluidmultiphase_ele_action.hpp"

#include <Teuchos_RCP.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace MAT
{
  class Material;
}
namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    class PoroFluidMultiPhaseEleParameter;

    namespace POROFLUIDMANAGER
    {
      template <int, int>
      class VariableManagerInterface;

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief general interface to variable manager for minimal access (non template)

      These classes manage all accesses to primary variables
       (in contrast to the phase manager, which manages values of secondary).

      \note The 'phi' vector denotes the 'generic' primary variable, which
            can be a pressure, saturation or a differential pressure. This
            does not (and needs not) to differentiate between the actual choice
            of primary variables (this is the responsibility of the phase manager).

      This class is a pure virtual non-template class, which provides access to
      the generic state phi and the scalar state (coming from a ScaTra coupling).

      \author vuong
      */
      class VariableManagerMinAccess
      {
       public:
        //! constructor
        VariableManagerMinAccess(){};

        //! destructor
        virtual ~VariableManagerMinAccess() = default;

        virtual const std::vector<double>* Phinp() const = 0;
        virtual const std::vector<double>* Scalarnp() const = 0;
      };


      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/

      /*!
      \brief general interface to variable manager (template)

      The idea is, that there are the methods extract_element_and_node_values(...)
      and EvaluateGPVariables(..), which need to be called before evaluation.
      extract_element_and_node_values(...) reads the node values associated with the element
      from the global state vector and EvaluateGPVariables(..) performs the interpolation
      to the gauss points.
      All other methods are (more or less) constant access methods.

      As fixed sized CORE::LINALG::Matrix is used for saving the values, almost
      all variables managers are templated by the number of space dimensions 'nsd'
      and the number of element nodes 'nen'.

      This class is a pure virtual template interface class.
      It implements the factory method, which builds the concrete variable manager
      based on the element action.

      \author vuong
      */
      template <int nsd, int nen>
      class VariableManagerInterface : public VariableManagerMinAccess
      {
       public:
        //! constructor
        VariableManagerInterface(){};

        //! factory method
        static Teuchos::RCP<VariableManagerInterface<nsd, nen>> create_variable_manager(
            const DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter& para,
            const POROFLUIDMULTIPHASE::Action& action, Teuchos::RCP<CORE::MAT::Material> mat,
            const int numdofpernode, const int numfluidphases);

        //! extract element and node values from the discretization
        //! dofsetnum is the number of the porofluid-dofset on the current element
        //! default is set to zero, if called from a porofluidmultiphase-element
        //! otherwise it has to be explicitly passed from the caller
        virtual void extract_element_and_node_values(const CORE::Elements::Element& ele,
            const DRT::Discretization& discretization, CORE::Elements::Element::LocationArray& la,
            CORE::LINALG::Matrix<nsd, nen>& xyze, const int dofsetnum = 0) = 0;

        //! evaluate variables at gauss point
        virtual void EvaluateGPVariables(
            const CORE::LINALG::Matrix<nen, 1>& funct,  //! array for shape functions
            const CORE::LINALG::Matrix<nsd, nen>&
                derxy  //! array for shape function derivatives w.r.t x,y,z
            ) = 0;

        //! check if EvaluateGPVariables was called
        virtual void CheckIsEvaluated() const = 0;
        //! check if extract_element_and_node_values was called
        virtual void CheckIsExtracted() const = 0;
        //! return number of DOFs pre node
        virtual int NumDofPerNode() const = 0;

        //! @name Access methods
        const std::vector<double>* Phinp() const override = 0;
        virtual const CORE::LINALG::Matrix<nen, 1>* ElementPhinp(const int k) const = 0;
        virtual bool element_has_valid_vol_frac_pressure(const int ivolfrac) const = 0;
        virtual bool element_has_valid_vol_frac_species(const int ivolfrac) const = 0;
        virtual const std::vector<CORE::LINALG::Matrix<nsd, 1>>* GradPhinp() const = 0;
        virtual const std::vector<double>* Phidtnp() const = 0;
        virtual const std::vector<double>* Hist() const = 0;
        virtual const CORE::LINALG::Matrix<nsd, 1>* ConVelnp() const = 0;
        virtual double DivConVelnp() const = 0;
        virtual const CORE::LINALG::Matrix<nsd, nen>* EConVelnp() const = 0;
        virtual const CORE::LINALG::Matrix<nsd, 1>* Dispnp() const = 0;
        const std::vector<double>* Scalarnp() const override = 0;
        virtual const std::vector<CORE::LINALG::Matrix<nsd, 1>>* GradScalarnp() const = 0;
        //@}
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/

      /*!
      \brief  base class for all core variable managers

      This class is the base class for all variable managers. It basically
      implements some safety checks.

      \author vuong
      */
      template <int nsd, int nen>
      class VariableManagerBase : public VariableManagerInterface<nsd, nen>
      {
       public:
        //! constructor
        VariableManagerBase(int numdofpernode)
            : VariableManagerInterface<nsd, nen>(),
              numdofpernode_(numdofpernode),
              isextracted_(false),
              isevaluated_(false){};

        //! check if EvaluateGPVariables has been called
        void CheckIsEvaluated() const override
        {
          if (not isevaluated_)
            FOUR_C_THROW("EvaluateGPVariables has not been called on variable manager!");
        };

        //! check if extract_element_and_node_values has been called
        void CheckIsExtracted() const override
        {
          if (not isextracted_)
            FOUR_C_THROW(
                "extract_element_and_node_values has not been called on variable manager!");
        };

        //! return number of DOFs per node
        int NumDofPerNode() const override { return numdofpernode_; }

        //! @name Access methods (throw error by default)
        const std::vector<double>* Phinp() const override
        {
          FOUR_C_THROW("Access method Phinp() not implemented! Wrong VariableManager?");
          return nullptr;
        };
        const CORE::LINALG::Matrix<nen, 1>* ElementPhinp(const int k) const override
        {
          FOUR_C_THROW("Access method ElementPhinp() not implemented! Wrong VariableManager?");
          return nullptr;
        };
        bool element_has_valid_vol_frac_pressure(const int ivolfrac) const override
        {
          FOUR_C_THROW(
              "Access method element_has_valid_vol_frac_pressure() not implemented! Wrong "
              "VariableManager?");
          return 0.0;
        };
        bool element_has_valid_vol_frac_species(const int ivolfrac) const override
        {
          FOUR_C_THROW(
              "Access method element_has_valid_vol_frac_species() not implemented! Wrong "
              "VariableManager?");
          return 0.0;
        };
        const std::vector<CORE::LINALG::Matrix<nsd, 1>>* GradPhinp() const override
        {
          FOUR_C_THROW("Access method GradPhinp() not implemented! Wrong VariableManager?");
          return nullptr;
        };
        const std::vector<double>* Phidtnp() const override
        {
          FOUR_C_THROW("Access method Phidtnp() not implemented! Wrong VariableManager?");
          return nullptr;
        };

        const std::vector<double>* Hist() const override
        {
          FOUR_C_THROW("Access method Hist() not implemented! Wrong VariableManager?");
          return nullptr;
        };
        const CORE::LINALG::Matrix<nsd, 1>* ConVelnp() const override
        {
          FOUR_C_THROW("Access method ConVelnp() not implemented! Wrong VariableManager?");
          return nullptr;
        };
        double DivConVelnp() const override
        {
          FOUR_C_THROW("Access method DivConVelnp() not implemented! Wrong VariableManager?");
          return 0.0;
        };
        const CORE::LINALG::Matrix<nsd, nen>* EConVelnp() const override
        {
          FOUR_C_THROW("Access method EConVelnp() not implemented! Wrong VariableManager?");
          return nullptr;
        };
        const CORE::LINALG::Matrix<nsd, 1>* Dispnp() const override
        {
          FOUR_C_THROW("Access method Dispnp() not implemented! Wrong VariableManager?");
          return nullptr;
        };
        const std::vector<double>* Scalarnp() const override
        {
          FOUR_C_THROW("Access method Salarnp() not implemented! Wrong VariableManager?");
          return nullptr;
        };
        const std::vector<CORE::LINALG::Matrix<nsd, 1>>* GradScalarnp() const override
        {
          FOUR_C_THROW("Access method GradScalarnp() not implemented! Wrong VariableManager?");
          return nullptr;
        };

        //@}

       protected:
        //! number of DOFs per node
        const int numdofpernode_;
        //! flag if ExtracElementAndNodeValues was called
        bool isextracted_;
        //! flag if EvaluateGPVariables was called
        bool isevaluated_;
      };


      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/

      /*!
      \brief   minimal variable manager only holding the generic state vector 'phi'

      This class is the minimal version of a concrete variable managers. It only
      holds the generic state vector 'phi' (this can be e.g. the pressure or
      the saturation).

      \author vuong
      */
      template <int nsd, int nen>
      class VariableManagerPhi : public VariableManagerBase<nsd, nen>
      {
       public:
        //! constructor
        VariableManagerPhi(int numdofpernode)
            : VariableManagerBase<nsd, nen>(numdofpernode),
              ephinp_(numdofpernode),
              phinp_(numdofpernode, 0.0){};

        //! extract node values related to the state vector 'phinp'
        //! dofsetnum is the number of the porofluid-dofset on the current element
        //! default is set to zero, if called from a porofluidmultiphase-element
        //! otherwise it has to be explicitly passed from the caller
        void extract_element_and_node_values(const CORE::Elements::Element& ele,
            const DRT::Discretization& discretization, CORE::Elements::Element::LocationArray& la,
            CORE::LINALG::Matrix<nsd, nen>& xyze, const int dofsetnum = 0) override;

        //! evaluate state vector at gauss point
        void EvaluateGPVariables(
            const CORE::LINALG::Matrix<nen, 1>& funct,  //! array for shape functions
            const CORE::LINALG::Matrix<nsd, nen>&
                derxy  //! array for shape function derivatives w.r.t x,y,z
            ) override;

        //! access method
        const std::vector<double>* Phinp() const override
        {
          this->CheckIsEvaluated();
          return &phinp_;
        }

       protected:
        const CORE::LINALG::Matrix<nen, 1>* ElementPhinp(const int k) const override
        {
          this->CheckIsExtracted();
          return &ephinp_[k];
        }

        //! state variables at t_(n+1) or t_(n+alpha_F)
        std::vector<CORE::LINALG::Matrix<nen, 1>> ephinp_;

        //! scalar at t_(n+1) or t_(n+alpha_F)
        std::vector<double> phinp_;
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/

      /*!
      \brief  variable manager only holding the state vector 'phi' and its gradient

      This class only holds the generic state vector 'phi' (this can be e.g. the pressure or
      the saturation) and its gradient.

      \author vuong
      */
      template <int nsd, int nen>
      class VariableManagerPhiGradPhi : public VariableManagerPhi<nsd, nen>
      {
       public:
        //! constructor
        VariableManagerPhiGradPhi(int numdofpernode)
            : VariableManagerPhi<nsd, nen>(numdofpernode), gradphi_(numdofpernode){};

        //! evaluate phi and its gradient at gauss point
        void EvaluateGPVariables(
            const CORE::LINALG::Matrix<nen, 1>& funct,  //! array for shape functions
            const CORE::LINALG::Matrix<nsd, nen>&
                derxy  //! array for shape function derivatives w.r.t x,y,z
            ) override;

        //! access method
        const std::vector<CORE::LINALG::Matrix<nsd, 1>>* GradPhinp() const override
        {
          this->CheckIsEvaluated();
          return &gradphi_;
        }

       private:
        //! spatial gradient of current scalar value
        std::vector<CORE::LINALG::Matrix<nsd, 1>> gradphi_;
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief  base class for all decorator/wrapper classes, i.e. extensions to core variable
      managers

      The idea is to use a core variable manager (e.g. VariableManagerPhi) and extend
      it when the evaluation to be performed demands it. For this a decorator
       pattern is used, i.e. this class wraps another variable manager, extending it if necessary.

      This is the base class for all decorators.

      \author vuong
      */
      template <int nsd, int nen>
      class VariableManagerDecorator : public VariableManagerInterface<nsd, nen>
      {
       public:
        //! constructor
        VariableManagerDecorator(Teuchos::RCP<VariableManagerInterface<nsd, nen>> varmanager)
            : VariableManagerInterface<nsd, nen>(), varmanager_(varmanager){};

        //! @name Access methods
        const std::vector<double>* Phinp() const override { return varmanager_->Phinp(); };
        const CORE::LINALG::Matrix<nen, 1>* ElementPhinp(const int k) const override
        {
          return varmanager_->ElementPhinp(k);
        };
        double DivConVelnp() const override { return varmanager_->DivConVelnp(); };
        const std::vector<CORE::LINALG::Matrix<nsd, 1>>* GradPhinp() const override
        {
          return varmanager_->GradPhinp();
        };
        const std::vector<double>* Phidtnp() const override { return varmanager_->Phidtnp(); };

        const std::vector<double>* Hist() const override { return varmanager_->Hist(); };
        const CORE::LINALG::Matrix<nsd, 1>* ConVelnp() const override
        {
          return varmanager_->ConVelnp();
        };
        bool element_has_valid_vol_frac_pressure(const int ivolfrac) const override
        {
          return varmanager_->element_has_valid_vol_frac_pressure(ivolfrac);
        };
        bool element_has_valid_vol_frac_species(const int ivolfrac) const override
        {
          return varmanager_->element_has_valid_vol_frac_species(ivolfrac);
        };
        const CORE::LINALG::Matrix<nsd, nen>* EConVelnp() const override
        {
          return varmanager_->EConVelnp();
        };
        const CORE::LINALG::Matrix<nsd, 1>* Dispnp() const override
        {
          return varmanager_->Dispnp();
        };
        const std::vector<double>* Scalarnp() const override { return varmanager_->Scalarnp(); };
        const std::vector<CORE::LINALG::Matrix<nsd, 1>>* GradScalarnp() const override
        {
          return varmanager_->GradScalarnp();
        };

        //@}

        //! check if EvaluateGPVariables was called
        void CheckIsEvaluated() const override { varmanager_->CheckIsEvaluated(); }

        //! check if extract_element_and_node_values was called
        void CheckIsExtracted() const override { varmanager_->CheckIsExtracted(); }

        //! return number of DOFs per node
        int NumDofPerNode() const override { return varmanager_->NumDofPerNode(); }

       protected:
        //! wrapped variable manager
        Teuchos::RCP<VariableManagerInterface<nsd, nen>> varmanager_;
      };


      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/

      /*!
      \brief  decorator/wrapper class for variable manager, extensions to time derivatives

      This class entends a wraped variable manager by the instationary variables
      (time derivatives, and 'history'-vector)

      \author vuong
      */
      template <int nsd, int nen>
      class VariableManagerInstat : public VariableManagerDecorator<nsd, nen>
      {
       public:
        //! constructor
        VariableManagerInstat(Teuchos::RCP<VariableManagerInterface<nsd, nen>> varmanager)
            : VariableManagerDecorator<nsd, nen>(varmanager),
              ephidtnp_(varmanager->NumDofPerNode()),
              ehist_(varmanager->NumDofPerNode()),
              phidtnp_(varmanager->NumDofPerNode(), 0.0),
              hist_(varmanager->NumDofPerNode(), 0.0){};

        //! extract node values related to time derivatives
        //! dofsetnum is the number of the porofluid-dofset on the current element
        //! default is set to zero, if called from a porofluidmultiphase-element
        //! otherwise it has to be explicitly passed from the caller
        void extract_element_and_node_values(const CORE::Elements::Element& ele,
            const DRT::Discretization& discretization, CORE::Elements::Element::LocationArray& la,
            CORE::LINALG::Matrix<nsd, nen>& xyze, const int dofsetnum = 0) override;

        //! evaluate variables at gauss point
        void EvaluateGPVariables(
            const CORE::LINALG::Matrix<nen, 1>& funct,  //! array for shape functions
            const CORE::LINALG::Matrix<nsd, nen>&
                derxy  //! array for shape function derivatives w.r.t x,y,z
            ) override;

        //! @name Access methods

        //! get time derivative of state phi
        const std::vector<double>* Phidtnp() const override
        {
          this->varmanager_->CheckIsEvaluated();
          return &phidtnp_;
        }
        //! get history vector
        const std::vector<double>* Hist() const override
        {
          this->varmanager_->CheckIsEvaluated();
          return &hist_;
        }
        //@}

       private:
        //! time derivatives of state variables at t_(n+1)
        std::vector<CORE::LINALG::Matrix<nen, 1>> ephidtnp_;
        //! history vector of transported scalars
        std::vector<CORE::LINALG::Matrix<nen, 1>> ehist_;

        //! time derivative of scalar at t_(n+1)
        std::vector<double> phidtnp_;

        //! history data (or acceleration)
        std::vector<double> hist_;
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/

      /*!
      \brief  decorator/wrapper class for variable manager, extensions to moving domain variables

      This class entends a wraped variable manager by the  variables needed for an ALE
      formulation, i.e. mesh displacements and convective velocity.

      \author vuong
      */
      template <int nsd, int nen>
      class VariableManagerStruct : public VariableManagerDecorator<nsd, nen>
      {
       public:
        //! constructor
        VariableManagerStruct(
            int ndsvel, int ndsdisp, Teuchos::RCP<VariableManagerInterface<nsd, nen>> varmanager)
            : VariableManagerDecorator<nsd, nen>(varmanager),
              ndsvel_(ndsvel),
              ndsdisp_(ndsdisp),
              econvelnp_(true),  // initialized to zero
              edispnp_(true),    // initialized to zero
              divconvelint_(0.0),
              convelint_(true),
              dispint_(true){};

        //! extract variables related to structure coupling
        //! dofsetnum is the number of the porofluid-dofset on the current element
        //! default is set to zero, if called from a porofluidmultiphase-element
        //! otherwise it has to be explicitly passed from the caller
        void extract_element_and_node_values(const CORE::Elements::Element& ele,
            const DRT::Discretization& discretization, CORE::Elements::Element::LocationArray& la,
            CORE::LINALG::Matrix<nsd, nen>& xyze, const int dofsetnum = 0) override;

        //! evaluate variables at gauss point
        void EvaluateGPVariables(
            const CORE::LINALG::Matrix<nen, 1>& funct,  //! array for shape functions
            const CORE::LINALG::Matrix<nsd, nen>&
                derxy  //! array for shape function derivatives w.r.t x,y,z
            ) override;

        //! @name Access methods
        const CORE::LINALG::Matrix<nsd, 1>* ConVelnp() const override
        {
          this->varmanager_->CheckIsEvaluated();
          return &convelint_;
        }
        const CORE::LINALG::Matrix<nsd, 1>* Dispnp() const override
        {
          this->varmanager_->CheckIsEvaluated();
          return &dispint_;
        }
        double DivConVelnp() const override
        {
          this->varmanager_->CheckIsEvaluated();
          return divconvelint_;
        }
        const CORE::LINALG::Matrix<nsd, nen>* EConVelnp() const override
        {
          this->varmanager_->CheckIsEvaluated();
          return &econvelnp_;
        }

        //@}

       private:
        // number of dofset associated with velocity related dofs
        const int ndsvel_;
        // number of dofset associated with displacement related dofs
        const int ndsdisp_;

        //! nodal velocity values at t_(n+1) or t_(n+alpha_F)
        CORE::LINALG::Matrix<nsd, nen> econvelnp_;
        //! nodal displacement values for ALE
        CORE::LINALG::Matrix<nsd, nen> edispnp_;

        // velocity divergence required for conservative form
        double divconvelint_;

        // structure velocity
        CORE::LINALG::Matrix<nsd, 1> convelint_;

        // gauss point displacements
        CORE::LINALG::Matrix<nsd, 1> dispint_;
      };


      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/

      /*!
      \brief  decorator/wrapper class for variable manager, extensions to ScaTra coupling

      This class entends a wraped variable manager by the  variables needed for a ScaTra
      coupling, i.e. the scalar values.

      \author vuong
      */
      template <int nsd, int nen>
      class VariableManagerScalar : public VariableManagerDecorator<nsd, nen>
      {
       public:
        //! constructor
        VariableManagerScalar(
            int ndsscalar, Teuchos::RCP<VariableManagerInterface<nsd, nen>> varmanager)
            : VariableManagerDecorator<nsd, nen>(varmanager),
              ndsscalar_(ndsscalar),
              escalarnp_(),
              scalarnp_(),
              gradscalarnp_(){};

        //! extract node values related to ScaTra coupling
        //! dofsetnum is the number of the porofluid-dofset on the current element
        //! default is set to zero, if called from a porofluidmultiphase-element
        //! otherwise it has to be explicitly passed from the caller
        void extract_element_and_node_values(const CORE::Elements::Element& ele,
            const DRT::Discretization& discretization, CORE::Elements::Element::LocationArray& la,
            CORE::LINALG::Matrix<nsd, nen>& xyze, const int dofsetnum = 0) override;

        void EvaluateGPVariables(
            const CORE::LINALG::Matrix<nen, 1>& funct,  //! array for shape functions
            const CORE::LINALG::Matrix<nsd, nen>&
                derxy  //! array for shape function derivatives w.r.t x,y,z
            ) override;

        //! access method
        const std::vector<double>* Scalarnp() const override
        {
          this->varmanager_->CheckIsEvaluated();
          return &scalarnp_;
        };

        //! access method
        const std::vector<CORE::LINALG::Matrix<nsd, 1>>* GradScalarnp() const override
        {
          this->varmanager_->CheckIsEvaluated();
          return &gradscalarnp_;
        }

       private:
        // number of dofset associated with scalar related dofs
        const int ndsscalar_;

        //! nodal scalar values for scatra coupling
        std::vector<CORE::LINALG::Matrix<nen, 1>> escalarnp_;

        //! scalar values
        std::vector<double> scalarnp_;

        //! spatial gradient of current scalar value
        std::vector<CORE::LINALG::Matrix<nsd, 1>> gradscalarnp_;
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/

      /*!
      \brief  decorator/wrapper class for variable manager, extension to maximum nodal volume
      fraction

      This class evaluates for each element the maximum value of a specific volume fraction at its
      nodes. This is necessary to decide if the calculation of volume fraction pressures inside an
      element makes sense returns a bool in element_has_valid_vol_frac_pressure to decide if volume
      fraction pressure has to be evaluated

      \author kremheller
      */
      template <int nsd, int nen>
      class VariableManagerMaximumNodalVolFracValue : public VariableManagerDecorator<nsd, nen>
      {
       public:
        //! constructor
        VariableManagerMaximumNodalVolFracValue(const int numvolfrac,
            Teuchos::RCP<VariableManagerInterface<nsd, nen>> varmanager,
            Teuchos::RCP<CORE::MAT::Material> multiphasemat)
            : VariableManagerDecorator<nsd, nen>(varmanager),
              numvolfrac_(numvolfrac),
              ele_has_valid_volfrac_press_(numvolfrac_, false),
              ele_has_valid_volfrac_spec_(numvolfrac_, false),
              multiphasemat_(multiphasemat){};

        //! extract node values related to time derivatives
        //! dofsetnum is the number of the porofluid-dofset on the current element
        //! default is set to zero, if called from a porofluidmultiphase-element
        //! otherwise it has to be explicitly passed from the caller
        void extract_element_and_node_values(const CORE::Elements::Element& ele,
            const DRT::Discretization& discretization, CORE::Elements::Element::LocationArray& la,
            CORE::LINALG::Matrix<nsd, nen>& xyze, const int dofsetnum = 0) override;

        //! evaluate variables at gauss point
        void EvaluateGPVariables(
            const CORE::LINALG::Matrix<nen, 1>& funct,  //! array for shape functions
            const CORE::LINALG::Matrix<nsd, nen>&
                derxy  //! array for shape function derivatives w.r.t x,y,z
            ) override;

        //! @name Access methods
        bool element_has_valid_vol_frac_pressure(const int ivolfrac) const override
        {
          this->varmanager_->CheckIsExtracted();
          if (ivolfrac >= numvolfrac_)
            FOUR_C_THROW(
                "%i is bigger than the number of volume fractions %i in the VariableManager",
                ivolfrac + 1, numvolfrac_);

          return ele_has_valid_volfrac_press_[ivolfrac];
        }

        bool element_has_valid_vol_frac_species(const int ivolfrac) const override
        {
          this->varmanager_->CheckIsExtracted();
          if (ivolfrac >= numvolfrac_)
            FOUR_C_THROW(
                "%i is bigger than the number of volume fractions %i in the VariableManager",
                ivolfrac + 1, numvolfrac_);

          return ele_has_valid_volfrac_spec_[ivolfrac];
        }

        //@}

       private:
        const int numvolfrac_;

        //! check if volume fraction pressure equation can be evaluated within this element
        std::vector<bool> ele_has_valid_volfrac_press_;

        //! check if volume fraction species equation can be evaluated within this element
        std::vector<bool> ele_has_valid_volfrac_spec_;

        Teuchos::RCP<CORE::MAT::Material> multiphasemat_;
      };

    }  // namespace POROFLUIDMANAGER

  }  // namespace ELEMENTS
}  // namespace DRT



FOUR_C_NAMESPACE_CLOSE

#endif
