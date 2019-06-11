/*----------------------------------------------------------------------*/
/*!
 \brief helper class for evaluation of the governing equation of multiphase porou flow

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUID_EVALUATOR_H_
#define SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUID_EVALUATOR_H_

#include "porofluidmultiphase_ele_action.H"

#include "../drt_lib/drt_dserror.H"

#include "../linalg/linalg_fixedsizematrix.H"

#include "Teuchos_RCP.hpp"

#include <vector>

class Epetra_SerialDenseMatrix;
class Epetra_SerialDenseVector;

namespace DRT
{
  namespace UTILS
  {
    // forward declaration
    class VariableExprFunction;
  }  // namespace UTILS
  namespace ELEMENTS
  {
    class PoroFluidMultiPhaseEleParameter;

    namespace POROFLUIDMANAGER
    {
      class PhaseManagerInterface;
      template <int, int>
      class VariableManagerInterface;
    }  // namespace POROFLUIDMANAGER

    namespace POROFLUIDEVALUATOR
    {
      template <int, int>
      class EvaluatorInterface;


      /*!
      \brief A helper class for element assembly of the porous multiphase flow equations

      The thing is, that in this formulation for the porous multiphase flow equations,
      one equation is somewhat special. That is, that one equation is actually the
      sum of all phases incorporating the phase constraint, i.e. that all saturations
      sum up to 1.

      For this propose these small classes handle in which DOF, i.e. in which row of
      the element matrix, the terms are assembled into. For now, there are three
      assembly classes:

         1. a general, pure virtual interface
         2. standard assemble, i.e. assemble into one phase, that is the phase that
            currently evaluated
         3. assemble into two phases. The current phase and the summed up phase

      \author vuong
      */

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      //! general interface to helper class for assembly into element matrix
      class AssembleInterface
      {
       public:
        //! constructor
        AssembleInterface(const bool inittimederiv) : inittimederiv_(inittimederiv){};

        //! destructor
        virtual ~AssembleInterface(){};

        virtual int NumPhasesToAssembleInto() const = 0;

        virtual int PhaseToAssembleInto(int iassemble, int numdofpernode) const = 0;

        bool CalcInitTimeDeriv() const { return inittimederiv_; };

       private:
        // do we calculate the inital time derivative?
        const bool inittimederiv_;
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      //! helper class for standard assembly into element matrix
      class AssembleStandard : public AssembleInterface
      {
       public:
        //! constructor
        AssembleStandard(int curphase, const bool inittimederiv)
            : AssembleInterface(inittimederiv), curphase_(curphase){};

        //! destructor
        virtual ~AssembleStandard(){};

        int NumPhasesToAssembleInto() const { return 1; };

        int PhaseToAssembleInto(int iassemble, int numdofpernode) const { return curphase_; };

       private:
        const int curphase_;
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      //! helper class for assembly into additional row in element matrix
      class AssembleAlsoIntoOtherPhase : public AssembleInterface
      {
       public:
        //! constructor
        AssembleAlsoIntoOtherPhase(int curphase, int otherphase, const bool inittimederiv)
            : AssembleInterface(inittimederiv), phasestoassemble_(2)
        {
          phasestoassemble_[0] = curphase;
          phasestoassemble_[1] = otherphase;
        };

        //! destructor
        virtual ~AssembleAlsoIntoOtherPhase(){};

        int NumPhasesToAssembleInto() const { return 2; };

        int PhaseToAssembleInto(int iassemble, int numdofpernode) const
        {
          return phasestoassemble_[iassemble];
        };

       private:
        std::vector<int> phasestoassemble_;
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/

      /*!
      \brief Evaluator class for the element matrices

      These classes do the actual work, they assemble the terms within the porous multiphase
      flow equations.

      The idea is simple. Each additive term in the equation has its own evaluator.
      This way, single summands can be turned on and off easily. The key methods are
      EvaluateMatrix(..) and EvaluateVector(..). One assembles the linearization, the
      other the RHS vector.

      This class is a general interface class for evaluation of the  element matrix and the RHS
      vector. It templated by the space dimensions 'nsd' and the number of nodes 'nen'. It comprises
      the pure virtual functions EvaluateMatrix(..) and EvaluateVector(..) and the factory method
      CreateEvaluator(). The factory method is the central place, where the terms are defined.

      \author vuong
      */
      template <int nsd, int nen>
      class EvaluatorInterface
      {
       public:
        //! constructor
        EvaluatorInterface(){};

        //! destructor
        virtual ~EvaluatorInterface(){};

        //! factory method
        static Teuchos::RCP<EvaluatorInterface<nsd, nen>> CreateEvaluator(
            const DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter& para,
            const POROFLUIDMULTIPHASE::Action& action, int numdofpernode, int numfluidphases,
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager);

        //! evaluate matrixes (stiffness)
        virtual void EvaluateMatrix(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor times time-integration factor
            double fac            //!< domain-integration factor
            ) = 0;

        //! evaluate vectors (RHS vector)
        virtual void EvaluateVector(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
            double fac      //!< domain-integration factor
            ) = 0;

        //! evaluate off-diagonal coupling matrix with structure
        virtual void EvaluateMatrixODStruct(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor times time-integration factor
            double fac,           //!< domain-integration factor
            double det) = 0;

        //! evaluate off-diagonal coupling matrix with scatra
        virtual void EvaluateMatrixODScatra(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor times time-integration factor
            double fac            //!< domain-integration factor
            ) = 0;
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief Evaluator class for evaluation of element matrix/vectors for multiple phases

      This class wraps multiple evaluators. For evaluation of matrix and vector, it just loops over
      all single evaluators.

      \author vuong
      */
      template <int nsd, int nen>
      class MultiEvaluator : public EvaluatorInterface<nsd, nen>
      {
       public:
        //! constructor
        MultiEvaluator(){};

        //! destructor
        virtual ~MultiEvaluator(){};

        //! evaluate matrixes (stiffness)
        void EvaluateMatrix(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        )
        {
          // loop over the evaluators
          typename std::vector<Teuchos::RCP<EvaluatorInterface<nsd, nen>>>::iterator it;
          for (it = evaluators_.begin(); it != evaluators_.end(); it++)
            (*it)->EvaluateMatrix(elemat, funct, derxy, numdofpernode, phasemanager,
                variablemanager, timefacfac, fac);
        };

        //! evaluate vectors (RHS vector)
        void EvaluateVector(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
            double fac      //!< domain-integration factor
        )
        {
          // loop over the evaluators
          typename std::vector<Teuchos::RCP<EvaluatorInterface<nsd, nen>>>::iterator it;
          for (it = evaluators_.begin(); it != evaluators_.end(); it++)
            (*it)->EvaluateVector(elevec, funct, derxy, xyze, numdofpernode, phasemanager,
                variablemanager, rhsfac, fac);
        };

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStruct(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det)
        {
          // loop over the evaluators
          typename std::vector<Teuchos::RCP<EvaluatorInterface<nsd, nen>>>::iterator it;
          for (it = evaluators_.begin(); it != evaluators_.end(); it++)
          {
            (*it)->EvaluateMatrixODStruct(elemat, funct, deriv, derxy, xjm, numdofpernode,
                phasemanager, variablemanager, timefacfac, fac, det);
          }
        };

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatra(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        )
        {
          // loop over the evaluators
          typename std::vector<Teuchos::RCP<EvaluatorInterface<nsd, nen>>>::iterator it;
          for (it = evaluators_.begin(); it != evaluators_.end(); it++)
          {
            (*it)->EvaluateMatrixODScatra(elemat, funct, derxy, numdofpernode, phasemanager,
                variablemanager, timefacfac, fac);
          }
        };

        //! add an evaluator to the list of evaluators
        void AddEvaluator(Teuchos::RCP<EvaluatorInterface<nsd, nen>> evaluator)
        {
          evaluators_.push_back(evaluator);
        };

       private:
        //! list of all evaluators
        std::vector<Teuchos::RCP<EvaluatorInterface<nsd, nen>>> evaluators_;
      };


      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief general base class for evaluation of one term

      This class  is the base class for single summand evaluators. It comprises an assembler object,
      defining in which rows the term is to be assembled into. The pure virtual
      methods EvaluateMatrixAndAssemble(...) and EvaluateVectorAndAssemble(...) defined
      the actual term and are to be implemented in derived classes.

      \author vuong
      */
      template <int nsd, int nen>
      class EvaluatorBase : public EvaluatorInterface<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorBase(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : assembler_(assembler), myphase_(curphase){};

        //! destructor
        virtual ~EvaluatorBase(){};

        //! evaluate matrixes (stiffness)
        void EvaluateMatrix(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        )
        {
          // the assembler class decides, where the terms are assembled into
          for (int iassemble = 0; iassemble < assembler_->NumPhasesToAssembleInto(); iassemble++)
          {
            // call the actual evaluation and assembly of the respective term (defined by derived
            // class)
            EvaluateMatrixAndAssemble(elemat, funct, derxy, myphase_,
                assembler_->PhaseToAssembleInto(iassemble, numdofpernode), numdofpernode,
                phasemanager, variablemanager, timefacfac, fac, assembler_->CalcInitTimeDeriv());
          }
        };

        //! evaluate vectors (RHS vector)
        void EvaluateVector(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
            double fac      //!< domain-integration factor
        )
        {
          // the assembler class decides, where the terms are assembled into
          for (int iassemble = 0; iassemble < assembler_->NumPhasesToAssembleInto(); iassemble++)
          {
            // call the actual evaluation and assembly of the respective term (defined by derived
            // class)
            EvaluateVectorAndAssemble(elevec, funct, derxy, xyze, myphase_,
                assembler_->PhaseToAssembleInto(iassemble, numdofpernode), numdofpernode,
                phasemanager, variablemanager, rhsfac, fac, assembler_->CalcInitTimeDeriv());
          }
        };

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStruct(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det)
        {
          // the assembler class decides, where the terms are assembled into
          for (int iassemble = 0; iassemble < assembler_->NumPhasesToAssembleInto(); iassemble++)
          {
            // call the actual evaluation and assembly of the respective term (defined by derived
            // class)
            EvaluateMatrixODStructAndAssemble(elemat, funct, deriv, derxy, xjm, myphase_,
                assembler_->PhaseToAssembleInto(iassemble, numdofpernode), numdofpernode,
                phasemanager, variablemanager, timefacfac, fac, det);
          }
        };

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatra(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        )
        {
          // the assembler class decides, where the terms are assembled into
          for (int iassemble = 0; iassemble < assembler_->NumPhasesToAssembleInto(); iassemble++)
          {
            // call the actual evaluation and assembly of the respective term (defined by derived
            // class)
            EvaluateMatrixODScatraAndAssemble(elemat, funct, derxy, myphase_,
                assembler_->PhaseToAssembleInto(iassemble, numdofpernode), numdofpernode,
                phasemanager, variablemanager, timefacfac, fac);
          }
        };

       protected:
        // actual evaluation and assembly of the respective term in the stiffness matrix (defined by
        // derived class)
        virtual void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of phase to add into
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
            ) = 0;

        // actual evaluation and assembly of the respective term in the RHS vector (defined by
        // derived class)
        virtual void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of phase to add into
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
            ) = 0;

        // actual evaluation and assembly of the respective term in the off-diagonal matrix (defined
        // by derived class)
        virtual void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of phase to add into
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det) = 0;

        // actual evaluation and assembly of the respective term in the off-diagonal matrix (defined
        // by derived class)
        virtual void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of phase to add into
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
            ) = 0;

        // OD-mesh linearization of diffusive term
        void CalcDiffODMesh(Epetra_SerialDenseMatrix& mymat,  //!< element matrix to be filled
            const LINALG::Matrix<nsd, nen>& deriv, const LINALG::Matrix<nsd, nen>& derxy,
            const LINALG::Matrix<nsd, nsd>& xjm, const LINALG::Matrix<nsd, 1>& diffflux,
            const LINALG::Matrix<nsd, 1>& refgrad, const LINALG::Matrix<nsd, 1>& grad,
            const double timefacfac, const double difffac, const int numdofpernode,
            const int phasetoadd);

        // OD-mesh linearization of fac (Jacobian)
        void CalcLinFacODMesh(Epetra_SerialDenseMatrix& mymat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,                //! array for shape functions
            const LINALG::Matrix<nsd, nen>& derxy, const double vrhs, const int numdofpernode,
            const int phasetoadd);

        // OD-mesh linearization of divergence term
        void CalcDivVelODMesh(Epetra_SerialDenseMatrix& mymat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,                //! array for shape functions
            const LINALG::Matrix<nsd, nen>& deriv, const LINALG::Matrix<nsd, nen>& derxy,
            const LINALG::Matrix<nsd, nsd>& xjm, const LINALG::Matrix<nsd, nsd>& gridvelderiv,
            const double timefacfac, const double fac, const double det, const int numdofpernode,
            const int phasetoadd);

        // linearization of a term scaled with saturation after fluid dofs
        void SaturationLinearizationFluid(Epetra_SerialDenseMatrix& mymat,
            const LINALG::Matrix<nen, 1>& funct, const double prefac, const int numdofpernode,
            const int numfluidphases, const int curphase, const int phasetoadd,
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager);

        // linearization of a term scaled with porosity after fluid dofs
        void PorosityLinearizationFluid(Epetra_SerialDenseMatrix& mymat,
            const LINALG::Matrix<nen, 1>& funct, const double prefac, const int numdofpernode,
            const int phasetoadd, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager);

       private:
        //! assemble strategy
        Teuchos::RCP<AssembleInterface> assembler_;
        //! phase the term is associated with
        const int myphase_;
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of convective term into the element matrix

      This class implements the convective term $(w, v \nabla \cdot S )$.

      \note this term is not used, since the equations are written in a Lagrangian description
            w.r.t. skeleton.

      \author vuong
      */
      template <int nsd, int nen>
      class EvaluatorConv : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorConv(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorConv(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of divergence of the (mesh) velocity field

      This class implements the term $(w, \nabla \cdot v^s )$.

      \author vuong
      */
      template <int nsd, int nen>
      class EvaluatorDivVel : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorDivVel(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorDivVel(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of divergence of the (mesh) velocity field, scaled by saturation

      This class implements the term $(w, S \nabla \cdot v^s )$.

      \author vuong
      */
      template <int nsd, int nen>
      class EvaluatorSatDivVel : public EvaluatorDivVel<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorSatDivVel(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorDivVel<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorSatDivVel(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of biot stabilization terms

      This class implements the term $(w, \tau R_{struct} )$.

      \author vuong
      */
      template <int nsd, int nen>
      class EvaluatorBiotStab : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorBiotStab(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorBiotStab(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };


      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of diffusive term into the element matrix

      This class implements the term $( \nabla w, K \nabla p )$.

      \author vuong
      */
      template <int nsd, int nen>
      class EvaluatorDiff : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorDiff(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorDiff(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of reactive term into the element matrix

      This class implements all kinds of reactive terms, defined by the phasemanager.

      \author vuong
      */
      template <int nsd, int nen>
      class EvaluatorReac : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorReac(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorReac(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of mass term (pressure) into the element matrix

      This class implements the term $( w,porosity S/K \frac{\partial p}{\partial t} )$.

      \author vuong
      */
      template <int nsd, int nen>
      class EvaluatorMassPressure : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorMassPressure(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorMassPressure(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );

        //! get transient term for rhs and OD
        double GetRhsTrans(int curphase,  //!< index of current phase
            int phasetoadd,               //!< index of current phase
            int numdofpernode,            //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
            double fac      //!< domain-integration factor);
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of mass term (solid pressure) into the element matrix

      This class implements the term $( w, \frac{(1-\porosity) }{K_s} \frac{\partial p^s}{\partial
      t} )$.

      \author vuong
      */
      template <int nsd, int nen>
      class EvaluatorMassSolidPressure : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorMassSolidPressure(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorMassSolidPressure(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );

        //! get transient term for rhs and OD
        double GetRhsTrans(int curphase,  //!< index of current phase
            int phasetoadd,               //!< index of current phase
            int numdofpernode,            //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
            double fac      //!< domain-integration factor);
        );
      };


      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of mass term (solid pressure), scaled with saturation into the
      element matrix

      This class implements the term $( w, S\frac{(1-\porosity) }{K_s} \frac{\partial p^s}{\partial
      t} )$.

      \author vuong
      */
      template <int nsd, int nen>
      class EvaluatorMassSolidPressureSat : public EvaluatorMassSolidPressure<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorMassSolidPressureSat(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorMassSolidPressure<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorMassSolidPressureSat(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of mass term (solid saturation) into the element matrix

      This class implements the term $( w, \porosity \frac{\partial S}{\partial t} )$.

      \author vuong
      */
      template <int nsd, int nen>
      class EvaluatorMassSaturation : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorMassSaturation(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorMassSaturation(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );

        //! get transient term for rhs and OD
        double GetRhsTrans(int curphase,  //!< index of current phase
            int phasetoadd,               //!< index of current phase
            int numdofpernode,            //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
            double fac      //!< domain-integration factor);
        );
      };


      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief helper class for evaluation of pressure and saturation

      This class implements the post processing of pressures and saturation at the nodes.

      \author vuong
      */
      template <int nsd, int nen>
      class EvaluatorPressureAndSaturation : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorPressureAndSaturation(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorPressureAndSaturation(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief helper class for evaluation of the solid pressure

      This class implements the post processing of the solid pressure at the nodes.

      \author vuong
      */
      template <int nsd, int nen>
      class EvaluatorSolidPressure : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorSolidPressure(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorSolidPressure(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief helper class for finding the volume fraction pressure dofs which
             do not have to be evaluated

      The DOFs where the volume fraction pressure is valid, i.e. has a physical meaning
      are those where the respective volume fraction is greater than a threshold (MIN_VOLFRAC)
      These are identified here be setting ones into the valid_volfracpress_dofs_-vector

      \author kremheller
      */
      template <int nsd, int nen>
      class EvaluatorValidVolFracPressures : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorValidVolFracPressures(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorValidVolFracPressures(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief helper class for evaluation of porosity

      This class implements the post processing of the porosity at the nodes.

      \author kremheller
      */
      template <int nsd, int nen>
      class EvaluatorPorosity : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorPorosity(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorPorosity(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief helper class for evaluation of domain integrals

      This class implements the evaluation of domain integrals which can be used for output

      \author kremheller
      */
      template <int nsd, int nen>
      class EvaluatorDomainIntegrals : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorDomainIntegrals(Teuchos::RCP<AssembleInterface> assembler, int curphase,
            std::vector<int> domainint_funct, int numscal)
            : EvaluatorBase<nsd, nen>(assembler, curphase),
              domainint_funct_(domainint_funct),
              numscal_(numscal){};

        //! destructor
        virtual ~EvaluatorDomainIntegrals(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );

        //! cast to VarExp-function
        inline DRT::UTILS::VariableExprFunction& Function(int functnum) const;

        //! vector holding the functions to be integrated
        std::vector<int> domainint_funct_;

        //! number of scalars in system
        int numscal_;
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief helper class for reconstruction of flux

      This class implements the linearization of the flux reconstruction matrix (L_2 projection).
      Only the matrix! For RHS see class ReconstructFluxRHS.

      \author vuong
      */
      template <int nsd, int nen>
      class ReconstructFluxLinearization : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        ReconstructFluxLinearization(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~ReconstructFluxLinearization(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief helper class for reconstruction of flux

      This class implements the RHS the flux reconstruction matrix (L_2 projection).
      Only the RHS! For the matrix see ReconstructFlux.

      \author vuong
      */
      template <int nsd, int nen>
      class ReconstructFluxRHS : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        ReconstructFluxRHS(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~ReconstructFluxRHS(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of additional volume fraction terms in fluid equations
      (instationary terms) into the element matrix

      This class implements the term $( w, \frac{-\sum^volfrac \phi_volfrac) }{K_s} \frac{\partial
      p^s}{\partial t}
                                          -\sum^volfrac \frac{\partial\phi_volfrac) }{\partial t})

      \author kremheller
      */
      template <int nsd, int nen>
      class EvaluatorVolFracAddInstatTerms : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorVolFracAddInstatTerms(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorVolFracAddInstatTerms(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );

        //! get transient term for rhs and OD
        double GetRhs(int curphase,  //!< index of current phase
            int phasetoadd,          //!< index of current phase
            int numdofpernode,       //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,  //!< time-integration factor for rhs times domain-integration factor
            double fac      //!< domain-integration factor);
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of additional term in fluid equation introduced by volume
      fractions:
             - divergence of the (mesh) velocity field times sum of volume fraction

      This class implements the term $(w, -\sum^volfrac \phi_volfrac \nabla \cdot v^s )$.

      \author kremheller
      */
      template <int nsd, int nen>
      class EvaluatorVolFracAddDivVelTerm : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorVolFracAddDivVelTerm(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorVolFracAddDivVelTerm(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of additional term in fluid equation introduced by volume
      fractions:
             - divergence of the (mesh) velocity field times sum of volume fraction scaled with
      saturation

      This class implements the term $(w, S* -\sum^volfrac \phi_volfrac \nabla \cdot v^s )$.

      \author kremheller
      */
      template <int nsd, int nen>
      class EvaluatorVolFracAddDivVelTermSat : public EvaluatorVolFracAddDivVelTerm<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorVolFracAddDivVelTermSat(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorVolFracAddDivVelTerm<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorVolFracAddDivVelTermSat(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of additional volume fraction terms in fluid equations
      (instationary terms) scaled with saturation into the element matrix

      This class implements the term $( w, S* ( \frac{-\sum^volfrac \phi_volfrac) }{K_s}
      \frac{\partial p^s}{\partial t}
                                          -\sum^volfrac \frac{\partial\phi_volfrac) }{\partial t}) )

      \author kremheller
      */
      template <int nsd, int nen>
      class EvaluatorVolFracAddInstatTermsSat : public EvaluatorVolFracAddInstatTerms<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorVolFracAddInstatTermsSat(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorVolFracAddInstatTerms<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorVolFracAddInstatTermsSat(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of instationary term for volume fractions

      This class implements the term $( w, rho \frac{\partial volfrac^i}{\partial t} )$.
      It is assembled into the equation for volume fractions and for volume fraction pressures

      \author kremheller
      */
      template <int nsd, int nen>
      class EvaluatorVolFracInstat : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorVolFracInstat(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorVolFracInstat(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of divergence of the (mesh) velocity field for volume fractions

      This class implements the term $(w, rho volfrac \nabla \cdot v^s )$.
      It is assembled into the equation for volume fractions and for volume fraction pressures

      \author kremheller
      */
      template <int nsd, int nen>
      class EvaluatorVolFracDivVel : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorVolFracDivVel(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorVolFracDivVel(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of diffusive term into the element matrix

      This class implements the term $( \nabla w, D \nabla volfrac )$.

      \author kremheller
      */
      template <int nsd, int nen>
      class EvaluatorVolFracDiff : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorVolFracDiff(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorVolFracDiff(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of reactive term of volume fractions into the element matrix

      This class implements the term $(  w, reac )$.

      \author kremheller
      */
      template <int nsd, int nen>
      class EvaluatorVolFracReac : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorVolFracReac(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorVolFracReac(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of additional flux depending on Scatra primary variable

      This class implements the term $( \nabla w, D \nabla phi_scatra )$.

      \author kremheller
      */
      template <int nsd, int nen>
      class EvaluatorVolFracAddFlux : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorVolFracAddFlux(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorVolFracAddFlux(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of diffusive volfrac pressure term into the element matrix

      This class implements the term $( \nabla w, k/\mu \nabla volfrac_pressure )$.

      \author kremheller
      */
      template <int nsd, int nen>
      class EvaluatorVolFracPressureDiff : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorVolFracPressureDiff(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorVolFracPressureDiff(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

      /*----------------------------------------------------------------------*
       * **********************************************************************
       *----------------------------------------------------------------------*/
      /*!
      \brief class for evaluation of reactive term of volume fraction pressures
             into the element matrix

      This class implements the term $(  w, reac )$.

      \author kremheller
      */
      template <int nsd, int nen>
      class EvaluatorVolFracPressureReac : public EvaluatorBase<nsd, nen>
      {
       public:
        //! constructor
        EvaluatorVolFracPressureReac(Teuchos::RCP<AssembleInterface> assembler, int curphase)
            : EvaluatorBase<nsd, nen>(assembler, curphase){};

        //! destructor
        virtual ~EvaluatorVolFracPressureReac(){};

       protected:
        //! evaluate element matrix
        void EvaluateMatrixAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            bool inittimederiv    //!< calculate only parts for initial time derivative
        );

        //! evaluate element RHS vector
        void EvaluateVectorAndAssemble(
            std::vector<Epetra_SerialDenseVector*>& elevec,  //!< element vector to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nen>& xyze,  //!< current element coordinates
            int curphase,                          //!< index of current phase
            int phasetoadd,                        //!< index of current phase
            int numdofpernode,                     //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double rhsfac,      //!< time-integration factor for rhs times domain-integration factor
            double fac,         //!< domain-integration factor
            bool inittimederiv  //!< calculate only parts for initial time derivative
        );

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODStructAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                deriv,  //! array for shape function derivatives w.r.t r,s,t
            const LINALG::Matrix<nsd, nen>&
                derxy,  //! array for shape function derivatives w.r.t x,y,z
            const LINALG::Matrix<nsd, nsd>& xjm,
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac,           //!< domain-integration factor times time-integration factor
            double det);

        //! evaluate off-diagonal coupling matrix with structure
        void EvaluateMatrixODScatraAndAssemble(
            std::vector<Epetra_SerialDenseMatrix*>& elemat,  //!< element matrix to be filled
            const LINALG::Matrix<nen, 1>& funct,             //! array for shape functions
            const LINALG::Matrix<nsd, nen>&
                derxy,          //! array for shape function derivatives w.r.t x,y,z
            int curphase,       //!< index of current phase
            int phasetoadd,     //!< index of current phase
            int numdofpernode,  //!< total number of DOFs/phases
            const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,  //!< phase manager
            const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>&
                variablemanager,  //!< variable manager
            double timefacfac,    //!< domain-integration factor
            double fac            //!< domain-integration factor times time-integration factor
        );
      };

    }  // namespace POROFLUIDEVALUATOR

  }  // namespace ELEMENTS
}  // namespace DRT


#endif /* SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUID_EVALUATOR_H_ */
