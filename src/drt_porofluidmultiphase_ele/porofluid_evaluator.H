/*----------------------------------------------------------------------*/
/*!
 \file porofluid_evaluator.H

 \brief helper class for evaluation of the governing equation of multiphase porou flow

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
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
namespace ELEMENTS
{
class PoroFluidMultiPhaseEleParameter;

namespace POROFLUIDMANAGER
{
  class PhaseManagerInterface;
  template <int,int> class VariableManagerInterface;
}

namespace POROFLUIDEVALUATOR
{
template<int,int> class EvaluatorInterface;


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
  AssembleInterface(const bool inittimederiv):
  inittimederiv_(inittimederiv)
  {};

  //! destructor
  virtual ~AssembleInterface(){};

  virtual int NumPhasesToAssembleInto() const =0;

  virtual int PhaseToAssembleInto(int iassemble, int numdofpernode) const =0;

  bool CalcInitTimeDeriv() const {return inittimederiv_;};

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
  AssembleStandard(int curphase, const bool inittimederiv):
  AssembleInterface(inittimederiv),
  curphase_(curphase)
  {};

  //! destructor
  virtual ~AssembleStandard(){};

  int NumPhasesToAssembleInto() const {return 1;};

  int PhaseToAssembleInto(int iassemble, int numdofpernode) const {return curphase_;};

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
  AssembleAlsoIntoOtherPhase(int curphase ,int otherphase, const bool inittimederiv):
    AssembleInterface(inittimederiv),
      phasestoassemble_(2)
  {
    phasestoassemble_[0] = curphase;
    phasestoassemble_[1] = otherphase;
  };

  //! destructor
  virtual ~AssembleAlsoIntoOtherPhase(){};

  int NumPhasesToAssembleInto() const {return 2;};

  int PhaseToAssembleInto(int iassemble, int numdofpernode) const {return phasestoassemble_[iassemble];};

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

This class is a general interface class for evaluation of the  element matrix and the RHS vector.
It templated by the space dimensions 'nsd' and the number of nodes 'nen'.
It comprises the pure virtual functions EvaluateMatrix(..) and EvaluateVector(..)
and the factory method CreateEvaluator(). The factory method is the central place,
 where the terms are defined.

\author vuong
*/
template<int nsd, int nen>
class EvaluatorInterface
{
public:
  //! constructor
  EvaluatorInterface(){};

  //! destructor
  virtual ~EvaluatorInterface(){};

  //! factory method
  static  Teuchos::RCP< EvaluatorInterface<nsd,nen> > CreateEvaluator(
      const DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter& para,
      const POROFLUIDMULTIPHASE::Action&                    action,
      int                                                   numdofpernode,
      const POROFLUIDMANAGER::PhaseManagerInterface&        phasemanager);

  //! evaluate matrixes (stiffness)
  virtual void EvaluateMatrix(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         iter,           //!< iteration index (gauss point or node)
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor times time-integration factor
      double                                                      fac             //!< domain-integration factor
    ) = 0;

  //! evaluate vectors (RHS vector)
  virtual void EvaluateVector(
      std::vector<Epetra_SerialDenseVector*>&                     elevec,           //!< element vector to be filled
      const LINALG::Matrix<nen,1>&                                funct,           //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,           //! array for shape function derivatives w.r.t x,y,z
      int                                                         iter,            //!< iteration index (gauss point or node)
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac              //!< domain-integration factor
      )=0;

  //! evaluate off-diagonal coupling matrix with structure
  virtual void EvaluateMatrixOD(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              deriv,          //! array for shape function derivatives w.r.t r,s,t
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      const LINALG::Matrix<nsd,nsd>&                              xjm,
      int                                                         iter,           //!< iteration index (gauss point or node)
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor times time-integration factor
      double                                                      fac,            //!< domain-integration factor
      double                                                      det
  ) = 0;

};

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*!
\brief Evaluator class for evaluation of element matrix/vectors for multiple phases

This class wraps multiple evaluators. For evaluation of matrix and vector, it just loops over all
single evaluators.

\author vuong
*/
template<int nsd, int nen>
class MultiEvaluator : public EvaluatorInterface<nsd,nen>
{
public:
  //! constructor
  MultiEvaluator(){};

  //! destructor
  virtual ~MultiEvaluator(){};

  //! evaluate matrixes (stiffness)
  void EvaluateMatrix(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         iter,            //!< iteration index (gauss point or node)
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac             //!< domain-integration factor times time-integration factor
    )
  {
    // loop over the evaluators
    typename std::vector<Teuchos::RCP<EvaluatorInterface<nsd,nen> > >::iterator it;
    for(it=evaluators_.begin();it!=evaluators_.end();it++)
      (*it)->EvaluateMatrix(
            elemat,
            funct,
            derxy,
            iter,
            numdofpernode,
            phasemanager,
            variablemanager,
            timefacfac,
            fac);
  };

  //! evaluate vectors (RHS vector)
  void EvaluateVector(
      std::vector<Epetra_SerialDenseVector*>&                     elevec,           //!< element vector to be filled
      const LINALG::Matrix<nen,1>&                                funct,           //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,           //! array for shape function derivatives w.r.t x,y,z
      int                                                         iter,            //!< iteration index (gauss point or node)
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac              //!< domain-integration factor
      )
  {
    // loop over the evaluators
    typename std::vector<Teuchos::RCP<EvaluatorInterface<nsd,nen> > >::iterator it;
    for(it=evaluators_.begin();it!=evaluators_.end();it++)
      (*it)->EvaluateVector(
          elevec,
          funct,
          derxy,
          iter,
          numdofpernode,
          phasemanager,
          variablemanager,
          rhsfac,
          fac);
  };

  //! evaluate off-diagonal coupling matrix with structure
  void EvaluateMatrixOD(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              deriv,          //! array for shape function derivatives w.r.t r,s,t
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      const LINALG::Matrix<nsd,nsd>&                              xjm,
      int                                                         iter,            //!< iteration index (gauss point or node)
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      double                                                      det
    )
  {
    // loop over the evaluators
    typename std::vector<Teuchos::RCP<EvaluatorInterface<nsd,nen> > >::iterator it;
    for(it=evaluators_.begin();it!=evaluators_.end();it++)
    {
      (*it)->EvaluateMatrixOD(
            elemat,
            funct,
            deriv,
            derxy,
            xjm,
            iter,
            numdofpernode,
            phasemanager,
            variablemanager,
            timefacfac,
            fac,
            det);
    }
  };

  //! add an evaluator to the list of evaluators
  void AddEvaluator(Teuchos::RCP<EvaluatorInterface<nsd,nen> > evaluator)
  {evaluators_.push_back(evaluator);};

private:
  //! list of all evaluators
  std::vector<Teuchos::RCP<EvaluatorInterface<nsd,nen> > > evaluators_;

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
template<int nsd, int nen>
class EvaluatorBase : public EvaluatorInterface<nsd,nen>
{
public:
  //! constructor
  EvaluatorBase(Teuchos::RCP<AssembleInterface> assembler,int curphase)
  :  assembler_(assembler),
     myphase_(curphase)
  {};

  //! destructor
  virtual ~EvaluatorBase(){};

  //! evaluate matrixes (stiffness)
  void EvaluateMatrix(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         iter,           //!< iteration index (gauss point or node)
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac            //!< domain-integration factor times time-integration factor
    )
  {
    // the assembler class decides, where the terms are assembled into
    for(int iassemble=0; iassemble<assembler_->NumPhasesToAssembleInto(); iassemble++)
    {
      // call the actual evaluation and assembly of the respective term (defined by derived class)
       EvaluateMatrixAndAssemble(
          elemat,
          funct,
          derxy,
          myphase_,
          assembler_->PhaseToAssembleInto(iassemble,numdofpernode),
          numdofpernode,
          phasemanager,
          variablemanager,
          timefacfac,
          fac,
          assembler_->CalcInitTimeDeriv()
        );
    }
  };

  //! evaluate vectors (RHS vector)
  void EvaluateVector(
      std::vector<Epetra_SerialDenseVector*>&                     elevec,           //!< element vector to be filled
      const LINALG::Matrix<nen,1>&                                funct,           //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,           //! array for shape function derivatives w.r.t x,y,z
      int                                                         iter,            //!< iteration index (gauss point or node)
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac             //!< domain-integration factor
      )
  {
    // the assembler class decides, where the terms are assembled into
    for(int iassemble=0; iassemble<assembler_->NumPhasesToAssembleInto(); iassemble++)
    {
      // call the actual evaluation and assembly of the respective term (defined by derived class)
      EvaluateVectorAndAssemble(
          elevec,
          funct,
          derxy,
          myphase_,
          assembler_->PhaseToAssembleInto(iassemble,numdofpernode),
          numdofpernode,
          phasemanager,
          variablemanager,
          rhsfac,
          fac,
          assembler_->CalcInitTimeDeriv()
        );
    }
  };

  //! evaluate off-diagonal coupling matrix with structure
  void EvaluateMatrixOD(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              deriv,          //! array for shape function derivatives w.r.t r,s,t
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      const LINALG::Matrix<nsd,nsd>&                              xjm,
      int                                                         iter,           //!< iteration index (gauss point or node)
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,           //!< domain-integration factor times time-integration factor
      double                                                      det
    )
  {
    // the assembler class decides, where the terms are assembled into
    for(int iassemble=0; iassemble<assembler_->NumPhasesToAssembleInto(); iassemble++)
    {
      // call the actual evaluation and assembly of the respective term (defined by derived class)
       EvaluateMatrixODAndAssemble(
          elemat,
          funct,
          deriv,
          derxy,
          xjm,
          myphase_,
          assembler_->PhaseToAssembleInto(iassemble,numdofpernode),
          numdofpernode,
          phasemanager,
          variablemanager,
          timefacfac,
          fac,
          det
        );
    }
  };

protected:

  // actual evaluation and assembly of the respective term in the stiffness matrix (defined by derived class)
  virtual void EvaluateMatrixAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of phase to add into
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
    )=0;

  // actual evaluation and assembly of the respective term in the RHS vector (defined by derived class)
  virtual void EvaluateVectorAndAssemble(
      std::vector<Epetra_SerialDenseVector*>&                     elevec,           //!< element vector to be filled
      const LINALG::Matrix<nen,1>&                                funct,           //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,           //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,        //!< index of current phase
      int                                                         phasetoadd,     //!< index of phase to add into
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac,             //!< domain-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
      )=0;

  // actual evaluation and assembly of the respective term in the off-diagonal matrix (defined by derived class)
  virtual void EvaluateMatrixODAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              deriv,          //! array for shape function derivatives w.r.t r,s,t
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      const LINALG::Matrix<nsd,nsd>&                              xjm,
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of phase to add into
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      double                                                      det
      )=0;

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
template<int nsd, int nen>
class EvaluatorConv : public EvaluatorBase<nsd,nen>
{
public:
  //! constructor
  EvaluatorConv(Teuchos::RCP<AssembleInterface> assembler,int curphase):
    EvaluatorBase<nsd,nen>(assembler,curphase)
  {};

  //! destructor
  virtual ~EvaluatorConv(){};

protected:

  //! evaluate element matrix
  void EvaluateMatrixAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
    );

  //! evaluate element RHS vector
  void EvaluateVectorAndAssemble(
      std::vector<Epetra_SerialDenseVector*>&                     elevec,           //!< element vector to be filled
      const LINALG::Matrix<nen,1>&                                funct,           //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,           //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,        //!< index of current phase
      int                                                         phasetoadd,      //!< index of current phase
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac,             //!< domain-integration factor
      bool                                                        inittimederiv    //!< calculate only parts for initial time derivative
      );

  //! evaluate off-diagonal coupling matrix with structure
  void EvaluateMatrixODAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              deriv,          //! array for shape function derivatives w.r.t r,s,t
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      const LINALG::Matrix<nsd,nsd>&                              xjm,
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      double                                                      det
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
template<int nsd, int nen>
class EvaluatorDivVel : public EvaluatorBase<nsd,nen>
{
public:
  //! constructor
  EvaluatorDivVel(Teuchos::RCP<AssembleInterface> assembler,int curphase):
    EvaluatorBase<nsd,nen>(assembler,curphase)
  {};

  //! destructor
  virtual ~EvaluatorDivVel(){};

protected:

  //! evaluate element matrix
  void EvaluateMatrixAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
    );

  //! evaluate element RHS vector
  void EvaluateVectorAndAssemble(
      std::vector<Epetra_SerialDenseVector*>&                     elevec,           //!< element vector to be filled
      const LINALG::Matrix<nen,1>&                                funct,           //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,           //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,        //!< index of current phase
      int                                                         phasetoadd,      //!< index of current phase
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac,             //!< domain-integration factor
      bool                                                        inittimederiv    //!< calculate only parts for initial time derivative
      );

  //! evaluate off-diagonal coupling matrix with structure
  void EvaluateMatrixODAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              deriv,          //! array for shape function derivatives w.r.t r,s,t
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      const LINALG::Matrix<nsd,nsd>&                              xjm,
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      double                                                      det
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
template<int nsd, int nen>
class EvaluatorSatDivVel : public EvaluatorDivVel<nsd,nen>
{
public:
  //! constructor
  EvaluatorSatDivVel(Teuchos::RCP<AssembleInterface> assembler,int curphase):
    EvaluatorDivVel<nsd,nen>(assembler,curphase)
  {};

  //! destructor
  virtual ~EvaluatorSatDivVel(){};

protected:

  //! evaluate element matrix
  void EvaluateMatrixAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
    );

  //! evaluate element RHS vector
  void EvaluateVectorAndAssemble(
      std::vector<Epetra_SerialDenseVector*>&                     elevec,           //!< element vector to be filled
      const LINALG::Matrix<nen,1>&                                funct,           //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,           //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,        //!< index of current phase
      int                                                         phasetoadd,      //!< index of current phase
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac,             //!< domain-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
      );

  //! evaluate off-diagonal coupling matrix with structure
  void EvaluateMatrixODAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              deriv,          //! array for shape function derivatives w.r.t r,s,t
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      const LINALG::Matrix<nsd,nsd>&                              xjm,
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      double                                                      det
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
template<int nsd, int nen>
class EvaluatorDiff : public EvaluatorBase<nsd,nen>
{
public:
  //! constructor
  EvaluatorDiff(Teuchos::RCP<AssembleInterface> assembler,int curphase):
    EvaluatorBase<nsd,nen>(assembler,curphase)
  {};

  //! destructor
  virtual ~EvaluatorDiff(){};

protected:

  //! evaluate element matrix
  void EvaluateMatrixAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
  );

  //! evaluate element RHS vector
  void EvaluateVectorAndAssemble(
      std::vector<Epetra_SerialDenseVector*>&                     elevec,           //!< element vector to be filled
      const LINALG::Matrix<nen,1>&                                funct,           //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,           //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,        //!< index of current phase
      int                                                         phasetoadd,      //!< index of current phase
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac,             //!< domain-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
      );

  //! evaluate off-diagonal coupling matrix with structure
  void EvaluateMatrixODAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              deriv,          //! array for shape function derivatives w.r.t r,s,t
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      const LINALG::Matrix<nsd,nsd>&                              xjm,
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      double                                                      det
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
template<int nsd, int nen>
class EvaluatorReac : public EvaluatorBase<nsd,nen>
{
public:
  //! constructor
  EvaluatorReac(Teuchos::RCP<AssembleInterface> assembler,int curphase):
    EvaluatorBase<nsd,nen>(assembler,curphase)
  {};

  //! destructor
  virtual ~EvaluatorReac(){};

protected:

  //! evaluate element matrix
  void EvaluateMatrixAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
    );

  //! evaluate element RHS vector
  void EvaluateVectorAndAssemble(
      std::vector<Epetra_SerialDenseVector*>&                     elevec,           //!< element vector to be filled
      const LINALG::Matrix<nen,1>&                                funct,           //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,           //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,        //!< index of current phase
      int                                                         phasetoadd,      //!< index of current phase
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac,             //!< domain-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
      );

  //! evaluate off-diagonal coupling matrix with structure
  void EvaluateMatrixODAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              deriv,          //! array for shape function derivatives w.r.t r,s,t
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      const LINALG::Matrix<nsd,nsd>&                              xjm,
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      double                                                      det
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
template<int nsd, int nen>
class EvaluatorMassPressure : public EvaluatorBase<nsd,nen>
{
public:
  //! constructor
  EvaluatorMassPressure(Teuchos::RCP<AssembleInterface> assembler,int curphase):
    EvaluatorBase<nsd,nen>(assembler,curphase)
  {};

  //! destructor
  virtual ~EvaluatorMassPressure(){};

protected:

  //! evaluate element matrix
  void EvaluateMatrixAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
  );

  //! evaluate element RHS vector
  void EvaluateVectorAndAssemble(
      std::vector<Epetra_SerialDenseVector*>&                     elevec,           //!< element vector to be filled
      const LINALG::Matrix<nen,1>&                                funct,           //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,           //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,        //!< index of current phase
      int                                                         phasetoadd,      //!< index of current phase
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac,             //!< domain-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
      );

  //! evaluate off-diagonal coupling matrix with structure
  void EvaluateMatrixODAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              deriv,          //! array for shape function derivatives w.r.t r,s,t
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      const LINALG::Matrix<nsd,nsd>&                              xjm,
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      double                                                      det
      );

  //! get transient term for rhs and OD
  double GetRhsTrans(
      int                                                         curphase,        //!< index of current phase
      int                                                         phasetoadd,      //!< index of current phase
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac              //!< domain-integration factor);
      );

};

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*!
\brief class for evaluation of mass term (solid pressure) into the element matrix

This class implements the term $( w, \frac{(1-\porosity) }{K_s} \frac{\partial p^s}{\partial t} )$.

\author vuong
*/
template<int nsd, int nen>
class EvaluatorMassSolidPressure : public EvaluatorBase<nsd,nen>
{
public:
  //! constructor
  EvaluatorMassSolidPressure(Teuchos::RCP<AssembleInterface> assembler,int curphase):
    EvaluatorBase<nsd,nen>(assembler,curphase)
  {};

  //! destructor
  virtual ~EvaluatorMassSolidPressure(){};

protected:

  //! evaluate element matrix
  void EvaluateMatrixAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
    );

  //! evaluate element RHS vector
  void EvaluateVectorAndAssemble(
      std::vector<Epetra_SerialDenseVector*>&                     elevec,           //!< element vector to be filled
      const LINALG::Matrix<nen,1>&                                funct,           //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,           //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,        //!< index of current phase
      int                                                         phasetoadd,      //!< index of current phase
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac,             //!< domain-integration factor
      bool                                                        inittimederiv    //!< calculate only parts for initial time derivative
      );

  //! evaluate off-diagonal coupling matrix with structure
  void EvaluateMatrixODAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              deriv,          //! array for shape function derivatives w.r.t r,s,t
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      const LINALG::Matrix<nsd,nsd>&                              xjm,
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      double                                                      det
      );

  //! get transient term for rhs and OD
  double GetRhsTrans(
      int                                                         curphase,        //!< index of current phase
      int                                                         phasetoadd,      //!< index of current phase
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac              //!< domain-integration factor);
      );

};


/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*!
\brief class for evaluation of mass term (solid pressure), scaled with saturation into the element matrix

This class implements the term $( w, S\frac{(1-\porosity) }{K_s} \frac{\partial p^s}{\partial t} )$.

\author vuong
*/
template<int nsd, int nen>
class EvaluatorMassSolidPressureSat : public EvaluatorMassSolidPressure<nsd,nen>
{
public:
  //! constructor
  EvaluatorMassSolidPressureSat(Teuchos::RCP<AssembleInterface> assembler,int curphase):
    EvaluatorMassSolidPressure<nsd,nen>(assembler,curphase)
  {};

  //! destructor
  virtual ~EvaluatorMassSolidPressureSat(){};

protected:

  //! evaluate element matrix
  void EvaluateMatrixAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
    );

  //! evaluate element RHS vector
  void EvaluateVectorAndAssemble(
      std::vector<Epetra_SerialDenseVector*>&                     elevec,           //!< element vector to be filled
      const LINALG::Matrix<nen,1>&                                funct,           //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,           //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,        //!< index of current phase
      int                                                         phasetoadd,      //!< index of current phase
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac,             //!< domain-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
      );

  //! evaluate off-diagonal coupling matrix with structure
  void EvaluateMatrixODAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              deriv,          //! array for shape function derivatives w.r.t r,s,t
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      const LINALG::Matrix<nsd,nsd>&                              xjm,
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      double                                                      det
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
template<int nsd, int nen>
class EvaluatorMassSaturation : public EvaluatorBase<nsd,nen>
{
public:
  //! constructor
  EvaluatorMassSaturation(Teuchos::RCP<AssembleInterface> assembler,int curphase):
    EvaluatorBase<nsd,nen>(assembler,curphase)
  {};

  //! destructor
  virtual ~EvaluatorMassSaturation(){};

protected:

  //! evaluate element matrix
  void EvaluateMatrixAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
    );

  //! evaluate element RHS vector
  void EvaluateVectorAndAssemble(
      std::vector<Epetra_SerialDenseVector*>&                     elevec,           //!< element vector to be filled
      const LINALG::Matrix<nen,1>&                                funct,           //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,           //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,        //!< index of current phase
      int                                                         phasetoadd,      //!< index of current phase
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac,             //!< domain-integration factor
      bool                                                        inittimederiv    //!< calculate only parts for initial time derivative
      );

  //! evaluate off-diagonal coupling matrix with structure
  void EvaluateMatrixODAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              deriv,          //! array for shape function derivatives w.r.t r,s,t
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      const LINALG::Matrix<nsd,nsd>&                              xjm,
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      double                                                      det
      );

  //! get transient term for rhs and OD
  double GetRhsTrans(
      int                                                         curphase,        //!< index of current phase
      int                                                         phasetoadd,      //!< index of current phase
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac              //!< domain-integration factor);
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
template<int nsd, int nen>
class EvaluatorPressureAndSaturation : public EvaluatorBase<nsd,nen>
{
public:
  //! constructor
  EvaluatorPressureAndSaturation(Teuchos::RCP<AssembleInterface> assembler,int curphase)
   : EvaluatorBase<nsd,nen>(assembler,curphase)
     {};

  //! destructor
  virtual ~EvaluatorPressureAndSaturation(){};

protected:

  //! evaluate element matrix
  void EvaluateMatrixAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
    );

  //! evaluate element RHS vector
  void EvaluateVectorAndAssemble(
      std::vector<Epetra_SerialDenseVector*>&                     elevec,           //!< element vector to be filled
      const LINALG::Matrix<nen,1>&                                funct,           //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,           //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,        //!< index of current phase
      int                                                         phasetoadd,      //!< index of current phase
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac,             //!< domain-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
      );

  //! evaluate off-diagonal coupling matrix with structure
  void EvaluateMatrixODAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              deriv,          //! array for shape function derivatives w.r.t r,s,t
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      const LINALG::Matrix<nsd,nsd>&                              xjm,
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      double                                                      det
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
template<int nsd, int nen>
class EvaluatorSolidPressure : public EvaluatorBase<nsd,nen>
{
public:
  //! constructor
  EvaluatorSolidPressure(Teuchos::RCP<AssembleInterface> assembler,int curphase)
  : EvaluatorBase<nsd,nen>(assembler,curphase)
    {};

  //! destructor
  virtual ~EvaluatorSolidPressure(){};

protected:

  //! evaluate element matrix
  void EvaluateMatrixAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
    );

  //! evaluate element RHS vector
  void EvaluateVectorAndAssemble(
      std::vector<Epetra_SerialDenseVector*>&                     elevec,           //!< element vector to be filled
      const LINALG::Matrix<nen,1>&                                funct,           //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,           //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,        //!< index of current phase
      int                                                         phasetoadd,      //!< index of current phase
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac,             //!< domain-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
      );

  //! evaluate off-diagonal coupling matrix with structure
  void EvaluateMatrixODAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              deriv,          //! array for shape function derivatives w.r.t r,s,t
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      const LINALG::Matrix<nsd,nsd>&                              xjm,
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      double                                                      det
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
template<int nsd, int nen>
class EvaluatorPorosity : public EvaluatorBase<nsd,nen>
{
public:
  //! constructor
  EvaluatorPorosity(Teuchos::RCP<AssembleInterface> assembler,int curphase)
  : EvaluatorBase<nsd,nen>(assembler,curphase)
    {};

  //! destructor
  virtual ~EvaluatorPorosity(){};

protected:

  //! evaluate element matrix
  void EvaluateMatrixAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
    );

  //! evaluate element RHS vector
  void EvaluateVectorAndAssemble(
      std::vector<Epetra_SerialDenseVector*>&                     elevec,           //!< element vector to be filled
      const LINALG::Matrix<nen,1>&                                funct,           //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,           //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,        //!< index of current phase
      int                                                         phasetoadd,      //!< index of current phase
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac,             //!< domain-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
      );

  //! evaluate off-diagonal coupling matrix with structure
  void EvaluateMatrixODAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              deriv,          //! array for shape function derivatives w.r.t r,s,t
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      const LINALG::Matrix<nsd,nsd>&                              xjm,
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      double                                                      det
      );

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
template<int nsd, int nen>
class ReconstructFluxLinearization : public EvaluatorBase<nsd,nen>
{
public:
  //! constructor
  ReconstructFluxLinearization(Teuchos::RCP<AssembleInterface> assembler,int curphase)
  : EvaluatorBase<nsd,nen>(assembler,curphase)
    {};

  //! destructor
  virtual ~ReconstructFluxLinearization(){};

protected:

  //! evaluate element matrix
  void EvaluateMatrixAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
    );

  //! evaluate element RHS vector
  void EvaluateVectorAndAssemble(
      std::vector<Epetra_SerialDenseVector*>&                     elevec,           //!< element vector to be filled
      const LINALG::Matrix<nen,1>&                                funct,           //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,           //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,        //!< index of current phase
      int                                                         phasetoadd,      //!< index of current phase
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac,             //!< domain-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
      );

  //! evaluate off-diagonal coupling matrix with structure
  void EvaluateMatrixODAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              deriv,          //! array for shape function derivatives w.r.t r,s,t
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      const LINALG::Matrix<nsd,nsd>&                              xjm,
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      double                                                      det
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
template<int nsd, int nen>
class ReconstructFluxRHS : public EvaluatorBase<nsd,nen>
{
public:
  //! constructor
  ReconstructFluxRHS(Teuchos::RCP<AssembleInterface> assembler,int curphase)
  : EvaluatorBase<nsd,nen>(assembler,curphase)
    {};

  //! destructor
  virtual ~ReconstructFluxRHS(){};

protected:

  //! evaluate element matrix
  void EvaluateMatrixAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
    );

  //! evaluate element RHS vector
  void EvaluateVectorAndAssemble(
      std::vector<Epetra_SerialDenseVector*>&                     elevec,           //!< element vector to be filled
      const LINALG::Matrix<nen,1>&                                funct,           //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,           //! array for shape function derivatives w.r.t x,y,z
      int                                                         curphase,        //!< index of current phase
      int                                                         phasetoadd,      //!< index of current phase
      int                                                         numdofpernode,   //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,    //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager, //!< variable manager
      double                                                      rhsfac,          //!< time-integration factor for rhs times domain-integration factor
      double                                                      fac,             //!< domain-integration factor
      bool                                                        inittimederiv   //!< calculate only parts for initial time derivative
      );

  //! evaluate off-diagonal coupling matrix with structure
  void EvaluateMatrixODAndAssemble(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              deriv,          //! array for shape function derivatives w.r.t r,s,t
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      const LINALG::Matrix<nsd,nsd>&                              xjm,
      int                                                         curphase,       //!< index of current phase
      int                                                         phasetoadd,     //!< index of current phase
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac,            //!< domain-integration factor times time-integration factor
      double                                                      det
      );
};



} //namespace POROFLUIDEVALUATOR

} // namespace ELEMENTS
} // namespace DRT


#endif /* SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUID_EVALUATOR_H_ */
