/*----------------------------------------------------------------------*/
/*!
 \file porofluid_evaluator.H

 \brief helper class for evaluation of the governing equation of multiphase porou flow

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
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

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
//! general interface to helper class for assembly into element matrix
class AssembleInterface
{
public:
  //! constructor
  AssembleInterface(){};

  //! destructor
  virtual ~AssembleInterface(){};

  virtual int NumPhasesToAssembleInto() const =0;

  virtual int PhaseToAssembleInto(int iassemble, int numdofpernode) const =0;

};

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
//! helper class for standard assembly into element matrix
class AssembleStandard : public AssembleInterface
{
public:
  //! constructor
  AssembleStandard(int curphase):
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
  AssembleAlsoIntoOtherPhase(int curphase ,int otherphase):
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
//! general interface class for evaluation of element matrix
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

  virtual void EvaluateMatrix(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         iter,           //!< iteration index (gauss point or node)
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac             //!< domain-integration factor times time-integration factor
    ) = 0;

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

};

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
//! helper class for evaluation of element matrix/vectors for multiple phases
template<int nsd, int nen>
class MultiEvaluator : public EvaluatorInterface<nsd,nen>
{
public:
  //! constructor
  MultiEvaluator(){};

  //! destructor
  virtual ~MultiEvaluator(){};

  //! evaluate matrixes
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

  //! evaluate vectors
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

  void AddEvaluator(Teuchos::RCP<EvaluatorInterface<nsd,nen> > evaluator)
  {evaluators_.push_back(evaluator);};

private:
  std::vector<Teuchos::RCP<EvaluatorInterface<nsd,nen> > > evaluators_;

};


/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
//! general base class for evaluation of one term
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


  void EvaluateMatrix(
      std::vector<Epetra_SerialDenseMatrix*>&                     elemat,         //!< element matrix to be filled
      const LINALG::Matrix<nen,1>&                                funct,          //! array for shape functions
      const LINALG::Matrix<nsd,nen>&                              derxy,          //! array for shape function derivatives w.r.t x,y,z
      int                                                         iter,           //!< iteration index (gauss point or node)
      int                                                         numdofpernode,  //!< total number of DOFs/phases
      const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,   //!< phase manager
      const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,//!< variable manager
      double                                                      timefacfac,     //!< domain-integration factor
      double                                                      fac             //!< domain-integration factor times time-integration factor
    )
  {
    for(int iassemble=0; iassemble<assembler_->NumPhasesToAssembleInto(); iassemble++)
    {
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
          fac
        );
    }
  };

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
    for(int iassemble=0; iassemble<assembler_->NumPhasesToAssembleInto(); iassemble++)
    {
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
          fac
        );
    }
  };

protected:

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
      double                                                      fac             //!< domain-integration factor times time-integration factor
    )=0;

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
      double                                                      fac              //!< domain-integration factor
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
//! class for evaluation of convective term into the element matrix
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
      double                                                      fac             //!< domain-integration factor times time-integration factor
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
      double                                                      fac              //!< domain-integration factor
      );

};

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
//! class for evaluation of divergence of the (mesh) velocity field
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
      double                                                      fac             //!< domain-integration factor times time-integration factor
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
      double                                                      fac              //!< domain-integration factor
      );

};

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
//! class for evaluation of divergence of the (mesh) velocity field, scaled by saturation
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
      double                                                      fac             //!< domain-integration factor times time-integration factor
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
      double                                                      fac              //!< domain-integration factor
      );

};

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
//! class for evaluation of diffusive term into the element matrix
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
      double                                                      fac             //!< domain-integration factor times time-integration factor
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
      double                                                      fac              //!< domain-integration factor
      );
};

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
//! class for evaluation of reactive term into the element matrix
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
      double                                                      fac             //!< domain-integration factor times time-integration factor
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
      double                                                      fac              //!< domain-integration factor
      );

};

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
//! class for evaluation of mass term (pressure) into the element matrix
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
      double                                                      fac             //!< domain-integration factor times time-integration factor
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
      double                                                      fac              //!< domain-integration factor
      );

};

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
//! class for evaluation of mass term (solid pressure) into the element matrix
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
      double                                                      fac             //!< domain-integration factor times time-integration factor
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
      double                                                      fac              //!< domain-integration factor
      );

};


/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
//! class for evaluation of mass term (solid pressure) into the element matrix
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
      double                                                      fac             //!< domain-integration factor times time-integration factor
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
      double                                                      fac              //!< domain-integration factor
      );

};

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
//! class for evaluation of mass term (solid saturation) into the element matrix
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
      double                                                      fac             //!< domain-integration factor times time-integration factor
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
      double                                                      fac              //!< domain-integration factor
      );

};


/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
//! helper class for evaluation of pressure and saturation
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
      double                                                      fac             //!< domain-integration factor times time-integration factor
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
      double                                                      fac              //!< domain-integration factor
      );

};

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
//! helper class for evaluation of pressure and saturation
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
      double                                                      fac             //!< domain-integration factor times time-integration factor
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
      double                                                      fac              //!< domain-integration factor
      );

};

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
//! helper class for reconstruction of flux
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
      double                                                      fac             //!< domain-integration factor times time-integration factor
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
      double                                                      fac              //!< domain-integration factor
      );

};

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
//! helper class for reconstruction of flux
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
      double                                                      fac             //!< domain-integration factor times time-integration factor
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
      double                                                      fac              //!< domain-integration factor
      );
};


} //namespace POROFLUIDEVALUATOR

} // namespace ELEMENTS
} // namespace DRT


#endif /* SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUID_EVALUATOR_H_ */
