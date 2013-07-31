/*!----------------------------------------------------------------------
\file mortar_abstract_strategy.H

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/
#ifndef MORTAR_STRATEGY_BASE_H
#define MORTAR_STRATEGY_BASE_H

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace LINALG
{
  class MapExtractor;
  class Solver;
  class SparseOperator;
  class SparseMatrix;
}

namespace IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}

namespace MORTAR
{

// forward declarations

/*!
\brief Abstract base class for mortar solution strategies

Every specific solution algorithm (e.g. mortar contact with Lagrange multipliers or
mortar meshtying with penalty method) has to be specified in a corresponding derived
subclass defining the concrete algorithmic steps.

\author popp (popp@lnm.mw.tum.de)
*/
class StrategyBase
{
  public:

    //! @name Enums and Friends
    enum QuantityType // can be called by StoreNodalQuantities() or StoreDMtoNodes()
    {
      lmcurrent,
      lmold,
      lmupdate,
      lmuzawa,
      activeold,
      jump,
      dm,
      pentrac,
      wear
    };
    //@}

    /*!
    \brief Standard Constructor

    Creates the strategy base object and initializes all global variables.

    \param probdiscret (in): Underlying problem discretization
    \param params (in): List of meshtying/contact parameters
    \param dim (in): Global problem dimension
    \param comm (in): A communicator object
    \param alphaf (in): Midpoint for Gen-alpha time integration
    \param maxdof (in): Highest dof number in global problem

    */
    StrategyBase(DRT::Discretization& probdiscret,
                 Teuchos::ParameterList params,
                 int dim, Teuchos::RCP<Epetra_Comm> comm,
                 double alphaf, int maxdof);

    /*!
    \brief Destructor

    */
    virtual ~StrategyBase() {};

    //! @name Access methods

    /*!
    \brief Get parameter list

    */
    Teuchos::ParameterList& Params() { return scontact_; }

    /*!
    \brief Get problem dimension

    */
    const int& Dim() { return dim_; }

    /*!
    \brief Get Epetra communicator

    */
    const Epetra_Comm& Comm() const { return *comm_; }

    /*!
    \brief Get the underlying problem discretization

    */
    const DRT::Discretization& ProblemDiscret() const { return probdiscret_; }

    /*!
    \brief Get the underlying problem dof row map

    */
    const Teuchos::RCP<Epetra_Map> ProblemDofs() { return probdofs_; }

    /*!
    \brief Get the underlying problem node row map

    */
    const Teuchos::RCP<Epetra_Map> ProblemNodes() { return probnodes_; }

    /*!
    \brief Get the underlying problem dof row map

    */
    const Teuchos::RCP<Epetra_Map> ProblemElements() { return probeles_; }

    //@}

    //! @name Purely virtual functions
    
    // All these functions are defined in one or more specific derived classes,
    // such as CONTACT::ContactLagrangeStrategy or CONTACT::MeshtyingPenaltyStrategy.
    // As the base class MORTAR::StrategyBase is always called from the control routine
    // (time integrator), these functions need to be defined purely virtual here.

    virtual const Teuchos::RCP<Epetra_Map> ActiveRowNodes() = 0;
    virtual bool ActiveSetConverged() = 0;
    virtual bool ActiveSetSemiSmoothConverged() = 0;
    virtual void ApplyForceStiffCmt(Teuchos::RCP<Epetra_Vector> dis, Teuchos::RCP<LINALG::SparseOperator>& kt, Teuchos::RCP<Epetra_Vector>& f, bool predictor = false) = 0;
    virtual void CollectMapsForPreconditioner(Teuchos::RCP<Epetra_Map>& MasterDofMap, Teuchos::RCP<Epetra_Map>& SlaveDofMap, Teuchos::RCP<Epetra_Map>& InnerDofMap, Teuchos::RCP<Epetra_Map>& ActiveDofMap ) = 0;
    virtual double ConstraintNorm() = 0;
    virtual Teuchos::RCP<Epetra_Vector> ContactNorStress() = 0;
    virtual Teuchos::RCP<Epetra_Vector> ContactTanStress() = 0;
    virtual Teuchos::RCP<Epetra_Vector> ContactWear() = 0;
    virtual Teuchos::RCP<LINALG::SparseMatrix> DMatrix() = 0;
    virtual void DoReadRestart(IO::DiscretizationReader& reader, Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual void DoWriteRestart(Teuchos::RCP<Epetra_Vector>& activetoggle,Teuchos::RCP<Epetra_Vector>& sliptoggle,Teuchos::RCP<Epetra_Vector>& weightedwear,Teuchos::RCP<Epetra_Vector>& realwear) = 0;
    virtual void Evaluate(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual void EvaluateMeshtying(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff, Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual Teuchos::RCP<LINALG::SparseMatrix> EvaluateNormals(Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual void EvaluateReferenceState(int step,const Teuchos::RCP<Epetra_Vector> vec) = 0;
    virtual void EvaluateRelMov() = 0;
    virtual void EvaluateRelMovPredict() = 0;
    virtual void ForceRefConfig() = 0;
    virtual bool Friction() = 0;
    virtual void InitEvalInterface() = 0;
    virtual void InitEvalMortar() = 0;
    virtual void Initialize() = 0;
    virtual void InitializeUzawa(Teuchos::RCP<LINALG::SparseOperator>& kteff, Teuchos::RCP<Epetra_Vector>& feff) = 0;
    virtual double InitialPenalty() = 0;
    virtual void InterfaceForces(bool output = false) = 0;
    virtual double Inttime() = 0;
    virtual void Inttime_init()=0;
    virtual bool IsInContact() = 0;
    virtual Teuchos::RCP<Epetra_Vector> LagrMult() = 0;
    virtual Teuchos::RCP<Epetra_Vector> LagrMultOld() = 0;
    virtual Teuchos::RCP<Epetra_Vector> ConstrRhs() = 0;
    virtual Teuchos::RCP<Epetra_Vector> LagrMultSolveIncr() = 0;
    virtual void MeshInitialization() = 0;
    virtual Teuchos::RCP<LINALG::SparseMatrix> MMatrix() = 0;
    virtual void MortarCoupling(const Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual const int NumberOfActiveNodes() = 0;
    virtual const int NumberOfSlipNodes() = 0;
    virtual void OutputStresses() = 0;
    virtual void OutputWear() = 0;
    virtual void Print(std::ostream& os) const = 0;
    virtual void PrintActiveSet() = 0;
    virtual void Recover(Teuchos::RCP<Epetra_Vector> disi) = 0;
    virtual void RedistributeContact(Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual void RedistributeMeshtying() = 0;
    virtual void ResetActiveSet() = 0;
    virtual void ResetPenalty() = 0;
    virtual void RestrictMeshtyingZone() = 0;
    virtual void SaddlePointSolve(LINALG::Solver& solver, LINALG::Solver& fallbacksolver, Teuchos::RCP<LINALG::SparseOperator> kdd, Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<Epetra_Vector> dirichtoggle, int numiter) = 0;
    virtual void SaveReferenceState(Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual void SetState(const std::string& statename, const Teuchos::RCP<Epetra_Vector> vec) = 0;
    virtual const Teuchos::RCP<Epetra_Map> SlipRowNodes() = 0;
    virtual void StoreDirichletStatus(Teuchos::RCP<LINALG::MapExtractor> dbcmaps) = 0;
    virtual void StoreNodalQuantities(MORTAR::StrategyBase::QuantityType type) = 0;
    virtual void Update(int istep, Teuchos::RCP<Epetra_Vector> dis) = 0;
    virtual void UpdateActiveSet() = 0;
    virtual void UpdateActiveSetSemiSmooth() = 0;
    virtual void UpdateAugmentedLagrange() = 0;
    virtual void UpdateConstraintNorm(int uzawaiter = 0) = 0;
    virtual void VisualizeGmsh(const int step, const int iter = 0) = 0;
    virtual bool WasInContact() = 0;
    virtual bool WasInContactLastTimeStep() = 0;
    virtual bool Wear() = 0;
 
    //@}

  protected:

    // don't want cctor (= operator impossible anyway for abstract class)
    StrategyBase(const StrategyBase& old);

    DRT::Discretization&        probdiscret_;      // underlying problem discretization (not only interfaces)
    Teuchos::RCP<Epetra_Map>    probdofs_;         // underlying problem dof row map (not only interfaces)
    Teuchos::RCP<Epetra_Map>    probnodes_;        // underlying problem node row map (not only interfaces)
    Teuchos::RCP<Epetra_Map>    probeles_;         // underlying problem element row map (not only interfaces)

    Teuchos::RCP<Epetra_Comm>   comm_;             // communicator
    Teuchos::ParameterList      scontact_;         // containing contact input parameters
    int                         dim_;              // dimension of problem (2D or 3D)
    double                      alphaf_;           // Generalized-alpha parameter (0.0 for statics)
    bool                        parredist_;        // flag indicating parallel redistribution status
    int                         maxdof_;           // highest dof number in problem discretization

}; // class StrategyBase
} // namespace MORTAR

#endif  // #ifndef MORTAR_STRATEGY_BASE_H
