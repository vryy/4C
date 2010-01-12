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
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef MORTAR_STRATEGY_BASE_H
#define MORTAR_STRATEGY_BASE_H

#include "Teuchos_RefCountPtr.hpp"
#include "mortar_interface.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_sparsematrix.H"
#include "../drt_io/io.H"

using namespace std;
using namespace Teuchos;

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
      pentrac
    };
    //@}

    /*!
    \brief Standard Constructor

    Creates the strategy base object and initializes all global variables.

    */
    explicit StrategyBase(RCP<Epetra_Map> problemrowmap,
                          Teuchos::ParameterList params,
                          int dim, RCP<Epetra_Comm> comm, double alphaf);

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
    \brief Return row map of global problem

    */
    const RCP<Epetra_Map> ProblemRowMap() { return problemrowmap_; }
        
    //@}

    //! @name Purely virtual functions
    
    // All these functions are defined in one or more specific derived classes,
    // such as CONTACT::ContactLagrangeStrategy or CONTACT::MeshtyingPenaltyStrategy.
    // As the base class MORTAR::StrategyBase is always called from the control routine
    // (time integrator), these functions need to be defined purely virtual here.
    
    // Access methods
    
    
    
   
        
    // Evaluation methods
    
    virtual double ConstraintNorm() = 0;
    virtual void ContactForces(RCP<Epetra_Vector> fresm) = 0;
    virtual void DoReadRestart(IO::DiscretizationReader& reader, RCP<Epetra_Vector> dis) = 0;
    virtual void Evaluate(RCP<LINALG::SparseOperator>& kteff, RCP<Epetra_Vector>& feff, RCP<Epetra_Vector> dis) = 0;
    virtual void EvaluateMeshtying(RCP<LINALG::SparseOperator>& kteff, RCP<Epetra_Vector>& feff, RCP<Epetra_Vector> dis) = 0;
    virtual void InitializeUzawa(RCP<LINALG::SparseOperator>& kteff, RCP<Epetra_Vector>& feff) = 0;
    virtual double InitialPenalty() = 0;
    virtual RCP<Epetra_Vector> LagrMultOld() = 0;  
    virtual void MeshInitialization() = 0;
    virtual void MortarCoupling(const RCP<Epetra_Vector> dis) = 0;
    virtual void Print(ostream& os) const = 0;
    virtual void PrintActiveSet() = 0;
    virtual void Recover(RCP<Epetra_Vector> disi) = 0;
    virtual void ResetPenalty() = 0;
    virtual void SetState(const string& statename, const RCP<Epetra_Vector> vec) = 0;
    virtual void StoreDirichletStatus(RCP<LINALG::MapExtractor> dbcmaps) = 0;
    virtual void StoreNodalQuantities(MORTAR::StrategyBase::QuantityType type) = 0;
    virtual void Update(int istep) = 0;
    virtual void UpdateAugmentedLagrange() = 0;
    virtual void UpdateConstraintNorm(int uzawaiter = 0) = 0;
    
    //@}

  protected:

    RCP<Epetra_Map> problemrowmap_;     // global problem dof row map (not only interfaces)
    Teuchos::ParameterList scontact_;   // containing contact input parameters
    int dim_;                           // dimension of problem (2D or 3D)
    RCP<Epetra_Comm> comm_;             // communicator
    double alphaf_;                     // Generalized-alpha parameter (0.0 for statics)

}; // class StrategyBase
} // namespace MORTAR

#endif  // #ifndef MORTAR_STRATEGY_BASE_H
#endif  // #ifdef CCADISCRET
