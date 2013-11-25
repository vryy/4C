/*----------------------------------------------------------------------*/
/*!
\file ad_str_windkessel_merged.H

\brief Adapter Layer for Structures with 0D Windkessel coupling

<pre>
Maintainer: Marc Hirschvogel
            hirschvogel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15270
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */


#ifndef AD_STR_WINDKESSEL_MERGED_H
#define AD_STR_WINDKESSEL_MERGED_H

/*----------------------------------------------------------------------*/
/* headers */
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>

//#include "ad_str_fsiwrapper.H"
#include "ad_str_wrapper.H"

namespace DRT
{
  class Condition;
}

// forward declarations
namespace STR
{
  namespace AUX
  {
    class MapExtractor;
  }
}

namespace LINALG
{
  class MapExtractor;
}

/*----------------------------------------------------------------------*/
/* adapting adapter */
namespace ADAPTER {

  /*====================================================================*/
  /*!
   * \brief Adapter to Windkessel structural time integration.
   * This class wraps one of the standard adapters for structural time
   * integration. The results are modified and/or merged to account for the
   * additional degrees of freedom of the Windkessel pressures.
   *
   * \date 10/13
   */
  class StructureWindkesselMerged : public StructureWrapper
  {
  public:

    /// Constructor
    StructureWindkesselMerged(Teuchos::RCP<Structure> stru);

    void SetPressure(Teuchos::RCP<Epetra_Vector> couppres);

    /// initial guess of Newton's method
    Teuchos::RCP<const Epetra_Vector> InitialGuess();

    /// right-hand-side of Newton's method
    Teuchos::RCP<const Epetra_Vector> RHS();

    /// unknown displacements at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> Dispnp();

    /// known displacements at \f$t_{n}\f$
    Teuchos::RCP<const Epetra_Vector> Dispn();

    /// dof map of vector of unknowns
    Teuchos::RCP<const Epetra_Map> DofRowMap();


    /// direct access to system matrix
    Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

    /// direct access to system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();

    /// External force at t_{n+1}
//    Teuchos::RCP<const Epetra_Vector> FExtn();

    /// update displacement and evaluate elements
    void Evaluate(
      Teuchos::RCP<const Epetra_Vector> dispstepinc  ///< solution increment between time step n and n+1
    );

    /// communication object at the interface
    virtual Teuchos::RCP<const STR::AUX::MapExtractor> Interface() const
    {
      return interface_;
    }

    //! Return MapExtractor for Dirichlet boundary conditions
    virtual const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
    {
      return structure_->GetDBCMapExtractor();
    };

    /// domain map of system matrix (do we really need this?)
    virtual const Epetra_Map& DomainMap();

    /// are there any Windkessel boundary conditions?
     virtual bool HaveWindkessel()
     {
       return structure_->HaveWindkessel();
     };

    /// Return bool indicating if Windkessel functions are defined
    virtual const Teuchos::RCP< UTILS::WindkesselManager> GetWindkesselManager()
    {
      return structure_->GetWindkesselManager();
    };

    virtual INPAR::STR::STC_Scale GetSTCAlgo()
    { return structure_->GetSTCAlgo(); };

    virtual Teuchos::RCP<LINALG::SparseMatrix> GetSTCMat()
    {
      dserror("FSI with merged structural constraints does not work in combination with STC!");
      return structure_->GetSTCMat();
    }

    //! Update iteration
    //! Add residual increment to pressures stored in Windkessel manager
    virtual void UpdateIterIncrWindkessel
    (
        Teuchos::RCP<Epetra_Vector> presincr ///< pressure increment
    )
    {
      structure_->UpdateIterIncrWindkessel(presincr);
    }

    /// @name Apply interface forces

    //@}

    /// Integrate from t1 to t2
    virtual int Integrate() { return structure_->Integrate(); }

    /// map between coupling ID and conditions on structure
    std::map<int,DRT::Condition*> coupcond_;

  private:

    /// the interface map setup for interface <-> full translation
    Teuchos::RCP<STR::AUX::MapExtractor> interface_;

    /// the Windkessel map setup for full <-> stuct+windkessel transition
    Teuchos::RCP<LINALG::MapExtractor> windkmerger_;

    /// the complete non-overlapping degree of freedom row map for structure and Windkessel pressures
    Teuchos::RCP<Epetra_Map> dofrowmap_;

    /// @name local copies of input parameters
    //{@
    Teuchos::RCP<DRT::Discretization> discret_;  ///< the discretisation
    Teuchos::RCP<Teuchos::ParameterList> ioparams_;  ///< I/O flags ... not sure of really needed
    Teuchos::RCP<Teuchos::ParameterList> sdynparams_;  ///< dynamic control flags ... used,
                                                       ///< but could/should be circumvented
    Teuchos::RCP<Teuchos::ParameterList> xparams_;  ///< eXtra input parameters
    Teuchos::RCP<LINALG::Solver> solver_;  ///< the solver
    Teuchos::RCP<IO::DiscretizationWriter> output_;  ///< the output writer

    //@}

  };  // class StructureWindkessel

}  // namespace ADAPTER

/*----------------------------------------------------------------------*/
#endif
