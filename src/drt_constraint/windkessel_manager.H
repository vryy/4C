/*!----------------------------------------------------------------------
\file windkessel_manager.H

\brief Class controlling Windkessel functions and containing the necessary data

<pre>
Maintainer: Marc Hirschvogel
            hirschvogel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15270
</pre>

*----------------------------------------------------------------------*/

#ifndef WINDKESSEL_MANAGER_H
#define WINDKESSEL_MANAGER_H

#include "../drt_adapter/adapter_algorithmbase.H"
#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Teuchos_ParameterList.hpp>

using namespace Teuchos;

namespace ADAPTER
{
  class StructureWindkesselMerged;
}
namespace DRT
{
  class Condition;
}

namespace DRT
{
  class Discretization;
}

namespace LINALG
{
	class SparseMatrix;
  class SparseOperator;
  class MultiMapExtractor;
}

namespace UTILS
{

  //forward declarations
  class Windkessel;
  class WindkesselDofSet;
  
  class WindkesselManager
  {
  public:

    /*!
      \brief Constructor of Windkessel manager, allocating the Windkessels
    */
    WindkesselManager
    (
      RCP<DRT::Discretization> disc,  ///< standard discretization
      RCP<Epetra_Vector> disp,  ///< current displacement
      ParameterList param  ///<  parameterlist from time integration algorithm
    );
    /*!
        \brief Destructor

     */
    ~WindkesselManager()
    {
        return;
    };

    /*!
      \brief Change stiffness matrix and force vector according to the Windkessel functions.
      Values of pressure are taken from intern variable.
      Difference between current and prescribed values is calculated and stored as well.
    */
    void StiffnessAndInternalForces
    (
      const double time,  ///< time at end of time step
      RCP<Epetra_Vector> displast,  ///< displacement at beginning of time step
      RCP<Epetra_Vector> disp,  ///< displacement at end of time step
      ParameterList scalelist
    );

    /*!
     \brief Return norm of difference between actual and Windkessel values
    */
    double GetWindkesselRHSNorm() const
    {
      double foo;
      windkesselrhsm_->Norm2(&foo);
      return foo;
    };

    /*!
         \brief Return number of Windkessel functions
    */
    int GetNumberOfWindkesselFunctions() const
    {
        return numWindkesselID_;
    };

    /*!
     \brief Scale all pressures by a double d
    */
    void ScalePres
    (
      double d   ///< scale factor
    )
    {
        presn_->Scale(d);
        return;
    };

    /*!
         \brief Update Windkessel variables
    */
    void UpdateTimeStep();

    /// Add a vector as residual increment to the vector of pressures
    void UpdatePres
    (
      RCP<Epetra_Vector> presincrement  ///< vector to add
    );

    ///
    void EvaluateNeumannWindkesselCoupling
    (
      RCP<Epetra_Vector> actpres
    );

    void SetPres
    (
      Teuchos::RCP<Epetra_Vector> actpres
    );

    /*!
         \brief Return differences between prescribed and actual value of Windkessel number i
    */
    double GetWindkesselRHS
    (
      int i ///< ID of Windkessel of interest
    ) const
    {
        return (*windkesselrhsm_)[i];
    }

    /// return vector of differences between prescribed and actual values
    RCP<Epetra_Vector> GetWindkesselRHS() const
    {
      return windkesselrhsm_;
    }

    /*!
     \brief Return EpetraMap that determined distribution of Windkessel functions and pressures over processors
    */
    RCP<Epetra_Map> GetWindkesselMap() const
    {
      return windkesselmap_;
    };

    //! Return the additional rectangular matrix, constructed for pressure evaluation
    RCP<LINALG::SparseOperator> GetCoupOffdiagVolD() //const
    {
      return coupoffdiag_vol_d_;
    };

    //! Return the additional rectangular matrix, constructed for pressure evaluation
    RCP<LINALG::SparseOperator> GetCoupOffdiagFextP() //const
    {
      return coupoffdiag_fext_p_;
    };

    //! Return the additional rectangular matrix, constructed for pressure evaluation
    RCP<LINALG::SparseMatrix> GetWindkesselStiffness() //const
    {
      return windkesselstiffness_;
    };

    /*!
      \brief Return pressure for Windkessel i
    */
    double GetPres
    (
      int i ///< ID of Windkessel of interest
    ) const
    {
      return (*presn_)[i];
    };

    /*!
      \brief Return pressure vector
    */
    RCP<Epetra_Vector> GetPresVector() const
    {
      return presn_;
    };

    /*!
      \brief Return pressure of last converged step
    */
    RCP<Epetra_Vector> GetPresVectorOld() const
    {
      return pres_;
    };

    /*!
     \brief Return if there are Windkessels
    */
    bool HaveWindkessel() const
    {
      return havewindkessel_;
    };

    /*!
     \brief Return current value
    */
    double GetActVolValue
    (
      int i ///< ID of Windkessel of interest
    ) const
    {
        return (*voln_)[i];
    };

    /// Reset pressures
    void SetPresVector
    (
      RCP<Epetra_Vector> newpres  ///< new pressures
    )
    {
      presn_->Update(1.0,*newpres,0.0);
      pres_->Update(1.0,*newpres,0.0);
      return;
    }

    void PrintPresFlux() const;

    /// Return Reference base values to write restart
    RCP<Epetra_Vector> GetRefVolValue() const
    {
      return vol_;
    }
    
    //! switch Windkessel matrix to block matrix
    void UseBlockMatrix(Teuchos::RCP<const LINALG::MultiMapExtractor> domainmaps,
                        Teuchos::RCP<const LINALG::MultiMapExtractor> rangemaps);


  private:

    // don't want = operator, cctor and destructor
    WindkesselManager operator = (const WindkesselManager& old);
    WindkesselManager(const WindkesselManager& old);

    // underlying structure
    Teuchos::RCP< ::ADAPTER::StructureWindkesselMerged> structure_;
    //Teuchos::RCP<Epetra_Vector> couppres_;

    /// map between coupling ID and conditions on structure
    std::map<int,DRT::Condition*> coupcond_;

    RCP<DRT::Discretization> actdisc_; ///< discretization where elements of Windkessel boundary live in
    RCP<WindkesselDofSet> windkesseldofset_; ///< degrees of freedom of pressures
    RCP<Epetra_Map> windkesselmap_;  ///< unique map of Windkessel values
    RCP<Epetra_Map> redwindkesselmap_;   ///< fully redundant map of Windkessel values
    RCP<Epetra_Export> windkimpo_;  ///< importer for fully redundant Windkessel vector into distributed one
    RCP<Epetra_Vector> actpres_;
    RCP<Epetra_Vector> pres_; ///< pres vector at t_{n}
    RCP<Epetra_Vector> presn_; ///< pres vector at t_{n+1}
    RCP<Epetra_Vector> presm_; ///< pres vector at mid-point
    RCP<Epetra_Vector> presrate_; ///< pres rate vector at t_{n}
    RCP<Epetra_Vector> presraten_; ///< pres rate vector at t_{n+1}
    RCP<Epetra_Vector> presratem_; ///< pres rate vector at mid-point
    RCP<Epetra_Vector> flux_; ///< flux vector at t_{n}
    RCP<Epetra_Vector> fluxn_; ///< flux vector at t_{n+1}
    RCP<Epetra_Vector> fluxm_; ///< flux vector at mid-point
    RCP<Epetra_Vector> vol_; ///< flux vector at t_{n}
    RCP<Epetra_Vector> voln_; ///< flux vector at t_{n+1}
    RCP<Epetra_Vector> volm_; ///< flux vector at mid-point
    RCP<Epetra_Vector> windkesselrhsm_; ///< Windkessel rhs vector at mid-point
    RCP<Epetra_Vector> windk_resi_rhs_; ///< Windkessel rhs contribution from resistance
    RCP<Epetra_Vector> windk_comp_rhs_; ///< Windkessel rhs contribution from compliance
    RCP<Epetra_Vector> presnprint_;
    RCP<Epetra_Vector> fluxnprint_;
    int WindkesselID_;  ///< smallest Windkessel bc id
    int offsetID_;  ///< smallest Windkessel bc id
    int numWindkesselID_;  ///< number of Windkessel bcs
    std::vector<int> currentID;
    RCP<LINALG::SparseOperator> coupoffdiag_vol_d_;  ///< additional rectangular matrix
    RCP<LINALG::SparseOperator> coupoffdiag_fext_p_;  ///< additional rectangular matrix
    RCP<LINALG::SparseMatrix> windkesselstiffness_;  ///< additional rectangular matrix
    bool havewindkessel_;  ///< are there Windkessel bcs at all?
    RCP<Windkessel> rc_;
    int myrank_;

  }; //class
}
#endif /*WINDKESSELMANAGER_H*/
