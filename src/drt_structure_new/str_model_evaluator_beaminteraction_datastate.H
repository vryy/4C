/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_beaminteraction_datastate.H

\brief Global state data container for the structural (time)
       integration

\maintainer Jonas Eichinger

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef STR_MODELEVALUATOR_BEAMINTERACTION_DATASTATE
#define STR_MODELEVALUATOR_BEAMINTERACTION_DATASTATE

#include "../drt_lib/drt_dserror.H"

#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_timestepping/timintmstep.H"

#include <Teuchos_RCP.hpp>

// forward declaration
class Epetra_FEVector;
class Epetra_Time;
class Epetra_Comm;
namespace LINALG {
  class SparseOperator;
  class SparseMatrix;
}
namespace DRT {
  class Discretization;
namespace ELEMENTS {
  class Beam3Base;
}
}
namespace STR {
  class ModelEvaluator;
namespace MODELEVALUATOR {

  /** \brief Global state data container for the beaminteraction model
   *
   * This data container holds everything that needs to be updated each
   * iteration step
   */
  class BeamInteractionDataState
  {
  public:
    /// constructor
    BeamInteractionDataState();

    /// destructor
    virtual ~BeamInteractionDataState() {};

    /// initialize class variables
    void Init();

    /// setup of the new class variables
    void Setup(Teuchos::RCP<DRT::Discretization> const& ia_discret);

  protected:
    inline const bool& IsInit() const
    { return isinit_; };

    inline const bool& IsSetup() const
    { return issetup_; };

    inline void CheckInitSetup() const
    {
      if (!IsInit() or !IsSetup())
        dserror("Call Init() and Setup() first!");
    }

    inline void CheckInit() const
    {
      if (!IsInit())
        dserror("Init() has not been called, yet!");
    }

  public:

    /// @name General purpose algorithm members
    ///@{

    /// ID of actual processor in parallel
    int const& GetMyRank() const
    {
      CheckInitSetup();
      return myrank_;
    };

    ///@}

    /// @name search/interaction related stuff
    ///@{

    /// get extended bin to ele map
    std::map<int, std::set<int> > const& GetExtBinToEleMap() const
    {
      CheckInitSetup();
      return extbintoelemap_;
    };

    /// get mutable extended bin to ele map
    std::map<int, std::set<int> >& GetMutableExtBinToEleMap()
    {
      CheckInitSetup();
      return extbintoelemap_;
    };

    /// get extended ele to bin map
    std::map<int, std::set<int> > const& GetExtEleToBinMap() const
    {
      CheckInitSetup();
      return exteletobinmap_;
    };

    /// get extended ele to bin map
    std::set<int> const& GetExtEleToBinSet(int const i)
    {
      CheckInitSetup();
      return exteletobinmap_[i];
    };

    /// get mutable extended ele to bin map
    std::map<int, std::set<int> >& GetMutableExtEleToBinMap()
    {
      CheckInitSetup();
      return exteletobinmap_;
    };
    ///@}

    /// @name Get state variables (read only access)
    ///@{

    /// Return displacements \f$D_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> GetDisNp() const
    {
      CheckInitSetup();
      return disnp_;
    }

    /// Return displacements \f$D_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> GetDisColNp() const
    {
      CheckInitSetup();
      return discolnp_;
    }

    /// Return displacements \f$D_{n}\f$
    Teuchos::RCP<const Epetra_Vector> GetDisN() const
    {
      CheckInitSetup();
      return (*dis_)(0);
    }


    /// Return internal force \f$fint_{n}\f$
    Teuchos::RCP<const Epetra_FEVector> GetForceN() const
    {
      CheckInitSetup();
      return forcen_;
    }

    /// Return internal force \f$fint_{n+1}\f$
    Teuchos::RCP<const Epetra_FEVector> GetForceNp() const
    {
      CheckInitSetup();
      return forcenp_;
    }

    /// @name Get system matrices (read only access)
    ///@{
    /// returns the entire structural jacobian
    Teuchos::RCP<const LINALG::SparseMatrix> GetStiff() const
    {
      CheckInitSetup();
      return stiff_;
    }

    ///@}

    /// @name Get mutable state variables (read and write access)
    ///@{

    /// Return displacements \f$D_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& GetMutableDisNp()
    {
      CheckInitSetup();
      return disnp_;
    }

    /// Return displacements \f$D_{n+1}\f$
    Teuchos::RCP<Epetra_Vector>& GetMutableDisColNp()
    {
      CheckInitSetup();
      return discolnp_;
    }

    /// Return displacements \f$D_{n}\f$
    Teuchos::RCP<Epetra_Vector> GetMutableDisN()
    {
      CheckInitSetup();
      return (*dis_)(0);
    }

    /// Return multi-displacement vector \f$D_{n}, D_{n-1}, ...\f$
    Teuchos::RCP< ::TIMINT::TimIntMStep<Epetra_Vector> > GetMutableMultiDis()
    {
      CheckInitSetup();
      return dis_;
    }

    /// Return internal force \f$fint_{n}\f$
    Teuchos::RCP<Epetra_FEVector>& GetMutableForceN()
    {
      CheckInitSetup();
      return forcen_;
    }

    /// Return internal force \f$fint_{n+1}\f$
    Teuchos::RCP<Epetra_FEVector>& GetMutableForceNp()
    {
      CheckInitSetup();
      return forcenp_;
    }

    ///@}

    /// @name Get mutable system matrices
    ///@{
    /// returns the entire structural jacobian
    Teuchos::RCP<LINALG::SparseMatrix>& GetMutableStiff()
    {
      CheckInitSetup();
      return stiff_;
    }

    ///@}

  protected:
    /// @name variables for internal use only
    ///@{
    /// flag indicating if Init() has been called
    bool isinit_;

    /// flag indicating if Setup() has been called
    bool issetup_;

    ///@}

  private:

    /// @name General purpose algorithm members
    ///@{

    /// ID of actual processor in parallel
    int myrank_;

    ///@}

    /// @name search/interaction related stuff
    ///@{
    //! extended bin to ele map
    std::map<int, std::set<int> > extbintoelemap_;

    //! extended ele to bin map
    std::map<int, std::set<int> > exteletobinmap_;

     ///@}

    /// @name Global state vectors
    ///@{

    /// global displacements \f${D}_{n}, D_{n-1}, ...\f$
    Teuchos::RCP< ::TIMINT::TimIntMStep<Epetra_Vector> > dis_;

      /// global displacements \f${D}_{n+1}\f$ at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> disnp_;

    /// global displacements \f${D}_{n+1}\f$ at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> discolnp_;

    /// global internal force vector at \f$t_{n}\f$
    Teuchos::RCP<Epetra_FEVector> forcen_;

    /// global internal force vector at \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_FEVector> forcenp_;

    ///@}
    /// @name System matrices
    ///@{
    /// supposed to hold the entire jacobian (saddle point system if desired)
    Teuchos::RCP<LINALG::SparseMatrix> stiff_;

    ///@}
  };
}
}


#endif
