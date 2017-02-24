/*----------------------------------------------------------------------------*/
/**
\file xstr_xstructure_state.H

\brief XStructure state handling ( only one single XFEM-discretization )

\maintainer Michael Hiermeier

\date Jun 28, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/
#ifndef SRC_DRT_STRUCTURE_XSTRUCTURE_XSTR_XSTRUCTURE_STATE_H_
#define SRC_DRT_STRUCTURE_XSTRUCTURE_XSTR_XSTRUCTURE_STATE_H_

#include "../drt_xfem/xfield_state.H"
#include "../drt_structure_new/str_timint_basedataglobalstate.H"

namespace XSTR {
  /// Known and supported XSTR state vectors
  enum StateVectorName
  {
    state_vec_disnp,
    state_vec_velnp,
    state_vec_accnp,
    state_vec_fintn,
    state_vec_fintnp,
    state_vec_fextn,
    state_vec_fextnp,
    state_vec_freactnp,
    state_vec_finertialnp,
    state_vec_fviscon,
    state_vec_fvisconp,
    state_vec_fstructold
  };

  inline std::string StateVectorName2String( enum StateVectorName name )
  {
    std::stringstream enum_id;
    enum_id << " ( enum = " << name << " )";

    switch ( name )
    {
      case state_vec_disnp:
        return "state_vec_disnp" + enum_id.str();
      case state_vec_velnp:
        return "state_vec_velnp" + enum_id.str();
      case state_vec_accnp:
        return "state_vec_accnp" + enum_id.str();
      case state_vec_fintn:
        return "state_vec_fintn" + enum_id.str();
      case state_vec_fintnp:
        return "state_vec_fintnp" + enum_id.str();
      case state_vec_fextn:
        return "state_vec_fextn" + enum_id.str();
      case state_vec_fextnp:
        return "state_vec_fextnp" + enum_id.str();
      case state_vec_freactnp:
        return "state_vec_freactnp" + enum_id.str();
      case state_vec_finertialnp:
        return "state_vec_finertialnp" + enum_id.str();
      case state_vec_fviscon:
        return "state_vec_fviscon" + enum_id.str();
      case state_vec_fvisconp:
        return "state_vec_fvisconp" + enum_id.str();
      case state_vec_fstructold:
        return "state_vec_fstructold" + enum_id.str();
      default:
        return "unknwon XSTR::StateVectorName" + enum_id.str();
    }
    exit( EXIT_FAILURE );
  }

  /// Known and supported XSTR state matrices
  enum StateMatrixName
  {
    state_mat_jac = 0,
    state_mat_stiff = 1,
    state_mat_mass = 2,
    state_mat_damp = 3
  };

  inline std::string StateMatrixName2String( enum StateMatrixName name )
  {
    std::stringstream enum_id;
    enum_id << " ( enum = " << name << " )";
    switch ( name )
    {
      case state_mat_jac:
        return "state_mat_jac" + enum_id.str();
      case state_mat_stiff:
        return "state_mat_stiff" + enum_id.str();
      case state_mat_mass:
        return "state_mat_mass" + enum_id.str();
      case state_mat_damp:
        return "state_mat_damp" + enum_id.str();
      default:
        return "unknwon XSTR::StateMatrixName";
    }
    exit( EXIT_FAILURE );
  }

  typedef std::map<enum XSTR::StateVectorName, Teuchos::RCP<Epetra_Vector> * > XStateVecMap;
  typedef std::map<enum XSTR::StateMatrixName, Teuchos::RCP<LINALG::SparseOperator> * > XStateMatMap;

  class XStructureState : public XFEM::XFieldState,
                          public STR::TIMINT::BaseDataGlobalState
  {
  public:
    /// constructor
    XStructureState();

    /// destructor
    virtual ~XStructureState() {};

    /// setup of member variables
    virtual void Setup();

    virtual void TransferToNewState(
        const DRT::DiscretizationInterface& new_discret,
        XFEM::XFieldState & new_state ) const;

    virtual void ResetNonStandardDofs( const DRT::DiscretizationInterface & full_discret );

    virtual bool Destroy();

    virtual void SetNewState( const XFEM::XFieldState & xstate );

    void SetNewState( const XSTR::XStructureState & xstate );

    inline XStateVecMap::iterator XStateVectorBegin()
    {
      return xstate_vectors_.begin();
    }

    inline XStateVecMap::const_iterator XStateVectorBegin() const
    {
      return xstate_vectors_.begin();
    }

    inline XStateVecMap::iterator XStateVectorEnd()
    {
      return xstate_vectors_.end();
    }

    inline XStateVecMap::const_iterator XStateVectorEnd() const
    {
      return xstate_vectors_.end();
    }

    inline XStateMatMap::iterator XStateMatrixBegin()
    {
      return xstate_matrices_.begin();
    }

    inline XStateMatMap::const_iterator XStateMatrixBegin() const
    {
      return xstate_matrices_.begin();
    }

    inline XStateMatMap::iterator XStateMatrixEnd()
    {
      return xstate_matrices_.end();
    }

    inline XStateMatMap::const_iterator XStateMatrixEnd() const
    {
      return xstate_matrices_.end();
    }

  protected:
    inline void CheckInit() const
    {
      XFEM::XFieldState::CheckInit();
      STR::TIMINT::BaseDataGlobalState::CheckInit();
    }

    inline void CheckInitSetup() const
    {
      XFEM::XFieldState::CheckInitSetup();
      STR::TIMINT::BaseDataGlobalState::CheckInitSetup();
    }

  private:

    void RegisterStateVectorsInMap();

    void RegisterStateMatricesInMap();

    Teuchos::RCP<Epetra_Map> BuildNewStandardDofRowMap(
        const DRT::DiscretizationInterface & new_discret ) const;

    void NaturalExtensionOfNonStandardDofValues(
        enum StateVectorName state_name,
        const DRT::DiscretizationInterface & new_discret,
        Epetra_Vector & new_xstate_vec ) const;

    void NaturalExtensionOfNonStandardDofValues(
        const DRT::DiscretizationInterface & new_discret,
        Epetra_Vector & new_xstate_vec ) const;

  private:
    XStateVecMap xstate_vectors_;
    XStateMatMap xstate_matrices_;
  };  // class XStructureState
} // XSTR


#endif /* SRC_DRT_STRUCTURE_XSTRUCTURE_XSTR_XSTRUCTURE_STATE_H_ */
