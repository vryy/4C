/*----------------------------------------------------------------------*/
/*!
\file thr_contact.H

<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                               mgit 06/11 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef THR_CONTACT_H
#define THR_CONTACT_H

/*----------------------------------------------------------------------*
 | headers                                                   mgit 06/11 |
 *----------------------------------------------------------------------*/
#include "thrtimint.H"
#include "thr_aux.H"

#include "thrtimint_mstep.H"
#include "../linalg/linalg_sparsematrix.H"

#include "../drt_mortar/mortar_manager_base.H"
#include "../drt_contact/meshtying_manager.H"
#include "../drt_contact/contact_manager.H"

// belongs to thermo namepsace          
namespace THR
{

  /*!
  \brief Class for all thermo contact specific issues

  \author gitterle (gitterle@lnm.mw.tum.de)
  */

  class ThermoContactMan 
  {

    public:
    
    /*!
    \brief Standard Constructor
  
    The constructor takes the structural contact manager and the 
    structural and thermal discretizations. 
    
    */
    ThermoContactMan(Teuchos::RCP<MORTAR::ManagerBase> cmtman,
                     Teuchos::RCP<DRT::Discretization> discretstruct,
                     Teuchos::RCP<DRT::Discretization> discretthermo);
  
    //! Destructor
    virtual ~ThermoContactMan() { ; }
    
    //! Return thermcondLM matrix
    Teuchos::RCP<LINALG::SparseMatrix> ThermCondLMMatrix() { return thermcondLM_; }
    
    //! Return thermal Lagrange multiplier
    Teuchos::RCP<Epetra_Vector> ThermLM() { return z_; }
    
    //! Apply thermal contact
    void ApplyThermoContact(Teuchos::RCP<LINALG::SparseMatrix>& tang,
                            Teuchos::RCP<Epetra_Vector>& feff,
                            Teuchos::RCP<Epetra_Vector>& temp,
                            double dt);
    
    // convert maps form structure discretization to thermo discretization
    void ConvertMaps(RCP<Epetra_Map>& slavedofs,
                     RCP<Epetra_Map>& activedofs,
                     RCP<Epetra_Map>& masterdofs);
    
    //! transform mortar matrices from structure dofs in thermo dofs
    // This needs
    // 1. Reduction of mortar matrices
    // 2. Transformation into thermo dofs
    void TransformDM(LINALG::SparseMatrix& dmatrix,
                     LINALG::SparseMatrix& mmatrix,
                     RCP<Epetra_Map>& slavedofs,
                     RCP<Epetra_Map>& masterdofs);
    
    //! transform matrix A
    void TransformA(LINALG::SparseMatrix& amatrix,
                    RCP<Epetra_Map>& slavedofs);

    //! apply B matrix for nodes in contact
    void AssembleB(LINALG::SparseMatrix& bmatrix);

    //! assemble mechanical dissipation for master nodes
    void AssembleMechDissMaster(Epetra_Vector& mechdissrate,double dt);

    //! assemble mechanical dissipation for slave nodes
    void AssembleMechDissSlave(Epetra_Vector& mechdissrate,double dt);

    //! assemble thermal contact condition
    void AssembleThermContCondition(LINALG::SparseMatrix& thermcontLM,
                                    LINALG::SparseMatrix& thermcontTEMP,
                                    Epetra_Vector& thermcontRHS,
                                    LINALG::SparseMatrix& dmatrix,
                                    LINALG::SparseMatrix& mmatrix,
                                    LINALG::SparseMatrix& amatrix,
                                    RCP<Epetra_Map> activedofs,
                                    RCP<Epetra_Map> masterdofs,
                                    Teuchos::RCP<Epetra_Vector>& temp,
                                    double dt);
    
    //! initialize thermal lagrange multipliers
    void InitializeThermLM(RCP<Epetra_Map> sthermdofs);
    
    //! recover thermal lagrange multipliers from temperature
    void RecoverThermLM(RCP<Epetra_Vector> tempi);
  
    protected:
  
    //! copy constructor is NOT wanted
    ThermoContactMan(const ThermoContactMan& old);
  
    Teuchos::RCP<MORTAR::ManagerBase> cmtman_;         // structural contact manager     
    Teuchos::RCP<DRT::Discretization> discretstruct_;  // structural discretisation
    Teuchos::RCP<DRT::Discretization> discretthermo_;  // thermal discretisation
    
    RCP<Epetra_Comm>                      comm_;       // an Epetra_comm
    
    // some RCP's needed for contact,
    // will be moved to seperate class soon!
    RCP<LINALG::SparseMatrix> thermcondLM_; // stiffness of thermal condition
                                            // with respect to Lagrange multipliers
    RCP<Epetra_Vector> z_;                  // thermal Lagrange multipliers
    RCP<Epetra_Vector> fs_;                 // slave side effective heat fluxes (needed for LM)
    RCP<LINALG::SparseMatrix> invd_;        // inverse of Mortar matrix D (needed for LM)
    RCP<LINALG::SparseMatrix> ksn_;         // stiffness block K_sn (needed for LM)
    RCP<LINALG::SparseMatrix> ksm_;         // stiffness block K_sm (needed for LM)
    RCP<LINALG::SparseMatrix> kss_;         // stiffness block K_ss (needed for LM)

    
  };  // class ThermoContactMan

}  // namespace THR

/*----------------------------------------------------------------------*/
#endif  // #ifndef THR_CONTACT_H
#endif  // #ifdef CCADISCRET
