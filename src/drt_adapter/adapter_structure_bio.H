/*----------------------------------------------------------------------*/
/*!
\file adapter_structure_bio.H

<pre>
Maintainer: Mirella Coroneo
            coroneo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef ADAPTER_STRUCTURE_BIO_H
#define ADAPTER_STRUCTURE_BIO_H

#include "../drt_adapter/adapter_structure.H"
#include "../drt_ale/ale.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "adapter_coupling.H"

#include "adapter_algorithmbase.H"

namespace ADAPTER
{
  class StructureBio : public AlgorithmBase,
        public ADAPTER::StructureBaseAlgorithm,
        public ADAPTER::ScaTraBaseAlgorithm,
        public ALE::AleBaseAlgorithm
  {
  public:

    /// constructor using a Epetra_Comm
    explicit StructureBio(
        Epetra_Comm& comm,              			///< communicator
        const Teuchos::ParameterList& prbdyn, 		///< problem-specific parameters
        bool isale,  								///< do we need an ALE formulation of the structure?
        const int disnum=0, 						///< scatra discretization number (default: 0)
        const std::string condname = "FSICoupling"	///< name of condition that defines structure-ale coupling
        );

    /// virtual destructor to support polymorph destruction
    virtual ~StructureBio();

    /// solve structure-ale
    virtual void StructAleSolve(Teuchos::RCP<Epetra_Vector> idisp);

    /// provide access to algorithm parameters
    virtual const Teuchos::ParameterList& AlgoParameters(){
        return params_;
    }

    /// apply the displacements at the interface
    void ApplyInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp);

    void ChangeConfig(Teuchos::RCP<Epetra_Vector> structdisp);

  protected:

    //! @name Transfer helpers

    /// field transform
    virtual Teuchos::RCP<Epetra_Vector> AleToStructField(Teuchos::RCP<Epetra_Vector> iv) const;

    /// field transform
    virtual Teuchos::RCP<Epetra_Vector> AleToStructField(Teuchos::RCP<const Epetra_Vector> iv) const;

    /// interface transform
    virtual Teuchos::RCP<Epetra_Vector> StructToAle(Teuchos::RCP<Epetra_Vector> iv) const;

    /// interface transform
    virtual Teuchos::RCP<Epetra_Vector> StructToAle(Teuchos::RCP<const Epetra_Vector> iv) const;

    Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structurebasealgorithm;
   // Structure& structure;

  private:

    /// problem-specific parameter list
    const Teuchos::ParameterList& params_;

    /// coupling of structure and ale (interface only)
    ADAPTER::Coupling icoupsa_;

    /// coupling of structure and ale (whole field)
    ADAPTER::Coupling coupsa_;

    RCP<DRT::Discretization> structuredis;

    ///  displacement at time t^{n+1} ???
    RCP<Epetra_Vector> dispnp_;

    Teuchos::RCP<Epetra_Vector> structdisp;

  };

} // namespace ADAPTER

#endif // ADAPTER_STRUCTURE_BIO_H
#endif // CCADISCRET
