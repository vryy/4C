/*----------------------------------------------------------------------*/
/*!
\file ad_str_bio.H

<pre>
Maintainer: Mirella Coroneo
            coroneo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/


#ifndef ADAPTER_STRUCTURE_BIO_H
#define ADAPTER_STRUCTURE_BIO_H

#include "ad_str_structure.H"
#include "../drt_ale/ale.H"
#include "adapter_scatra_base_algorithm.H"
#include "adapter_algorithmbase.H"

namespace ADAPTER
{
  // forward declarations

  class Coupling;
  class FSIStructureWrapper;


  class StructureBio
	  : public AlgorithmBase
  {
  public:

    /// constructor using a Epetra_Comm
    explicit StructureBio(
        const Epetra_Comm& comm,              			///< communicator
        const Teuchos::ParameterList& prbdyn, 		///< problem-specific parameters
        bool isale,  								///< do we need an ALE formulation of the structure?
        const std::string disname="scatra", 						///< scatra discretization name (default: "scatra")
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

    virtual void ReadRestart(int step);

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


  private:

    /// problem-specific parameter list
    const Teuchos::ParameterList& params_;

    /// coupling of structure and ale (interface only)
    Teuchos::RCP<ADAPTER::Coupling> icoupsa_;

    /// coupling of structure and ale (whole field)
    Teuchos::RCP<ADAPTER::Coupling> coupsa_;

    ///  displacement at time t^{n+1} ???
    Teuchos::RCP<Epetra_Vector> dispnp_;

    Teuchos::RCP<Epetra_Vector> structdisp;

    Teuchos::RCP<FSIStructureWrapper> structure_;

    Teuchos::RCP<ALE::AleBaseAlgorithm> ale;
  };

} // namespace ADAPTER

#endif // ADAPTER_STRUCTURE_BIO_H
