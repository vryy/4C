/*!----------------------------------------------------------------------
\file ssi_partitioned_2wc_protrusionformation.H

\brief specialization of ssi2wc, which includes "structale"-surface growth

<pre>
\maintainer Andreas Rauch
            rauch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#ifndef SSI_PARTITIONED2_WC_PROTRUSIONFORMATION_H_
#define SSI_PARTITIONED2_WC_PROTRUSIONFORMATION_H_

#include "../drt_ssi/ssi_partitioned_2wc.H"
#include "../drt_inpar/inpar_ssi.H"

namespace ADAPTER
{
  class FSIStructureWrapper;
  class AleFluidWrapper;
  class Coupling;
}

namespace DRT
{
  class Element;
}

namespace SSI
{
  //! specialized class to deal with partitioned 2WC SSI during protrusion formation in cell migration.
  class SSI_Part2WC_PROTRUSIONFORMATION : public SSI_Part2WC
  {

  public:

    //!explicit constructor
    explicit SSI_Part2WC_PROTRUSIONFORMATION(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams,
        const Teuchos::ParameterList& structparams,
        const std::string struct_disname,
        const std::string scatra_disname);

    //!virtual destructor
    virtual ~SSI_Part2WC_PROTRUSIONFORMATION(){};

    //! access to ale field
    const Teuchos::RCP<ADAPTER::AleFluidWrapper> AleField(){return ale_;};

    //! access to scatra field
    const Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm>&  ScaTraField(){return scatra_;};

  protected:

    //! perform iteration step of ale field
    void DoAleStep(Teuchos::RCP<Epetra_Vector> growthincrement);

    virtual void DoStructStep();

    //! access to scatra field
    const Teuchos::RCP< ::ADAPTER::FSIStructureWrapper>&  SpecStructureField(){return specialized_structure_;};

    //! return coupling object
    Teuchos::RCP< ::ADAPTER::Coupling> AleStruCoupling(){return coupalestru_; };

    /*!
    \brief parameter space mapping between configurations

    \param Xtarget           (out) : new material coordinate of considered node (spatialtomaterial=true), or
                                     new spatial  coordinate of considered node (spatialtomaterial=false)
    \param Xsource           (in)  : current mesh coordinate of configuration that is NOT supposed to be transported
    \param ElementPtr        (in)  : pointer to elements adjacent to considered node
    \param numelements       (in)  : number of elements adjacent to considered node
    \param spatialtomaterial (in)  : true if considered node already has correct spatial coordinate
                                     and corresponding material coordinate is to be determined
     */
    void AdvectionMap(
        double* Xtarget,            // out
        double* Xsource,            // in
        DRT::Element** ElementPtr,  // in
        int numelements,            // in
        bool spatialtomaterial);    // in

    //! evaluate growth law for cell
    void EvaluateGrowth();

    //! update material configuration
    void UpdateMatConf();

    //! update spatial configuration
    void UpdateSpatConf();

    //! transform from structure to ale dofs
    virtual Teuchos::RCP<Epetra_Vector> StructureToAle(Teuchos::RCP<Epetra_Vector> vec);

    //! transform from structure to ale dofs
    virtual Teuchos::RCP<Epetra_Vector> StructureToAle(Teuchos::RCP<const Epetra_Vector> vec);

    //! bettina: transform from ale to structure dofs
    virtual Teuchos::RCP<Epetra_Vector> AleToStructure(Teuchos::RCP<Epetra_Vector> vec);

    //! bettina: transform from ale to structure dofs
    virtual Teuchos::RCP<Epetra_Vector> AleToStructure(Teuchos::RCP<const Epetra_Vector> vec);

    //! setup discretizations and dofsets
    virtual void SetupDiscretizations(const Epetra_Comm& comm, const std::string struct_disname, const std::string scatra_disname);

    //! convergence check of outer loop
    virtual bool ConvergenceCheck(int itnum);

  private:

    //! communicator
    int myrank_;

    //! growth step increment (increment containing growth from t_n -> t_{n+1})
    Teuchos::RCP<Epetra_Vector> growth_;

    //! growth iteration increment of the growth loop (increment containing growth from iteration i -> i+1)
    Teuchos::RCP<Epetra_Vector> growthincrement_;

    //! yet unknown vector
    Teuchos::RCP<Epetra_Vector> delta_ale_;

    // coupling objects declaration
    Teuchos::RCP< ::ADAPTER::Coupling>               coupalestru_;   // ale struct coupling
    Teuchos::RCP< ::ADAPTER::Coupling>               coupstrualei_;  // ale struct coupling on ale interface

    //! underlying ALE
    Teuchos::RCP<ADAPTER::AleFluidWrapper> ale_;

    /// underlying structure of the SSI problem
    Teuchos::RCP< ::ADAPTER::FSIStructureWrapper>    specialized_structure_;


  }; // class
} // namespace SSI

#endif /* SSI_PARTITIONED2_WC__PROTRUSIONFORMATION_H_ */
