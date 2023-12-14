/*----------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problems using a Dirichlet-Neumann partitioned approach
       with volume coupling

\level 3

*/
/*----------------------------------------------------------------------*/



#ifndef BACI_FSI_DIRICHLETNEUMANN_VOLCOUPL_H
#define BACI_FSI_DIRICHLETNEUMANN_VOLCOUPL_H

#include "baci_coupling_adapter.H"
#include "baci_coupling_adapter_volmortar.H"
#include "baci_fsi_dirichletneumann_disp.H"

BACI_NAMESPACE_OPEN

namespace CORE::GEO
{
  class SearchTree;
}

namespace FSI
{
  class InterfaceCorrector;

  /// Dirichlet-Neumann VolCoupled system
  class DirichletNeumannVolCoupl : public DirichletNeumannDisp
  {
    friend class DirichletNeumannFactory;

   protected:
    /**
     *  \brief constructor
     *
     * You will have to use the FSI::DirichletNeumannFactory to create an instance of this class
     */
    explicit DirichletNeumannVolCoupl(const Epetra_Comm& comm);

   public:
    /// setup this object
    void Setup() override;

   protected:
    /// setup
    void SetupCouplingStructAle(const Teuchos::ParameterList& fsidyn, const Epetra_Comm& comm);

    /// setup
    void SetupInterfaceCorrector(const Teuchos::ParameterList& fsidyn, const Epetra_Comm& comm);

    /** \brief interface fluid operator
     *
     * Solve the fluid field problem.  Since the fluid field is the Dirichlet partition, the
     * interface displacement is prescribed as a Dirichlet boundary condition.
     *
     * \param[in] idisp interface displacement
     * \param[in] fillFlag Type of evaluation in computeF() (cf. NOX documentation for details)
     *
     * \returns interface force
     */
    Teuchos::RCP<Epetra_Vector> FluidOp(
        Teuchos::RCP<Epetra_Vector> idisp, const FillType fillFlag) final;


    void ExtractPreviousInterfaceSolution() override;

    /// structure to ale mapping
    Teuchos::RCP<Epetra_Vector> StuctureToAle(Teuchos::RCP<Epetra_Vector> iv) const;

    /// structure to ale mapping
    Teuchos::RCP<Epetra_Vector> StructureToAle(Teuchos::RCP<const Epetra_Vector> iv) const;

    /// ale to structure mapping
    Teuchos::RCP<Epetra_Vector> AleToStructure(Teuchos::RCP<Epetra_Vector> iv) const;

    /// ale to structure
    Teuchos::RCP<Epetra_Vector> AleToStructure(Teuchos::RCP<const Epetra_Vector> iv) const;

    /// coupling of structure and ale at the interface
    Teuchos::RCP<CORE::ADAPTER::MortarVolCoupl> coupsa_;

    /// coupling of structure and ale at the interface
    Teuchos::RCP<InterfaceCorrector> icorrector_;
  };

  class VolCorrector;

  class InterfaceCorrector
  {
   public:
    /// constructor
    InterfaceCorrector() : idisp_(Teuchos::null){};

    /// destructor
    virtual ~InterfaceCorrector() = default;

    virtual void Setup(Teuchos::RCP<ADAPTER::FluidAle> fluidale);

    void SetInterfaceDisplacements(
        Teuchos::RCP<Epetra_Vector>& idisp_struct, CORE::ADAPTER::Coupling& icoupfs);

    virtual void CorrectInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp_fluid,
        Teuchos::RCP<FLD::UTILS::MapExtractor> const& finterface);

   private:
    Teuchos::RCP<const Epetra_Vector> idisp_;
    Teuchos::RCP<CORE::ADAPTER::Coupling> icoupfs_;

    Teuchos::RCP<Epetra_Vector> deltadisp_;
    Teuchos::RCP<ADAPTER::FluidAle> fluidale_;

    Teuchos::RCP<VolCorrector> volcorrector_;
  };


  class VolCorrector
  {
   public:
    /// constructor
    VolCorrector() : dim_(-1){};

    /// destructor
    virtual ~VolCorrector() = default;

    virtual void Setup(const int dim, Teuchos::RCP<ADAPTER::FluidAle> fluidale);

    virtual void CorrectVolDisplacements(Teuchos::RCP<ADAPTER::FluidAle> fluidale,
        Teuchos::RCP<Epetra_Vector> deltadisp, Teuchos::RCP<Epetra_Vector> idisp_fluid,
        Teuchos::RCP<FLD::UTILS::MapExtractor> const& finterface);

   private:
    virtual void CorrectVolDisplacementsParaSpace(Teuchos::RCP<ADAPTER::FluidAle> fluidale,
        Teuchos::RCP<Epetra_Vector> deltadisp, Teuchos::RCP<Epetra_Vector> idisp_fluid,
        Teuchos::RCP<FLD::UTILS::MapExtractor> const& finterface);

    virtual void CorrectVolDisplacementsPhysSpace(Teuchos::RCP<ADAPTER::FluidAle> fluidale,
        Teuchos::RCP<Epetra_Vector> deltadisp, Teuchos::RCP<Epetra_Vector> idisp_fluid,
        Teuchos::RCP<FLD::UTILS::MapExtractor> const& finterface);

    void InitDopNormals();

    std::map<int, CORE::LINALG::Matrix<9, 2>> CalcBackgroundDops(
        Teuchos::RCP<DRT::Discretization> searchdis);

    CORE::LINALG::Matrix<9, 2> CalcDop(DRT::Element& ele);

    std::vector<int> Search(
        DRT::Element& ele, std::map<int, CORE::LINALG::Matrix<9, 2>>& currentKDOPs);

    //! Spatial dimension of the problem
    int dim_;

    //! Searchtree for mortar evaluations
    Teuchos::RCP<CORE::GEO::SearchTree> searchTree_;

    //! Dop normals for search algorithm
    CORE::LINALG::Matrix<9, 3> dopnormals_;

    std::map<int, std::vector<int>> fluidaleelemap_;

    std::map<int, std::vector<int>> fluidalenodemap_;
    std::map<int, std::vector<int>> fluidalenodeFSImap_;
  };


}  // namespace FSI

BACI_NAMESPACE_CLOSE

#endif
