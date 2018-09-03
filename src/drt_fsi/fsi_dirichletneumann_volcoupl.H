/*----------------------------------------------------------------------*/
/*!
\file fsi_dirichletneumann_volcoupl.H

\brief ...

\level 3

<pre>
\maintainer  Philip Farah
             farah@lnm.mw.tum.de
             089 - 289-10362
</pre>
*/
/*----------------------------------------------------------------------*/



#ifndef FSI_DIRICHLETNEUMANN_VOLCOUPL_H
#define FSI_DIRICHLETNEUMANN_VOLCOUPL_H

//#include "fsi_partitioned.H"
#include "fsi_dirichletneumann.H"
#include "../drt_adapter/adapter_coupling_volmortar.H"

namespace GEO
{
  class SearchTree;
}

namespace FSI
{
  class InterfaceCorrector;

  /// Dirichlet Neumann VolCoupled system
  /*!
     documentation

    \author u.kue
    \date 06/07
   */
  class DirichletNeumannVolCoupl : public DirichletNeumann
  {
   public:
    // constructor
    explicit DirichletNeumannVolCoupl(const Epetra_Comm& comm);

    /// setup this object
    virtual void Setup();

   protected:
    // setup
    void SetupCouplingStructAle(const Teuchos::ParameterList& fsidyn, const Epetra_Comm& comm);

    // setup
    void SetupInterfaceCorrector(const Teuchos::ParameterList& fsidyn, const Epetra_Comm& comm);

    /// interface fluid operator
    virtual Teuchos::RCP<Epetra_Vector> FluidOp(
        Teuchos::RCP<Epetra_Vector> idisp, const FillType fillFlag);


    void ExtractPreviousInterfaceSolution();

    /// structure to ale mapping
    Teuchos::RCP<Epetra_Vector> StuctureToAle(Teuchos::RCP<Epetra_Vector> iv) const;

    /// structure to ale mapping
    Teuchos::RCP<Epetra_Vector> StructureToAle(Teuchos::RCP<const Epetra_Vector> iv) const;

    /// ale to structure mapping
    Teuchos::RCP<Epetra_Vector> AleToStructure(Teuchos::RCP<Epetra_Vector> iv) const;

    /// ale to structure
    Teuchos::RCP<Epetra_Vector> AleToStructure(Teuchos::RCP<const Epetra_Vector> iv) const;

   private:
    //! we couple based on displacements
    bool displacementcoupling_;

    /// coupling of structure and ale at the interface
    Teuchos::RCP<ADAPTER::MortarVolCoupl> coupsa_;

    /// coupling of structure and ale at the interface
    Teuchos::RCP<InterfaceCorrector> icorrector_;
  };

  class VolCorrector;

  class InterfaceCorrector
  {
   public:
    // constructor
    InterfaceCorrector() : idisp_(Teuchos::null){};

    // destructor
    virtual ~InterfaceCorrector(){};

    virtual void Setup(Teuchos::RCP<ADAPTER::FluidAle> fluidale);

    void SetInterfaceDisplacements(
        Teuchos::RCP<Epetra_Vector>& idisp_struct, ADAPTER::Coupling& icoupfs);

    virtual void CorrectInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp_fluid,
        Teuchos::RCP<FLD::UTILS::MapExtractor> const& finterface);

   private:
    Teuchos::RCP<const Epetra_Vector> idisp_;
    Teuchos::RCP<ADAPTER::Coupling> icoupfs_;

    Teuchos::RCP<Epetra_Vector> deltadisp_;
    Teuchos::RCP<ADAPTER::FluidAle> fluidale_;

    Teuchos::RCP<VolCorrector> volcorrector_;
  };


  class VolCorrector
  {
   public:
    // constructor
    VolCorrector() : dim_(-1){};

    // destructor
    virtual ~VolCorrector(){};

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

    std::map<int, LINALG::Matrix<9, 2>> CalcBackgroundDops(
        Teuchos::RCP<DRT::Discretization> searchdis);

    LINALG::Matrix<9, 2> CalcDop(DRT::Element& ele);

    std::vector<int> Search(DRT::Element& ele, std::map<int, LINALG::Matrix<9, 2>>& currentKDOPs);

    int dim_;
    Teuchos::RCP<GEO::SearchTree> searchTree_;
    LINALG::Matrix<9, 3> dopnormals_;  /// dop normals for seach algorithm

    std::map<int, std::vector<int>> fluidaleelemap_;

    std::map<int, std::vector<int>> fluidalenodemap_;
    std::map<int, std::vector<int>> fluidalenodeFSImap_;
  };


}  // namespace FSI

#endif
