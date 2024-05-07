/*----------------------------------------------------------------------*/
/*! \file
\brief A finite element for simulation transport phenomena

\level 1

*/
/*----------------------------------------------------------------------*/
#include "4C_scatra_ele.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_discretization_fem_general_utils_gausspoints.hpp"
#include "4C_discretization_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_fluid_ele_nullspace.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_elasthyper.hpp"
#include "4C_mat_elchmat.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_list_chemoreac.hpp"
#include "4C_mat_list_chemotaxis.hpp"
#include "4C_mat_list_reactions.hpp"
#include "4C_mat_myocard.hpp"
#include "4C_mat_scatra_chemotaxis.hpp"
#include "4C_mat_scatra_reaction.hpp"
#include "4C_scatra_ele_calc_utils.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


DRT::ELEMENTS::TransportType DRT::ELEMENTS::TransportType::instance_;

DRT::ELEMENTS::TransportType& DRT::ELEMENTS::TransportType::Instance() { return instance_; }

CORE::COMM::ParObject* DRT::ELEMENTS::TransportType::Create(const std::vector<char>& data)
{
  DRT::ELEMENTS::Transport* object = new DRT::ELEMENTS::Transport(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::TransportType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "TRANSP" or eletype == "CONDIF2" or eletype == "CONDIF3")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Transport(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::TransportType::Create(const int id, const int owner)
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Transport(id, owner));
  return ele;
}


void DRT::ELEMENTS::TransportType::NodalBlockInformation(
    DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = dwele->NumDofPerNode(*(dwele->Nodes()[0]));
  dimns = numdf;
  nv = numdf;

  if (GLOBAL::Problem::Instance(0)->GetProblemType() == GLOBAL::ProblemType::elch)
  {
    if (nv > 1)  // only when we have more than 1 dof per node!
    {
      nv -= 1;  // ion concentrations
      np = 1;   // electric potential
    }
  }
}

CORE::LINALG::SerialDenseMatrix DRT::ELEMENTS::TransportType::ComputeNullSpace(
    DRT::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return FLD::ComputeFluidNullSpace(node, numdof, dimnsp);
}

void DRT::ELEMENTS::TransportType::SetupElementDefinition(
    std::map<std::string, std::map<std::string, INPUT::LineDefinition>>& definitions)
{
  std::map<std::string, INPUT::LineDefinition>& defs = definitions["TRANSP"];

  defs["HEX8"] = INPUT::LineDefinition::Builder()
                     .AddIntVector("HEX8", 8)
                     .AddNamedInt("MAT")
                     .AddNamedString("TYPE")
                     .AddOptionalNamedDoubleVector("FIBER1", 3)
                     .Build();

  defs["HEX20"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("HEX20", 20)
                      .AddNamedInt("MAT")
                      .AddNamedString("TYPE")
                      .AddOptionalNamedDoubleVector("FIBER1", 3)
                      .Build();

  defs["HEX27"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("HEX27", 27)
                      .AddNamedInt("MAT")
                      .AddNamedString("TYPE")
                      .AddOptionalNamedDoubleVector("FIBER1", 3)
                      .Build();

  defs["NURBS27"] = INPUT::LineDefinition::Builder()
                        .AddIntVector("NURBS27", 27)
                        .AddNamedInt("MAT")
                        .AddNamedString("TYPE")
                        .AddOptionalNamedDoubleVector("FIBER1", 3)
                        .Build();

  defs["NURBS8"] = INPUT::LineDefinition::Builder()
                       .AddIntVector("NURBS8", 8)
                       .AddNamedInt("MAT")
                       .AddNamedString("TYPE")
                       .AddOptionalNamedDoubleVector("FIBER1", 3)
                       .Build();

  defs["TET4"] = INPUT::LineDefinition::Builder()
                     .AddIntVector("TET4", 4)
                     .AddNamedInt("MAT")
                     .AddNamedString("TYPE")
                     .AddOptionalNamedDoubleVector("FIBER1", 3)
                     .Build();

  defs["TET10"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("TET10", 10)
                      .AddNamedInt("MAT")
                      .AddNamedString("TYPE")
                      .AddOptionalNamedDoubleVector("FIBER1", 3)
                      .Build();

  defs["WEDGE6"] = INPUT::LineDefinition::Builder()
                       .AddIntVector("WEDGE6", 6)
                       .AddNamedInt("MAT")
                       .AddNamedString("TYPE")
                       .AddOptionalNamedDoubleVector("FIBER1", 3)
                       .Build();

  defs["WEDGE15"] = INPUT::LineDefinition::Builder()
                        .AddIntVector("WEDGE15", 15)
                        .AddNamedInt("MAT")
                        .AddNamedString("TYPE")
                        .AddOptionalNamedDoubleVector("FIBER1", 3)
                        .Build();

  defs["PYRAMID5"] = INPUT::LineDefinition::Builder()
                         .AddIntVector("PYRAMID5", 5)
                         .AddNamedInt("MAT")
                         .AddNamedString("TYPE")
                         .AddOptionalNamedDoubleVector("FIBER1", 3)
                         .Build();

  defs["QUAD4"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("QUAD4", 4)
                      .AddNamedInt("MAT")
                      .AddNamedString("TYPE")
                      .AddOptionalNamedDoubleVector("FIBER1", 3)
                      .Build();

  defs["QUAD8"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("QUAD8", 8)
                      .AddNamedInt("MAT")
                      .AddNamedString("TYPE")
                      .AddOptionalNamedDoubleVector("FIBER1", 3)
                      .Build();

  defs["QUAD9"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("QUAD9", 9)
                      .AddNamedInt("MAT")
                      .AddNamedString("TYPE")
                      .AddOptionalNamedDoubleVector("FIBER1", 3)
                      .Build();

  defs["TRI3"] = INPUT::LineDefinition::Builder()
                     .AddIntVector("TRI3", 3)
                     .AddNamedInt("MAT")
                     .AddNamedString("TYPE")
                     .AddOptionalNamedDoubleVector("FIBER1", 3)
                     .Build();

  defs["TRI6"] = INPUT::LineDefinition::Builder()
                     .AddIntVector("TRI6", 6)
                     .AddNamedInt("MAT")
                     .AddNamedString("TYPE")
                     .AddOptionalNamedDoubleVector("FIBER1", 3)
                     .Build();

  defs["NURBS4"] = INPUT::LineDefinition::Builder()
                       .AddIntVector("NURBS4", 4)
                       .AddNamedInt("MAT")
                       .AddNamedString("TYPE")
                       .AddOptionalNamedDoubleVector("FIBER1", 3)
                       .Build();

  defs["NURBS9"] = INPUT::LineDefinition::Builder()
                       .AddIntVector("NURBS9", 9)
                       .AddNamedInt("MAT")
                       .AddNamedString("TYPE")
                       .AddOptionalNamedDoubleVector("FIBER1", 3)
                       .Build();

  defs["LINE2"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("LINE2", 2)
                      .AddNamedInt("MAT")
                      .AddNamedString("TYPE")
                      .AddOptionalNamedDoubleVector("FIBER1", 3)
                      .Build();

  defs["LINE3"] = INPUT::LineDefinition::Builder()
                      .AddIntVector("LINE3", 3)
                      .AddNamedInt("MAT")
                      .AddNamedString("TYPE")
                      .AddOptionalNamedDoubleVector("FIBER1", 3)
                      .Build();

  defs["NURBS2"] = INPUT::LineDefinition::Builder()
                       .AddIntVector("NURBS2", 2)
                       .AddNamedInt("MAT")
                       .AddNamedString("TYPE")
                       .AddOptionalNamedDoubleVector("FIBER1", 3)
                       .Build();

  defs["NURBS3"] = INPUT::LineDefinition::Builder()
                       .AddIntVector("NURBS3", 3)
                       .AddNamedInt("MAT")
                       .AddNamedString("TYPE")
                       .AddOptionalNamedDoubleVector("FIBER1", 3)
                       .Build();
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TransportType::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::Transport* actele = dynamic_cast<DRT::ELEMENTS::Transport*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to Transport element failed");
    actele->Initialize();
  }
  return 0;
}


DRT::ELEMENTS::TransportBoundaryType DRT::ELEMENTS::TransportBoundaryType::instance_;

DRT::ELEMENTS::TransportBoundaryType& DRT::ELEMENTS::TransportBoundaryType::Instance()
{
  return instance_;
}

Teuchos::RCP<DRT::Element> DRT::ELEMENTS::TransportBoundaryType::Create(
    const int id, const int owner)
{
  // return Teuchos::rcp( new TransportBoundary( id, owner ) );
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                             gjb 05/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Transport::Transport(int id, int owner)
    : DRT::Element(id, owner),
      distype_(CORE::FE::CellType::dis_none),
      name_(),
      vis_map_(),
      numdofpernode_(-1),
      impltype_(INPAR::SCATRA::impltype_undefined)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        gjb 05/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Transport::Transport(const DRT::ELEMENTS::Transport& old)
    : DRT::Element(old),
      distype_(old.distype_),
      name_(old.name_),
      vis_map_(old.vis_map_),
      numdofpernode_(old.numdofpernode_),
      impltype_(old.impltype_)
{
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Transport and return pointer to it (public) |
 |                                                            gjb 05/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Transport::Clone() const
{
  DRT::ELEMENTS::Transport* newelement = new DRT::ELEMENTS::Transport(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |  create material class (public)                            gjb 07/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Transport::SetMaterial(const int index, Teuchos::RCP<CORE::MAT::Material> mat)
{
  // the standard part:
  DRT::Element::SetMaterial(index, mat);

  if (mat->MaterialType() == CORE::Materials::m_scatra or
      mat->MaterialType() == CORE::Materials::m_scatra_aniso or
      mat->MaterialType() == CORE::Materials::m_scatra_multiscale or
      mat->MaterialType() == CORE::Materials::m_myocard or
      mat->MaterialType() == CORE::Materials::m_mixfrac or
      mat->MaterialType() == CORE::Materials::m_sutherland or
      mat->MaterialType() == CORE::Materials::m_tempdepwater or
      mat->MaterialType() == CORE::Materials::m_arrhenius_pv or
      mat->MaterialType() == CORE::Materials::m_ferech_pv or
      mat->MaterialType() == CORE::Materials::m_ion or
      mat->MaterialType() == CORE::Materials::m_th_fourier_iso or
      mat->MaterialType() == CORE::Materials::m_thermostvenant or
      mat->MaterialType() == CORE::Materials::m_yoghurt or
      mat->MaterialType() == CORE::Materials::m_soret or
      mat->MaterialType() == CORE::Materials::m_scatra_multiporo_fluid or
      mat->MaterialType() == CORE::Materials::m_scatra_multiporo_volfrac or
      mat->MaterialType() == CORE::Materials::m_scatra_multiporo_solid or
      mat->MaterialType() == CORE::Materials::m_scatra_multiporo_temperature or
      (mat->MaterialType() == CORE::Materials::m_electrode and
          impltype_ == INPAR::SCATRA::impltype_std))
    numdofpernode_ = 1;  // we only have a single scalar
  else if (mat->MaterialType() == CORE::Materials::m_electrode)
    numdofpernode_ = 2;  // concentration and electric potential
  else if (mat->MaterialType() == CORE::Materials::m_matlist)  // we have a system of scalars
  {
    const MAT::MatList* actmat = static_cast<const MAT::MatList*>(mat.get());
    numdofpernode_ = actmat->NumMat();

    // for problem type ELCH we have one additional degree of freedom per node
    // for the electric potential
    if (GLOBAL::Problem::Instance()->GetProblemType() == GLOBAL::ProblemType::elch)
    {
      for (int ii = 0; ii < numdofpernode_; ++ii)
      {
        // In the context of ELCH the only valid material combination is m_matlist and m_ion
        if (actmat->MaterialById(actmat->MatID(ii))->MaterialType() != CORE::Materials::m_ion)
          FOUR_C_THROW(
              "In the context of ELCH the material Mat_matlist can be only used in combination "
              "with Mat_ion");
      }
      numdofpernode_ += 1;
    }
    // for problem type LOMA, only combination of Arrhenius-type species (first)
    // and temperature (last) equation possible in this specific order
    // in case of matlist
    else if (GLOBAL::Problem::Instance()->GetProblemType() == GLOBAL::ProblemType::loma)
    {
      // only two-equation systems, for the time being: check!
      if (numdofpernode_ > 2)
        FOUR_C_THROW(
            "Only two-equation systems (one species and one temperature equation for "
            "Arrhenius-type systems, for the time being!");

      // check that first equations are species equations and that temperature
      // equation is last equation
      for (int ii = 0; ii < (numdofpernode_ - 1); ++ii)
      {
        if (actmat->MaterialById(actmat->MatID(ii))->MaterialType() !=
            CORE::Materials::m_arrhenius_spec)
          FOUR_C_THROW(
              "For problem type LOMA, only combination of Arrhenius-type species (first equations) "
              "and temperature (last equation) possible in this specific order in case of matlist: "
              "one of the first equations is not a species equation!");
      }
      if (actmat->MaterialById(actmat->MatID(numdofpernode_ - 1))->MaterialType() !=
          CORE::Materials::m_arrhenius_temp)
        FOUR_C_THROW(
            "For problem type LOMA, only combination of Arrhenius-type species (first equations) "
            "and temperature (last equation) possible in this specific order in case of matlist: "
            "last equation is not a temperature equation!");
    }
  }
  else if (mat->MaterialType() ==
           CORE::Materials::m_matlist_reactions)  // we have a system of reactive scalars
  {
    // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are in
    // a diamond inheritance structure
    const MAT::MatListReactions* actmat = dynamic_cast<const MAT::MatListReactions*>(mat.get());
    numdofpernode_ = actmat->NumMat();

    for (int ii = 0; ii < numdofpernode_; ++ii)
    {
      // In the context of reactions the only valid material combination is m_matlist and m_scatra
      if (actmat->MaterialById(actmat->MatID(ii))->MaterialType() != CORE::Materials::m_scatra and
          actmat->MaterialById(actmat->MatID(ii))->MaterialType() !=
              CORE::Materials::m_scatra_multiporo_fluid and
          actmat->MaterialById(actmat->MatID(ii))->MaterialType() !=
              CORE::Materials::m_scatra_multiporo_volfrac and
          actmat->MaterialById(actmat->MatID(ii))->MaterialType() !=
              CORE::Materials::m_scatra_multiporo_temperature and
          actmat->MaterialById(actmat->MatID(ii))->MaterialType() !=
              CORE::Materials::m_scatra_multiporo_solid)
        FOUR_C_THROW(
            "The material Mat_matlist_reaction only supports MAT_scatra and MAT_scatra_multiporo "
            "as valid main Material");
    }

    int numreac = actmat->NumReac();
    for (int jj = 0; jj < numreac; ++jj)
    {
      // In the context of reactions the only valid material combination is m_matlist and
      // m_scatra_reaction
      if (actmat->MaterialById(actmat->ReacID(jj))->MaterialType() !=
              CORE::Materials::m_scatra_reaction and
          actmat->MaterialById(actmat->ReacID(jj))->MaterialType() !=
              CORE::Materials::m_scatra_reaction_poroECM)
        FOUR_C_THROW(
            "The material MAT_matlist_reaction only supports MAT_scatra_reaction and "
            "MAT_scatra_reaction_poro as valid reaction Material");

      // some safty check for the MAT_scatra_reaction materials
      const Teuchos::RCP<const MAT::ScatraReactionMat>& reacmat =
          Teuchos::rcp_static_cast<const MAT::ScatraReactionMat>(
              actmat->MaterialById(actmat->ReacID(jj)));
      const int stoichlength = reacmat->NumScal();
      if (stoichlength != numdofpernode_)
        FOUR_C_THROW(
            "The number of scalars in your MAT_scatra_reaction material with ID %i does not fit to "
            "the number of scalars!",
            actmat->ReacID(jj));
    }
  }
  else if (mat->MaterialType() ==
           CORE::Materials::m_matlist_chemotaxis)  // we have a system of chemotactic scalars
  {
    // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are in
    // a diamond inheritance structure
    const MAT::MatListChemotaxis* actmat = dynamic_cast<const MAT::MatListChemotaxis*>(mat.get());
    numdofpernode_ = actmat->NumMat();

    for (int ii = 0; ii < numdofpernode_; ++ii)
    {
      // In the context of chemotaxis the only valid material combination is m_matlist and m_scatra
      if (actmat->MaterialById(actmat->MatID(ii))->MaterialType() != CORE::Materials::m_scatra)
        FOUR_C_THROW(
            "The material Mat_matlist_chemotaxis only supports MAT_scatra as valid main Material");
    }

    int numpair = actmat->NumPair();
    for (int jj = 0; jj < numpair; ++jj)
    {
      // In the context of chemotaxis the only valid material combination is m_matlist and
      // m_scatra_chemotaxis
      if (actmat->MaterialById(actmat->PairID(jj))->MaterialType() !=
          CORE::Materials::m_scatra_chemotaxis)
        FOUR_C_THROW(
            "The material MAT_matlist_chemotaxis only supports MAT_scatra_chemotaxis as valid "
            "reaction Material");

      // some safty check for the MAT_scatra_chemotaxis materials
      const Teuchos::RCP<const MAT::ScatraChemotaxisMat>& reacmat =
          Teuchos::rcp_static_cast<const MAT::ScatraChemotaxisMat>(
              actmat->MaterialById(actmat->PairID(jj)));
      const int pairlength = reacmat->Pair()->size();
      if (pairlength != numdofpernode_)
        FOUR_C_THROW(
            "The number of scalars in your MAT_scatra_chemotaxis material with ID %i does not fit "
            "to the number of scalars!",
            actmat->PairID(jj));
    }
  }
  else if (mat->MaterialType() ==
           CORE::Materials::m_matlist_chemoreac)  // we have a system of chemotactic scalars
  {
    // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are in
    // a diamond inheritance structure
    const MAT::MatListChemoReac* actmat = dynamic_cast<const MAT::MatListChemoReac*>(mat.get());
    numdofpernode_ = actmat->NumMat();

    for (int ii = 0; ii < numdofpernode_; ++ii)
    {
      // In the context of reactions/chemotaxis the only valid material combination is m_matlist and
      // m_scatra
      if (actmat->MaterialById(actmat->MatID(ii))->MaterialType() != CORE::Materials::m_scatra)
        FOUR_C_THROW(
            "The material Mat_matlist_chemoreac only supports MAT_scatra as valid main Material");
    }

    int numreac = actmat->NumReac();
    for (int jj = 0; jj < numreac; ++jj)
    {
      // In the context of reactions the only valid material combination is m_matlist and
      // m_scatra_reaction
      if (actmat->MaterialById(actmat->ReacID(jj))->MaterialType() !=
          CORE::Materials::m_scatra_reaction)
        FOUR_C_THROW(
            "The material MAT_matlist_reaction only supports MAT_scatra_reaction as valid reaction "
            "Material");

      // some safty check for the MAT_scatra_reaction materials
      const Teuchos::RCP<const MAT::ScatraReactionMat>& reacmat =
          Teuchos::rcp_static_cast<const MAT::ScatraReactionMat>(
              actmat->MaterialById(actmat->ReacID(jj)));
      const int stoichlength = reacmat->NumScal();
      if (stoichlength != numdofpernode_)
        FOUR_C_THROW(
            "The number of scalars in your MAT_scatra_reaction material with ID %i does not fit to "
            "the number of scalars!",
            actmat->ReacID(jj));
    }

    int numpair = actmat->NumPair();
    for (int jj = 0; jj < numpair; ++jj)
    {
      // In the context of chemotaxis the only valid material combination is m_matlist and
      // m_scatra_chemotaxis
      if (actmat->MaterialById(actmat->PairID(jj))->MaterialType() !=
          CORE::Materials::m_scatra_chemotaxis)
        FOUR_C_THROW(
            "The material MAT_matlist_chemotaxis only supports MAT_scatra_chemotaxis as valid "
            "reaction Material");

      // some safty check for the MAT_scatra_chemotaxis materials
      const Teuchos::RCP<const MAT::ScatraChemotaxisMat>& reacmat =
          Teuchos::rcp_static_cast<const MAT::ScatraChemotaxisMat>(
              actmat->MaterialById(actmat->PairID(jj)));
      const int pairlength = reacmat->Pair()->size();
      if (pairlength != numdofpernode_)
        FOUR_C_THROW(
            "The number of scalars in your MAT_scatra_chemotaxis material with ID %i does not fit "
            "to the number of scalars!",
            actmat->PairID(jj));
    }
  }
  else if (mat->MaterialType() == CORE::Materials::m_elchmat)
  {
    const MAT::ElchMat* actmat = static_cast<const MAT::ElchMat*>(mat.get());

    numdofpernode_ = actmat->NumDOF();
  }
  else
    FOUR_C_THROW("Transport element got unsupported material type %d", mat->MaterialType());

  return;
}

/*----------------------------------------------------------------------*
 |  create material class (public)                            gjb 07/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Transport::SetMaterial(int matnum, DRT::Element* oldele)
{
  SetMaterial(0, MAT::Factory(matnum));

  Teuchos::RCP<CORE::MAT::Material> mat = Material();

  if (mat->MaterialType() == CORE::Materials::m_myocard)
  {
    Teuchos::RCP<MAT::Myocard> actmat = Teuchos::rcp_dynamic_cast<MAT::Myocard>(mat);

    Teuchos::RCP<MAT::ElastHyper> somat =
        Teuchos::rcp_dynamic_cast<MAT::ElastHyper>(oldele->Material());
    if (somat == Teuchos::null) FOUR_C_THROW("cast to ElastHyper failed");

    // copy fiber information from solid material to scatra material (for now, only one fiber
    // vector)
    std::vector<CORE::LINALG::Matrix<3, 1>> fibervecs(0);
    somat->GetFiberVecs(fibervecs);
    actmat->Setup(fibervecs[0]);
  }
}

/*----------------------------------------------------------------------*
 |  Return the shape of a Transport element                      (public) |
 |                                                            gjb 05/08 |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::Transport::Shape() const { return distype_; }

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gjb 05/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Transport::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // add base class Element
  Element::Pack(data);

  // add internal data
  AddtoPack(data, name_);
  AddtoPack(data, vis_map_);
  AddtoPack(data, numdofpernode_);
  AddtoPack(data, distype_);
  AddtoPack(data, impltype_);
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gjb 05/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Transport::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);

  // extract internal data
  ExtractfromPack(position, data, name_);
  ExtractfromPack(position, data, vis_map_);
  ExtractfromPack(position, data, numdofpernode_);
  distype_ = static_cast<CORE::FE::CellType>(ExtractInt(position, data));
  impltype_ = static_cast<INPAR::SCATRA::ImplType>(ExtractInt(position, data));

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

/*----------------------------------------------------------------------*
 |  Return number of lines of this element (public)           gjb 07/08 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Transport::NumLine() const
{
  return CORE::FE::getNumberOfElementLines(distype_);
}


/*----------------------------------------------------------------------*
 |  Return number of surfaces of this element (public)        gjb 07/08 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Transport::NumSurface() const
{
  return CORE::FE::getNumberOfElementSurfaces(distype_);
}


/*----------------------------------------------------------------------*
 | Return number of volumes of this element (public)          gjb 07/08 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Transport::NumVolume() const
{
  return CORE::FE::getNumberOfElementVolumes(distype_);
}



/*----------------------------------------------------------------------*
 |  print this element (public)                               gjb 05/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Transport::Print(std::ostream& os) const
{
  os << "Transport element";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  " << CORE::FE::CellTypeToString(distype_) << std::endl;
  std::cout << std::endl;
  std::cout << "Number DOF per Node: " << numdofpernode_ << std::endl;
  std::cout << std::endl;
  std::cout << "Type of scalar transport: " << SCATRA::ImplTypeToString(impltype_) << std::endl;
  std::cout << std::endl;
}


/*----------------------------------------------------------------------*
 |  get vector of lines            (public)                  g.bau 03/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Transport::Lines()
{
  return CORE::COMM::GetElementLines<TransportBoundary, Transport>(*this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          g.bau 03/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Transport::Surfaces()
{
  return CORE::COMM::GetElementSurfaces<TransportBoundary, Transport>(*this);
}

/*----------------------------------------------------------------------*
 | set implementation type                                   fang 02/15 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Transport ::SetImplType(const INPAR::SCATRA::ImplType impltype)
{
  // set implementation type
  impltype_ = impltype;
}

/*----------------------------------------------------------------------*
 |  init the element                                        vuong08/16 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Transport::Initialize()
{
  Teuchos::RCP<CORE::MAT::Material> mat = Material();
  // for now, we only need to do something in case of reactions (for the initialization of functions
  // in case of reactions by function)
  if (mat->MaterialType() == CORE::Materials::m_matlist_reactions or
      mat->MaterialType() == CORE::Materials::m_matlist_chemoreac)
  {
    // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are in
    // a diamond inheritance structure
    Teuchos::RCP<MAT::MatListReactions> actmat =
        Teuchos::rcp_dynamic_cast<MAT::MatListReactions>(mat);
    actmat->Initialize();
  }
  else if (mat->MaterialType() == CORE::Materials::m_myocard)
  {
    Teuchos::RCP<MAT::Myocard> actmat = Teuchos::rcp_dynamic_cast<MAT::Myocard>(mat);
    int deg = 0;
    if (this->Degree() == 1)
      deg = 4 * this->Degree();
    else
      deg = 3 * this->Degree();
    Teuchos::RCP<CORE::FE::GaussPoints> quadrature(
        CORE::FE::GaussPointCache::Instance().Create(this->Shape(), deg));
    int gp = quadrature->NumPoints();
    if (actmat->Parameter() != nullptr and
        !actmat->MyocardMat())  // in case we are not in post-process mode
    {
      actmat->SetGP(gp);
      actmat->Initialize();
    }
  }

  return 0;
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*----------------------------------------------------------------------*
 |  ctor (public)                                             gjb 01/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TransportBoundary::TransportBoundary(int id, int owner, int nnode,
    const int* nodeids, DRT::Node** nodes, DRT::ELEMENTS::Transport* parent, const int lsurface)
    : DRT::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  SetParentMasterElement(parent, lsurface);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        gjb 01/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TransportBoundary::TransportBoundary(const DRT::ELEMENTS::TransportBoundary& old)
    : DRT::FaceElement(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it     (public) gjb 01/09 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::TransportBoundary::Clone() const
{
  DRT::ELEMENTS::TransportBoundary* newelement = new DRT::ELEMENTS::TransportBoundary(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Return shape of this element                    (public)  gjb 01/09 |
 *----------------------------------------------------------------------*/
CORE::FE::CellType DRT::ELEMENTS::TransportBoundary::Shape() const
{
  return CORE::FE::getShapeOfBoundaryElement(NumNode(), ParentElement()->Shape());
}

/*----------------------------------------------------------------------*
 |  Pack data (public)                                        gjb 01/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TransportBoundary::Pack(CORE::COMM::PackBuffer& data) const
{
  FOUR_C_THROW("This TransportBoundary element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data (public)                                      gjb 01/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TransportBoundary::Unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("This TransportBoundary element does not support communication");
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                               gjb 01/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TransportBoundary::Print(std::ostream& os) const
{
  os << "TransportBoundary element";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  " << CORE::FE::CellTypeToString(Shape()) << std::endl;
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 | Return number of lines of boundary element (public)        gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TransportBoundary::NumLine() const
{
  return CORE::FE::getNumberOfElementLines(Shape());
}

/*----------------------------------------------------------------------*
 |  Return number of surfaces of boundary element (public)    gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TransportBoundary::NumSurface() const
{
  return CORE::FE::getNumberOfElementSurfaces(Shape());
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                              gjb 01/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::TransportBoundary::Lines()
{
  FOUR_C_THROW("Lines of TransportBoundary not implemented");
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                              gjb 01/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::TransportBoundary::Surfaces()
{
  FOUR_C_THROW("Surfaces of TransportBoundary not implemented");
}

FOUR_C_NAMESPACE_CLOSE
