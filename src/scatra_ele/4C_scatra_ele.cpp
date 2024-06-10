/*----------------------------------------------------------------------*/
/*! \file
\brief A finite element for simulation transport phenomena

\level 1

*/
/*----------------------------------------------------------------------*/
#include "4C_scatra_ele.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_fluid_ele_nullspace.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
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


Discret::ELEMENTS::TransportType Discret::ELEMENTS::TransportType::instance_;

Discret::ELEMENTS::TransportType& Discret::ELEMENTS::TransportType::Instance() { return instance_; }

Core::Communication::ParObject* Discret::ELEMENTS::TransportType::Create(
    const std::vector<char>& data)
{
  Discret::ELEMENTS::Transport* object = new Discret::ELEMENTS::Transport(-1, -1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::TransportType::Create(
    const std::string eletype, const std::string eledistype, const int id, const int owner)
{
  if (eletype == "TRANSP" or eletype == "CONDIF2" or eletype == "CONDIF3")
  {
    Teuchos::RCP<Core::Elements::Element> ele =
        Teuchos::rcp(new Discret::ELEMENTS::Transport(id, owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::TransportType::Create(
    const int id, const int owner)
{
  Teuchos::RCP<Core::Elements::Element> ele =
      Teuchos::rcp(new Discret::ELEMENTS::Transport(id, owner));
  return ele;
}


void Discret::ELEMENTS::TransportType::nodal_block_information(
    Core::Elements::Element* dwele, int& numdf, int& dimns, int& nv, int& np)
{
  numdf = dwele->NumDofPerNode(*(dwele->Nodes()[0]));
  dimns = numdf;
  nv = numdf;

  if (Global::Problem::Instance(0)->GetProblemType() == Core::ProblemType::elch)
  {
    if (nv > 1)  // only when we have more than 1 dof per node!
    {
      nv -= 1;  // ion concentrations
      np = 1;   // electric potential
    }
  }
}

Core::LinAlg::SerialDenseMatrix Discret::ELEMENTS::TransportType::ComputeNullSpace(
    Core::Nodes::Node& node, const double* x0, const int numdof, const int dimnsp)
{
  return FLD::ComputeFluidNullSpace(node, numdof, dimnsp);
}

void Discret::ELEMENTS::TransportType::setup_element_definition(
    std::map<std::string, std::map<std::string, Input::LineDefinition>>& definitions)
{
  std::map<std::string, Input::LineDefinition>& defs = definitions["TRANSP"];

  defs["HEX8"] = Input::LineDefinition::Builder()
                     .AddIntVector("HEX8", 8)
                     .AddNamedInt("MAT")
                     .AddNamedString("TYPE")
                     .add_optional_named_double_vector("FIBER1", 3)
                     .Build();

  defs["HEX20"] = Input::LineDefinition::Builder()
                      .AddIntVector("HEX20", 20)
                      .AddNamedInt("MAT")
                      .AddNamedString("TYPE")
                      .add_optional_named_double_vector("FIBER1", 3)
                      .Build();

  defs["HEX27"] = Input::LineDefinition::Builder()
                      .AddIntVector("HEX27", 27)
                      .AddNamedInt("MAT")
                      .AddNamedString("TYPE")
                      .add_optional_named_double_vector("FIBER1", 3)
                      .Build();

  defs["NURBS27"] = Input::LineDefinition::Builder()
                        .AddIntVector("NURBS27", 27)
                        .AddNamedInt("MAT")
                        .AddNamedString("TYPE")
                        .add_optional_named_double_vector("FIBER1", 3)
                        .Build();

  defs["NURBS8"] = Input::LineDefinition::Builder()
                       .AddIntVector("NURBS8", 8)
                       .AddNamedInt("MAT")
                       .AddNamedString("TYPE")
                       .add_optional_named_double_vector("FIBER1", 3)
                       .Build();

  defs["TET4"] = Input::LineDefinition::Builder()
                     .AddIntVector("TET4", 4)
                     .AddNamedInt("MAT")
                     .AddNamedString("TYPE")
                     .add_optional_named_double_vector("FIBER1", 3)
                     .Build();

  defs["TET10"] = Input::LineDefinition::Builder()
                      .AddIntVector("TET10", 10)
                      .AddNamedInt("MAT")
                      .AddNamedString("TYPE")
                      .add_optional_named_double_vector("FIBER1", 3)
                      .Build();

  defs["WEDGE6"] = Input::LineDefinition::Builder()
                       .AddIntVector("WEDGE6", 6)
                       .AddNamedInt("MAT")
                       .AddNamedString("TYPE")
                       .add_optional_named_double_vector("FIBER1", 3)
                       .Build();

  defs["WEDGE15"] = Input::LineDefinition::Builder()
                        .AddIntVector("WEDGE15", 15)
                        .AddNamedInt("MAT")
                        .AddNamedString("TYPE")
                        .add_optional_named_double_vector("FIBER1", 3)
                        .Build();

  defs["PYRAMID5"] = Input::LineDefinition::Builder()
                         .AddIntVector("PYRAMID5", 5)
                         .AddNamedInt("MAT")
                         .AddNamedString("TYPE")
                         .add_optional_named_double_vector("FIBER1", 3)
                         .Build();

  defs["QUAD4"] = Input::LineDefinition::Builder()
                      .AddIntVector("QUAD4", 4)
                      .AddNamedInt("MAT")
                      .AddNamedString("TYPE")
                      .add_optional_named_double_vector("FIBER1", 3)
                      .Build();

  defs["QUAD8"] = Input::LineDefinition::Builder()
                      .AddIntVector("QUAD8", 8)
                      .AddNamedInt("MAT")
                      .AddNamedString("TYPE")
                      .add_optional_named_double_vector("FIBER1", 3)
                      .Build();

  defs["QUAD9"] = Input::LineDefinition::Builder()
                      .AddIntVector("QUAD9", 9)
                      .AddNamedInt("MAT")
                      .AddNamedString("TYPE")
                      .add_optional_named_double_vector("FIBER1", 3)
                      .Build();

  defs["TRI3"] = Input::LineDefinition::Builder()
                     .AddIntVector("TRI3", 3)
                     .AddNamedInt("MAT")
                     .AddNamedString("TYPE")
                     .add_optional_named_double_vector("FIBER1", 3)
                     .Build();

  defs["TRI6"] = Input::LineDefinition::Builder()
                     .AddIntVector("TRI6", 6)
                     .AddNamedInt("MAT")
                     .AddNamedString("TYPE")
                     .add_optional_named_double_vector("FIBER1", 3)
                     .Build();

  defs["NURBS4"] = Input::LineDefinition::Builder()
                       .AddIntVector("NURBS4", 4)
                       .AddNamedInt("MAT")
                       .AddNamedString("TYPE")
                       .add_optional_named_double_vector("FIBER1", 3)
                       .Build();

  defs["NURBS9"] = Input::LineDefinition::Builder()
                       .AddIntVector("NURBS9", 9)
                       .AddNamedInt("MAT")
                       .AddNamedString("TYPE")
                       .add_optional_named_double_vector("FIBER1", 3)
                       .Build();

  defs["LINE2"] = Input::LineDefinition::Builder()
                      .AddIntVector("LINE2", 2)
                      .AddNamedInt("MAT")
                      .AddNamedString("TYPE")
                      .add_optional_named_double_vector("FIBER1", 3)
                      .Build();

  defs["LINE3"] = Input::LineDefinition::Builder()
                      .AddIntVector("LINE3", 3)
                      .AddNamedInt("MAT")
                      .AddNamedString("TYPE")
                      .add_optional_named_double_vector("FIBER1", 3)
                      .Build();

  defs["NURBS2"] = Input::LineDefinition::Builder()
                       .AddIntVector("NURBS2", 2)
                       .AddNamedInt("MAT")
                       .AddNamedString("TYPE")
                       .add_optional_named_double_vector("FIBER1", 3)
                       .Build();

  defs["NURBS3"] = Input::LineDefinition::Builder()
                       .AddIntVector("NURBS3", 3)
                       .AddNamedInt("MAT")
                       .AddNamedString("TYPE")
                       .add_optional_named_double_vector("FIBER1", 3)
                       .Build();
}

/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::TransportType::Initialize(Core::FE::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    Discret::ELEMENTS::Transport* actele =
        dynamic_cast<Discret::ELEMENTS::Transport*>(dis.lColElement(i));
    if (!actele) FOUR_C_THROW("cast to Transport element failed");
    actele->initialize();
  }
  return 0;
}


Discret::ELEMENTS::TransportBoundaryType Discret::ELEMENTS::TransportBoundaryType::instance_;

Discret::ELEMENTS::TransportBoundaryType& Discret::ELEMENTS::TransportBoundaryType::Instance()
{
  return instance_;
}

Teuchos::RCP<Core::Elements::Element> Discret::ELEMENTS::TransportBoundaryType::Create(
    const int id, const int owner)
{
  // return Teuchos::rcp( new TransportBoundary( id, owner ) );
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                             gjb 05/08 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Transport::Transport(int id, int owner)
    : Core::Elements::Element(id, owner),
      distype_(Core::FE::CellType::dis_none),
      name_(),
      vis_map_(),
      numdofpernode_(-1),
      impltype_(Inpar::ScaTra::impltype_undefined)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        gjb 05/08 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::Transport::Transport(const Discret::ELEMENTS::Transport& old)
    : Core::Elements::Element(old),
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
Core::Elements::Element* Discret::ELEMENTS::Transport::Clone() const
{
  Discret::ELEMENTS::Transport* newelement = new Discret::ELEMENTS::Transport(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |  create material class (public)                            gjb 07/08 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Transport::SetMaterial(
    const int index, Teuchos::RCP<Core::Mat::Material> mat)
{
  // the standard part:
  Core::Elements::Element::SetMaterial(index, mat);

  if (mat->MaterialType() == Core::Materials::m_scatra or
      mat->MaterialType() == Core::Materials::m_scatra_multiscale or
      mat->MaterialType() == Core::Materials::m_myocard or
      mat->MaterialType() == Core::Materials::m_sutherland or
      mat->MaterialType() == Core::Materials::m_ion or
      mat->MaterialType() == Core::Materials::m_th_fourier_iso or
      mat->MaterialType() == Core::Materials::m_thermostvenant or
      mat->MaterialType() == Core::Materials::m_soret or
      mat->MaterialType() == Core::Materials::m_scatra_multiporo_fluid or
      mat->MaterialType() == Core::Materials::m_scatra_multiporo_volfrac or
      mat->MaterialType() == Core::Materials::m_scatra_multiporo_solid or
      mat->MaterialType() == Core::Materials::m_scatra_multiporo_temperature or
      (mat->MaterialType() == Core::Materials::m_electrode and
          impltype_ == Inpar::ScaTra::impltype_std))
    numdofpernode_ = 1;  // we only have a single scalar
  else if (mat->MaterialType() == Core::Materials::m_electrode)
    numdofpernode_ = 2;  // concentration and electric potential
  else if (mat->MaterialType() == Core::Materials::m_matlist)  // we have a system of scalars
  {
    const Mat::MatList* actmat = static_cast<const Mat::MatList*>(mat.get());
    numdofpernode_ = actmat->NumMat();

    // for problem type ELCH we have one additional degree of freedom per node
    // for the electric potential
    if (Global::Problem::Instance()->GetProblemType() == Core::ProblemType::elch)
    {
      for (int ii = 0; ii < numdofpernode_; ++ii)
      {
        // In the context of ELCH the only valid material combination is m_matlist and m_ion
        if (actmat->MaterialById(actmat->MatID(ii))->MaterialType() != Core::Materials::m_ion)
          FOUR_C_THROW(
              "In the context of ELCH the material Mat_matlist can be only used in combination "
              "with Mat_ion");
      }
      numdofpernode_ += 1;
    }
  }
  else if (mat->MaterialType() ==
           Core::Materials::m_matlist_reactions)  // we have a system of reactive scalars
  {
    // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are in
    // a diamond inheritance structure
    const Mat::MatListReactions* actmat = dynamic_cast<const Mat::MatListReactions*>(mat.get());
    numdofpernode_ = actmat->NumMat();

    for (int ii = 0; ii < numdofpernode_; ++ii)
    {
      // In the context of reactions the only valid material combination is m_matlist and m_scatra
      if (actmat->MaterialById(actmat->MatID(ii))->MaterialType() != Core::Materials::m_scatra and
          actmat->MaterialById(actmat->MatID(ii))->MaterialType() !=
              Core::Materials::m_scatra_multiporo_fluid and
          actmat->MaterialById(actmat->MatID(ii))->MaterialType() !=
              Core::Materials::m_scatra_multiporo_volfrac and
          actmat->MaterialById(actmat->MatID(ii))->MaterialType() !=
              Core::Materials::m_scatra_multiporo_temperature and
          actmat->MaterialById(actmat->MatID(ii))->MaterialType() !=
              Core::Materials::m_scatra_multiporo_solid)
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
              Core::Materials::m_scatra_reaction and
          actmat->MaterialById(actmat->ReacID(jj))->MaterialType() !=
              Core::Materials::m_scatra_reaction_poroECM)
        FOUR_C_THROW(
            "The material MAT_matlist_reaction only supports MAT_scatra_reaction and "
            "MAT_scatra_reaction_poro as valid reaction Material");

      // some safty check for the MAT_scatra_reaction materials
      const Teuchos::RCP<const Mat::ScatraReactionMat>& reacmat =
          Teuchos::rcp_static_cast<const Mat::ScatraReactionMat>(
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
           Core::Materials::m_matlist_chemotaxis)  // we have a system of chemotactic scalars
  {
    // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are in
    // a diamond inheritance structure
    const Mat::MatListChemotaxis* actmat = dynamic_cast<const Mat::MatListChemotaxis*>(mat.get());
    numdofpernode_ = actmat->NumMat();

    for (int ii = 0; ii < numdofpernode_; ++ii)
    {
      // In the context of chemotaxis the only valid material combination is m_matlist and m_scatra
      if (actmat->MaterialById(actmat->MatID(ii))->MaterialType() != Core::Materials::m_scatra)
        FOUR_C_THROW(
            "The material Mat_matlist_chemotaxis only supports MAT_scatra as valid main Material");
    }

    int numpair = actmat->NumPair();
    for (int jj = 0; jj < numpair; ++jj)
    {
      // In the context of chemotaxis the only valid material combination is m_matlist and
      // m_scatra_chemotaxis
      if (actmat->MaterialById(actmat->PairID(jj))->MaterialType() !=
          Core::Materials::m_scatra_chemotaxis)
        FOUR_C_THROW(
            "The material MAT_matlist_chemotaxis only supports MAT_scatra_chemotaxis as valid "
            "reaction Material");

      // some safty check for the MAT_scatra_chemotaxis materials
      const Teuchos::RCP<const Mat::ScatraChemotaxisMat>& reacmat =
          Teuchos::rcp_static_cast<const Mat::ScatraChemotaxisMat>(
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
           Core::Materials::m_matlist_chemoreac)  // we have a system of chemotactic scalars
  {
    // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are in
    // a diamond inheritance structure
    const Mat::MatListChemoReac* actmat = dynamic_cast<const Mat::MatListChemoReac*>(mat.get());
    numdofpernode_ = actmat->NumMat();

    for (int ii = 0; ii < numdofpernode_; ++ii)
    {
      // In the context of reactions/chemotaxis the only valid material combination is m_matlist and
      // m_scatra
      if (actmat->MaterialById(actmat->MatID(ii))->MaterialType() != Core::Materials::m_scatra)
        FOUR_C_THROW(
            "The material Mat_matlist_chemoreac only supports MAT_scatra as valid main Material");
    }

    int numreac = actmat->NumReac();
    for (int jj = 0; jj < numreac; ++jj)
    {
      // In the context of reactions the only valid material combination is m_matlist and
      // m_scatra_reaction
      if (actmat->MaterialById(actmat->ReacID(jj))->MaterialType() !=
          Core::Materials::m_scatra_reaction)
        FOUR_C_THROW(
            "The material MAT_matlist_reaction only supports MAT_scatra_reaction as valid reaction "
            "Material");

      // some safty check for the MAT_scatra_reaction materials
      const Teuchos::RCP<const Mat::ScatraReactionMat>& reacmat =
          Teuchos::rcp_static_cast<const Mat::ScatraReactionMat>(
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
          Core::Materials::m_scatra_chemotaxis)
        FOUR_C_THROW(
            "The material MAT_matlist_chemotaxis only supports MAT_scatra_chemotaxis as valid "
            "reaction Material");

      // some safty check for the MAT_scatra_chemotaxis materials
      const Teuchos::RCP<const Mat::ScatraChemotaxisMat>& reacmat =
          Teuchos::rcp_static_cast<const Mat::ScatraChemotaxisMat>(
              actmat->MaterialById(actmat->PairID(jj)));
      const int pairlength = reacmat->Pair()->size();
      if (pairlength != numdofpernode_)
        FOUR_C_THROW(
            "The number of scalars in your MAT_scatra_chemotaxis material with ID %i does not fit "
            "to the number of scalars!",
            actmat->PairID(jj));
    }
  }
  else if (mat->MaterialType() == Core::Materials::m_elchmat)
  {
    const Mat::ElchMat* actmat = static_cast<const Mat::ElchMat*>(mat.get());

    numdofpernode_ = actmat->NumDOF();
  }
  else
    FOUR_C_THROW("Transport element got unsupported material type %d", mat->MaterialType());

  return;
}

/*----------------------------------------------------------------------*
 |  create material class (public)                            gjb 07/08 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Transport::SetMaterial(int matnum, Core::Elements::Element* oldele)
{
  SetMaterial(0, Mat::Factory(matnum));

  Teuchos::RCP<Core::Mat::Material> mat = Material();

  if (mat->MaterialType() == Core::Materials::m_myocard)
  {
    Teuchos::RCP<Mat::Myocard> actmat = Teuchos::rcp_dynamic_cast<Mat::Myocard>(mat);

    Teuchos::RCP<Mat::ElastHyper> somat =
        Teuchos::rcp_dynamic_cast<Mat::ElastHyper>(oldele->Material());
    if (somat == Teuchos::null) FOUR_C_THROW("cast to ElastHyper failed");

    // copy fiber information from solid material to scatra material (for now, only one fiber
    // vector)
    std::vector<Core::LinAlg::Matrix<3, 1>> fibervecs(0);
    somat->GetFiberVecs(fibervecs);
    actmat->Setup(fibervecs[0]);
  }
}

/*----------------------------------------------------------------------*
 |  Return the shape of a Transport element                      (public) |
 |                                                            gjb 05/08 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::Transport::Shape() const { return distype_; }

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gjb 05/08 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Transport::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
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
void Discret::ELEMENTS::Transport::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Element::Unpack(basedata);

  // extract internal data
  ExtractfromPack(position, data, name_);
  ExtractfromPack(position, data, vis_map_);
  ExtractfromPack(position, data, numdofpernode_);
  distype_ = static_cast<Core::FE::CellType>(ExtractInt(position, data));
  impltype_ = static_cast<Inpar::ScaTra::ImplType>(ExtractInt(position, data));

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", (int)data.size(), position);
}

/*----------------------------------------------------------------------*
 |  Return number of lines of this element (public)           gjb 07/08 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::Transport::NumLine() const
{
  return Core::FE::getNumberOfElementLines(distype_);
}


/*----------------------------------------------------------------------*
 |  Return number of surfaces of this element (public)        gjb 07/08 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::Transport::NumSurface() const
{
  return Core::FE::getNumberOfElementSurfaces(distype_);
}


/*----------------------------------------------------------------------*
 | Return number of volumes of this element (public)          gjb 07/08 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::Transport::NumVolume() const
{
  return Core::FE::getNumberOfElementVolumes(distype_);
}



/*----------------------------------------------------------------------*
 |  print this element (public)                               gjb 05/08 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Transport::Print(std::ostream& os) const
{
  os << "Transport element";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  " << Core::FE::CellTypeToString(distype_) << std::endl;
  std::cout << std::endl;
  std::cout << "Number DOF per Node: " << numdofpernode_ << std::endl;
  std::cout << std::endl;
  std::cout << "Type of scalar transport: " << ScaTra::ImplTypeToString(impltype_) << std::endl;
  std::cout << std::endl;
}


/*----------------------------------------------------------------------*
 |  get vector of lines            (public)                  g.bau 03/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Transport::Lines()
{
  return Core::Communication::GetElementLines<TransportBoundary, Transport>(*this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          g.bau 03/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Transport::Surfaces()
{
  return Core::Communication::GetElementSurfaces<TransportBoundary, Transport>(*this);
}

/*----------------------------------------------------------------------*
 | set implementation type                                   fang 02/15 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::Transport ::SetImplType(const Inpar::ScaTra::ImplType impltype)
{
  // set implementation type
  impltype_ = impltype;
}

/*----------------------------------------------------------------------*
 |  init the element                                        vuong08/16 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::Transport::initialize()
{
  Teuchos::RCP<Core::Mat::Material> mat = Material();
  // for now, we only need to do something in case of reactions (for the initialization of functions
  // in case of reactions by function)
  if (mat->MaterialType() == Core::Materials::m_matlist_reactions or
      mat->MaterialType() == Core::Materials::m_matlist_chemoreac)
  {
    // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction are in
    // a diamond inheritance structure
    Teuchos::RCP<Mat::MatListReactions> actmat =
        Teuchos::rcp_dynamic_cast<Mat::MatListReactions>(mat);
    actmat->Initialize();
  }
  else if (mat->MaterialType() == Core::Materials::m_myocard)
  {
    Teuchos::RCP<Mat::Myocard> actmat = Teuchos::rcp_dynamic_cast<Mat::Myocard>(mat);
    int deg = 0;
    if (this->Degree() == 1)
      deg = 4 * this->Degree();
    else
      deg = 3 * this->Degree();
    Teuchos::RCP<Core::FE::GaussPoints> quadrature(
        Core::FE::GaussPointCache::Instance().Create(this->Shape(), deg));
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
Discret::ELEMENTS::TransportBoundary::TransportBoundary(int id, int owner, int nnode,
    const int* nodeids, Core::Nodes::Node** nodes, Discret::ELEMENTS::Transport* parent,
    const int lsurface)
    : Core::Elements::FaceElement(id, owner)
{
  SetNodeIds(nnode, nodeids);
  BuildNodalPointers(nodes);
  set_parent_master_element(parent, lsurface);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        gjb 01/09 |
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::TransportBoundary::TransportBoundary(
    const Discret::ELEMENTS::TransportBoundary& old)
    : Core::Elements::FaceElement(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it     (public) gjb 01/09 |
 *----------------------------------------------------------------------*/
Core::Elements::Element* Discret::ELEMENTS::TransportBoundary::Clone() const
{
  Discret::ELEMENTS::TransportBoundary* newelement =
      new Discret::ELEMENTS::TransportBoundary(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Return shape of this element                    (public)  gjb 01/09 |
 *----------------------------------------------------------------------*/
Core::FE::CellType Discret::ELEMENTS::TransportBoundary::Shape() const
{
  return Core::FE::getShapeOfBoundaryElement(num_node(), parent_element()->Shape());
}

/*----------------------------------------------------------------------*
 |  Pack data (public)                                        gjb 01/09 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::TransportBoundary::Pack(Core::Communication::PackBuffer& data) const
{
  FOUR_C_THROW("This TransportBoundary element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data (public)                                      gjb 01/09 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::TransportBoundary::Unpack(const std::vector<char>& data)
{
  FOUR_C_THROW("This TransportBoundary element does not support communication");
  return;
}



/*----------------------------------------------------------------------*
 |  print this element (public)                               gjb 01/09 |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::TransportBoundary::Print(std::ostream& os) const
{
  os << "TransportBoundary element";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << "DiscretizationType:  " << Core::FE::CellTypeToString(Shape()) << std::endl;
  std::cout << std::endl;
  return;
}

/*----------------------------------------------------------------------*
 | Return number of lines of boundary element (public)        gjb 01/09 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::TransportBoundary::NumLine() const
{
  return Core::FE::getNumberOfElementLines(Shape());
}

/*----------------------------------------------------------------------*
 |  Return number of surfaces of boundary element (public)    gjb 01/09 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::TransportBoundary::NumSurface() const
{
  return Core::FE::getNumberOfElementSurfaces(Shape());
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                              gjb 01/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::TransportBoundary::Lines()
{
  FOUR_C_THROW("Lines of TransportBoundary not implemented");
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                              gjb 01/09 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::TransportBoundary::Surfaces()
{
  FOUR_C_THROW("Surfaces of TransportBoundary not implemented");
}

FOUR_C_NAMESPACE_CLOSE
