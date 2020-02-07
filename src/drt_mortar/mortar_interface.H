/*----------------------------------------------------------------------*/
/*! \file
\brief One mortar coupling interface

\level 1

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/
#ifndef MORTAR_INTERFACE_H
#define MORTAR_INTERFACE_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

#include "../drt_inpar/inpar_mortar.H"
#include "../drt_lib/drt_dserror.H"
#include "../headers/pairedvector.H"
#include "mortar_paramsinterface.H"

#include <boost/algorithm/string/predicate.hpp>  // case insensitive string compare


// forward declarations
namespace BINSTRATEGY
{
  class BinningStrategy;
}  // namespace BINSTRATEGY

namespace DRT
{
  class Discretization;
  class Element;
  class Exporter;
}  // namespace DRT

namespace LINALG
{
  class SparseMatrix;
}  // namespace LINALG

namespace MORTAR
{
  // forward declarations
  class MortarInterface;
  class MortarNode;
  class MortarElement;
  class IntElement;
  class BinaryTree;
  class IntCell;
  class ParamsInterface;

  /*! \brief State type for the set state routine
   *
   * \author hiermeier */
  enum StateType
  {
    state_old_displacement,            //!< old displacement state
    state_new_displacement,            //!< new displacement state
    state_lagrange_multiplier,         //!< lagrange multiplier (only necessary for poro)
    state_fvelocity,                   //!< fvelocity (only necessary for poro)
    state_svelocity,                   //!< svelocity (only necessary for poro)
    state_fpressure,                   //!< fpressure (only necessary for poro)
    state_temperature,                 //!< temperature (only necessary for TSI)
    state_thermo_lagrange_multiplier,  //!< thermo lagrange multiplier (only necessary for TSI)
    state_vague                        //!< undefined state type
  };

  /*! \brief Map state type enum to std::string
   *
   *  \author hiermeier */
  static inline std::string StateType2String(const enum StateType& type)
  {
    switch (type)
    {
      case state_old_displacement:
        return "displacement";
        break;
      case state_new_displacement:
        return "displacement new";
        break;
      case state_lagrange_multiplier:
        return "lm";
        break;
      case state_fvelocity:
        return "fveclocity";
        break;
      case state_svelocity:
        return "svelocity";
        break;
      case state_fpressure:
        return "fpressure";
        break;
      case state_temperature:
        return "temperature";
        break;
      case state_thermo_lagrange_multiplier:
        return "thermo_lm";
        break;
      default:
        return "unknown MORTAR::StateType";
        break;
    }
    return "";
  };

  /*! \brief Map std::string to state type enum
   *
   *  \author hiermeier */
  static inline enum MORTAR::StateType String2StateType(const std::string& name)
  {
    MORTAR::StateType type = state_vague;
    if (boost::iequals(name, "displacement"))
      type = state_new_displacement;
    else if (boost::iequals(name, "olddisplacement"))
      type = state_old_displacement;
    else if (boost::iequals(name, "lm"))
      type = state_lagrange_multiplier;
    else if (boost::iequals(name, "fvelocity"))
      type = state_fvelocity;
    else if (boost::iequals(name, "svelocity"))
      type = state_svelocity;
    else if (boost::iequals(name, "fpressure"))
      type = state_fpressure;
    else if (boost::iequals(name, "temperature"))
      type = state_temperature;
    else if (boost::iequals(name, "thermo_lm"))
      type = state_thermo_lagrange_multiplier;
    else
    {
      dserror("Unknown state type name. (name = %s)", name.c_str());
    }
    return type;
  };

  /*----------------------------------------------------------------------------*/
  /** \brief Mortar interface data container
   *
   *  This class is supposed to contain all relevant members for the mortar
   *  contact/mesh-tying interfaces. The external storage in this object, instead
   *  of the actual interface class itself, makes it possible to share the interface
   *  data between different interface objects w/o the need of copying them.
   *
   *  \author hiermeier \date 03/17 */
  class InterfaceDataContainer
  {
   public:
    /// constructor
    InterfaceDataContainer();

    /// destructor
    virtual ~InterfaceDataContainer() = default;

    /// @name accessors
    /// @{

    inline int& Id() { return id_; }

    inline int Id() const { return id_; }

    inline Teuchos::RCP<const Epetra_Comm>& CommPtr() { return comm_; }

    inline const Epetra_Comm& Comm() const { return *comm_; }

    inline Teuchos::RCP<Epetra_Comm>& lComm() { return lcomm_; }

    inline Teuchos::RCP<const Epetra_Comm> lComm() const { return lcomm_; }

    inline std::map<int, int>& ProcMap() { return procmap_; }

    inline const std::map<int, int>& ProcMap() const { return procmap_; }

    inline bool& IsRedistributed() { return redistributed_; }

    inline bool IsRedistributed() const { return redistributed_; }

    inline Teuchos::RCP<DRT::Discretization>& IDiscret() { return idiscret_; }

    inline Teuchos::RCP<const DRT::Discretization> IDiscret() const { return idiscret_; }

    inline int& Dim() { return dim_; }

    inline int Dim() const { return dim_; }

    inline Teuchos::ParameterList& IMortar() { return imortar_; }

    inline const Teuchos::ParameterList& IMortar() const { return imortar_; }

    inline enum INPAR::MORTAR::ShapeFcn& ShapeFcn() { return shapefcn_; }

    inline enum INPAR::MORTAR::ShapeFcn ShapeFcn() const { return shapefcn_; }

    inline bool& IsQuadSlave() { return quadslave_; }

    inline bool IsQuadSlave() const { return quadslave_; }

    inline enum INPAR::MORTAR::RedundantStorage& RedundantStorage() { return redundant_; }

    inline enum INPAR::MORTAR::RedundantStorage RedundantStorage() const { return redundant_; }

    inline Teuchos::RCP<Epetra_Map>& OldNodeColMap() { return oldnodecolmap_; }

    inline Teuchos::RCP<const Epetra_Map> OldNodeColMap() const { return oldnodecolmap_; }

    inline Teuchos::RCP<Epetra_Map>& OldEleColMap() { return oldelecolmap_; }

    inline Teuchos::RCP<const Epetra_Map> OldEleColMap() const { return oldelecolmap_; }

    inline Teuchos::RCP<Epetra_Map>& SNodeRowMap() { return snoderowmap_; }

    inline Teuchos::RCP<const Epetra_Map> SNodeRowMap() const { return snoderowmap_; }

    inline Teuchos::RCP<Epetra_Map>& SNodeColMap() { return snodecolmap_; }

    inline Teuchos::RCP<const Epetra_Map> SNodeColMap() const { return snodecolmap_; }

    inline Teuchos::RCP<Epetra_Map>& MNodeRowMap() { return mnoderowmap_; }

    inline Teuchos::RCP<const Epetra_Map> MNodeRowMap() const { return mnoderowmap_; }

    inline Teuchos::RCP<Epetra_Map>& PSNodeRowMap() { return psnoderowmap_; }

    inline Teuchos::RCP<const Epetra_Map> PSNodeRowMap() const { return psnoderowmap_; }

    inline Teuchos::RCP<Epetra_Map>& PMNodeRowMap() { return pmnoderowmap_; }

    inline Teuchos::RCP<const Epetra_Map> PMNodeRowMap() const { return pmnoderowmap_; }

    inline Teuchos::RCP<Epetra_Map>& MNodeColMap() { return mnodecolmap_; }

    inline Teuchos::RCP<const Epetra_Map> MNodeColMap() const { return mnodecolmap_; }

    inline Teuchos::RCP<Epetra_Map>& SNodeRowMapBound() { return snoderowmapbound_; }

    inline Teuchos::RCP<const Epetra_Map> SNodeRowMapBound() const { return snoderowmapbound_; }


    inline Teuchos::RCP<Epetra_Map>& SNodeColMapBound() { return snodecolmapbound_; }

    inline Teuchos::RCP<const Epetra_Map> SNodeColMapBound() const { return snodecolmapbound_; }

    inline Teuchos::RCP<Epetra_Map>& MNodeRowMapNoBound() { return mnoderowmapnobound_; }

    inline Teuchos::RCP<const Epetra_Map> MNodeRowMapNoBound() const { return mnoderowmapnobound_; }

    inline Teuchos::RCP<Epetra_Map>& MNodeColMapNoBound() { return mnodecolmapnobound_; }

    inline Teuchos::RCP<const Epetra_Map> MNodeColMapNoBound() const { return mnodecolmapnobound_; }

    inline Teuchos::RCP<Epetra_Map>& SEleRowMap() { return selerowmap_; }

    inline Teuchos::RCP<const Epetra_Map> SEleRowMap() const { return selerowmap_; }

    inline Teuchos::RCP<Epetra_Map>& SEleColMap() { return selecolmap_; }

    inline Teuchos::RCP<const Epetra_Map> SEleColMap() const { return selecolmap_; }

    inline Teuchos::RCP<Epetra_Map>& MEleRowMap() { return melerowmap_; }

    inline Teuchos::RCP<const Epetra_Map> MEleRowMap() const { return melerowmap_; }

    inline Teuchos::RCP<Epetra_Map>& MEleColMap() { return melecolmap_; }

    inline Teuchos::RCP<const Epetra_Map> MEleColMap() const { return melecolmap_; }

    inline Teuchos::RCP<Epetra_Map>& SDofRowMap() { return sdofrowmap_; }

    inline Teuchos::RCP<const Epetra_Map> SDofRowMap() const { return sdofrowmap_; }

    inline Teuchos::RCP<Epetra_Map>& SDofColMap() { return sdofcolmap_; }

    inline Teuchos::RCP<const Epetra_Map> SDofColMap() const { return sdofcolmap_; }

    inline Teuchos::RCP<Epetra_Map>& MDofRowMap() { return mdofrowmap_; }

    inline Teuchos::RCP<const Epetra_Map> MDofRowMap() const { return mdofrowmap_; }

    inline Teuchos::RCP<Epetra_Map>& MDofColMap() { return mdofcolmap_; }

    inline Teuchos::RCP<const Epetra_Map> MDofColMap() const { return mdofcolmap_; }

    inline Teuchos::RCP<Epetra_Map>& PSDofRowMap() { return psdofrowmap_; }

    inline Teuchos::RCP<Epetra_Map>& PMDofRowMap() { return pmdofrowmap_; }

    inline Teuchos::RCP<const Epetra_Map> PMDofRowMap() const { return pmdofrowmap_; }

    inline Teuchos::RCP<const Epetra_Map> PSDofRowMap() const { return psdofrowmap_; }

    inline Teuchos::RCP<Epetra_Map>& PLmDofRowMap() { return plmdofmap_; }

    inline Teuchos::RCP<const Epetra_Map> PLmDofRowMap() const { return plmdofmap_; }

    inline Teuchos::RCP<Epetra_Map>& LmDofRowMap() { return lmdofmap_; }

    inline Teuchos::RCP<const Epetra_Map> LmDofRowMap() const { return lmdofmap_; }

    inline int& MaxDofGlobal() { return maxdofglobal_; }

    inline int MaxDofGlobal() const { return maxdofglobal_; }

    inline INPAR::MORTAR::SearchAlgorithm& SearchAlgorithm() { return searchalgo_; }

    inline enum INPAR::MORTAR::SearchAlgorithm SearchAlgorithm() const { return searchalgo_; }

    inline Teuchos::RCP<MORTAR::BinaryTree>& BinaryTree() { return binarytree_; }

    inline Teuchos::RCP<const MORTAR::BinaryTree> BinaryTree() const { return binarytree_; }

    inline double& SearchParam() { return searchparam_; }

    inline double SearchParam() const { return searchparam_; }

    inline bool& SearchUseAuxPos() { return searchuseauxpos_; }

    inline bool SearchUseAuxPos() const { return searchuseauxpos_; }

    inline double& IntTimeInterface() { return inttime_interface_; }

    inline double IntTimeInterface() const { return inttime_interface_; }

    inline bool& IsNurbs() { return nurbs_; }

    inline bool IsNurbs() const { return nurbs_; }

    inline bool& IsPoro() { return poro_; }

    inline bool IsPoro() const { return poro_; }

    inline Teuchos::RCP<DRT::Exporter>& SlExporterPtr() { return sl_exporter_; }

    inline DRT::Exporter& Exporter()
    {
      if (sl_exporter_.is_null()) dserror("The exporter has not been initialized.");

      return *sl_exporter_;
    }

    void SetMaSharingRefInterfacePtr(const MortarInterface* interface_ptr)
    {
      masharingrefinterface_ = interface_ptr;
    }

    const MortarInterface* GetMaSharingRefInterfacePtr() const { return masharingrefinterface_; }

    const MortarInterface& GetMaSharingRefInterface() const
    {
      if (masharingrefinterface_ == NULL) dserror("NULL pointer");
      return *masharingrefinterface_;
    }

    inline bool& IsEhl() { return ehl_; }

    inline bool IsEhl() const { return ehl_; }

    /// @}

    /// @name DataContainer specific functions
    /// @{

    inline bool IsInit() const { return isinit_; }

    inline void SetIsInit(bool isinit) { isinit_ = isinit; }

    /// @}

   private:
    //! unique interface id
    int id_;

    //! @name Communication and parallelism
    //! @{

    //! communicator
    Teuchos::RCP<const Epetra_Comm> comm_;

    //! interface local communicator
    Teuchos::RCP<Epetra_Comm> lcomm_;

    //! mapping global -> local communicator PIDs
    std::map<int, int> procmap_;

    //! @}

    //! redistribution for this time step?
    bool redistributed_;

    //! master sharing reference interface ptr
    const MortarInterface* masharingrefinterface_ = NULL;

    //! the discretization of the mortar interface
    Teuchos::RCP<DRT::Discretization> idiscret_;

    //! Spatial dimension of problem (2D or 3D)
    int dim_;

    //! containing contact input parameters of interface
    Teuchos::ParameterList imortar_;

    //! employed type of shape function set
    INPAR::MORTAR::ShapeFcn shapefcn_;

    //! flag indicating quadratic 2d/3d slave elements
    bool quadslave_;

    //! employed type of redundancy in storage of interface
    INPAR::MORTAR::RedundantStorage redundant_;

    //! @name Maps
    //! @{

    //! column map of all interface nodes (overlap=1)
    Teuchos::RCP<Epetra_Map> oldnodecolmap_;

    //! column map of all interface elements (overlap=1)
    Teuchos::RCP<Epetra_Map> oldelecolmap_;

    //! row map of all slave nodes
    Teuchos::RCP<Epetra_Map> snoderowmap_;

    //! column map of all slave nodes
    Teuchos::RCP<Epetra_Map> snodecolmap_;

    //! row map of all master nodes
    Teuchos::RCP<Epetra_Map> mnoderowmap_;

    //! column map of all master nodes
    Teuchos::RCP<Epetra_Map> mnodecolmap_;

    //! row map of slave nodes (+ boundary nodes)
    Teuchos::RCP<Epetra_Map> snoderowmapbound_;

    //! col map of slave nodes (+ boundary nodes)
    Teuchos::RCP<Epetra_Map> snodecolmapbound_;

    //! row map of master nodes (- boundary nodes)
    Teuchos::RCP<Epetra_Map> mnoderowmapnobound_;

    //! col map of master nodes (- boundary nodes)
    Teuchos::RCP<Epetra_Map> mnodecolmapnobound_;

    //! row map of all slave elements
    Teuchos::RCP<Epetra_Map> selerowmap_;

    //! column map of all slave elements
    Teuchos::RCP<Epetra_Map> selecolmap_;

    //! row map of all master elements
    Teuchos::RCP<Epetra_Map> melerowmap_;

    //! column map of all master elements
    Teuchos::RCP<Epetra_Map> melecolmap_;

    //! row map of all slave dofs
    Teuchos::RCP<Epetra_Map> sdofrowmap_;

    //! column map of all slave dofs
    Teuchos::RCP<Epetra_Map> sdofcolmap_;

    //! row map of all master dofs
    Teuchos::RCP<Epetra_Map> mdofrowmap_;

    //! column map of all master dofs
    Teuchos::RCP<Epetra_Map> mdofcolmap_;

    //! row map of all slave dofs before any redistribution took place
    Teuchos::RCP<Epetra_Map> psdofrowmap_;

    //! row map of all master dofs before any redistribution took place
    Teuchos::RCP<Epetra_Map> pmdofrowmap_;

    //! row map of all Lagrange multiplier dofs before any redistribution took place
    Teuchos::RCP<Epetra_Map> plmdofmap_;

    //! row map of all slave nodes before any redistribution took place
    Teuchos::RCP<Epetra_Map> psnoderowmap_;

    //! row map of all master nodes before any redistribution took place
    Teuchos::RCP<Epetra_Map> pmnoderowmap_;

    //! row map of all Lagrange multiplier dofs
    Teuchos::RCP<Epetra_Map> lmdofmap_;

    //! @}

    //! maximum dof ID in global discretization
    int maxdofglobal_;

    //! @name Search algorithm
    //! @{

    //! type of search algorithm
    INPAR::MORTAR::SearchAlgorithm searchalgo_;

    //! binary searchtree
    Teuchos::RCP<MORTAR::BinaryTree> binarytree_;

    //! search parameter
    double searchparam_;

    //! use auxiliary position when computing dops
    bool searchuseauxpos_;

    //! @}

    //! integration time
    double inttime_interface_;

    //! flag for nurbs shape functions
    bool nurbs_;

    //! flag if poro contact problem!
    bool poro_;

    //! flag if ehl contact problem!
    bool ehl_;

    Teuchos::RCP<DRT::Exporter> sl_exporter_;

    //! Is this data container object initialized?
    bool isinit_;

  };  // class InterfaceDataContainer

  /*----------------------------------------------------------------------------*/
  /*!
  \brief One mortar coupling interface

  \author gee (gee@lnm.mw.tum.de), popp (popp@lnm.mw.tum.de)
  */
  class MortarInterface
  {
   protected:
    /*!
    \brief Constructor (only for derived classes)

    @param[in] interfaceData_ptr Interface data container
    */
    MortarInterface(const Teuchos::RCP<MORTAR::InterfaceDataContainer>& interfaceData_ptr);

   public:
    /** \brief Create a new Mortar interface object
     *
     *  This method creates first a new interface data object and subsequently
     *  a new interface object.
     *
     *  \param id (in) : unique interface ID
     *  \param comm (in) : communicator object
     *  \param spatialDim (in) : spatial dimension of the problem
     *  \param imortar (in) : global contact/mesh-tying parameter-list
     *  \param redundant (in) : type of redundancy in interface storage
     *
     *  \author hiermeier \date 03/17 */
    static Teuchos::RCP<MortarInterface> Create(const int id, const Epetra_Comm& comm,
        const int spatialDim, const Teuchos::ParameterList& imortar,
        INPAR::MORTAR::RedundantStorage redundant);

    /*!
    \brief Standard constructor creating empty mortar interface

    \param interfaceData_ptr (in): Interface data container
    \param id (in): Unique interface id
    \param comm (in): A communicator object
    \param spatialDim (in): Global problem dimension
    \param imortar (in): Global contact parameter list
    \param redundant (in): type of redundancy in interface storage
    */
    MortarInterface(const Teuchos::RCP<InterfaceDataContainer>& interfaceData_ptr, const int id,
        const Epetra_Comm& comm, const int spatialDim, const Teuchos::ParameterList& imortar,
        INPAR::MORTAR::RedundantStorage redundant);

    /*!
    \brief Destructor
    */
    virtual ~MortarInterface() = default;

    /*!
    \brief Get unique ID of this interface

    @return Unique interface ID
    */
    inline int Id() const { return id_; }

    /*!
    \brief Print this Interface

    \param[in] os Output stream used for printing
    */
    virtual void Print(std::ostream& os) const;

    /*!
    \brief Check whether interface was called FillComplete

    @return Boolean flag to indicate if FillComplete has been called
    */
    bool Filled() const;

    /*!
    \brief Print parallel distribution of this Interface
    */
    void PrintParallelDistribution() const;

    /*!
    \brief Get communicator
    */
    const Epetra_Comm& Comm() const { return *comm_; }

    /*!
    \brief Get local interface communicator (Filled()==true is prerequisite)
           This Epetra_Comm only contains the participating processors of the interface
    */
    const Epetra_Comm* lComm() const { return lcomm_.get(); }

    /*!
    \brief Get global -> local processor map
    */
    std::map<int, int>& Procmap() { return procmap_; }

    //! @name Access methods

    /*!
    \brief Get discretization of this interface
    */
    DRT::Discretization& Discret() const { return *idiscret_; }

    /*!
    \brief Get problem dimension

    Note that only 2D and 3D are possible here as this refers to the global
    problem dimension. On interface level this corresponds to 1D interfaces
    (dim_==2) and 2D interfaces (dim_==3)!
    */
    int Dim() const { return dim_; };

    /*!
    \brief Get interface contact parameter list
    */
    Teuchos::ParameterList& InterfaceParams() { return imortar_; };
    const Teuchos::ParameterList& InterfaceParams() const { return imortar_; };

    /*!
    \brief Get quadratic 2d/3d slave element flag

    Returns TRUE if at least one higher-order 2d/3d slave element in interface.
    */
    const bool& Quadslave() { return quadslave_; };

    /*!
    \brief Get flag indicating nurbs shape functions
    */
    bool IsNurbs() { return nurbs_; };

    /*!
    \brief Get redundant storage flag

    Returns status of redundancy in interface storage (nodes, elements).
    */
    INPAR::MORTAR::RedundantStorage Redundant() const { return redundant_; };

    /*!
    \brief Get type of search algorithm
    */
    INPAR::MORTAR::SearchAlgorithm SearchAlg() const { return searchalgo_; };

    /*!
    \brief Get search algorithm parameter
    */
    double SearchParam() const { return searchparam_; };

    /*!
    \brief bool whether auxiliary position is used when computing dops
    */
    bool SearchUseAuxPos() const { return searchuseauxpos_; };

    /*!
    \brief Get column map of all interface nodes (Filled()==true is prerequisite)
    */
    const Teuchos::RCP<Epetra_Map>& OldColNodes() const
    {
      if (Filled())
        return oldnodecolmap_;
      else
        dserror("MORTAR::MortarInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get column map of all interface elements (Filled()==true is prerequisite)
    */
    const Teuchos::RCP<Epetra_Map>& OldColElements() const
    {
      if (Filled())
        return oldelecolmap_;
      else
        dserror("MORTAR::MortarInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of slave nodes (Filled()==true is prerequisite)
    */
    const Teuchos::RCP<Epetra_Map>& SlaveRowNodes() const
    {
      if (Filled())
        return snoderowmap_;
      else
        dserror("MORTAR::MortarInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of slave nodes (Filled()==true is prerequisite)

    \note This is the slave row map before any parallel redistribution took place.
    */
    const Teuchos::RCP<Epetra_Map>& PSlaveRowNodes() const
    {
      if (not Filled()) dserror("MORTAR::MortarInterface::FillComplete was not called");

      return interfaceData_->PSNodeRowMap();
    }

    /*!
    \brief Get row map of master nodes (Filled()==true is prerequisite)
    */
    const Teuchos::RCP<Epetra_Map>& MasterRowNodes() const
    {
      if (Filled())
        return mnoderowmap_;
      else
        dserror("MORTAR::MortarInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of master nodes (Filled()==true is prerequisite)

    \note This is the master row map before any parallel redistribution took place.
    */
    const Teuchos::RCP<Epetra_Map>& PMasterRowNodes() const
    {
      if (not Filled()) dserror("MORTAR::MortarInterface::FillComplete was not called");

      return interfaceData_->PMNodeRowMap();
    }

    /*!
    \brief Get column map of slave nodes (Filled()==true is prerequisite)
    */
    const Teuchos::RCP<Epetra_Map>& SlaveColNodes() const
    {
      if (Filled())
        return snodecolmap_;
      else
        dserror("MORTAR::MortarInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get column map of master nodes (Filled()==true is prerequisite)
    */
    const Teuchos::RCP<Epetra_Map>& MasterColNodes() const
    {
      if (Filled())
        return mnodecolmap_;
      else
        dserror("MORTAR::MortarInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of slave nodes + boundary nodes (Filled()==true is prerequisite)
    */
    const Teuchos::RCP<Epetra_Map>& SlaveRowNodesBound() const
    {
      if (Filled())
        return snoderowmapbound_;
      else
        dserror("MORTAR::MortarInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get column map of slave nodes + boundary nodes (Filled()==true is prerequisite)
    */
    const Teuchos::RCP<Epetra_Map>& SlaveColNodesBound() const
    {
      if (Filled())
        return snodecolmapbound_;
      else
        dserror("MORTAR::MortarInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of master nodes - boundary nodes (Filled()==true is prerequisite)
    */
    const Teuchos::RCP<Epetra_Map>& MasterRowNodesNoBound() const
    {
      if (Filled())
        return mnoderowmapnobound_;
      else
        dserror("MORTAR::MortarInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get column map of master nodes - boundary nodes (Filled()==true is prerequisite)
    */
    const Teuchos::RCP<Epetra_Map>& MasterColNodesNoBound() const
    {
      if (Filled())
        return mnodecolmapnobound_;
      else
        dserror("MORTAR::MortarInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of slave elements (Filled()==true is prerequisite)
    */
    const Teuchos::RCP<Epetra_Map>& SlaveRowElements() const
    {
      if (Filled())
        return selerowmap_;
      else
        dserror("MORTAR::MortarInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of master elements (Filled()==true is prerequisite)
    */
    const Teuchos::RCP<Epetra_Map>& MasterRowElements() const
    {
      if (Filled())
        return melerowmap_;
      else
        dserror("MORTAR::MortarInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get column map of slave elements (Filled()==true is prerequisite)
    */
    const Teuchos::RCP<Epetra_Map>& SlaveColElements() const
    {
      if (Filled())
        return selecolmap_;
      else
        dserror("MORTAR::MortarInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get column map of master elements (Filled()==true is prerequisite)
    */
    const Teuchos::RCP<Epetra_Map>& MasterColElements() const
    {
      if (Filled())
        return melecolmap_;
      else
        dserror("MORTAR::MortarInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of slave dofs (Filled()==true is prerequisite)
    */
    const Teuchos::RCP<Epetra_Map>& SlaveRowDofs() const
    {
      if (Filled())
        return sdofrowmap_;
      else
        dserror("MORTAR::MortarInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get column map of slave dofs (Filled()==true is prerequisite)
    */
    const Teuchos::RCP<Epetra_Map>& SlaveColDofs() const
    {
      if (Filled())
        return sdofcolmap_;
      else
        dserror("MORTAR::MortarInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of slave dofs (Filled()==true is prerequisite)

    \note This is the slave dof map before any parallel redistribution took place.
    */
    const Teuchos::RCP<Epetra_Map>& PSlaveRowDofs() const
    {
      if (not Filled()) dserror("MORTAR::MortarInterface::FillComplete was not called");

      return interfaceData_->PSDofRowMap();
    }

    /*!
    \brief Get row map of master dofs (Filled()==true is prerequisite)
    */
    const Teuchos::RCP<Epetra_Map>& MasterRowDofs() const
    {
      if (Filled())
        return mdofrowmap_;
      else
        dserror("MORTAR::MortarInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get column map of master dofs (Filled()==true is prerequisite)
    */
    const Teuchos::RCP<Epetra_Map>& MasterColDofs() const
    {
      if (Filled())
        return mdofcolmap_;
      else
        dserror("MORTAR::MortarInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get row map of master dofs (Filled()==true is prerequisite)

    \note This is the master dof map before any parallel redistribution took place.
    */
    const Teuchos::RCP<Epetra_Map>& PMasterRowDofs() const
    {
      if (not Filled()) dserror("MORTAR::MortarInterface::FillComplete was not called");

      return interfaceData_->PMDofRowMap();
    }

    /*!
    \brief Get map of Lagrange multiplier dofs (Filled()==true is prerequisite)
    */
    const Teuchos::RCP<Epetra_Map>& LagMultDofs() const
    {
      if (Filled())
        return lmdofmap_;
      else
        dserror("MORTAR::MortarInterface::FillComplete was not called");
      exit(EXIT_FAILURE);  // calm down the compiler
    }

    /*!
    \brief Get maximum dof ID in global discretization
    */
    int MaxDofGlobal() const
    {
      if (maxdofglobal_ < 0) dserror("MaxDofGlobal not yet initialized");
      return maxdofglobal_;
    }

    //@}

    //! @name Evlauation methods

    /*!
    \brief Add a MORTAR::MortarNode to the interface (Filled()==true NOT prerequisite)

    \param mrtrnode (in): Teuchos::rcp to a mortar node

    \return Filled()==false
    */
    void AddMortarNode(Teuchos::RCP<MORTAR::MortarNode> mrtrnode);

    /*!
    \brief Add a MORTAR::MortarElement to the interface (Filled()==true is prerequisite)

    \param mrtrele (in): Teuchos::rcp to a mortar element

    \return Filled()==false

    */
    void AddMortarElement(Teuchos::RCP<MORTAR::MortarElement> mrtrele);

    /*!
    \brief delete all elements/nodes on one side of the interface

    \param slaveside (in): bool to choose side of interface that will be deleted

    */
    void RemoveSingleInterfaceSide(bool slaveside);

    /*!
    \brief Finalize construction of mortar interface

    The methods completes construction phase of a mortar interface. It creates all row/column maps
    of the mortar interface discretization.

    If we have arrived at the final parallel distribution, we have to ask the underlying
    DRT::Discretization to assign degrees of freedom. Since this is very expensive, let's do this
    only if requested by the user/algorithm.

    \sa ExtendInterfaceGhosting()

    \param isFinalParallelDistribution (in): Is this the final parallel distribution?
    \param maxdof (in): maximum dof ID in global discretization
                        (if default=0, then no lmdofrowmap is created)
    \param meanVelocity (in): mean velocity of this interface (for contact with ghosting by binning)

    \return Filled()==true
    */
    void FillComplete(const bool isFinalParallelDistribution, const int maxdof = 0,
        const double meanVelocity = 0.0);

    /*!
    \brief Redistribute interface among all procs

    When first creating a mortar interface, its parallel distribution
    is simply copied from the underlying problem discretization. This,
    of course, is not the optimal parallel distribution for evaluating
    the mortar coupling terms, as the interface ownership might be
    restricted to only very few processors. Moreover, no parallel
    scalability can be achieved with this procedure, because adding
    processors to the problem discretization does not automatically
    mean adding processors to the interface discretization.

    Thus, an independent parallel distribution of the interface is
    desirable, which divides the interface among all available
    processors. Redistribute() is the method to achieve this.
    Internally, we call ZOLTAN to re-partition both slave and
    master side of the interface independently. This results in new
    "optimal" node/element maps of the interface discretization.
    Note that after Redistribute(), we must call FillComplete() again.

    */
    virtual void Redistribute();

    /*!
    \brief Create binary search tree

    The methods creates a binary tree object for efficient search.

    */
    virtual void CreateSearchTree();

    /*!
      \brief Restrict slave sets to actual meshtying zone

      This update is done ONCE in the initialization phase, when the
      slave surface only partially overlaps with the master surface.
      It restricts the slave sets (nodes, dofs, but NOT the elements)
      to the actual meshtying zone for the whole simulation.

      */
    void RestrictSlaveSets();

    /*!
    \brief Update interface Lagrange multiplier sets

    This update is usually only done ONCE in the initialization phase
    and sets up the Lagrange mutliplier set (only dofs) for the whole
    simulation. Yet, in the case of self contact the sets need to be
    updated again and again during simulation time, as the slave/master
    status and thus the LM set is assigned dynamically.

    */
    inline void UpdateLagMultSets(int offset_if) { UpdateLagMultSets(offset_if, false); };
    inline void UpdateLagMultSets(int offset_if, const bool& redistributed)
    {
      lmdofmap_ = UpdateLagMultSets(offset_if, redistributed, *sdofrowmap_);
    };
    Teuchos::RCP<Epetra_Map> UpdateLagMultSets(
        int offset_if, const bool& redistributed, const Epetra_Map& ref_map) const;

    /*! \brief Store the unredistributed local slave and lagrange multiplier maps
     *
     *  \author hiermeier
     *  \date 07/2016 */
    void StoreUnredistributedMaps();

    /*! \brief Redistribute the lagrange multiplier sets in a deterministic way
     *  in correlation to the redistributed slave dofs.
     *
     *  Please note, that this function becomes unnecessary as soon as we use
     *  dofsets for the Lagrange and displacement dofs!
     *
     *  \author hiermeier
     *  \date 07/2016 */
    Teuchos::RCP<Epetra_Map> RedistributeLagMultSets() const;

    /*!
    \brief Initialize / reset mortar interface

    */
    virtual void Initialize();

    /*!
    \brief Set current deformation state

    \param[in] Enum to encode type of state
    \param[in] Vector with state data
    */
    void SetState(const enum StateType& statetype, const Epetra_Vector& vec);

    /*!
    \brief Create integration cells for interface

    */
    void EvaluateGeometry(std::vector<Teuchos::RCP<MORTAR::IntCell>>& intcells);

    /*! @name Evaluate mortar interface
     *
     * This is the main routine of the MortarInterface class, where nodal normals
     * are computed, search is performed, mortar segments are set up and the entries
     * of the mortar matrices D and M are integrated. If the boolean flag "nonlinear"
     * is set to true (only for contact), then nonlinear mortar coupling is performed
     * including evaluation of the weighted gap g and linearizations of all mortar quantities. */
    /// @{

    /** \brief Evaluate mortar interface
     *
     *  No input necessary. */
    void Evaluate() { Evaluate(0, 0, 0, Teuchos::null); };

    /** \brief Evaluate mortar interface
     *
     *  \param rriter      (in)     : round robin iteration */
    void Evaluate(int rriter) { Evaluate(rriter, 0, 0, Teuchos::null); };

    /** \brief Evaluate mortar interface
     *
     *  \param rriter      (in)     : round robin iteration
     *  \param step        (in)     : current step counter */
    void Evaluate(int rriter, const int& step) { Evaluate(rriter, step, 0, Teuchos::null); };

    /** \brief Evaluate mortar interface
     *
     *  \param rriter      (in)     : round robin iteration
     *  \param step        (in)     : current step counter
     *  \param iter        (in)     : current iteration number */
    void Evaluate(int rriter, const int& step, const int& iter)
    {
      Evaluate(rriter, step, iter, Teuchos::null);
    };

    /** \brief Evaluate mortar interface
     *
     *  \param mparams_ptr (in/out) : mortar parameter interface pointer */
    void Evaluate(Teuchos::RCP<MORTAR::ParamsInterface> mparams_ptr)
    {
      Evaluate(0, mparams_ptr->GetStepNp(), mparams_ptr->GetNlnIter(), mparams_ptr);
    };

    /** \brief Evaluate mortar interface
     *
     *  \param rriter      (in)     : round robin iteration
     *  \param mparams_ptr (in/out) : mortar parameter interface pointer */
    void Evaluate(int rriter, Teuchos::RCP<MORTAR::ParamsInterface> mparams_ptr)
    {
      Evaluate(rriter, mparams_ptr->GetStepNp(), mparams_ptr->GetNlnIter(), mparams_ptr);
    };

    /** \brief Evaluate mortar interface
     *
     *  \param rriter      (in)     : round robin iteration
     *  \param step        (in)     : current step counter
     *  \param iter        (in)     : current iteration number
     *  \param mparams_ptr (in/out) : mortar parameter interface pointer */
    void Evaluate(int rriter, const int& step, const int& iter,
        Teuchos::RCP<MORTAR::ParamsInterface> mparams_ptr);

    /// @}

    /*!
    \brief Evaluate nodal normals

    */
    virtual void EvaluateNodalNormals() const;

    /*!
    \brief Evaluate nodal normals

    This methods computes the nodal normal and can be used to extract the nodal normals
    via the passed map. The nodal normals are stored in the map as vectors with the nodal gid
    as key. This function was only added because we want to use the NodalNormal computation
    of the mortar interface in UQ simulations.
    */
    virtual void EvaluateNodalNormals(std::map<int, std::vector<double>>& mynormals);

    /*!
    \brief Integrate Mortar matrices D and M and gap g on slave/master overlaps

    */
    virtual bool MortarCoupling(MORTAR::MortarElement* sele,
        std::vector<MORTAR::MortarElement*> mele,
        const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

    /*!
    \brief Assemble lagrange multipliers into global z vector (penalty strategy)

    */
    void AssembleLM(Epetra_Vector& zglobal);

    /*!
    \brief Assemble Mortar matrices D and M

    */
    void AssembleDM(LINALG::SparseMatrix& dglobal, LINALG::SparseMatrix& mglobal);

    /*!
    \brief Assemble Mortar matrix D

    */
    void AssembleD(LINALG::SparseMatrix& dglobal);

    /*!
    \brief Assemble Mortar matrix M

    */
    void AssembleM(LINALG::SparseMatrix& mglobal);

    /*!
    \brief Assemble matrix of normals N

    */
    void AssembleNormals(LINALG::SparseMatrix& nglobal);

    /*!
    \brief Assemble transformation matrices T and T^(-1)

    These matrices need to be applied to the slave displacements
    in the cases of dual LM interpolation for tet10/hex20 meshes
    in 3D. Here, the displacement basis functions have been modified
    in order to assure positivity of the D matrix entries and at
    the same time biorthogonality. Thus, to scale back the modified
    discrete displacements \hat{d} to the nodal discrete displacements
    {d}, we have to apply the transformation matrix T and vice
    versa with the transformation matrix T^(-1).
    */
    void AssembleTrafo(
        LINALG::SparseMatrix& trafo, LINALG::SparseMatrix& invtrafo, std::set<int>& donebefore);

    /*!
    \brief Detect actual meshtying zone (node by node)
    */
    void DetectTiedSlaveNodes(int& founduntied);

    /*!
    \brief Ghost underlying volume elements
    */
    void CreateVolumeGhosting();

    /*!
    \brief Binary tree search algorithm for potentially coupling
           slave / master pairs (element-based algorithm)
    */
    virtual bool EvaluateSearchBinarytree();

    /*!
    \brief Brute force search algorithm for potentially coupling
           slave / master pairs (element-based algorithm)
    */
    void EvaluateSearchBruteForce(const double& eps);

    /*!
    \brief find meles for one snode
    */
    void FindMEles(MortarNode& mrtrnode, std::vector<MortarElement*>& meles) const;

    /*!
    \brief return integration time for current interface
    */
    double Inttime() { return inttime_interface_; };

    /*!
    \brief return bool if interface is redistributed
    */
    bool& IsRedistributed() { return redistributed_; };

    //! @}

    //! @name Visualization and Debugging methods

    /*!
    \brief Visualize mortar stuff with gmsh (may be overloaded)
    */
    virtual void VisualizeGmsh(const int step, const int iter);

    /*!
    \brief Print shape function type (enum)
    */
    void PrintShapeFcn() { std::cout << shapefcn_ << std::endl; };

    void SetPoroFlag(bool poro) { poro_ = poro; }
    void SetEhlFlag(bool ehl) { ehl_ = ehl; }

    //@}

    bool HasMaSharingRefInterface() const;
    void AddMaSharingRefInterface(const MortarInterface* ref_interface);

    //! @name Output
    //! @{

    /*!
    \brief Output results for postprocessing/visualization of this interface

    Use this interface's discretization write to write output of this interface.

    \note This is purely for visualization. Writing/reading restart data is still dealt with by the
    structure discretization.

    \param[in] ouputParams Parameter list with stuff required by interfaces to write output
    */
    virtual void PostprocessQuantities(const Teuchos::ParameterList& outputParams);

    //! @}

   protected:
    //! @name Output
    //! @{

    /*!
    \brief Check validity of parameter list for restart/output writing

    Writing restart data and results is controlled via parameters in the
    \c outParams parameter list. We check whether required entries exist and
    maybe throw an error.

    To add a entry to the list of required entries, just add a new line with
    \code
    requiredEntries.push_back("<new entry name>");
    \endcode

    \param[in] outParams Parameter list with output configuration and auxiliary output data
    \param[in] requiredEntries List of required parameter list entries

    \sa OutputState()
     */
    bool CheckOutputList(const Teuchos::ParameterList& outParams,
        const std::vector<std::string>& requiredEntries) const;

    //! @}

    const MortarInterface* GetMaSharingRefInterfacePtr() const;

    /*!
    \brief Redistribute the master side of the interface

    \param[out] rownodes Redistributed node row map of master side after
    \param[out] colnodes Redistributed node column map of master side
    \param[in] roweles Element row map of master side
    \param[in] comm Communicator object
    \param[in] parts Number of desired subdomains after redistribution
    \param[in] imbalance Max. relative imbalance of subdomain size
    */
    void RedistributeMasterSide(Teuchos::RCP<Epetra_Map>& rownodes,
        Teuchos::RCP<Epetra_Map>& colnodes, const Teuchos::RCP<Epetra_Map>& roweles,
        const Teuchos::RCP<Epetra_Comm>& comm, const int parts, const double imbalance) const;

    /** \brief function called at the beginning of each MortarCoupling call
     *
     *  \param sele        (in): pointer to the current slave element
     *  \param mele        (in): pointer to the current master element
     *  \param mparams_ptr (in): mortar parameter interface pointer
     *
     *  \author hiermeier
     *  \date 11/16 */
    virtual void PreMortarCoupling(const MORTAR::MortarElement* sele,
        const std::vector<MORTAR::MortarElement*> mele,
        const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr) const {
        /* does nothing in the default case */};

    /** \brief function called at the end of each MortarCoupling call
     *
     *  \param sele        (in): pointer to the current slave element
     *  \param mele        (in): pointer to the current master element
     *  \param mparams_ptr (in): mortar parameter interface pointer
     *
     *  \author hiermeier
     *  \date 11/16 */
    virtual void PostMortarCoupling(const MORTAR::MortarElement* sele,
        const std::vector<MORTAR::MortarElement*> mele,
        const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr) const {
        /* does nothing in the default case */};

    /*!
    \brief Evaluate segment-to-segment coupling (mortar...)

    */
    virtual void EvaluateSTS(
        const Epetra_Map& selecolmap, const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

    /*!
    \brief Evaluate node-to-segment coupling

    */
    virtual void EvaluateNTS();

    /*!
    \brief Evaluate line-to-segment coupling

    */
    virtual void EvaluateLTS();

    /*!
    \brief Evaluate line-to-segment coupling

    */
    virtual void EvaluateLTL();

    /*!
    \brief Evaluate segment-to-line coupling

    */
    virtual void EvaluateSTL();

    /*!
    \brief Export nodal normals

    This method exports / communicates the nodal normal vector from row to
    column map layout (needed for coupling evaluation). The reason behind this
    is that only the row nodes can correctly evaluate the nodal normal based on
    averaging of the adjacent element normals. A column node might not know of
    all adjacent elements and thus would compute wrong nodal normals by itself.

    */
    virtual void ExportNodalNormals() const;

    /*!
    \brief build nodal normals and other stuff before real coupling
           starts

    */
    virtual void PreEvaluate(const int& step, const int& iter);

    /*!
    \brief Evaluate mortar interface

    Here, we decide if we evaluate STS (mortar), NTS, GPTS... contact.
    Thus, it can be interpreted as switch between constraints discretizations

    */
    virtual void EvaluateCoupling(const Epetra_Map& selecolmap, const Epetra_Map* snoderowmap,
        const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

    /*!
    \brief do scaling and other operations after real coupling

    */
    virtual void PostEvaluate(const int step = 0, const int iter = 0);

    /*!
    \brief find master nodes for one snode

    */
    void FindMNodes(
        MortarNode& mrtrnode, std::vector<MortarElement*>& meles, std::vector<MortarNode*>& mnodes);

    /*!
    \brief Initialize data container for nodes and elements

    */
    virtual void InitializeDataContainer();

    /*!
    \brief check and init possible corner/edge modification

    */
    virtual void InitializeCornerEdge();

    /*!
    \brief Set element areas

    */
    virtual void SetElementAreas();

    /*!
    \brief Split MortarElements into IntElements for 3D quadratic coupling

    */
    bool SplitIntElements(
        MORTAR::MortarElement& ele, std::vector<Teuchos::RCP<MORTAR::IntElement>>& auxele);

    /*!
    \brief Update interface master and slave sets

    This update is usually only done ONCE in the initialization phase
    and sets up the slave and master sets (elements, nodes, dofs) for
    the whole simulation. Yet, in the case of self contact the sets
    need to be updated again and again during simulation time, as the
    slave/master status is assigned dynamically.

    */
    virtual void UpdateMasterSlaveSets();

   private:
    /*!
    \brief Extend interface ghosting

    To guarantee that the search algorithms can find _every_ possible slave/master pair
    -- even if entities reside on different processors --
    we afford the luxury to ghost some or all nodes
    on all processors in the general mortar coupling framework.

    Technically, we compute extended node/element column maps and the export the discretization
    based on these extended column maps. The actual extension strategy is defined in the input file.

    The tasks to be performed depend on the status of the parallel distriution: Some tasks are only
    necessary, if the parallel distribution is final and won't change anymore, before using the
    interface discretization.

    \note We'll do ghosting _NOT ONLY_ on procs that do own or ghost any of the nodes in the natural
    distribution of the interface discretization #idiscret_, but really on _ALL_ procs. This makes
    dynamic redistribution easier.

    \note In some cases (self contact, sliding ALE mortar coupling), we still
    need the SLAVE nodes and elements in fully overlapping column layout,
    too. In the case of self contact, this is due to the fact that contact
    search is then based on the contact interface as a whole without
    initially distinguishing between slave and master sides. In general,
    however we do not need (or even want) slave redundancy. Redundancy of the
    master side is controlled via the input file.

    \post The interface discretization is _NOT_ FillComplete.

    \param[in] isFinalParallelDistribution Is this parallel distribution final?
    \param[in] meanVelocity Mean velocity of this interface (needed for contact with binning)
    */
    void ExtendInterfaceGhosting(const bool isFinalParallelDistribution, const double meanVelocity);

    /*!
    \brief Update master and slave interface maps for DOFs

    Update #sdofrowmap_, #sdofcolmap_, #mdofrowmap_, #mdofcolmap_.

    \pre A valid node column map of interface discretization exists.
    \post DOF row/column maps for slave/master side exist.
    */
    void UpdateMasterSlaveDofMaps();

    /*!
    \brief Update master and slave interface maps for elements based on current maps in the
    interface discretization

    \pre A valid node column map of interface discretization exists.

    \sa UpdateMasterSlaveElementMaps(const Epetra_Map&, const Epetra_Map&)
    */
    void UpdateMasterSlaveElementMaps();

    /*!
    \brief Update master and slave interface maps for elements using given maps

    \post Element row/column maps for slave/master side -- i.e. #selerowmap_, #selecolmap_,
    #melerowmap_, #melecolmap_ -- exist and have been updated.

    \param[in] elementRowMap Row map of interface elements
    \param[in] elementColumnMap Column map of interface elements

    \warning Use the row map only for checking ownership of an element. Do not loop the row map,
    because we sometimes pass a column map to the row map argument.
    */
    void UpdateMasterSlaveElementMaps(
        const Epetra_Map& elementRowMap, const Epetra_Map& elementColumnMap);

    /*!
    \brief Update master and slave interface maps for nodes based on current maps in the interface
    discretization

    \pre A valid node column map of interface discretization exists.

    \sa UpdateMasterSlaveNodeMaps(const Epetra_Map&, const Epetra_Map&)
    */
    void UpdateMasterSlaveNodeMaps();

    /*!
    \brief Update master and slave interface maps for nodes

    \post Node row/column maps for slave/master side -- i.e. #snoderowmap_, #snodecolmap_,
    #mnoderowmap_, #mnodecolmap_ -- exist and have been updated.

    \param[in] nodeRowMap Row map of interface nodes
    \param[in] nodeColumnMap Column map of interface nodes

    \warning Use the row map only for checking ownership of a node. Do not loop the row map, because
    we sometimes pass a column map to the row map argument.
    */
    void UpdateMasterSlaveNodeMaps(const Epetra_Map& nodeRowMap, const Epetra_Map& nodeColumnMap);

    /*!
    \brief Setup the binning strategy for geometrically motivated extended ghosting

    \note The average velocity is only important, if binning is performed in the deformed
    configuration, e.g. as in contact problems. For meshtying, where binning occurs during mortar
    coupling in the reference configuration, it suffices to set the velocity to 0.0.

    @param[in] meanVelocity Current absolute value of the mean velocity of this interface
    @return Binning strategy object ready to be used
    */
    Teuchos::RCP<BINSTRATEGY::BinningStrategy> SetupBinningStrategy(const double meanVelocity);

    /*!
    \brief Create interface local communicator

    Find all procs that have business on this interface (own or ghost nodes/elements). Then, build a
    Epetra_Comm that contains only those procs. This intra-communicator will be used to handle most
    stuff on this interface so the interface will not block all other procs

    ToDo (mayr.mt) Move this to MORTAR::UTILS::INTERFACE. This does not have to be a function in
    this class.
    */
    void CreateInterfaceLocalCommunicator();

    /*!
    \brief check and initialize possible cross point modification

    */
    void InitializeCrossPoints();

    /*!
    \brief check and initialize for linear Lagrange multiplier interpolation for 2nd order FE

    */
    void InitializeLagMultLin();

    /*!
    \brief check and initialize for constant Lagrange muliplier interpolation for 2nd order FE

    */
    void InitializeLagMultConst();

    /*!
    \brief Setup interface discretization

    Creates either a regular DRT::Discretization or a DRT::NurbsDiscretization.
    */
    void CreateInterfaceDiscretization();

    /*!
    \brief Set type of shape functions

    Read type of shape function form input file parameter list and store it in #shapefcn_.
    */
    void SetShapeFunctionType();

    /// pointer to the interface data object
    Teuchos::RCP<InterfaceDataContainer> interfaceData_;

   protected:
    //! don't want = operator
    MortarInterface operator=(const MortarInterface& old) = delete;

    //! don't want copy constructor
    MortarInterface(const MortarInterface& old) = delete;

    /** @name References to the interface data container content
     *
     * \remark Please add no new member variables to this class and use the
     *  corresponding data container, instead! If you have any questions
     *  concerning this, do not hesitate and ask me.
     *                                                          hiermeier 03/17 */
    /// @{

    int& id_;                                ///< ref. to unique interface id
    Teuchos::RCP<const Epetra_Comm>& comm_;  ///< ref. to communicator
    Teuchos::RCP<Epetra_Comm>& lcomm_;       ///< ref. to interface local communicator
    std::map<int, int>& procmap_;            ///< ref. to mapping global -> local communicator PIDs
    bool& redistributed_;                    ///< ref. to redistribution for this time step?

    Teuchos::RCP<DRT::Discretization>&
        idiscret_;                     ///< ref. to the discretization of the mortar interface
    int& dim_;                         ///< ref. to dimension of problem (2D or 3D)
    Teuchos::ParameterList& imortar_;  ///< ref. to containing contact input parameters of interface
    INPAR::MORTAR::ShapeFcn& shapefcn_;  ///< ref. to employed type of shape function set
    bool& quadslave_;                    ///< ref. to flag indicating quadratic 2d/3d slave elements
    INPAR::MORTAR::RedundantStorage&
        redundant_;  ///< ref. to employed type of redundancy in storage of interface

    Teuchos::RCP<Epetra_Map>&
        oldnodecolmap_;  ///< ref. to column map of all interface nodes (overlap=1)
    Teuchos::RCP<Epetra_Map>&
        oldelecolmap_;  ///< ref. to column map of all interface elements (overlap=1)

    Teuchos::RCP<Epetra_Map>& snoderowmap_;  ///< ref. to row map of all slave nodes
    Teuchos::RCP<Epetra_Map>& snodecolmap_;  ///< ref. to column map of all slave nodes
    Teuchos::RCP<Epetra_Map>& mnoderowmap_;  ///< ref. to row map of all master nodes
    Teuchos::RCP<Epetra_Map>& mnodecolmap_;  ///< ref. to column map of all master nodes

    Teuchos::RCP<Epetra_Map>&
        snoderowmapbound_;  ///< ref. to row map of slave nodes (+ boundary nodes)
    Teuchos::RCP<Epetra_Map>&
        snodecolmapbound_;  ///< ref. to col map of slave nodes (+ boundary nodes)
    Teuchos::RCP<Epetra_Map>&
        mnoderowmapnobound_;  ///< ref. to row map of master nodes (- boundary nodes)
    Teuchos::RCP<Epetra_Map>&
        mnodecolmapnobound_;  ///< ref. to col map of master nodes (- boundary nodes)

    Teuchos::RCP<Epetra_Map>& selerowmap_;  ///< ref. to row map of all slave elements
    Teuchos::RCP<Epetra_Map>& selecolmap_;  ///< ref. to column map of all slave elements
    Teuchos::RCP<Epetra_Map>& melerowmap_;  ///< ref. to row map of all master elements
    Teuchos::RCP<Epetra_Map>& melecolmap_;  ///< ref. to column map of all master elements

    Teuchos::RCP<Epetra_Map>& sdofrowmap_;  ///< ref. to row map of all slave dofs
    Teuchos::RCP<Epetra_Map>& sdofcolmap_;  ///< ref. to column map of all slave dofs
    Teuchos::RCP<Epetra_Map>& mdofrowmap_;  ///< ref. to row map of all master dofs
    Teuchos::RCP<Epetra_Map>& mdofcolmap_;  ///< ref. to column map of all master dofs

    Teuchos::RCP<Epetra_Map>&
        psdofrowmap_;  ///< ref. to row map of all slave dofs before any redistribution took place
    Teuchos::RCP<Epetra_Map>&
        plmdofmap_;  ///< ref. to row map of all lm dofs before any redistribution took place

    Teuchos::RCP<Epetra_Map>& lmdofmap_;  ///< ref. to row map of all Lagrange multiplier dofs
    int& maxdofglobal_;                   ///< ref. to maximum dof ID in global discretization

    INPAR::MORTAR::SearchAlgorithm& searchalgo_;    ///< ref. to type of search algorithm
    Teuchos::RCP<MORTAR::BinaryTree>& binarytree_;  ///< ref. to binary searchtree
    double& searchparam_;                           ///< ref. to search parameter
    bool& searchuseauxpos_;      ///< ref. to use auxiliary position when computing dops
    double& inttime_interface_;  ///< ref. to integration time

    bool& nurbs_;  ///< ref. to flag for nurbs shape functions
    bool& poro_;   ///< ref. to flag if poro contact problem!
    bool& ehl_;    ///< ref. to flag if ehl contact problem!

    /// @}
  };  // class MortarInterface
}  // namespace MORTAR


// << operator
std::ostream& operator<<(std::ostream& os, const MORTAR::MortarInterface& interface);


#endif  // #ifndef MORTAR_INTERFACE_H
