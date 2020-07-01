/*----------------------------------------------------------------------------*/
/** \file

\brief Implementation of a multi discretization wrapper for eXtended structural
       problems


\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_STRUCTURE_XSTRUCTURE_XSTR_MULTI_DISCRETIZATION_WRAPPER_H_
#define SRC_DRT_STRUCTURE_XSTRUCTURE_XSTR_MULTI_DISCRETIZATION_WRAPPER_H_

#include "../drt_lib/drt_discret_interface.H"
// necessary due to enum definition
#include "../drt_xfem/xfem_enum_lists.H"
#include "../drt_io/io.H"

namespace IO
{
  class DiscretizationWriter;
}  // namespace IO
namespace XFEM
{
  class MultiFieldMapExtractor;
}  // namespace XFEM
namespace XSTR
{
  namespace IO
  {
    class MultiDiscretizationWriter;
  }  // namespace IO

  /** \brief Implementation of a multi discretization wrapper for eXtended structural
   *  problems
   *
   *  This class wraps a couple of discretizations with shared joint interfaces and
   *  makes it possible to access the row information which is necessary to handle
   *  system vectors and system matrices without any further effort.
   *
   *  Nevertheless, it is not possible to access any column information. Here you
   *  are forced to access each single wrapped discretization separately, since
   *  the merged information is not unique.
   *
   *  \author hiermeier \date 10/16 */
  class MultiDiscretizationWrapper : public DRT::DiscretizationInterface
  {
   public:
    typedef std::map<enum XFEM::FieldName, Teuchos::RCP<DRT::DiscretizationInterface>> XDisMap;
    typedef std::pair<enum XFEM::FieldName, Teuchos::RCP<DRT::DiscretizationInterface>> XDisPair;
    typedef std::pair<enum XFEM::FieldName, Teuchos::RCP<const DRT::DiscretizationInterface>>
        cXDisPair;

   private:
    typedef std::vector<Teuchos::RCP<DRT::DiscretizationInterface>> XDisVec;
    typedef std::vector<Teuchos::RCP<const DRT::DiscretizationInterface>> cXDisVec;

    typedef std::map<enum XFEM::FieldName, std::vector<Teuchos::RCP<LINALG::SparseOperator>>>
        XSysMatMap;
    typedef std::pair<enum XFEM::FieldName, std::vector<Teuchos::RCP<LINALG::SparseOperator>>>
        XSysMatPair;

   public:
    /// Overriding base class overloaded functions
    using DiscretizationInterface::Dof;
    using DiscretizationInterface::EvaluateCondition;
    using DiscretizationInterface::EvaluateNeumann;
    using DiscretizationInterface::GetState;
    using DiscretizationInterface::NumDof;
    using DiscretizationInterface::SetState;

   public:
    /// constructor
    MultiDiscretizationWrapper();

    /// destructor
    virtual ~MultiDiscretizationWrapper(){};

    /** \brief Initialize class member variables
     *
     *  \param name                           (in): name of the multi discretization wrapper
     *  \param discret                        (in): discretization to wrap
     *  \param max_num_reserved_dofs_per_node (in): Maximum number of dofs per enriched node
     *
     *  \author hiermeier \date 10/16 */
    void Init(const std::string& name,
        const Teuchos::RCP<DRT::DiscretizationInterface>& discret_ptr,
        const int& max_num_reserved_dofs_per_node);

    /** \brief Initialize class member variables
     *
     *  \param name                           (in): name of the multi discretization wrapper
     *  \param discret_vec                    (in): vector of discretizations to wrap
     *  \param max_num_reserved_dofs_per_node (in): Maximum number of dofs per enriched node
     *
     *  \author hiermeier \date 10/16 */
    void Init(const std::string& name, const XDisVec& discret_vec,
        const int& max_num_reserved_dofs_per_node);

    /** \brief Setup class new class member variables
     *
     *  \author hiermeier \date 10/16 */
    virtual void Setup();

    /** \brief Add a new discretization
     *
     *  Add a new discretization to the wrapped discretization map. If a new
     *  discretization is added, the isfilled_ flag is set to FALSE!
     *
     *  \param discret_ptr     (in): pointer to the discretization to add.
     *  \param check_insertion (in): Check, if the discretization was already added
     *                               and throw an error if a discretization
     *                               with the same key has been found (default: TRUE).
     *                               If the check is switched off, we replace the
     *                               existing pointer by a new pointer.
     *
     *  \author hiermeier \date 10/16 */
    void AddDiscretization(const Teuchos::RCP<DRT::DiscretizationInterface>& discret_ptr)
    {
      AddDiscretization(discret_ptr, true);
    }
    void AddDiscretization(
        const Teuchos::RCP<DRT::DiscretizationInterface>& discret_ptr, const bool& check_insertion);

    void ExtractVector(const Epetra_MultiVector& full, enum XFEM::FieldName field,
        Epetra_MultiVector& partial, enum ::IO::VectorType vt = ::IO::dofvector) const;

    void InsertVector(const Epetra_MultiVector& partial, enum XFEM::FieldName field,
        Epetra_MultiVector& full, enum ::IO::VectorType vt = ::IO::dofvector) const;

    void AddVector(const Epetra_MultiVector& partial, enum XFEM::FieldName field,
        Epetra_MultiVector& full, double scale, enum ::IO::VectorType vt = ::IO::dofvector) const;

    /** \brief Access all the wrapped discretizations as std::map
     *
     *  \author hiermeier \date 10/16 */
    const XDisMap& DiscretMap() const
    {
      CheckInit();
      return discret_map_;
    }

    /** \brief Access one of the wrapped discretizations (as reference)
     *
     *  \author hiermeier \date 10/16 */
    const DRT::DiscretizationInterface& Discret(enum XFEM::FieldName field) const
    {
      CheckInit();
      return *(GetDiscretPtr(field));
    }
    DRT::DiscretizationInterface& Discret(enum XFEM::FieldName field)
    {
      CheckInit();
      return *(GetDiscretPtr(field));
    }

    /** \brief Access one of the wrapped discretizations (as rcp)
     *
     *  \author hiermeier \date 10/16 */
    Teuchos::RCP<const DRT::DiscretizationInterface> DiscretPtr(const XFEM::FieldName& field) const
    {
      CheckInit();
      return GetDiscretPtr(field).getConst();
    }
    const Teuchos::RCP<DRT::DiscretizationInterface>& DiscretPtr(const XFEM::FieldName& field)
    {
      CheckInit();
      return GetDiscretPtr(field);
    }

    /** \brief get the maximum number of reserved DoF's per node (set by Init())
     *
     *  \author hiermeier \date 10/16 */
    const int& MaxNumReservedDofsPerNode() const
    {
      CheckInit();
      return max_num_reserved_dofs_per_node_;
    }

    /** \brief Access the communicator object [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual const Epetra_Comm& Comm() const
    {
      CheckInit();
      return discret_map_.begin()->second->Comm();
    }

    /** \brief Get communicator associated with this class
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual Teuchos::RCP<::IO::DiscretizationWriter> Writer() const;

    /** \brief Get flag indicating whether FillComplete() has been called [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual inline bool Filled() const
    {
      bool is_filled = true;
      XDisMap::const_iterator cit;
      for (cit = discret_map_.begin(); cit != discret_map_.end(); ++cit)
      {
        if (not cit->second->Filled()) is_filled = false;
      }
      return (is_filled and isfilled_);
    }

    /** \brief Get name of the discretization [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual inline const std::string& Name() const
    {
      CheckInit();
      return name_;
    }

    /** \brief Print this discretization to os (Filled()==true NOT prerequisite)
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual void Print(std::ostream& os) const;

    /** \brief Complete construction of a discretization
     *         (Filled()==true NOT prerequisite) [derived]
     *
     *  This method calls the class specific version, setting the
     *  new additional parameter \c setupmapextractor to TRUE.
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual int FillComplete(bool assigndegreesoffreedom = true, bool initelements = true,
        bool doboundaryconditions = true)
    {
      return FillComplete(assigndegreesoffreedom, initelements, doboundaryconditions, true, true);
    }
    virtual int FillComplete(bool assigndegreesoffreedom = true, bool initelements = true,
        bool doboundaryconditions = true, bool buildsystemmaps = true,
        bool setupmapextractor = true);

    /** \brief Get map associated with the distribution of the ownership of nodes [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual const Epetra_Map* NodeRowMap() const;

    /** \brief Get map associated with the distribution of the ownership of elements [derive]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual const Epetra_Map* ElementRowMap() const;

    /** \brief Get processor local number of nodes owned by this processor [derived]
     *         (Filled()==true prerequisite)
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual int NumMyRowNodes() const;

    /** \brief Query whether a Node with global id gid is stored on this proc
     *         (Filled()==true NOT prerequisite) [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual bool HaveGlobalNode(int gid) const;

    /** \brief Get the node with global row id gid (Filled()==true NOT prerequisite) [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual DRT::Node* gNode(int gid) const;

    /** \brief Get the node with local row id lid (Filled()==true prerequisite) [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual DRT::Node* lRowNode(int lid) const;

    /** \brief Get processor local number of elements owned by this processor [derived]
     *         (Filled()==true prerequisite)
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual int NumMyRowElements() const;

    /** \brief Get the element with local row id lid (Filled()==true prerequisite) [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual DRT::Element* lRowElement(int lid) const;

    /** \brief Get the element with global id gid (Filled()==true NOT prerequisite) [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual DRT::Element* gElement(int gid) const;

    /** \brief Get degree of freedom row map (Filled()==true prerequisite) [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual const Epetra_Map* DofRowMap(unsigned nds = 0) const;

    /** \brief True if degrees of freedom where assigned [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual bool HaveDofs() const;

    /** \brief Get the number of DoF sets [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual int NumDofSets() const;

    /** \brief Get number of dofs for given node. [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual int NumDof(unsigned dofset_id, const DRT::Node* node) const;

    /** \brief Get number of standard (w/o enrichment) DoF's for given node.
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual int NumStandardDof(const unsigned& dofset_id, const DRT::Node* node) const;

    /** \brief Get the gid of all dofs of a node. [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual std::vector<int> Dof(unsigned dofset_id, const DRT::Node* node) const;

    /** \brief Get the gid of a dof for given node. [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual int Dof(unsigned dofset_id, const DRT::Node* node, const int ldof) const;

    /** \brief Get the gid of all dofs of a node. [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual void Dof(unsigned dofset_id, const DRT::Node* node, std::vector<int>& lm) const;

    /*! \brief Get the gid of all dofs of a node.
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual void Dof(std::vector<int>& dof, const DRT::Node* node, unsigned dofset_id,
        unsigned nodaldofset_id, const DRT::Element* element = NULL) const;

    /** \brief Get all conditions with a certain name [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual void GetCondition(const std::string& name, std::vector<DRT::Condition*>& out) const;
    virtual void GetCondition(
        const std::string& name, std::vector<Teuchos::RCP<DRT::Condition>>& out) const;

    /** \brief Get a condition with a certain name [derived]
     *
     *  \remark The implemented variant returns  the first meaningful candidate!
     *  Please note: Maybe it makes more sense to merge the first condition of
     *  each discretization to a new condition and to return this merged
     *  condition.
     *  Anyway, if you think, that the condition can be found in more than one
     *  discretization, it is always a better idea to use the GetCondition
     *  routines, which will return a std::vector of all conditions with the
     *  given name. In that way you can be sure, that you do not miss any.
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual DRT::Condition* GetCondition(const std::string& name) const;

    /** \brief Set a reference to a data vector [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual void SetState(
        unsigned nds, const std::string& name, Teuchos::RCP<const Epetra_Vector> state);

    /** \brief Get a reference to a data vector [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual Teuchos::RCP<const Epetra_Vector> GetState(unsigned nds, const std::string& name) const;

    /** \brief Clear solution state references [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual void ClearState(bool clearalldofsets = false);

    /** \brief Evaluate Neumann boundary conditions [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual void EvaluateNeumann(Teuchos::ParameterList& params, Epetra_Vector& full_systemvector,
        LINALG::SparseOperator* full_systemmatrix = NULL);

    /** \brief Call elements to evaluate [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual void Evaluate(Teuchos::ParameterList& params,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
        Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
        Teuchos::RCP<Epetra_Vector> systemvector3);

    /** \brief Evaluate Dirichlet boundary conditions [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual void EvaluateDirichlet(Teuchos::ParameterList& params,
        Teuchos::RCP<Epetra_Vector> systemvector, Teuchos::RCP<Epetra_Vector> systemvectord,
        Teuchos::RCP<Epetra_Vector> systemvectordd, Teuchos::RCP<Epetra_Vector> toggle,
        Teuchos::RCP<LINALG::MapExtractor> dbcmapextractor = Teuchos::null);

    /** \brief Evaluate a specific condition [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual void EvaluateCondition(Teuchos::ParameterList& params,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
        Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
        Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
        Teuchos::RCP<Epetra_Vector> systemvector3, const std::string& condstring,
        const int condid = -1);

    /** \brief Assemble scalar quantities across elements [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual void EvaluateScalars(
        Teuchos::ParameterList& params, Teuchos::RCP<Epetra_SerialDenseVector> scalars);

    /** \brief Compute the nullspace of the discretization [derived]
     *
     *  \note See the base class for more information.
     *
     *  \author hiermeier \date 10/16 */
    virtual void ComputeNullSpaceIfNecessary(
        Teuchos::ParameterList& solveparams, bool recompute = false);

   protected:
    inline bool IsInit() const { return isinit_; };

    inline bool IsSetup() const { return issetup_; };

    inline void CheckInitSetup() const
    {
      if (not IsInit() or not IsSetup())
        dserror("Call XSTR::MultiDiscretizationWrapper::Init() and Setup() first!");
    };

    inline void CheckInit() const
    {
      if (not IsInit())
        dserror("XSTR::MultiDiscretizationWrapper::Init() has not been called, yet!");
    };

    bool IsXFEM(enum XFEM::FieldName field) const;

    void CreateNewMap(const std::vector<int>& my_entries, Teuchos::RCP<Epetra_Map>& new_map) const;

   private:
    void Reset();

    void BuildMergedNodes();

    void BuildMergedNodeRowMap();

    void BuildMergedNodeRowPtr();

    void BuildMergedElements();

    void BuildMergedElementRowMap();

    void BuildMergedElementRowPtr();

    const Teuchos::RCP<LINALG::SparseOperator>& XFieldSysMat(
        enum XFEM::FieldName field, unsigned mat_id);

    inline const Teuchos::RCP<DRT::DiscretizationInterface>& GetDiscretPtr(
        enum XFEM::FieldName field) const
    {
      XDisMap::const_iterator cit = discret_map_.find(field);
      if (cit == discret_map_.end())
        dserror(
            "The given field name \"%s\" could not be found in the wrapped "
            "discretization map!",
            XFEM::FieldName2String(field).c_str());

      return cit->second;
    }

   protected:
    bool isinit_;

    bool issetup_;

    bool isfilled_;

    std::string name_;

    int max_num_reserved_dofs_per_node_;

    XDisMap discret_map_;

    XSysMatMap partial_system_matrices_;

    Teuchos::RCP<XFEM::MultiFieldMapExtractor> map_extractor_;

    Teuchos::RCP<XSTR::IO::MultiDiscretizationWriter> writers_wrapper_;

    std::map<int, Teuchos::RCP<DRT::Node>> merged_node_;

    Teuchos::RCP<const Epetra_Map> merged_node_row_map_;

    std::vector<DRT::Node*> merged_node_row_ptr_;

    std::map<int, Teuchos::RCP<DRT::Element>> merged_element_;

    Teuchos::RCP<Epetra_Map> merged_ele_row_map_;

    std::vector<DRT::Element*> merged_element_row_ptr_;
  };  // class MultiDiscretizationWrapper
}  // namespace XSTR



#endif /* SRC_DRT_STRUCTURE_XSTRUCTURE_XSTR_MULTI_DISCRETIZATION_WRAPPER_H_ */
