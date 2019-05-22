/*----------------------------------------------------------------------------*/
/*!

\brief Find connectivities in a discretization.

\maintainer Martin Kronbichler

\level 3
*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_LIB_CONNECTIVITY_H_
#define SRC_DRT_LIB_CONNECTIVITY_H_

#include <unordered_map>
#include <Teuchos_RCP.hpp>


class Epetra_Map;
class Epetra_BlockMap;
class Epetra_Comm;
namespace LINALG
{
  class SerialDenseVector;
}
namespace DRT
{
  class DiscretizationInterface;
  class Node;

  /*--------------------------------------------------------------------------*/
  /** @brief Find connectivities
   *
   *  This class enables e.g. the splitting of a discretization object into
   *  distinct non-interconnected subsets. Feel free to add additional related
   *  functionality.
   *
   *  @author hiermeier @date 03/18 */
  class Connectivity
  {
    typedef std::set<const DRT::Node*> node_set;
    typedef std::vector<std::set<const DRT::Node*>> node_sets;

    /**
     * 1. index = LID of the body on the sending proc
     * 2. index = LID/GID of the sending proc
     * 3. index = GIDs of the nodal ghosted values on the sending proc */
    typedef std::vector<std::unordered_map<int, std::vector<int>>> received_gid_set;

    // forward declaration
    struct Connection;

   public:
    enum class MapFormat
    {
      row,
      column
    };

    /// constructors needed
    Connectivity() = delete;
    Connectivity(const Connectivity& c) = delete;
    ~Connectivity() = delete;

    /** @brief Split a discretization into a set of non-interconnected
     *  node maps
     *
     *  Each of the node maps represents one body. The only current
     *  prerequisite is that a layer of ghosted nodes is used in parallel.
     *
     *  @param(in) full_discret  Filled full discretization object
     *  @param(in) format        Choose the desired format of the body node maps
     *  @param(out) body_node_maps  Vector containing the distinct body node maps
     *
     *  \author hiermeier \date 03/18 */
    static void splitDiscretizationIntoDistinctNodeMaps(
        const DRT::DiscretizationInterface& full_discret, const MapFormat format,
        std::vector<Teuchos::RCP<Epetra_Map>>& body_node_maps);

   private:
    static void fillMyConnectedNodeSet(
        node_set& connected_nset, const DRT::DiscretizationInterface& full_discret);

    static void receiveAnyOfMyGids(received_gid_set& all_received_ngids,
        const std::vector<node_set*>& my_new_body_nsets, const Epetra_Comm& comm);

    static Teuchos::RCP<Epetra_Map> createNodeMapFromBody(
        const Epetra_Comm& comm, const node_set& my_node_set, const MapFormat format);

    static void createBodyNodeMaps(const MapFormat format, const Epetra_Comm& comm,
        const Connection& connections, const node_sets& my_body_node_sets,
        std::vector<Teuchos::RCP<Epetra_Map>>& body_node_maps);

   private:
    /** \brief Handle the connection among different bodies over multiple procs
     *
     *  \author hiermeier \date 03/18 */
    struct Connection
    {
      Connection(int my_pid) : my_pid_(my_pid){};

      void add(int my_id, int other_id, int other_pid)
      {
        my_body_ids_.push_back(my_id);
        other_body_ids_.push_back(other_id);
        other_pids_.push_back(other_pid);
      }

      void print(std::ostream& os) const
      {
        for (unsigned i = 0; i < my_body_ids_.size(); ++i)
          os << "My body #" << my_body_ids_[i] << " on proc #" << my_pid_
             << " is connected to body #" << other_body_ids_[i] << " on proc #" << other_pids_[i]
             << std::endl;
      }

      void getOthers(int bid, std::vector<std::pair<int, int>>& others,
          const std::set<int>* ignore = NULL) const
      {
        others.clear();
        for (unsigned i = 0; i < my_body_ids_.size(); ++i)
        {
          if (my_body_ids_[i] == bid)
          {
            if (ignore and ignore->find(other_pids_[i]) != ignore->end()) continue;

            others.push_back(std::make_pair(other_body_ids_[i], other_pids_[i]));
          }
        }
      }

      int my_pid_ = -1;
      std::vector<int> my_body_ids_;
      std::vector<int> other_body_ids_;
      std::vector<int> other_pids_;
    };
  };  // class Connectivity
}  // namespace DRT



#endif /* SRC_DRT_LIB_CONNECTIVITY_H_ */
