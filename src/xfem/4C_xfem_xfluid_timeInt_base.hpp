/*----------------------------------------------------------------------*/
/*! \file

\brief provides the basic time integration classes "TimeInt", "TimeIntStd", "TimeIntEnr"

\level 2


*/
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_XFEM_XFLUID_TIMEINT_BASE_HPP
#define FOUR_C_XFEM_XFLUID_TIMEINT_BASE_HPP

#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_inpar_xfem.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Map.h>

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  class CutWizard;

  namespace Cut
  {
    class SideHandle;
    class VolumeCell;
  }  // namespace Cut
}  // namespace Core::Geo

namespace XFEM
{
  class XFEMDofSet;

  /*!
  \brief this class is the basic TIMEINT class for the projection, adaption or
         something else in XFEM-problems between consecutive time steps
   */
  class XfluidTimeintBase
  {
   public:
    /*!
    \brief This subclass holds data for a node that will get some new values.
           For the moment, only Semi-Lagrangean algorithms use the data.
     */
    class TimeIntData
    {
     public:
      /*========================================================================*/
      //! @name enumerations used for nodes for which values have to be constructed
      /*========================================================================*/

      //! status for both the used special algorithm and the state within the algorithm
      enum State
      {
        basicStd_,
        currSL_,
        nextSL_,
        failedSL_,
        doneStd_,
      };

      //! basic computation type due to FGI and FRS
      enum Type
      {
        predictor_ = 0,
        standard_ = 1
      };

      //! does to closest projection of a node lie on a side, a line or on a point?
      enum Projection
      {
        onSide_,
        onLine_,
        onPoint_,
        failed_
      };

      /*========================================================================*/
      //! @name constructor's for different time-integration algorithms, destructor
      /*========================================================================*/

      //! constructor for basic data for standard computation in Semi-Lagrangean algorithm
      TimeIntData(Core::Nodes::Node& node,  //! node for which SL-algorithm is called
          int nds_np,  //! nds (nodal dofset) number w.r.t new interface position, for which SL-algo
                       //! is called
          Core::LinAlg::Matrix<3, 1> vel,  //! velocity at point x (=x_Lagr(t^n+1))
          std::vector<Core::LinAlg::Matrix<3, 3>>
              velDeriv,  //! velocity gradient at point x (=x_Lagr(t^n+1))
          std::vector<Core::LinAlg::Matrix<1, 3>>
              presDeriv,                      //! pressure gradient at point x (=x_Lagr(t^n+1))
          Core::LinAlg::Matrix<3, 1> dispnp,  //! displacement at point x (=x_Lagr(t^n+1))
          Core::LinAlg::Matrix<3, 1>
              startpoint,  //! current startpoint (Lagrangean origin) approximation
          int searchedProcs, int counter, double dMin,
          Type newtype  //! basic computation type due to FGI and FRS
          )
          : node_(node),
            nds_np_(nds_np),
            vel_(vel),
            velDeriv_(velDeriv),
            presDeriv_(presDeriv),
            dispnp_(dispnp),
            startpoint_(startpoint),
            last_valid_ele_(-1),
            last_valid_vc_(nullptr),
            changedside_(false),
            accepted_(false),
            initial_eid_(-1),
            initial_ele_owner_(-1),
            proj_(failed_),
            searchedProcs_(searchedProcs),
            counter_(counter),
            dMin_(dMin),
            state_(basicStd_),
            type_(newtype)
      {
        return;
      }

      //! constructor for current data in Semi-Lagrangean algorithm (used for receiving data from
      //! other procs)
      TimeIntData(Core::Nodes::Node& node,
          int nds_np,  //! nds (nodal dofset) number w.r.t new interface position, for which SL-algo
                       //! is called
          Core::LinAlg::Matrix<3, 1>& vel, std::vector<Core::LinAlg::Matrix<3, 3>>& velDeriv,
          std::vector<Core::LinAlg::Matrix<1, 3>>& presDeriv,
          Core::LinAlg::Matrix<3, 1> dispnp,  //! displacement at point x (=x_Lagr(t^n+1))
          Core::LinAlg::Matrix<3, 1>& initialpoint, int initial_eid, int initial_ele_owner,
          Core::LinAlg::Matrix<3, 1>& startpoint, int searchedProcs, int counter, Type newtype)
          : node_(node),
            nds_np_(nds_np),
            vel_(vel),
            velDeriv_(velDeriv),
            presDeriv_(presDeriv),
            dispnp_(dispnp),
            startpoint_(startpoint),
            last_valid_ele_(-1),
            last_valid_vc_(nullptr),
            changedside_(false),
            accepted_(false),
            initialpoint_(initialpoint),
            initial_eid_(initial_eid),
            initial_ele_owner_(initial_ele_owner),
            proj_(failed_),
            searchedProcs_(searchedProcs),
            counter_(counter),
            dMin_(INFINITY),
            state_(currSL_),
            type_(newtype)
      {
        return;
      }

      // TODO: cleanup these constructors!

      //! constructor for failed data in Semi-Lagrange, standard computation
      TimeIntData(Core::Nodes::Node& node,
          int nds_np,  //! nds (nodal dofset) number w.r.t new interface position, for which SL-algo
                       //! is called
          Core::LinAlg::Matrix<3, 1>& vel, std::vector<Core::LinAlg::Matrix<3, 3>>& velDeriv,
          std::vector<Core::LinAlg::Matrix<1, 3>>& presDeriv,
          Core::LinAlg::Matrix<3, 1> dispnp,  //! displacement at point x (=x_Lagr(t^n+1))
          Core::LinAlg::Matrix<3, 1>& initialpoint, int initial_eid, int initial_ele_owner,
          Type newtype)
          : node_(node),
            nds_np_(nds_np),
            vel_(vel),
            velDeriv_(velDeriv),
            presDeriv_(presDeriv),
            dispnp_(dispnp),
            initialpoint_(initialpoint),
            initial_eid_(initial_eid),
            initial_ele_owner_(initial_ele_owner),
            proj_(failed_),
            state_(failedSL_),
            type_(newtype)
      {
        return;
      }

      //! constructor for done data in standard computation (at the end of Semi-Lagrangean
      //! algorithm)
      TimeIntData(Core::Nodes::Node& node,
          int nds_np,  //! nds (nodal dofset) number w.r.t new interface position, for which SL-algo
                       //! is called
          Core::LinAlg::Matrix<3, 1> dispnp,  //! displacement at point x (=x_Lagr(t^n+1))
          Core::LinAlg::Matrix<3, 1>& startpoint,
          std::vector<Core::LinAlg::Matrix<3, 1>>& velValues, std::vector<double>& presValues,
          Type newtype)
          : node_(node),
            nds_np_(nds_np),
            dispnp_(dispnp),
            startpoint_(startpoint),
            velValues_(velValues),
            presValues_(presValues),
            state_(doneStd_),
            type_(newtype)
      {
        return;
      }

      //@}



      /*========================================================================*/
      //! @name to-String methods
      /*========================================================================*/

      //! std::string representation of the data's state
      std::string stateToString()
      {
        std::string output;
        switch (state_)
        {
          case basicStd_:
            output = "basic standard";
            break;
          case currSL_:
            output = "current Semi-Lagrange";
            break;
          case nextSL_:
            output = "next Semi-Lagrange";
            break;
          case failedSL_:
            output = "failed Semi-Lagrange";
            break;
          case doneStd_:
            output = "done Standard";
            break;
          default:
            FOUR_C_THROW("unknown status in TimeIntData class");
            break;
        }
        return output;
      };

      //! std::string representation of the data's type
      std::string typeToString()
      {
        std::string output;
        switch (type_)
        {
          case standard_:
            output = "standard type";
            break;
          case predictor_:
            output = "predictor";
            break;
          default:
            FOUR_C_THROW("unknown type in TimeIntData class");
            break;
        }
        return output;
      }

      //@}


      /*========================================================================*/
      //! @name data for the start-point subclass
      /*========================================================================*/

      //------------------------------------------
      //! @name data for the node and it's dofset for which the XFEM-timeintegration algorithm is
      //! called
      //------------------------------------------
      Core::Nodes::Node node_;  //! node for which SL-algorithm is called, no pointer!
      int nds_np_;  //! nds (nodal dofset) number w.r.t new interface position, for which SL-algo is
                    //! called

      //------------------------------------------
      //! @name data for the quantities at the node
      //------------------------------------------
      Core::LinAlg::Matrix<3, 1> vel_;  //! velocity at point x (=x_Lagr(t^n+1))=node, only for the
                                        //! real velocity vector (transport velocity)
      std::vector<Core::LinAlg::Matrix<3, 3>>
          velDeriv_;  //! velocity gradient at point x (=x_Lagr(t^n+1))=node, averaged around
                      //! elements, for all vectors to be reconstructed
      std::vector<Core::LinAlg::Matrix<1, 3>>
          presDeriv_;  //! pressure gradient at point x (=x_Lagr(t^n+1))=node, averaged around
                       //! elements, for all vectors to be reconstructed

      Core::LinAlg::Matrix<3, 1> dispnp_;  //! ale displacement at point x (=x_Lagr(t^n+1))=node

      //------------------------------------------
      //! @name data for current startpoint (Lagrangean origin) approximation
      //------------------------------------------
      Core::LinAlg::Matrix<3, 1>
          startpoint_;        //! current startpoint (Lagrangean origin) approximation
      std::vector<int> nds_;  //! current dofset number for startpoint approximation
      std::vector<int>
          last_valid_nds_;  //! last valid nds vector (used, when iteration of Lagrangean origin
                            //! moves inside the structure, or changes the side
      int last_valid_ele_;  //! last valid element id (used, when iteration of Lagrangean origin
                            //! moves inside the structure, or changes the side
      Core::Geo::Cut::VolumeCell*
          last_valid_vc_;  //! last valid fluid volumecell (used, when iteration of Lagrangean
                           //! origin moves inside the structure, or changes the side
      bool changedside_;   //! changed side compared to initial point
      bool accepted_;      //! has the current approximation been accepted as Lagrangean origin

      //------------------------------------------
      //! @name data for initial point (the point where we start the search for the Lagrangean
      //! origin)
      //------------------------------------------
      Core::LinAlg::Matrix<3, 1>
          initialpoint_;       //! first startpoint approximation at the right side
      int initial_eid_;        //! element Id of element where the initialpoint lies in
      int initial_ele_owner_;  //! owner of the element where the initialpoint lies in
      Projection proj_;  //! projection of node on structural movement is on side/line/point/failed

      //------------------------------------------
      //! @name staff used for parallelization of the Semi-Lagrangean Newton-loop
      //------------------------------------------

      // TODO: cleanup these variables
      // are these variables used?
      int searchedProcs_;  // searched procnumber if point lies in
      int counter_;        // newton iteration counters for the points
      double dMin_;        // minimal distance in initialization

      //------------------------------------------
      //! @name computed/reconstructed pressure and velocity values
      //------------------------------------------
      std::vector<Core::LinAlg::Matrix<3, 1>> velValues_;  // computed velocity for a node
      std::vector<double> presValues_;                     // computed pressure for a node

      //------------------------------------------
      //! @name state, type of time-integration approach
      //------------------------------------------
      State
          state_;  //! status for both the used special algorithm and the state within the algorithm
      Type type_;  //! basic computation type due to FGI and FRS

     protected:
      explicit TimeIntData();  // don't want default constructor

    };  // end class TimeIntData


   public:
    //! constructor
    explicit XfluidTimeintBase(
        const Teuchos::RCP<Discret::Discretization> discret,      /// background discretization
        const Teuchos::RCP<Discret::Discretization> boundarydis,  /// cut discretization
        Teuchos::RCP<Core::Geo::CutWizard> wizard_old,  /// cut wizard w.r.t. old interface position
        Teuchos::RCP<Core::Geo::CutWizard> wizard_new,  /// cut wizard w.r.t. new interface position
        Teuchos::RCP<XFEM::XFEMDofSet> dofset_old,  /// XFEM dofset w.r.t. old interface position
        Teuchos::RCP<XFEM::XFEMDofSet> dofset_new,  /// XFEM dofset w.r.t. new interface position
        std::vector<Teuchos::RCP<Epetra_Vector>>
            oldVectors,                      /// vector of col-vectors w.r.t. old interface position
        Teuchos::RCP<Epetra_Vector> dispn,   /// displacment n
        Teuchos::RCP<Epetra_Vector> dispnp,  /// displacment n +1
        const Epetra_Map& olddofcolmap,      /// dofcolmap w.r.t. old interface position
        const Epetra_Map& newdofrowmap,      /// dofcolmap w.r.t. new interface position
        const Teuchos::RCP<std::map<int, std::vector<int>>>
            pbcmap  /// map of periodic boundary conditions
    );

    //! destructor
    virtual ~XfluidTimeintBase() = default;
    //! perform the computation
    virtual void compute(std::vector<Teuchos::RCP<Epetra_Vector>> newRowVectorsn,
        std::vector<Teuchos::RCP<Epetra_Vector>> newRowVectorsnp)
    {
      FOUR_C_THROW("Unused function! Use a function of the derived classes");
    };

    //! set computation type due to iteration counter
    void type(int iter, int iterMax);

    /*========================================================================*/
    //! @name  enumerations used for XFEM time-integration
    /*========================================================================*/

    //! computation type due to input data (FGI,FRS)
    enum FGIType
    {
      FRS1FGI1_,     /// first FRS iteration of first FGI
      FRS1FGINot1_,  /// first FRS iteration of first FGI
      FRSNot1_       /// first FGI
    };

    FGIType FGIType_;  // computation type due to input data (FGI,FRS)


    /*========================================================================*/
    //! @name access methods
    /*========================================================================*/

    //! time-integration for how many nodes called?
    int numNodes() { return timeIntData_->size(); };


   protected:
    //! initialize data to be set in every computation
    void handle_vectors(std::vector<Teuchos::RCP<Epetra_Vector>>& newRowVectorsn);

    //! return the number of Epetra vectors which shall get new values for a given node with
    //! according data
    size_t vector_size(TimeIntData* data) const { return vector_size(data->type_); };

    //! return the number of Epetra vectors which shall get new values for a given type
    size_t vector_size(TimeIntData::Type currtype) const
    {
      if (newVectors_.size() == 0)
        FOUR_C_THROW("call this function only when setting the final values");

      size_t size = 0;  // vector size
      switch (currtype)
      {
        case TimeIntData::predictor_:
          size = newVectors_.size();
          break;  // old solutions and new initialization of new solutions
        case TimeIntData::standard_:
          size = newVectors_.size() / 2;
          break;  // only first half with old solutions
        default:
          FOUR_C_THROW("undefined type");
          break;
      }
      return size;
    };

    //! overwrite an old state with a new state
    void reset_state(TimeIntData::State oldState, TimeIntData::State newState) const;

    //! clear all data having some state
    void clear_state(TimeIntData::State state  /// state of time int to clear
    ) const;


    /*========================================================================*/
    //! @name  interface-side compare method based on geometrical information
    /*========================================================================*/

    //! check if the current point x2 at time t^n changed the side compared to x1
    bool changed_side(Core::Elements::Element* ele1,  /// first element where x1 lies in
        Core::LinAlg::Matrix<3, 1>& x1,               /// global coordinates of point x1
        bool newTimeStep1,                            /// new/old timestep for point x1
        Core::Elements::Element* ele2,                /// second element where x2 lies in
        Core::LinAlg::Matrix<3, 1>& x2,               /// global coordinates of point x2
        bool newTimeStep2                             /// new/old timestep for point x2
    ) const
    {
      if (newTimeStep1 != newTimeStep2)
      {
        FOUR_C_THROW("how to check changing side at two different time steps?");
      }
      else
      {
        return changed_side_same_time(newTimeStep1, ele1, x1, ele2, x2);
      }

      return true;
    }


    //! check if the current point x2 changed the side compared to x1
    bool changed_side_same_time(
        const bool newTimeStep,          /// new/old timestep for both points x1 and x2
        Core::Elements::Element* ele1,   /// first element where x1 lies in
        Core::LinAlg::Matrix<3, 1>& x1,  /// global coordinates of point x1
        Core::Elements::Element* ele2,   /// second element where x2 lies in
        Core::LinAlg::Matrix<3, 1>& x2   /// global coordinates of point x2
    ) const;


    //! check if both element are neighbors sharing at least one common node and get the common
    //! nodes
    bool neighbors(Core::Elements::Element* ele1, Core::Elements::Element* ele2,
        std::set<int>& common_nodes) const;


    //! check if edge between x1 and x2 cuts the side
    bool call_side_edge_intersection(Core::Geo::Cut::SideHandle* sh,  /// side handle
        int sid,                                                      /// side id
        Core::LinAlg::Matrix<3, 1>& x1,  /// coordinates of edge's start point
        Core::LinAlg::Matrix<3, 1>& x2   /// coordinates of edge's end point
    ) const;


    //! check if edge between x1 and x2 cuts the side
    template <Core::FE::CellType sidetype>
    bool call_side_edge_intersection_t(Core::Geo::Cut::SideHandle* sh,  /// side handle
        int sid,                                                        /// side id
        Core::LinAlg::Matrix<3, 1>& x1,  /// coordinates of edge's start point
        Core::LinAlg::Matrix<3, 1>& x2   /// coordinates of edge's end point
    ) const;


    /*========================================================================*/
    //! @name periodic boundary conditions
    /*========================================================================*/

    //! add adjacebt elements for a periodic boundary node
    void add_pb_celements(
        const Core::Nodes::Node* node, std::vector<const Core::Elements::Element*>& eles) const;

    //! find the PBC node
    void find_pbc_node(
        const Core::Nodes::Node* node, Core::Nodes::Node*& pbcnode, bool& pbcnodefound) const;

    //@}

    /*========================================================================*/
    //! @name parallel communication
    /*========================================================================*/

    //! basic function for parallel sending of data
    void send_data(Core::Communication::PackBuffer& dataSend, int& dest, int& source,
        std::vector<char>& dataRecv) const;

    //! packing a node
    void pack_node(Core::Communication::PackBuffer& dataSend, Core::Nodes::Node& node) const;

    //! unpacking a node
    void unpack_node(std::vector<char>::size_type& posinData, std::vector<char>& dataRecv,
        Core::Nodes::Node& node) const;


    /*========================================================================*/
    //! @name local-global transformations
    /*========================================================================*/

    //! transformation of point x to local xi-coordinates of an integration cell
    void call_x_to_xi_coords(const Core::Elements::Element* ele,  /// pointer to element
        Core::LinAlg::Matrix<3, 1>& x,                            /// global coordinates of element
        Core::LinAlg::Matrix<3, 1>& xi,  /// determined local coordinates w.r.t ele
        const std::string state,         ///< state n or np?
        bool& pointInDomain              /// lies point in element ?
    ) const;

    //! call the computation of local coordinates for a polytop with corners given by the
    //! coordinates
    void call_x_to_xi_coords(
        Core::LinAlg::SerialDenseMatrix& nodecoords,  /// node coordinates of element
        Core::FE::CellType DISTYPE,                   /// discretization type
        Core::LinAlg::Matrix<3, 1>& x,                /// global coordinates of element
        Core::LinAlg::Matrix<3, 1>& xi,               /// determined local coordinates w.r.t ele
        bool& pointInDomain                           /// lies point in element ?
    ) const;

    //! compute local element coordinates and check whether the according point is inside the
    //! element
    template <Core::FE::CellType DISTYPE>
    void x_to_xi_coords(Core::LinAlg::SerialDenseMatrix& xyz,  /// node coordinates of element
        Core::LinAlg::Matrix<3, 1>& x,                         /// global coordinates of point
        Core::LinAlg::Matrix<3, 1>& xi,  /// determined local coordinates w.r.t ele
        bool& pointInCell                /// lies point in element?
    ) const;

    //! data at an arbitrary point lying in an element
    template <const int numnode, Core::FE::CellType DISTYPE>
    void eval_shape_and_deriv(Core::Elements::Element* element,  /// pointer to element
        Core::LinAlg::Matrix<3, 1>& xi,              /// local coordinates of point w.r.t element
        Core::LinAlg::Matrix<3, 3>& xji,             /// inverse of jacobian
        Core::LinAlg::Matrix<numnode, 1>& shapeFcn,  /// shape functions at point
        Core::LinAlg::Matrix<3, numnode>&
            shapeFcnDerivXY,  /// derivatives of shape function w.r.t global coordiantes xyz
        bool compute_deriv    /// shall derivatives and jacobian be computed
    ) const;

    /*========================================================================*/
    //! @name  data w.r.t old and new interface position
    /*========================================================================*/

    Teuchos::RCP<Discret::Discretization> discret_;      //! background discretization
    Teuchos::RCP<Discret::Discretization> boundarydis_;  //! cut discretization

    Teuchos::RCP<Core::Geo::CutWizard> wizard_old_;  //! cut wizard w.r.t. old interface position
    Teuchos::RCP<Core::Geo::CutWizard> wizard_new_;  //! cut wizard w.r.t. new interface position

    Teuchos::RCP<XFEM::XFEMDofSet> dofset_old_;  //! XFEM dofset w.r.t. old interface position
    Teuchos::RCP<XFEM::XFEMDofSet> dofset_new_;  //! XFEM dofset w.r.t. new interface position

    const Epetra_Map olddofcolmap_;  //! dofcolmap w.r.t. old interface position
    const Epetra_Map newdofrowmap_;  //! dofcolmap w.r.t. new interface position

    const std::vector<Teuchos::RCP<Epetra_Vector>>
        oldVectors_;  //! vector of col!-vectors w.r.t. old interface position
    std::vector<Teuchos::RCP<Epetra_Vector>>
        newVectors_;  //! vector of row!-vectors w.r.t. new interface position (overwritten with new
                      //! information for non-predictor case and filled otherwise)

    Teuchos::RCP<Epetra_Vector>
        dispn_;  //! col!-displacement vector for timestep n w.r.t. old interface position
    Teuchos::RCP<Epetra_Vector>
        dispnp_;  //! col!-displacement vector for timestep n + 1 w.r.t. old interface position

    Teuchos::RCP<std::vector<TimeIntData>>
        timeIntData_;  //! data-vector containing all data for computation with help of sub-class
                       //! TimeIntData

    //@}

    /*========================================================================*/
    //! @name stuff for parallel communication
    /*========================================================================*/

    const Teuchos::RCP<std::map<int, std::vector<int>>>
        pbcmap_;         //! map for master and slave elements of nodes
    const int myrank_;   //! current processor id
    const int numproc_;  //! number of processors

    //@}

    /*========================================================================*/
    //! constants
    /*========================================================================*/

    const int
        newton_max_iter_;      //! maximal number of newton iterations for Semi-Lagrangean algorithm
    const double limits_tol_;  //! newton tolerance for Semi-Lagrangean algorithm
    const double TOL_dist_;  //! tolerance to find the shortest distance of point to its projection
                             //! on the surface dis

    //@}


  };  // class TimeInt



  /*!
  \brief this class provides the basic functionality used for the reconstruction of standard degrees
         of freedom in XFEM-problems between consecutive time steps
   */
  class XfluidStd : public XfluidTimeintBase
  {
   public:
    /*========================================================================*/
    //! constructor/destructor
    /*========================================================================*/

    //! constructor
    explicit XfluidStd(XFEM::XfluidTimeintBase& timeInt,  ///< time integration base class object
        const std::map<int, std::vector<Inpar::XFEM::XFluidTimeInt>>&
            reconstr_method,                      ///< reconstruction map for nodes and its dofsets
        Inpar::XFEM::XFluidTimeInt& timeIntType,  ///< type of time integration
        const Teuchos::RCP<Epetra_Vector> veln,   ///< velocity at time t^n
        const double& dt,                         ///< time step size
        const bool initialize                     ///< is initialization?
    );

    /*========================================================================*/
    //! compute routines
    /*========================================================================*/

    //! perform the computation
    virtual void compute(std::vector<Teuchos::RCP<Epetra_Vector>>& newRowVectorsn);

    /*========================================================================*/
    //! element based routines
    /*========================================================================*/

    //! identify an element containing a point and additional data
    void elementSearch(
        Core::Elements::Element*& ele,   ///< pointer to element if point lies in a found element
        Core::LinAlg::Matrix<3, 1>& x,   ///< global coordiantes of point
        Core::LinAlg::Matrix<3, 1>& xi,  ///< determined local coordinates w.r.t ele
        bool& found                      ///< is element found?
    ) const;

    //! interpolate velocity and derivatives for a point in an element
    void getGPValues(Core::Elements::Element* ele,  ///< pointer to element
        Core::LinAlg::Matrix<3, 1>& xi,             ///< local coordinates of point w.r.t element
        std::vector<int>& nds,                      ///< nodal dofset of point for elemental nodes
        XFEM::XFEMDofSet& dofset,                   ///< XFEM dofset
        Core::LinAlg::Matrix<3, 1>& vel,            ///< determine velocity at point
        Core::LinAlg::Matrix<3, 3>& vel_deriv,      ///< determine velocity derivatives at point
        double& pres,                               ///< pressure
        Core::LinAlg::Matrix<1, 3>& pres_deriv,     ///< pressure gradient
        Teuchos::RCP<const Epetra_Vector> vel_vec,  ///< vector used for interpolating at gp
        bool compute_deriv = true                   ///< shall derivatives be computed?
    ) const;

    //! interpolate velocity and derivatives for a point in an element
    template <Core::FE::CellType DISTYPE>
    void getGPValuesT(Core::Elements::Element* ele,  ///< pointer to element
        Core::LinAlg::Matrix<3, 1>& xi,              ///< local coordinates of point w.r.t element
        std::vector<int>& nds,                       ///< nodal dofset of point for elemental nodes
        XFEM::XFEMDofSet& dofset,                    ///< xfem dofset
        Core::LinAlg::Matrix<3, 1>& vel,             ///< determine velocity at point
        Core::LinAlg::Matrix<3, 3>& vel_deriv,       ///< determine velocity derivatives at point
        double& pres,                                ///< pressure
        Core::LinAlg::Matrix<1, 3>& pres_deriv,      ///< pressure gradient
        Teuchos::RCP<const Epetra_Vector> vel_vec,   ///< vector used for interpolating at gp
        bool compute_deriv = true                    ///< shall derivatives be computed?
    ) const;

    //!  back-tracking of data at final Lagrangian origin of a point                       schott
    //!  04/12 *
    template <const int numnode>
    void extract_nodal_values_from_vector(
        Core::LinAlg::Matrix<3, numnode>& evel,  ///< element velocities
        Core::LinAlg::Matrix<numnode, 1>& epre,  ///< element pressure
        Teuchos::RCP<Epetra_Vector> vel_vec,     ///< global vector
        std::vector<int>& lm                     ///< local map
    )
    {
      const int nsd = 3;
      const int numdofpernode = nsd + 1;

      evel.Clear();
      epre.Clear();

      if (vel_vec == Teuchos::null) FOUR_C_THROW("vector is null");

      // extract local values of the global vectors
      std::vector<double> mymatrix(lm.size());
      Core::FE::ExtractMyValues(*vel_vec, mymatrix, lm);

      for (int inode = 0; inode < numnode; ++inode)  // number of nodes
      {
        for (int idim = 0; idim < nsd; ++idim)  // number of dimensions
        {
          (evel)(idim, inode) = mymatrix[idim + (inode * numdofpernode)];
        }
        (epre)(inode, 0) = mymatrix[nsd + (inode * numdofpernode)];
      }

      return;
    }

    /*========================================================================*/
    //! parallel routines
    /*========================================================================*/

    //! export data to startpoint processor when Semi-Lagrange algorithm failed
    void exportStartData();

    //! export final data to the proc where the according node is
    void exportFinalData();



   protected:
    /*========================================================================*/
    //! functions for computation of a reasonable startpoint to find the Lagrangean origin
    /*========================================================================*/

    //! determine a first starting point approximation for nodes  marked for Semi-Lagrangean
    //! algorithm
    void startpoints();

    //! project point on structural surface and track the projected point back at t^n
    void project_and_trackback(TimeIntData& data);

    //! Project the current point onto the structural surface given by point and side Ids
    bool project_to_surface(Core::LinAlg::Matrix<3, 1>& x,  ///< coords of point to be projected
        Core::LinAlg::Matrix<3, 1>& proj_x,                 ///< coords of projected point
        const std::string state,                            ///< state n or np?
        std::set<int>& points,  ///< node Ids of surface for that distance has to be computed
        std::set<int>& sides  ///< side Ids of surface for that and it's lines the distances have to
                              ///< be computed
    );

    //! find the nearest surface point, return if successful
    bool find_nearest_surf_point(
        Core::LinAlg::Matrix<3, 1>& x,       ///< coords of point to be projected
        Core::LinAlg::Matrix<3, 1>& proj_x,  ///< coords of projected point
        Core::Geo::Cut::VolumeCell* vc,      ///< volumcell on that's cut-surface we want to project
        const std::string state              ///< state n or np?
    );

    //! call and prepare the projection of point to side
    void call_project_on_side(
        Core::Elements::Element* side,  ///< pointer to structural side element
        const std::string state,        ///< state n or np?
        Core::LinAlg::Matrix<3, 1>&
            newNodeCoords,  ///< node coordinates of point that has to be projected
        int& proj_sid,      ///< id of side when projection lies on side
        double& min_dist,   ///< minimal distance, potentially updated
        Core::LinAlg::Matrix<3, 1>& proj_x_np,  ///< projection of point on this side
        Core::LinAlg::Matrix<2, 1>&
            proj_xi_side,  ///< local coordinates of projection of point w.r.t to this side
        TimeIntData::Projection& proj  ///< projection type
    );

    //! call and prepare the projection of point to line
    void call_project_on_line(
        Core::Elements::Element* side,  ///< pointer to structural side element
        Core::Elements::Element* line,  ///< pointer to structural line of side element
        int line_count,                 ///< local line id w.r.t side element
        const std::string state,        ///< state n or np?
        Core::LinAlg::Matrix<3, 1>&
            newNodeCoords,  ///< node coordinates of point that has to be projected
        double& min_dist,   ///< minimal distance, potentially updated
        Core::LinAlg::Matrix<3, 1>& proj_x_np,  ///< projection of point on this side
        std::map<std::vector<int>, std::vector<double>>&
            proj_xi_line,  ///< std::map<sorted nids, local line coordinates of projection of point
                           ///< w.r.t sides >
        std::map<std::vector<int>, std::vector<int>>&
            proj_lineid,  ///< std::map<sorted nids, local line id w.r.t sides>
        std::map<std::vector<int>, std::vector<int>>&
            proj_nid_line,             ///< std::map<sorted nids, side Ids>
        int& proj_sid,                 ///< id of side that contains the projected point
        TimeIntData::Projection& proj  ///< projection type
    );

    //! call and prepare the projection of point to point (distance computation)
    void call_project_on_point(Core::Nodes::Node* node,  ///< pointer to node
        const std::string state,                         ///< state n or np?
        Core::LinAlg::Matrix<3, 1>&
            newNodeCoords,  ///< node coordinates of point that has to be projected
        double& min_dist,   ///< minimal distance, potentially updated
        int& proj_nid_np,   ///< nid id of projected point on surface
        Core::LinAlg::Matrix<3, 1>& proj_x_np,  ///< projection of point on this side
        int& proj_sid,                          ///< id of side that contains the projected point
        std::map<std::vector<int>, std::vector<double>>
            proj_xi_line,  ///< std::map<side ID, local coordinates of projection of point w.r.t to
                           ///< this line>
        std::map<std::vector<int>, std::vector<int>>
            proj_lineid,  ///< std::map<side ID, local line id>
        std::map<std::vector<int>, std::vector<int>>
            proj_nid_line,             ///< std::map<side ID, vec<line Ids>>
        TimeIntData::Projection& proj  ///< projection type
    );

    //! project point from in normal direction onto corresponding side
    template <Core::FE::CellType side_distype, const int numdof>
    bool project_on_side(Core::LinAlg::SerialDenseMatrix& side_xyze,  ///< side's node coordinates
        const std::vector<int>& lm,                                   ///< local map
        const std::string state,                                      ///< state n or np
        Core::LinAlg::Matrix<3, 1>&
            x_gp_lin,  ///< global coordinates of point that has to be projected
        Core::LinAlg::Matrix<3, 1>& x_side,   ///< projected point on side
        Core::LinAlg::Matrix<2, 1>& xi_side,  ///< local coordinates of projected point w.r.t side
        double& dist                          ///< distance from point to its projection
    );

    //! project point in normal direction onto corresponding line
    template <Core::FE::CellType line_distype, const int numdof>
    bool project_on_line(Core::LinAlg::SerialDenseMatrix& line_xyze,  ///< line's node coordinates
        const std::vector<int>& lm,                                   ///< local map
        const std::string state,                                      ///< state n or np
        Core::LinAlg::Matrix<3, 1>&
            x_point_np,  ///< global coordinates of point that has to be projected
        Core::LinAlg::Matrix<3, 1>& x_line,  ///< projected point on line
        double& xi_line,                     ///< local coordinates of projected point w.r.t line
        double& dist                         ///< distance from point to its projection
    );

    //! compute distance (project) between two points
    void project_on_point(
        Core::LinAlg::SerialDenseMatrix& point_xyze,  ///< point's node coordinates
        const std::vector<int>& lm,                   ///< local map
        const std::string state,                      ///< state n or np
        Core::LinAlg::Matrix<3, 1>&
            x_point_np,  ///< global coordinates of point that has to be projected
        Core::LinAlg::Matrix<3, 1>& x_point,  ///< global coordinates of point on surface
        double& dist                          ///< distance from point to its projection
    );

    //! check if local coordinates are within limits
    template <Core::FE::CellType elementtype>
    bool within_limits(Core::LinAlg::Matrix<3, 1>& xsi_, const double TOL);

    //! compute reasonable start point for finding the Lagrangean origin, when projected point lies
    //! on a side
    void compute_start_point_side(Core::Elements::Element* side,  ///< pointer to side element
        Core::LinAlg::SerialDenseMatrix& side_xyze,               ///< side's node coordinates
        const std::vector<int>& lm,                               ///< local map
        Core::LinAlg::Matrix<2, 1>& xi_side,   ///< local coordinates of projected point w.r.t side
        double& dist,                          ///< distance from point to its projection
        Core::LinAlg::Matrix<3, 1>& proj_x_n,  ///< projected point at t^n
        Core::LinAlg::Matrix<3, 1>& start_point  ///< final start point
    );

    //! compute reasonable start point for finding the Lagrangean origin, when projected point lies
    //! on a line
    void compute_start_point_line(Core::Elements::Element* side1,  ///< pointer to side element
        Core::LinAlg::SerialDenseMatrix& side1_xyze,               ///< side's node coordinates
        Core::Elements::Element* side2,                            ///< pointer to side element
        Core::LinAlg::SerialDenseMatrix& side2_xyze,               ///< side's node coordinates
        const std::vector<int>& lm1,                               ///< local map
        const std::vector<int>& lm2,                               ///< local map
        double& dist,                            ///< distance from point to its projection
        Core::LinAlg::Matrix<3, 1>& proj_x_n,    ///< projected point at t^n
        Core::LinAlg::Matrix<3, 1>& start_point  ///< final start point
    );

    //! compute reasonable start point for finding the Lagrangean origin, when projected point lies
    //! on a point
    void compute_start_point_avg(
        const std::vector<Core::Elements::Element*>& sides,        ///< pointer to side element
        std::vector<Core::LinAlg::SerialDenseMatrix>& sides_xyze,  ///< side's node coordinates
        const std::vector<std::vector<int>>& sides_lm,             ///< local map
        double& dist,                            ///< distance from point to its projection
        Core::LinAlg::Matrix<3, 1>& proj_x_n,    ///< projected point at t^n
        Core::LinAlg::Matrix<3, 1>& start_point  ///< final start point
    );

    //! compute the normal vector on side at time t^n
    void callget_normal_side_tn(Core::Elements::Element* side,  ///< pointer to side element
        Core::LinAlg::Matrix<3, 1>& normal,                     ///< normal vector w.r.t side
        Core::LinAlg::SerialDenseMatrix& side_xyze,             ///< side's node coordinates
        const std::vector<int>& lm,                             ///< local map
        Core::LinAlg::Matrix<3, 1>& proj_x_n,                   ///< projected point on side
        Core::LinAlg::Matrix<2, 1>& xi_side  ///< local coordinates of projected point w.r.t side
    );

    //! compute the normal vector on side at time t^n
    template <Core::FE::CellType side_distype, const int numdof>
    void get_normal_side_tn(Core::LinAlg::Matrix<3, 1>& normal,  ///< normal vector w.r.t side
        Core::LinAlg::SerialDenseMatrix& side_xyze,              ///< side's node coordinates
        const std::vector<int>& lm,                              ///< local map
        Core::LinAlg::Matrix<3, 1>& proj_x_n,                    ///< projected point on side
        Core::LinAlg::Matrix<2, 1>& xi_side  ///< local coordinates of projected point w.r.t side
    );

    void call_get_projxn_line(
        Core::Elements::Element* side,         ///< pointer to structural side element
        Core::Elements::Element* line,         ///< pointer to structural line of side element
        Core::LinAlg::Matrix<3, 1>& proj_x_n,  ///< projected point on line
        double& xi_line                        ///< local coordinates of projected point w.r.t line
    );

    template <Core::FE::CellType line_distype, const int numdof>
    void get_projxn_line(Core::LinAlg::SerialDenseMatrix& line_xyze,  ///< line's node coordinates
        const std::vector<int>& lm,                                   ///< local map
        Core::LinAlg::Matrix<3, 1>& proj_x_n,                         ///< projected point on line
        double& xi_line  ///< local coordinates of projected point w.r.t line
    );

    //! add side's or line's interface displacements and set current node coordinates
    template <Core::FE::CellType distype, const int numdof>
    void addeidisp(Core::LinAlg::SerialDenseMatrix& xyze,  ///< node coordinates of side or line
        const Discret::Discretization& cutdis,             ///< cut discretization
        const std::string state,                           ///< state
        const std::vector<int>& lm                         ///< local map
    );

    //! set the final startvalues with respect to a node for finding the Lagrangean origin
    void set_final_data();

    /*========================================================================*/
    //! time-integration variables
    /*========================================================================*/

    Inpar::XFEM::XFluidTimeInt
        timeIntType_;  //! which computation/reconstruction algorithm for standard dofs

    Teuchos::RCP<Epetra_Vector> veln_;  //! velocity w.r.t old interface position in column map

    const double dt_;  //! time step size

  };  // class XFLUID_TIMEINT_STD

}  // namespace XFEM

FOUR_C_NAMESPACE_CLOSE

#endif
