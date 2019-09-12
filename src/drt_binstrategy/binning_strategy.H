/*----------------------------------------------------------------------*/
/*! \file

\brief Binning strategy for neighborhood search

\level 2

\maintainer Jonas Eichinger
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              ghamm 09/12 |
 *----------------------------------------------------------------------*/
#ifndef BINNING_STRATEGY_H
#define BINNING_STRATEGY_H

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 09/12 |
 *----------------------------------------------------------------------*/
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_particle_old.H"
#include "../drt_binstrategy/binning_strategy_utils.H"


#include "Epetra_MpiComm.h"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_CrsGraph.h>
#include <vector>

/*----------------------------------------------------------------------*
 | forward declarations                                     ghamm 09/12 |
 *----------------------------------------------------------------------*/
namespace DRT
{
  class Discretization;
  class Element;
  class Node;
}  // namespace DRT
namespace GEO
{
  namespace MESHFREE
  {
    class BoundingBox;
  }
}  // namespace GEO

/*----------------------------------------------------------------------*
 | binning strategy                                         ghamm 11/13 |
 *----------------------------------------------------------------------*/
namespace BINSTRATEGY
{
  // forward declaration
  class Less;

  class BinningStrategy
  {
   public:
    /// standard constructor
    BinningStrategy();

    /// constructor
    BinningStrategy(const Epetra_Comm& lcomm,  ///< local epetra communicator
        double cutoff_radius,                  ///< radius of influence
        LINALG::Matrix<3, 2> XAABB             ///< extended axis aligned bounding box
    );

    /// constructor
    BinningStrategy(const Epetra_Comm& lcomm  ///< local epetra communicator
    );

    /// destructor
    virtual ~BinningStrategy(){};

    /// initialize binning strategy
    virtual void Init(std::vector<Teuchos::RCP<DRT::Discretization>> const discret,
        std::vector<Teuchos::RCP<Epetra_Vector>> disnp =
            std::vector<Teuchos::RCP<Epetra_Vector>>());

    /// initialize binning strategy
    virtual void Init(const Epetra_Comm& lcomm  ///< local epetra communicator
    );

    /// setup binning strategy
    virtual void Setup(Teuchos::RCP<DRT::Discretization>& bindis,
        Teuchos::RCP<GEO::MESHFREE::BoundingBox> const pbb);

    /// get all bin centers needed for repartitioning
    virtual void GetAllBinCenters(
        Teuchos::RCP<Epetra_Map>& binrowmap, Teuchos::RCP<Epetra_MultiVector>& bincenters) const;

    /// distribute bins via recursive coordinate bisection
    virtual void DistributeBinsRecursCoordBisection(Teuchos::RCP<Epetra_Map>& binrowmap,
        Teuchos::RCP<Epetra_MultiVector>& bincenters,
        Teuchos::RCP<Epetra_MultiVector>& binweights) const;

    /// fill bins into bin discretization
    virtual void FillBinsIntoBinDiscretization(Teuchos::RCP<Epetra_Map> const& rowbins);

    /// get binning discretization pointer (nodes = particles, elements = bins)
    virtual inline Teuchos::RCP<DRT::Discretization>& BinDiscret() { return bindis_; }

    /// get const binning discretization pointer (nodes = particles, elements = bins)
    virtual inline Teuchos::RCP<DRT::Discretization> const& BinDiscret() const { return bindis_; }

    /// get list of boundary row bins
    virtual inline std::list<DRT::Element*>& BoundaryRowBins() { return boundaryrowbins_; }

    /// get const list with boundary row bins
    virtual inline std::list<DRT::Element*> const& BoundaryRowBins() const
    {
      return boundaryrowbins_;
    }

    /// get list of boundary col bins
    virtual inline std::set<int>& BoundaryColBinsIds() { return boundarycolbins_; }

    /// get const list with boundary col bins
    virtual inline std::set<int> const& BoundaryColBinsIds() const { return boundarycolbins_; }

    /// get dimension of particle problem
    virtual inline INPAR::PARTICLEOLD::ParticleDim const& ParticleDim() const
    {
      return particle_dim_;
    }

    /// get cutoff radius
    virtual inline double const& CutoffRadius() const { return cutoff_radius_; }

    /// set cutoff radius
    virtual inline void SetCutoffRadius(double cutoff_radius) { cutoff_radius_ = cutoff_radius; }

    /// get bin size
    virtual inline double const* BinSize() const { return &bin_size_[0]; }

    /// get bin per dir
    virtual inline int const* BinPerDir() const { return &bin_per_dir_[0]; }

    /// return flag whether pbcs are applied to the problem
    virtual bool HavePBC(const int dim = -1) const
    {
      if (dim == -1)
        return havepbc_;
      else
        return pbconoff_[dim];
    };

    /// return delta for pbc in x, y, or z direction
    virtual double PBCDelta(const int dim) const { return pbcdeltas_[dim]; };

    /// get bounding box dimensions
    virtual inline LINALG::Matrix<3, 2>& XAABB() { return XAABB_; }
    virtual inline LINALG::Matrix<3, 2> const& XAABB() const { return XAABB_; }

    /// get inverse bin size
    virtual inline double const* InvBinSize() const { return &inv_bin_size_[0]; }

    /// elements are assigned to bins (either slave or master side of mortar interface)
    virtual void DistributeElesToBins(
        const DRT::Discretization& mortardis,     ///< mortar interface discretization
        std::map<int, std::set<int>>& binelemap,  ///< map of bins and assigned elements
        bool isslave  ///< decide whether slave or master side is processed
        ) const;

    /// distribute elements to bins exploiting axis aligned bounding box idea
    virtual void DistributeRowElementsToBinsUsingEleXAABB(
        Teuchos::RCP<DRT::Discretization> const&
            discret,  ///< discret of which the elements are distributed to bins
        std::map<int, std::set<int>>& bintorowelemap,  ///< map of bins and assigned row elements
        Teuchos::RCP<Epetra_Vector> disnp =
            Teuchos::null  ///< statevector holding current displacements of discret's nodes
        ) const;

    /// distribute element to bins using axis aligned bounding box idea
    virtual void DistributeElementToBinsUsingEleXAABB(
        Teuchos::RCP<DRT::Discretization> const& discret, DRT::Element* eleptr,
        std::vector<int>& binIds, Teuchos::RCP<Epetra_Vector> const& disnp) const;

    /// transfer nodes and elements of input discret
    virtual void TransferNodesAndElements(Teuchos::RCP<DRT::Discretization>& discret,
        Teuchos::RCP<Epetra_Vector> disnp, std::map<int, std::set<int>>& bintorowelemap);

    /// communicate elements for element transfer
    virtual void CommunicateElements(Teuchos::RCP<DRT::Discretization>& discret,
        std::map<int, std::vector<DRT::Element*>> const& toranktosendeles) const;

    /// number of bins in ijk range
    virtual int GetNumberOfBinsInijkRange(int const ijk_range[6]) const;

    /// get ijk of node with a current position
    virtual void DistributeNodeToijk(Teuchos::RCP<DRT::Discretization> const& discret,
        DRT::Node const* const node, Teuchos::RCP<Epetra_Vector> const& disnp, int ijk[3]) const;

    /// add ijk to a given axis aligned ijk range
    virtual void AddijkToAxisAlignedijkRangeOfElement(int const ijk[3], int ijk_range[6]) const;

    /// add ijk to a given axis aligned ijk range
    virtual void AddijkToAxisAlignedijkRangeOfBeamElement(int const ijk[3], int ijk_range[6]) const;

    /// build axis aligned element bounding box for one noded rigid sphere element
    virtual void BuildAxisAlignedijkRangeForRigidSphere(
        Teuchos::RCP<DRT::Discretization> const& discret, DRT::Element const* const sphereele,
        Teuchos::RCP<Epetra_Vector> const& disnp, int ijk[3], int ijk_range[6]) const;

    /// elements from the underlying discretization are assigned to bins
    virtual void AssignElesToBins(Teuchos::RCP<DRT::Discretization> discret,
        std::map<int, std::set<int>> const& extendedfieldghosting) const;

    /// get elements of type bincontent in given bins
    void GetBinContent(std::set<DRT::Element*>& eles,
        std::vector<BINSTRATEGY::UTILS::BinContentType> bincontent, std::vector<int>& binIds,
        bool roweles = false);

    /// get elements of type bincontent in given bin
    void GetBinContent(DRT::Element* binptr, std::set<DRT::Element*>& eles,
        std::vector<BINSTRATEGY::UTILS::BinContentType> bincontent, bool roweles);

    /// remove elements of specific kinds from bins
    virtual void RemoveSpecificElesFromBins(BINSTRATEGY::UTILS::BinContentType bincontent);

    /// remove elements of all kinds from bins
    virtual void RemoveAllElesFromBins();

    /// nodes are assigned to bins
    virtual void DistributeNodesToBins(Teuchos::RCP<DRT::Discretization> discret,
        std::map<int, std::vector<int>>& nodesinbin,
        Teuchos::RCP<Epetra_Vector> = Teuchos::null) const;

    /// elements/nodes from all provided discr. are distributed according to bin distrib.
    virtual Teuchos::RCP<Epetra_Map> WeightedPartitioning(
        std::vector<Teuchos::RCP<DRT::Discretization>> discret,
        std::vector<Teuchos::RCP<Epetra_Map>>& stdelecolmap,
        std::vector<Teuchos::RCP<Epetra_Map>>& stdnodecolmap);

    /// weighted distribution of bins to procs according to number of nodes they contain
    virtual Teuchos::RCP<Epetra_Map> WeightedDistributionOfBinsToProcs(
        std::vector<Teuchos::RCP<DRT::Discretization>>& discret,
        std::vector<Teuchos::RCP<Epetra_Vector>>& disnp,
        std::vector<std::map<int, std::vector<int>>>& nodesinbin, double const& weight,
        bool repartition = false) const;

    /// determine boundary row bins
    virtual void DetermineBoundaryRowBins();

    /// determine boundary col bins
    virtual void DetermineBoundaryColBinsIds();

    /*!
    \brief Extend ghosthing of mortar master interface discretization

    @param mortardis Discretization of mortar interface
    @param initial_elecolmap Initial element col map for standard ghosting
    @param slavebinelemap Map of bins and assigned slave elements
    @param masterbinelemenap Map of bins and assigned master elements
    @return Element map with extended ghosting for master side of mortar problem
    */
    virtual Teuchos::RCP<Epetra_Map> ExtendGhosting(DRT::Discretization& mortardis,
        Teuchos::RCP<Epetra_Map> initial_elecolmap, std::map<int, std::set<int>>& slavebinelemap,
        std::map<int, std::set<int>>& masterbinelemap) const;

    /// extend ghosting according to scatra elements
    virtual void ExtendGhosting(
        Teuchos::RCP<DRT::Discretization> scatradis,  ///< scatra discretization
        std::map<int, std::set<int>>&
            escapedpartelemap,  ///< map of escaped particles and corresponding scatra elements
        std::map<int, std::set<int>>&
            myescapedpartelemap  ///< map of escaped particles and corresponding scatra elements
                                 ///< including particles from other procs
        ) const;

    /// extend ghosting according to bin distribution
    virtual Teuchos::RCP<Epetra_Map> ExtendGhosting(std::map<int, std::set<int>>& binelemap,
        std::map<int, std::set<int>>& ext_bintoele_ghosting,
        Teuchos::RCP<Epetra_Map> bincolmap) const;

    /// extend ghosting according to bin distribution
    virtual Teuchos::RCP<Epetra_Map> ExtendGhosting(const Epetra_Map* initial_elecolmap,
        std::map<int, std::set<int>>& binelemap,
        std::map<int, std::set<int>>& ext_bintoele_ghosting,
        Teuchos::RCP<Epetra_Map> rowbins = Teuchos::null,
        Teuchos::RCP<Epetra_Map> bincolmap = Teuchos::null) const;

    /// extend ghosting according to bin distribution
    virtual void ExtendBinGhosting(Teuchos::RCP<Epetra_Map> rowbins, std::set<int> const& colbins,
        bool assigndegreesoffreedom = true);

    /// do standard ghosting on problem discretization
    virtual void StandardDiscretizationGhosting(Teuchos::RCP<DRT::Discretization>& discret,
        Teuchos::RCP<Epetra_Map> const& rowbins, Teuchos::RCP<Epetra_Vector>& disnp,
        Teuchos::RCP<Epetra_Map>& stdelecolmap, Teuchos::RCP<Epetra_Map>& stdnodecolmap) const;

    virtual void ExtendGhosting(std::vector<Teuchos::RCP<DRT::Discretization>> dis) const;

    /// extend ghosting according to bin distribution
    virtual void ExtendEleGhosting(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<Epetra_Map> initial_elecolmap, Teuchos::RCP<Epetra_Map> bincolmap,
        bool assigndegreesoffreedom, bool initelements, bool doboundaryconditions) const;

    /// collect information of content in bins from other procs via Round Robin loop
    virtual void CollectInformation(Teuchos::RCP<Epetra_Map> rowbins,
        std::map<int, std::vector<int>>& nodesinbin,
        std::map<int, std::vector<int>>& nodesinmybins) const;

    /// convert i,j,k into bin id
    virtual int ConvertijkToGid(int* ijk  ///< i,j,k to be converted into a bin id
        ) const;

    /// get all bin ids for given range of ijk
    virtual void GidsInijkRange(const int* ijk_bound,  ///< given range of ijk
        std::set<int>& binIds,                         ///< all bin ids in the specified range
        bool checkexistence  ///< check can be added whether the gids are on this proc
        ) const;

    /// get all bin ids for given range of ijk
    virtual void GidsInijkRange(const int* ijk_bound,  ///< given range of ijk
        std::vector<int>& binIds,                      ///< all bin ids in the specified range
        bool checkexistence  ///< check can be added whether the gids are on this proc
        ) const;

    /// convert bin id into i,j,k
    virtual void ConvertGidToijk(int gid,  ///< bin id to be converted into i,j,k
        int* ijk                           ///< resulting i,j,k
        ) const;

    /// convert position into bin id
    virtual int ConvertPosToGid(
        const double* pos  ///< position of which the corresponding bin id is asked for
        ) const;

    /// convert position into i,j,k
    virtual void ConvertPosToijk(const double* pos,  ///< position to be converted into i,j,k
        int* ijk                                     ///< resulting i,j,k
        ) const;

    /// convert position into i,j,k
    virtual void ConvertPosToijk(
        const LINALG::Matrix<3, 1>& pos,  ///< position to be converted into i,j,k
        int* ijk                          ///< resulting i,j,k
        ) const;

    /// convert position into bin id
    virtual int ConvertPosToGid(const LINALG::Matrix<3, 1>&
            pos  ///< position of which the corresponding bin id is asked for
        ) const;

    /// get 26 neighboring bin ids to binId (if existing)
    virtual void GetNeighborBinIds(const int binId,  ///< bin id whose connectivity is asked for
        std::vector<int>& binIds                     ///< all neighboring bins on axes
        ) const;

    /// get 27 neighboring bin ids to binId and myself
    virtual void GetNeighborAndOwnBinIds(
        const int binId,          ///< bin id whose connectivity is asked for
        std::vector<int>& binIds  ///< all neighboring bins on axes
        ) const;

    /// get nodal coordinates of bin with given id
    virtual void GetBinCorners(const int binId,        ///< bin id of which corners are calculated
        std::vector<LINALG::Matrix<3, 1>>& bincorners  ///< corner position of given bin
        ) const;

    /// centroid position for given bin id is returned
    virtual LINALG::Matrix<3, 1> GetBinCentroid(
        const int binId  ///< bin id of which centroid is calculated
        ) const;

    /// get minimum bin size
    virtual double GetMinBinSize() const;

    /// get maximum bin size
    virtual double GetMaxBinSize() const;

    virtual void RevertExtendedGhosting(std::vector<Teuchos::RCP<DRT::Discretization>> dis,
        std::vector<Teuchos::RCP<Epetra_Map>>& stdelecolmap,
        std::vector<Teuchos::RCP<Epetra_Map>>& stdnodecolmap) const;

    // create XAABB for discretizations and compute cutoff with largest element in discret
    virtual void ComputeMinXAABBContainingAllElementsOfInputDiscrets(
        std::vector<Teuchos::RCP<DRT::Discretization>> discret,
        std::vector<Teuchos::RCP<Epetra_Vector>> disnp, LINALG::Matrix<3, 2>& XAABB,
        bool setmincutoff);

    // compute max cutoff as largest element in discret
    virtual double ComputeMinCutoffAsMaxEdgeLengthOfXAABBOfLargestEle(
        std::vector<Teuchos::RCP<DRT::Discretization>> discret,
        std::vector<Teuchos::RCP<Epetra_Vector>> disnp);

    /// bins are created based on XAABB and cutoff radius
    virtual void CreateBinsBasedOnCutoffAndXAABB(
        Teuchos::RCP<DRT::Discretization> dis =
            Teuchos::null  ///< underlying discret which can be used to compute bounding box
    );

    /// bins are created based on XAABB and cutoff radius
    virtual void CreateBins2D();

    /// bins are created based on XAABB and cutoff radius
    virtual void CreateBinsScatra(Teuchos::RCP<DRT::Discretization>
            dis  ///< underlying discret which can be used to compute bounding box
    );

    virtual void WriteBinOutput(int const step, double const time);

    /// shift current particle positions in case periodic boundary conditions are applied
    virtual void PeriodicBoundaryShift3D(
        LINALG::Matrix<3, 1>& d, LINALG::Matrix<3, 1> const X = LINALG::Matrix<3, 1>(true)) const;

    // create XAABB for discretization
    virtual void CreateXAABB(Teuchos::RCP<DRT::Discretization> dis);

    // create XAABB for discretization and compute cutoff with largest element in discret
    virtual void CreateXAABB(Teuchos::RCP<DRT::Discretization> discret,
        Teuchos::RCP<Epetra_Vector> disnp, LINALG::Matrix<3, 2>& XAABB, bool setcutoff = false);

    /// build periodic boundary conditions
    virtual void BuildPeriodicBC();

    /// create linear map with bin ids
    Teuchos::RCP<Epetra_Map> CreateLinearMapForNumbin(const Epetra_Comm& comm) const;

    ///
    void BinSizeSafetyCheck(
        double half_interaction_distance, Teuchos::RCP<const Epetra_Vector> dis_increment) const;

   private:
    /// get current position of node
    void GetCurrentNodePos(Teuchos::RCP<const DRT::Discretization> const discret,
        DRT::Node const* node, Teuchos::RCP<const Epetra_Vector> const disnp,
        double* currpos) const;

    /// communicate distribution of transferred elements to bins
    void CommunicateDistributionOfTransferredElementsToBins(
        Teuchos::RCP<DRT::Discretization>& discret,
        std::map<int, std::vector<std::pair<int, std::vector<int>>>> const& toranktosendbinids,
        std::map<int, std::set<int>>& bintorowelemap) const;

   private:
    Teuchos::RCP<DRT::Discretization> bindis_;  ///< particle discretization with bins

    Teuchos::RCP<DRT::Discretization> visbindis_;  ///< visualization discretization for bins

    std::list<DRT::Element*> boundaryrowbins_;  ///< list of boundary row bins

    std::set<int> boundarycolbins_;  ///< list of boundary col bins (one full layer)

    double cutoff_radius_;  ///< maximum influence radius

    LINALG::Matrix<3, 2> XAABB_;  ///< extendedAxisAlignedBoundingBox of bins

    //! If simulation domain is deforming, e.g. under shear
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> deforming_simulation_domain_handler;

    int writebinstype_;  ///< type of bin visualization

    double bin_size_[3];  ///< size of each bin in Cartesian coordinates

    double inv_bin_size_[3];  ///< inverse of size of each bin in Cartesian coordinates

    int bin_per_dir_[3];  ///< number of bins per direction

    int id_calc_bin_per_dir_[3];  ///< number of bins per direction for bin id calculation

    int id_calc_exp_bin_per_dir_[3];  ///< exponent 2^x = number of bins per direction for bin id
                                      ///< calculation

    bool havepbc_;  ///< flag whether periodic boundary conditions are specified for particles

    bool pbconoff_[3];  ///< flags for existence of pbcs in x, y, z direction

    double pbcdeltas_[3];  ///< deltas for pbcs in x, y, z direction

    const INPAR::PARTICLEOLD::ParticleDim particle_dim_;  ///< dimension of problem

    //! @name miscellaneous stuff
    //@{

    int myrank_;  ///< myrank

    //@}
  };

  /*!
  \brief Class for comparing Teuchos::RCP<DRT::Node> in a std::set

  @tparam ELEMENT
  @param first
  @param second
  */
  class Less
  {
   public:
    template <typename ELEMENT>
    bool operator()(const Teuchos::RCP<ELEMENT>& first, const Teuchos::RCP<ELEMENT>& second) const
    {
      return first->Id() < second->Id();
    }
  };

}  // namespace BINSTRATEGY


/*----------------------------------------------------------------------*/
#endif
