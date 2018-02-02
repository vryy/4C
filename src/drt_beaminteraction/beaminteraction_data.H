/*-----------------------------------------------------------*/
/*!
\file beaminteraction_data.H

\brief small data container for beam interaction framework

\maintainer Jonas Eichinger

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BEAMINTERACTION_DATA_H_
#define BEAMINTERACTION_DATA_H_

#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_beaminteraction.H"
#include <map>
#include <set>


// forward declarations
namespace DRT {
  class ParObject;
  class PackBuffer;
}
// forward declaration
namespace STR {
namespace TIMINT {
  class BaseDataGlobalState;
}
}
namespace BEAMINTERACTION {
  /*!
   * data container for input file parameters for submodel crosslinking in beam interaction
   * author eichinger*/
  class BeamInteractionParams
  {
  public:
    //! constructor
    BeamInteractionParams();

    //! destructor
    virtual ~BeamInteractionParams() {};

    //! initialize with the stuff coming from input file
    void Init();

    //! setup member variables
    void Setup();

    //! returns the isinit_ flag
    inline const bool& IsInit() const
    { return isinit_; };

    //! returns the issetup_ flag
    inline const bool& IsSetup() const
    { return issetup_; };

    //! Checks the init and setup status
    inline void CheckInitSetup() const
    {
      if (!IsInit() or !IsSetup())
        dserror("Call Init() and Setup() first!");
    }

    //! Checks the init status
    inline void CheckInit() const
    {
      if (!IsInit())
        dserror("Init() has not been called, yet!");
    }

    /// number of crosslinkers per type
    INPAR::BEAMINTERACTION::RepartitionStrategy GetRepartitionStrategy() const
    {
      CheckInitSetup();
      return rep_strategy_;
    };


  private:

    bool isinit_;

    bool issetup_;

    /// number of crosslinker that are initially set
    INPAR::BEAMINTERACTION::RepartitionStrategy rep_strategy_;

  };



namespace DATA
{

  //! struct to store crosslinker data and enable parallel redistribution
  struct CrosslinkerData
  {

  public:

    //! constructor
    CrosslinkerData();

    //! destructor
    virtual ~CrosslinkerData() {};


    //!@name data access functions
    //! @{

    int GetId() const { return id_; };

    LINALG::Matrix< 3, 1 > const & GetPosition() const { return pos_; };

    int GetNumberOfBonds() const { return numbond_; };

    std::vector< std::pair< int, int > > const & GetBSpots() const { return bspots_; };

    void SetId( int id ) { id_ = id; };

    void SetPosition( LINALG::Matrix< 3, 1 > const & clpos) { pos_ = clpos; };

    void SetNumberOfBonds( int clnumbond ) { numbond_ = clnumbond; };

    void SetBSpots( std::vector< std::pair< int, int > > const & clbspots ) { bspots_ = clbspots; };
    void SetBspot ( int bspotid, std::pair< int, int > bpsotone ) { bspots_[bspotid] = bpsotone; };

    //! @}

    //!@name methods for enabling parallel use of data container
    //! @{

    /*!
    \brief Pack this class so it can be communicated
    */
    void Pack( DRT::PackBuffer & data ) const;

    /*!
    \brief Unpack data from a char vector into this container
    */
    void Unpack( std::vector<char> const & data );

    //! @}

  private:

    //!@name linker member variables
    //! @{

    /// linker gid
    int id_;

    /// current position of crosslinker
    LINALG::Matrix< 3, 1 > pos_;

    /// number of active bonds
    int numbond_;

    /// gid of element to local number of bspot, [0] and [1] first and second bspot
    std::vector< std::pair< int, int > > bspots_;


    //! @}
  };

  //! struct to store crosslinker data and enable parallel redistribution
  struct BeamData
  {

  public:

    //! constructor
    BeamData();

    //! destructor
    virtual ~BeamData() {};


    //!@name data access functions
    //! @{

    int GetId() const { return id_; };

    std::map< int, LINALG::Matrix< 3, 1 > > const & GetBSpotPositions() const { return bspotpos_; };
    LINALG::Matrix< 3, 1 > const & GetBSpotPosition( int bspotid ) const { return bspotpos_.at(bspotid); };

    std::map< int, LINALG::Matrix< 3, 3 > > const & GetBSpotTriads() const { return bspottriad_; };
    LINALG::Matrix< 3, 3 > const & GetBSpotTriad( int bspotid ) const { return bspottriad_.at(bspotid); };

    std::map<int, std::set<int> > const & GetBSpotStatus() const { return bspotstatus_; };

    std::set<int> const & GetBSpotStatusAt( int bspotid ) const { return bspotstatus_.at(bspotid); };

    int GetNumberOfBindingSpots() const { return bspotstatus_.size();};


    void SetId( int id ) { id_ = id; };

    void SetBSpotPositions( std::map< int, LINALG::Matrix< 3, 1 > > const & bspotpos ) { bspotpos_ = bspotpos; };
    void SetBSpotPosition( int bspotid, LINALG::Matrix< 3, 1 > const & bspotpos ) { bspotpos_[bspotid] = bspotpos; };

    void SetBSpotTriads( std::map< int, LINALG::Matrix< 3, 3 > > const & bspottriad ) { bspottriad_ = bspottriad; };
    void SetBSpotTriad( int bspotid, LINALG::Matrix< 3, 3 > const & bspottriad ) { bspottriad_[bspotid] = bspottriad; };

    void SetBSpotStatus( std::map<int, std::set<int> > const & bspotstatus ) { bspotstatus_ = bspotstatus; };
    void SetBSpotStatus( int bspotid, std::set<int> clgids ) { bspotstatus_[bspotid] = clgids; };

    void EraseBondFromBindingSpot( int locbspotid, int clgid )
    { bspotstatus_[locbspotid].erase( clgid ); }
    void AddBondToBindingSpot( int locbspotid, int clgid )
    { bspotstatus_[locbspotid].insert( clgid ); }


    //! @}

    //!@name methods for enabling parallel use of data container
    //! @{

    /*!
    \brief Pack this class so it can be communicated
    */
    void Pack( DRT::PackBuffer & data ) const;

    /*!
    \brief Unpack data from a char vector into this container
    */
    void Unpack( std::vector<char> const & data );

    //! @}

  private:

    //!@name linker member variables
    //! @{

    /// beam gid
    int id_;

    /// current position at bindingspots (xi) (key is local number of binding spot)
    std::map< int, LINALG::Matrix< 3, 1 > > bspotpos_;

    /// current triad at bindingspots (xi) (key is local number of binding spot)
    std::map< int, LINALG::Matrix< 3, 3 > > bspottriad_;

    /// key is locn of bspot, holds gid of crosslinker to which it is bonded
    std::map< int, std::set<int> > bspotstatus_;

    //! @}
  };

  struct BindEventData
  {

    //! constructor
    BindEventData();

    //! destructor
    virtual ~BindEventData() {};

    //! Init
    void Init( int clgid, int elegid, int bspotlocn, int requestproc, int permission );


    //!@name data access functions
    //! @{

    int GetClId() const{ return clgid_; };

    int GetEleId() const { return elegid_; };

    int GetBSpotLocN() const { return bspotlocn_; };

    int GetRequestProc() const { return requestproc_; };

    int GetPermission() const { return permission_; };

    void SetClId( int clgid ) { clgid_ = clgid; };

    void SetEleId( int elegid) { elegid_ = elegid; };

    void SetBSpotLocN( int bspotlocn) { bspotlocn_ = bspotlocn; };

    void SetRequestProc( int requestproc) { requestproc_ = requestproc; };

    void SetPermission( int permission) { permission_ = permission; };

    //! @}

    //!@name methods for enabling parallel use of data container
    //! @{

    /*!
    \brief Pack this class so it can be communicated
    */
    void Pack( DRT::PackBuffer & data ) const;

    /*!
    \brief Unpack data from a char vector into this container
    */
    void Unpack( std::vector<char> const & data );

    //! @}

  private:

    //!@name member variables
    //! @{

    // gid of crosslinker
    int clgid_;

    // ele gid crosslinker wants to bind to
    int elegid_;

    // loc number of bspot on ele cl wants to bind to
    int bspotlocn_;

    // myrank, processor that is requesting
    int requestproc_;

    // permission/veto, if crosslinker is allowed to bind
    int permission_;

    //! @}

  };

  struct UnBindEventData
  {

    //! constructor
    UnBindEventData();

    //! destructor
    virtual ~UnBindEventData() {};

    //!@name data access functions
    //! @{

    int GetClId() const{ return clgid_; };

    std::pair< int, int > const & GetEleToUpdate() const { return eletoupdate_; };


    void SetClId( int clgid ) { clgid_ = clgid; };

    void SetEleToUpdate( std::pair< int, int > eletoupdate ) { eletoupdate_ = eletoupdate; };


    //! @}


    //!@name methods for enabling parallel use of data container
    //! @{

    /*!
    \brief Pack this class so it can be communicated
    */
    void Pack( DRT::PackBuffer & data ) const;

    /*!
    \brief Unpack data from a char vector into this container
    */
    void Unpack( std::vector<char> const & data );

    //! @}

  private:

    //!@name member variables
    //! @{

    // crosslinker (gid) that is no longer bonded to above binding spot
    int clgid_;

    // element gid (first) that needs to be updated at local binding (second)
    std::pair< int, int > eletoupdate_;


    //! @}

  };

  /// create respective data container from char vector
  template< typename T >
  T * CreateDataContainer( std::vector<char> const & data )
  {
    T * new_container = new T();
    new_container->Unpack(data);
    return new_container;
  };

}
}

#endif
