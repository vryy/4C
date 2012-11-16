/*!----------------------------------------------------------------------
\file drt_turbulent_inflow.H

\brief Methods to transfer a turbulent inflow profile a from a master
boundary to a slave boundary on a seperate domain.
The slave boundary must have an additional Dirichlet condition, the
master boundary will usually be a periodic boundary (but is not required
to)

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef TRANSFERTURBINF_H
#define TRANSFERTURBINF_H
#include "Teuchos_RCP.hpp"
#include "../linalg/linalg_utils.H"

/*

*/

//forward declarations
namespace DRT
{
  class Discretization;
}

namespace FLD
{

class TransferTurbulentInflowCondition
{
public:
  /*!
  \brief Standard Constructor

  */
  TransferTurbulentInflowCondition(
    RCP<DRT::Discretization>  dis,
    RCP<LINALG::MapExtractor> dbcmaps
    );

  /*!
  \brief Destructor

  */
  virtual ~TransferTurbulentInflowCondition();

  /*!
  \brief Transfer process copying values from master boundary to slave
         boundary (slave must be of Dirichlet type, otherwise this
         operation doesn't make to mauch sense)

         Intended to be called after ApplyDirichlet, overwriting the
         dummy Dirichlet values on the slave boundary by the values
         of the last time step on the master boundary
  */
  void Transfer(const Teuchos::RCP<Epetra_Vector> veln ,
                Teuchos::RCP<Epetra_Vector>       velnp,
                const double                      time);

private:

  //! there are two types of tranfer conditions. values are transferred
  //! from master to slave conditions
  enum ToggleType
  {
    none,
    master,
    slave
  };

  //! get all data from condtion
  void GetData(int & id,int & direction,ToggleType & type,const DRT::Condition*);

#ifdef PARALLEL
  //! receive a block in the round robin communication pattern
  void ReceiveBlock(
    vector<char>  & rblock  ,
    DRT::Exporter & exporter,
    MPI_Request   & request );

  //! send a block in the round robin communication pattern
  void SendBlock(
    vector<char>  & sblock  ,
    DRT::Exporter & exporter,
    MPI_Request   & request );
#endif

  //! unpack all master values contained in receive block
  void UnpackLocalMasterValues(
    vector<int>             & mymasters    ,
    vector<std::vector<double> > & mymasters_vel,
    vector<char>            & rblock
    );

  //! pack all master values into a send block
  void PackLocalMasterValues(
    vector<int>             & mymasters    ,
    vector<std::vector<double> > & mymasters_vel,
    DRT::PackBuffer         & sblock
    );

  //! for all values avaible on the processor, do the final setting of the value
  void SetValuesAvailableOnThisProc(
    vector<int>                 & mymasters,
    vector<std::vector<double> >     & mymasters_vel,
    Teuchos::RCP<Epetra_Vector>   velnp);

  //! flag active boundary condition (may be used to switch off everything)
  bool                             active_;

  //! the discretisation
  RCP<DRT::Discretization> dis_;

  //! info on DIirchlet boundary
  RCP<LINALG::MapExtractor> dbcmaps_;

  //! the connectivity of the boundary condition
  std::map<int,std::vector<int> >       midtosid_;

  //! time curve number
  int curve_;

};


} // end namespace FLD

#endif  // #ifndef TRANSFERTURBINF_H
