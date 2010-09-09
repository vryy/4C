
#ifdef STKADAPTIVE

#include "stk_discret_state.H"
#include "stk_discret.H"
#include "../stk_refine/stk_mesh.H"

#include <Isorropia_EpetraOrderer.hpp>
#include <Epetra_CrsMatrix.h>
#include <Epetra_FECrsMatrix.h>
#include <EpetraExt_Permutation.h>
#include <Epetra_Import.h>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_condition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
STK::DiscretizationState::DiscretizationState( STK::Discretization& dis )
  : dis_( dis )
{
  CreateEpetraMaps();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::DiscretizationState::Setup( const std::vector<stk::mesh::FieldBase*>& fields )
{
  dofset_.Setup( dis_, fields );

  // build dirichlet map

  STK::Mesh & mesh = dis_.GetMesh();
  stk::mesh::MetaData & meta = mesh.MetaData();
  stk::mesh::BulkData & bulk = mesh.BulkData();

  const Epetra_Map & dofrowmap = DofRowMap();
  const Epetra_Map & dofcolmap = DofColMap();
  Epetra_IntVector flag( dofcolmap );

  // Loop dirichlet conditions from higher to lower dimension since the
  // lower dimensions overwrite higher ones.

  static DRT::Condition::ConditionType dirichlet_types[] = {
    DRT::Condition::VolumeDirichlet,
    DRT::Condition::SurfaceDirichlet,
    DRT::Condition::LineDirichlet,
    DRT::Condition::PointDirichlet,
    DRT::Condition::none
  };

  for ( DRT::Condition::ConditionType * condtype = dirichlet_types;
        *condtype!=DRT::Condition::none;
        ++condtype )
  {
    const std::map<int, Teuchos::ParameterList> * dirichlet = dis_.Condition( *condtype );
    if ( dirichlet!=NULL )
    {
      for ( std::map<int, Teuchos::ParameterList>::const_iterator i=dirichlet->begin();
            i!=dirichlet->end();
            ++i )
      {
        //int id = i->first;
        const Teuchos::ParameterList & condflags = i->second;

        std::string name = condflags.get<std::string>( "Name" );
        const std::vector<int> & onoff = * condflags.get<Teuchos::RCP<std::vector<int> > >( "onoff" );

        stk::mesh::Part * condpart = meta.get_part( name );
        if ( condpart==NULL )
          dserror( "No condition part '%s'", name.c_str() );

        stk::mesh::Selector selector = *condpart & mesh.ActivePart();

        const std::vector<stk::mesh::Bucket*> & nodes = bulk.buckets( stk::mesh::Node );
        for ( std::vector<stk::mesh::Bucket*>::const_iterator j=nodes.begin();
              j!=nodes.end();
              ++j )
        {
          stk::mesh::Bucket & bucket = **j;
          if ( selector( bucket ) )
          {
            for ( stk::mesh::Bucket::iterator k=bucket.begin(); k!=bucket.end(); ++k )
            {
              stk::mesh::Entity & n = *k;
              std::vector<int> dofids;
              Dof( n.key(), dofids );
              if ( dofids.size() > onoff.size() )
                dserror( "too few Dirichlet flags" );
#if 0
              if ( condflags.get<DRT::Condition::ConditionType>( "Type" )==DRT::Condition::PointDirichlet )
              {
                std::copy( dofids.begin(), dofids.end(), std::ostream_iterator<int>( std::cout, " " ) );
                std::cout << " - ";
                std::copy( onoff.begin(), onoff.end(), std::ostream_iterator<int>( std::cout, " " ) );
                std::cout << "\n";
              }
#endif
              for ( unsigned k=0; k<dofids.size(); ++k )
              {
                int lid = dofcolmap.LID( dofids[k] );
                if ( lid < 0 )
                  dserror( "illegal lid" );
                flag[lid] = onoff[k];
              }
            }
          }
        }
      }
    }
  }

  std::vector<int> rowdbc;
  std::vector<int> coldbc;
  int num = dofcolmap.NumMyElements();
  for ( int i=0; i<num; ++i )
  {
    if ( flag[i]!=0 )
    {
      int gid = dofcolmap.GID( i );
      coldbc.push_back( gid );
      if ( dofrowmap.MyGID( gid ) )
      {
        rowdbc.push_back( gid );
      }
    }
  }

  rowdbcmap_ = Teuchos::rcp( new Epetra_Map( -1, rowdbc.size(), &rowdbc[0], 0, dis_.Comm() ) );
  coldbcmap_ = Teuchos::rcp( new Epetra_Map( -1, coldbc.size(), &coldbc[0], 0, dis_.Comm() ) );
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::DiscretizationState::EvaluateDirichlet( double time,
                                                  const std::vector<stk::mesh::FieldBase*> * v,
                                                  const std::vector<stk::mesh::FieldBase*> * dv,
                                                  const std::vector<stk::mesh::FieldBase*> * ddv )
{
  STK::Mesh & mesh = dis_.GetMesh();
  stk::mesh::MetaData & meta = mesh.MetaData();
  stk::mesh::BulkData & bulk = mesh.BulkData();

  //const Epetra_Map & dofrowmap = DofRowMap();

  unsigned deg = 0;  // highest degree of requested time derivative
  if ( dv!=NULL )
  {
    if ( v->size()!=dv->size() )
      dserror( "existing fields must match" );
    if ( ddv!=NULL )
    {
      if ( v->size()!=ddv->size() )
        dserror( "existing fields must match" );
      deg = 2;
    }
    else
    {
      deg = 1;
    }
  }

  static DRT::Condition::ConditionType dirichlet_types[] = {
    DRT::Condition::PointDirichlet,
    DRT::Condition::LineDirichlet,
    DRT::Condition::SurfaceDirichlet,
    DRT::Condition::VolumeDirichlet,
    DRT::Condition::none
  };

  stk::mesh::PartVector done;

  for ( DRT::Condition::ConditionType * condtype = dirichlet_types;
        *condtype!=DRT::Condition::none;
        ++condtype )
  {
    const std::map<int, Teuchos::ParameterList> * dirichlet = dis_.Condition( *condtype );
    if ( dirichlet!=NULL )
    {
      for ( std::map<int, Teuchos::ParameterList>::const_iterator i=dirichlet->begin();
            i!=dirichlet->end();
            ++i )
      {
        //int id = i->first;
        const Teuchos::ParameterList & condflags = i->second;

        //if ( condflags.get<DRT::Condition::ConditionType>( "Type" )==*condtype )
        {
          // find the part, find which node is covered by which part

          std::string name = condflags.get<std::string>( "Name" );
          const std::vector<int> & onoff = * condflags.get<Teuchos::RCP<std::vector<int> > >( "onoff" );
          const std::vector<int> & curve = * condflags.get<Teuchos::RCP<std::vector<int> > >( "curve" );
          const std::vector<int> & funct = * condflags.get<Teuchos::RCP<std::vector<int> > >( "funct" );
          const std::vector<double> & val = * condflags.get<Teuchos::RCP<std::vector<double> > >( "val" );

          stk::mesh::Part * condpart = meta.get_part( name );
          if ( condpart==NULL )
            dserror( "No condition part '%s'", name.c_str() );

          // iterate all active nodes including ghosted nodes
          stk::mesh::Selector selector = *condpart & mesh.ActivePart();

          const std::vector<stk::mesh::Bucket*> & nodes = bulk.buckets( stk::mesh::Node );
          for ( std::vector<stk::mesh::Bucket*>::const_iterator j=nodes.begin();
                j!=nodes.end();
                ++j )
          {
            stk::mesh::Bucket & bucket = **j;
            if ( selector( bucket ) )
            {
              bool isdone = false;
              for ( stk::mesh::PartVector::iterator k=done.begin(); k!=done.end(); ++k )
              {
                if ( has_superset( bucket, **k ) )
                {
                  isdone = true;
                  break;
                }
              }

              if ( not isdone )
              {
                for ( stk::mesh::Bucket::iterator k=bucket.begin(); k!=bucket.end(); ++k )
                {
                  stk::mesh::Entity & n = *k;

                  unsigned ifield = 0;
                  int field_base = 0;

                  int vfield_data_size = stk::mesh::field_data_size( *( *v )[ifield], bucket ) / sizeof( double );
                  double * vdata = reinterpret_cast<double*>( field_data( *( *v )[ifield], n ) );

                  // apply Dirichlet condition to node dofs

                  std::vector<int> dofs;
                  dis_.Dof( n.key(), dofs );
                  const int numdf = dofs.size();
                  for ( int j=0; j<numdf; ++j )
                  {
                    if ( j >= field_base+vfield_data_size )
                    {
                      field_base += vfield_data_size;
                      ifield += 1;
                      if ( ifield >= v->size() )
                        dserror( "to few fields for Dirichlet evaluation" );
                      vfield_data_size = stk::mesh::field_data_size( *( *v )[ifield], bucket ) / sizeof( double );
                      vdata = reinterpret_cast<double*>( field_data( *( *v )[ifield], n ) );
                    }

                    if ( onoff[j]!=0 )
                    {
                      //const int gid = dofs[j];
                      std::vector<double> value( deg+1, val[j] );

                      // factor given by time curve
                      std::vector<double> curvefac( deg+1, 1.0 );
                      int curvenum = curve[j];
                      if ( curvenum>=0 and time>=0.0 )
                      {
                        curvefac = DRT::Problem::Instance()->Curve(curvenum).FctDer( time, deg );
                      }
                      else
                      {
                        for (unsigned i=1; i<(deg+1); ++i)
                          curvefac[i] = 0.0;
                      }

                      // factor given by spatial function
                      double functfac = 1.0;
                      int funct_num = funct[j];
                      if (funct_num>0)
                      {
                        stk::mesh::VectorField & coords = dis_.GetMesh().Coordinates();
                        double * x = stk::mesh::field_data( coords , n );
                        functfac = DRT::Problem::Instance()->Funct(funct_num-1).Evaluate(j,x,time,NULL);
                      }

                      // apply factors to Dirichlet value
                      for (unsigned i=0; i<deg+1; ++i)
                      {
                        value[i] *= functfac * curvefac[i];
                      }

                      vdata[j - field_base] = value[0];
                      if ( deg > 0 )
                      {
                        double * dvdata = reinterpret_cast<double*>( field_data( *( *v )[ifield], n ) );
                        dvdata[j - field_base] = value[1];
                        if ( deg > 1 )
                        {
                          double * ddvdata = reinterpret_cast<double*>( field_data( *( *v )[ifield], n ) );
                          ddvdata[j - field_base] = value[2];
                        }
                      }
                    }
                  }
                }
              }
            }
          }

          done.push_back( condpart );
        }
      }
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::DiscretizationState::CreateEpetraMaps()
{
  {
    std::vector<int> vnids;
    stk::mesh::Selector selector = dis_.GetMesh().OwnedPart() & dis_.GetMesh().ActivePart();
    CreateEntityMap( selector, stk::mesh::Node, vnids );
    rownodemap_ = Teuchos::rcp( new Epetra_Map( -1,vnids.size(),&vnids[0],0,dis_.Comm() ) );
  }

  {
    std::vector<int> vnids;
    stk::mesh::Selector selector = dis_.GetMesh().ActivePart();
    CreateEntityMap( selector, stk::mesh::Node, vnids );
    colnodemap_ = Teuchos::rcp( new Epetra_Map( -1,vnids.size(),&vnids[0],0,dis_.Comm() ) );
  }

#if 0
  {
    std::vector<int> vnids;
    stk::mesh::Selector selector = dis_.GetMesh().OwnedPart() & dis_.GetMesh().ActivePart();
    CreateEntityMap( selector, stk::mesh::Element, vnids );
    rowelementmap_ = Teuchos::rcp( new Epetra_Map( -1,vnids.size(),&vnids[0],0,dis_.Comm() ) );
  }

  {
    std::vector<int> vnids;
    stk::mesh::Selector selector = dis_.GetMesh().ActivePart();
    CreateEntityMap( selector, stk::mesh::Element, vnids );
    colelementmap_ = Teuchos::rcp( new Epetra_Map( -1,vnids.size(),&vnids[0],0,dis_.Comm() ) );
  }
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::DiscretizationState::CreateEntityMap( const stk::mesh::Selector & selector, stk::mesh::EntityRank rank, std::vector<int> & vnids )
{
  // Use a set to gather nodal ids. This changes the order.
  std::set<int> nids;

  const std::vector<stk::mesh::Bucket*> & buckets = dis_.GetMesh().BulkData().buckets( rank );

  // loop buckets and entries
  for (std::vector<stk::mesh::Bucket*>::const_iterator i=buckets.begin();
       i!=buckets.end();
       ++i)
  {
    stk::mesh::Bucket & bucket = **i;
    if ( selector( bucket ) )
    {
      for (stk::mesh::Bucket::iterator ib=bucket.begin();
           ib!=bucket.end();
           ++ib)
      {
        stk::mesh::Entity & e = *ib;
        int id = static_cast<int>( e.identifier() );
        nids.insert(id);
      }
    }
  }

  vnids.reserve(nids.size());
  vnids.assign(nids.begin(),nids.end());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STK::DiscretizationState::LocationVector( stk::mesh::Entity & e,
                                               std::vector<int> & lm,
                                               std::vector<int> & lmowner )
{
  //stk::mesh::Bucket & bucket = e.bucket();
  //std::vector<stk::mesh::FieldBase*> fields;
  //algo.collect_unknowns( fields );

  lm.clear();
  lmowner.clear();

  for ( stk::mesh::PairIterRelation nodes = e.relations( stk::mesh::Node );
        not nodes.empty();
        ++nodes )
  {
    stk::mesh::Entity & n = * nodes->entity();

    std::vector<int> dof;
    Dof( n.key(), dof );
    std::copy( dof.begin(), dof.end(), std::back_inserter( lm ) );
    int owner = n.owner_rank();
    lmowner.resize( lm.size(), owner );
  }
}


#endif
