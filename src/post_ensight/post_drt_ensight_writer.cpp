/*!
  \file post_drt_ensight_writer.cpp

  \brief Ensight filter basis class

  \level 1

  \maintainer Martin Kronbichler
*/



#include "post_drt_ensight_writer.H"
#include "../post_drt_common/post_drt_common.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_fluid/fluid_rotsym_periodicbc_utils.H"
#include "../linalg/linalg_utils.H"

#include <string>

#include "../pss_full/pss_cpp.h"
extern "C" {
#include "../pss_full/pss_table_iter.h"
}

//! 6 Surfaces of a Hex27 element with 9 nodes per surface
const int Hex20_BaciToEnsightGold[20] =
        { 0,  1,  2,  3,
          4,  5,  6,  7,
          8,  9, 10, 11,
         16, 17, 18, 19,
         12, 13, 14, 15};

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EnsightWriter::EnsightWriter(PostField* field,
                             const std::string &filename)
:
    PostWriterBase(field, filename),
    nodeidgiven_(true),
    writecp_(false)
{
  // initialize proc0map_ correctly
  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();
  const Epetra_Map* noderowmap = dis->NodeRowMap();
  proc0map_ = LINALG::AllreduceEMap(*noderowmap,0);

  // sort proc0map_ so that we can loop it and get nodes in ascending order.
  std::vector<int> sortmap;
  sortmap.reserve(proc0map_->NumMyElements());
  sortmap.assign(proc0map_->MyGlobalElements(), proc0map_->MyGlobalElements()+proc0map_->NumMyElements());
  std::sort(sortmap.begin(), sortmap.end());
  proc0map_ = Teuchos::rcp(new Epetra_Map(-1, sortmap.size(), &sortmap[0], 0, proc0map_->Comm()));

  // get the number of elements for each distype (global numbers)
  numElePerDisType_ = GetNumElePerDisType(dis);

  // get the global ids of elements for each distype (global numbers)
  eleGidPerDisType_ = GetEleGidPerDisType(dis, numElePerDisType_);

  // map between distypes in BACI and existing Ensight strings
  // it includes only strings for cell types known in ensight
  // you need to manually switch to other types distypes before querying this map
  distype2ensightstring_.clear();
  distype2ensightstring_[DRT::Element::point1] = "point";
  distype2ensightstring_[DRT::Element::line2] = "bar2";
  distype2ensightstring_[DRT::Element::line3] = "bar2"; //"bar3";
  distype2ensightstring_[DRT::Element::hex8] = "hexa8";
  distype2ensightstring_[DRT::Element::hex20] = "hexa20";
  distype2ensightstring_[DRT::Element::tet4] = "tetra4";
  distype2ensightstring_[DRT::Element::tet10] = "tetra10";
  distype2ensightstring_[DRT::Element::nurbs8] = "hexa8";
  distype2ensightstring_[DRT::Element::nurbs27] = "hexa8";
  distype2ensightstring_[DRT::Element::nurbs4] = "quad4";
  distype2ensightstring_[DRT::Element::nurbs9] = "quad4";
  distype2ensightstring_[DRT::Element::quad4] = "quad4";
  distype2ensightstring_[DRT::Element::quad8] = "quad8";
  distype2ensightstring_[DRT::Element::tri3] = "tria3";
  distype2ensightstring_[DRT::Element::tri6] = "tria6";
  distype2ensightstring_[DRT::Element::wedge6] = "penta6";
  distype2ensightstring_[DRT::Element::wedge15]= "penta15";
  distype2ensightstring_[DRT::Element::pyramid5]= "pyramid5";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteFiles(PostFilterBase &filter)
{
  PostResult result(field_);

  // timesteps when the solution is written
  const std::vector<double> soltime = result.get_result_times(field_->name());
  const unsigned int numsoltimes = soltime.size();

  // When spatial approximation is based on NURBS, we want to write
  // real geometry data and control point information. Therefore,
  // we perform one standard writing step and one additional step
  // for the control points. Here, a new .case file is created which
  // ends with "_cp".
  int iter = 1;
  if(field_->problem()->SpatialApproximation()=="Nurbs")
    iter++;

  // For none-NURBS cases, this loop is just passed through once!
  for(int i=0;i<iter;++i)
  {
    // for control point output:
    // define aux string for control point (cp) files and set writecp_ to true
    // if necessary
    std::string aux = "";
    writecp_ = false;
    if(i==1)
    {
      writecp_ = true;
      aux = aux + "_cp";

      // should be cleared for control point output
      filesetmap_.clear();
    }

    ///////////////////////////////////
    //  write geometry file          //
    ///////////////////////////////////
    const std::string geofilename = filename_ + "_"+ field_->name() + aux + ".geo";
    const size_t found_path = geofilename.find_last_of("/\\");
    const std::string geofilename_nopath = geofilename.substr(found_path+1);
    WriteGeoFile(geofilename);
    std::vector<int> filesteps;
    filesteps.push_back(1);
    filesetmap_["geo"] = filesteps;
    std::vector<double> timesteps;
    if (soltime.size()>0)
      timesteps.push_back(soltime[0]);
    else
      timesteps.push_back(0.0);
    timesetmap_["geo"] = timesteps;
    // at the moment, we can only print out the first step -> to be changed

    ///////////////////////////////////
    //  write solution fields files  //
    ///////////////////////////////////
    filter.WriteAllResults(field_);

    // prepare the time sets and file sets for case file creation
    int setcounter = 0;
    int allresulttimeset = 0;
    for (std::map<std::string,std::vector<double> >::const_iterator entry = timesetmap_.begin(); entry != timesetmap_.end(); ++entry)
    {
      std::string key = entry->first;
      if ((entry->second).size()== numsoltimes)
      {
        if (allresulttimeset == 0)
        {
          setcounter++;
          allresulttimeset = setcounter;
        }
        timesetnumbermap_[key] = allresulttimeset; // reuse the default result time set, when possible
      }
      else
      {
        setcounter++;
        timesetnumbermap_[key] = setcounter; // a new time set number is needed
      }
    }

    // Paraview wants the geo file to be fileset number one
    setcounter = 1;
    for (std::map<std::string,std::vector<int> >::const_iterator entry = filesetmap_.begin(); entry != filesetmap_.end(); ++entry)
    {
      std::string key = entry->first;
      if (entry->first=="geo")
      {
        filesetnumbermap_[key] = 1;
      }
      else
      {
        setcounter++;
        filesetnumbermap_[key] = setcounter;
      }
    }

    ///////////////////////////////////
    //  now write the case file      //
    ///////////////////////////////////
    if (myrank_ == 0)
    {
      const std::string casefilename = filename_ + "_"+ field_->name() + aux + ".case";
      std::ofstream casefile;
      casefile.open(casefilename.c_str());
      casefile << "# created using post_drt_ensight\n"<< "FORMAT\n\n"<< "type:\tensight gold\n";

      casefile << "\nGEOMETRY\n\n";
      casefile << "model:\t"<<timesetnumbermap_["geo"]<<"\t"<<filesetnumbermap_["geo"]<<"\t"<< geofilename_nopath<< "\n";

      casefile << "\nVARIABLE\n\n";
      casefile << GetVariableSection(filesetmap_, variablenumdfmap_, variablefilenamemap_);

      casefile << "\nTIME\n\n";
      casefile << GetTimeSectionStringFromTimesets(timesetmap_);

      casefile << "\nFILE\n\n";
      casefile << GetFileSectionStringFromFilesets(filesetmap_);

      casefile.close();
    }
  } // end of loop

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteFilesChangingGeom(PostFilterBase &filter)
{
  // times and steps when the solution is written
  std::vector<double> soltime;
  std::vector<int> solstep;
  {
    PostResult result = PostResult(field_);
    result.get_result_timesandsteps(field_->name(), soltime, solstep);
  }
  const unsigned int numsoltimes = soltime.size();

  // prepare geometry file
  const std::string geofilename = filename_ + "_"+ field_->name() + ".geo";
  const size_t found_path = geofilename.find_last_of("/\\");
  const std::string geofilename_nopath = geofilename.substr(found_path+1);

  // open file
  std::ofstream geofile;
  if (myrank_ == 0)
  {
    geofile.open(geofilename.c_str());
    if (!geofile)
      dserror("failed to open file: %s", geofilename.c_str());
  }

  // header for geo file
  Write(geofile, "C Binary");
  if (geofile.is_open())
    geofile.close();

  std::map<std::string, std::vector<std::ofstream::pos_type> > resultfilepos;
  PostResult result = PostResult(field_);

  // loop over number of time steps
  for(unsigned int istep=0; istep<numsoltimes; istep++)
  {

    ///////////////////////////////////
    //  read in correct mesh         //
    ///////////////////////////////////

    int fieldpos = field_->field_pos();
    std::string fieldname = field_->name();
    field_->problem()->re_read_mesh(fieldpos, fieldname, solstep[istep]);

    ///////////////////////////////////
    //  write geometry file          //
    ///////////////////////////////////

    // open file
    if (myrank_ == 0)
    {
      geofile.open(geofilename.c_str(), std::ofstream::app | std::ofstream::binary);
      if (!geofile)
        dserror("failed to open file: %s", geofilename.c_str());
    }

    // print out one time step
    {
      WriteGeoFileOneTimeStep(geofile, resultfilepos, "geo");
    }

    if (geofile.is_open())
      geofile.close();

    timesetmap_["geo"].push_back(soltime[istep]);

    ///////////////////////////////////
    //  write solution fields files  //
    ///////////////////////////////////

    result.next_result();

    bool firststep = false;
    if( istep == 0 )
      firststep = true;

    bool laststep = false;
    if(istep+1 == numsoltimes)
      laststep = true;

    filter.WriteAllResultsOneTimeStep(result, firststep, laststep);
  }

  // append index table for geo file
  if (myrank_ == 0)
  {
    geofile.open(geofilename.c_str(), std::ofstream::app | std::ofstream::binary);
    if (!geofile)
      dserror("failed to open file: %s", geofilename.c_str());
  }

  WriteIndexTable(geofile, resultfilepos["geo"]);

  if (geofile.is_open())
    geofile.close();

  std::vector<int> filesteps;
  filesteps.push_back(numsoltimes);
  filesetmap_["geo"] = filesteps;

  // prepare the time sets and file sets for case file creation
  int setcounter = 0;
  int allresulttimeset = 0;
  for (std::map<std::string,std::vector<double> >::const_iterator entry = timesetmap_.begin(); entry != timesetmap_.end(); ++entry)
  {
    std::string key = entry->first;
    if ((entry->second).size()== numsoltimes)
    {
      if (allresulttimeset == 0)
      {
        setcounter++;
        allresulttimeset = setcounter;
      }
      timesetnumbermap_[key] = allresulttimeset; // reuse the default result time set, when possible
    }
    else
    {
      setcounter++;
      timesetnumbermap_[key] = setcounter; // a new time set number is needed
    }
  }

  // Paraview wants the geo file to be fileset number one
  setcounter = 1;
  for (std::map<std::string,std::vector<int> >::const_iterator entry = filesetmap_.begin(); entry != filesetmap_.end(); ++entry)
  {
    std::string key = entry->first;
    if (entry->first=="geo")
    {
      filesetnumbermap_[key] = 1;
    }
    else
    {
      setcounter++;
      filesetnumbermap_[key] = setcounter;
    }
  }

  ///////////////////////////////////
  //  now write the case file      //
  ///////////////////////////////////
  if (myrank_ == 0)
  {
    const std::string casefilename = filename_ + "_"+ field_->name() + ".case";
    std::ofstream casefile;
    casefile.open(casefilename.c_str());
    casefile << "# created using post_drt_ensight\n"<< "FORMAT\n\n"<< "type:\tensight gold\n";

    casefile << "\nGEOMETRY\n\n";
    casefile << "model:\t"<<timesetnumbermap_["geo"]<<"\t"<<filesetnumbermap_["geo"]<<"\t"<< geofilename_nopath<< "\n";

    casefile << "\nVARIABLE\n\n";
    casefile << GetVariableSection(filesetmap_, variablenumdfmap_, variablefilenamemap_);

    casefile << "\nTIME\n\n";
    casefile << GetTimeSectionStringFromTimesets(timesetmap_);

    casefile << "\nFILE\n\n";
    casefile << GetFileSectionStringFromFilesets(filesetmap_);

    casefile.close();
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteGeoFile(const std::string& geofilename)
{
  // open file
  std::ofstream geofile;
  if (myrank_ == 0)
  {
    geofile.open(geofilename.c_str());
    if (!geofile)
      dserror("failed to open file: %s", geofilename.c_str());
  }

  // header
  Write(geofile, "C Binary");

  // print out one
  // if more are needed, this has to go into a loop
  std::map<std::string, std::vector<std::ofstream::pos_type> > resultfilepos;

  {
    WriteGeoFileOneTimeStep(geofile, resultfilepos, "geo");
  }

  // append index table
  // TODO: ens_checker complains if this is turned!!!! but I can't see, whats wrong here a.ger 11/07
  // it is also correct to ommit WriteIndexTable, however the EnsightGold Format manual says,
  // it would improve performance to have it on...
  // Writing the index for the result fields is fine. Complains only for the geometry-file  gb 02/10
  // WriteIndexTable(geofile, resultfilepos["geo"]);

  if (geofile.is_open())
    geofile.close();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteGeoFileOneTimeStep(
  std::ofstream& file,
  std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
  const std::string name)
{
  std::vector<std::ofstream::pos_type>& filepos = resultfilepos[name];
  Write(file, "BEGIN TIME STEP");
  filepos.push_back(file.tellp());

  Write(file, field_->name() + " geometry");
  Write(file, "Comment");

  //nodeidgiven_ is set to true inside the class constructor
  if (nodeidgiven_)
    Write(file,"node id given");
  else
    Write(file, "node id assign");

  Write(file, "element id off");

  // part + partnumber + comment
  Write(file, "part");
  Write(file, field_->field_pos()+1);
  Write(file, field_->name() + " field");


  //switch between nurbs an others
  if(field_->problem()->SpatialApproximation()=="Nurbs" && !writecp_)
  {
    // cast dis to NurbsDiscretisation
    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*(field_->discretization())));

    if(nurbsdis==NULL)
    {
      dserror("This probably isn't a NurbsDiscretization\n");
    }

    // get number of patches
    int npatches = (nurbsdis->GetKnotVector())->ReturnNP();

    int totalnumvisp=0;

    // loop all patches
    for(int np=0;np<npatches;++np)
    {
      // get nurbs dis' knotvector sizes
      std::vector<int> nele_x_mele_x_lele(nurbsdis->Return_nele_x_mele_x_lele(np));

      int numvisp=1;

      for(unsigned rr=0;rr<nele_x_mele_x_lele.size();++rr)
      {
        numvisp*=2*(nele_x_mele_x_lele[rr])+1;
      }
      totalnumvisp+=numvisp;
    }
    if(myrank_==0)
    {
      std::cout << "Writing coordinates for " << totalnumvisp << " visualisation points\n";
    }
    Write(file, "coordinates");
    Write(file, totalnumvisp);
  }
  else
  {
    Write(file, "coordinates");
    Write(file, field_->num_nodes());
  }

  // write the grid information
  Teuchos::RCP<Epetra_Map> proc0map = WriteCoordinates(file, field_->discretization());
  proc0map_=proc0map; // update the internal map
  WriteCells(file, field_->discretization(), proc0map);

  Write(file, "END TIME STEP");
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> EnsightWriter::WriteCoordinates(
  std::ofstream& geofile,
  const Teuchos::RCP<DRT::Discretization> dis
  )
{

  if (myrank_==0)
  {
    std::cout << "(computing) coordinates for a ";
    std::cout << field_->problem()->SpatialApproximation();
    std::cout << " approximation\n";
  }

  // map for all visualisation points after they have been
  // communicated to proc 0
  Teuchos::RCP<Epetra_Map> proc0map;

  if(field_->problem()->SpatialApproximation()=="Polynomial"
     or field_->problem()->SpatialApproximation()=="Meshfree"
     or field_->problem()->SpatialApproximation()=="HDG")
  {
    WriteCoordinatesForPolynomialShapefunctions(
      geofile,dis,proc0map);
  }
  else if(field_->problem()->SpatialApproximation()=="Nurbs")
  {
    // write real geometry coordinates
    if(!writecp_)
      WriteCoordinatesForNurbsShapefunctions(
          geofile,dis,proc0map);
    // write control point coordinates
    else
      WriteCoordinatesForPolynomialShapefunctions(
        geofile,dis,proc0map);
  }
  else
  {
    dserror("Unknown spatial approximation.");
  }

  return proc0map;
}


/*----------------------------------------------------------------------*
  | write node connectivity for every element                  gjb 12/07 |
  *----------------------------------------------------------------------*/
void EnsightWriter::WriteCells(
  std::ofstream& geofile,
  const Teuchos::RCP<DRT::Discretization> dis,
  const Teuchos::RCP<Epetra_Map>& proc0map
  ) const
{
  const Epetra_Map* elementmap = dis->ElementRowMap();

  std::vector<int> nodevector;
  if (myrank_>0)
  {
    //reserve sufficient memory for storing the node connectivity
    //(ghosted nodes included)
    nodevector.reserve(dis->NumMyColNodes());
  }

  // for each found distype write block of the same typed elements
  NumElePerDisType::const_iterator iter;
  for (iter=numElePerDisType_.begin(); iter != numElePerDisType_.end(); ++iter)
  {
    const DRT::Element::DiscretizationType distypeiter = iter->first;
    const int ne = GetNumEleOutput(distypeiter, iter->second);
    const std::string ensightCellType = GetEnsightString(distypeiter);

    if (myrank_ == 0)
    {
      std::cout << "writing "<< iter->second<< " "<< DRT::DistypeToString(distypeiter) << " element(s) as "
           << ne << " " << ensightCellType << " ensight cell(s)..." << std::endl;
      Write(geofile, ensightCellType);
      Write(geofile,ne);
    }

    nodevector.clear();

    // loop all available elements
    for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
    {
      DRT::Element* const actele = dis->gElement(elementmap->GID(iele));
      if (actele->Shape() == distypeiter)
      {
        DRT::Node** const nodes = actele->Nodes();
        switch (actele->Shape())
        {
        case DRT::Element::point1:
        case DRT::Element::line2:
     // case DRT::Element::line3: // Ensight format supports line3, Paraview does not.
        case DRT::Element::hex8:
        case DRT::Element::quad4:
        case DRT::Element::quad8:
        case DRT::Element::tet4:
        case DRT::Element::tet10:
        case DRT::Element::tri3:
        case DRT::Element::tri6:
        case DRT::Element::wedge6:
        case DRT::Element::wedge15:
        case DRT::Element::pyramid5:
        {
          // standard case with direct support
          const int numnp = actele->NumNode();
          for (int inode=0; inode<numnp; ++inode)
          {
            if (myrank_==0) // proc0 can write its elements immediately
              Write(geofile, proc0map->LID(nodes[inode]->Id())+1);
            else // elements on other procs have to store their global node ids
              nodevector.push_back(nodes[inode]->Id());
          }
          break;
        }
        case DRT::Element::hex20:
        {
          // standard case with direct support
          const int numnp = actele->NumNode();
          for (int inode=0; inode<numnp; ++inode)
          {
            if (myrank_==0) // proc0 can write its elements immediately
              Write(geofile, proc0map->LID(nodes[Hex20_BaciToEnsightGold[inode]]->Id())+1);
            else // elements on other procs have to store their global node ids
              nodevector.push_back(nodes[Hex20_BaciToEnsightGold[inode]]->Id());
          }
          break;
        }
        case DRT::Element::hex16:
        {
          // write subelements
          for (int isubele=0; isubele<GetNumSubEle(DRT::Element::hex16); ++isubele)
            for (int isubnode=0; isubnode<8; ++isubnode)
              if (myrank_==0) // proc0 can write its elements immidiately
                Write(geofile, proc0map->LID(nodes[subhex16map[isubele][isubnode]]->Id())
                      +1);
              else // elements on other procs have to store their global node ids
                nodevector.push_back(nodes[subhex16map[isubele][isubnode]]->Id());
          break;
        }
        case DRT::Element::hex18:
        {
          // write subelements
          for (int isubele=0; isubele<GetNumSubEle(DRT::Element::hex18); ++isubele)
            for (int isubnode=0; isubnode<8; ++isubnode)
              if (myrank_==0) // proc0 can write its elements immidiately
                Write(geofile, proc0map->LID(nodes[subhex18map[isubele][isubnode]]->Id())
                      +1);
              else // elements on other procs have to store their global node ids
                nodevector.push_back(nodes[subhex18map[isubele][isubnode]]->Id());
          break;
        }
        case DRT::Element::hex27:
        {
          // write subelements
          for (int isubele=0; isubele<GetNumSubEle(DRT::Element::hex27); ++isubele)
            for (int isubnode=0; isubnode<8; ++isubnode)
              if (myrank_==0) // proc0 can write its elements immidiately
                Write(geofile, proc0map->LID(nodes[subhexmap[isubele][isubnode]]->Id())
                      +1);
              else // elements on other procs have to store their global node ids
                nodevector.push_back(nodes[subhexmap[isubele][isubnode]]->Id());
          break;
        }
        case DRT::Element::quad9:
        {
          // write subelements
          for (int isubele=0; isubele<GetNumSubEle(DRT::Element::quad9); ++isubele)
            for (int isubnode=0; isubnode<4; ++isubnode)
              if (myrank_==0) // proc0 can write its elements immidiately
                Write(geofile, proc0map->LID(nodes[subquadmap[isubele][isubnode]]->Id())
                      +1);
              else // elements on other procs have to store their global node ids
                nodevector.push_back(nodes[subquadmap[isubele][isubnode]]->Id());
          break;
        }
        case DRT::Element::line3:
        {
          // write subelements
          for (int isubele=0; isubele<GetNumSubEle(DRT::Element::line3); ++isubele)
            for (int isubnode=0; isubnode<2; ++isubnode)
              if (myrank_==0) // proc0 can write its elements immidiately
                Write(geofile, proc0map->LID(nodes[sublinemap[isubele][isubnode]]->Id())
                      +1);
              else // elements on other procs have to store their global node ids
                nodevector.push_back(nodes[sublinemap[isubele][isubnode]]->Id());
          break;
        }
        case DRT::Element::nurbs4:
        {
          if(!writecp_)
            WriteNurbsCell(
                actele->Shape(),
                actele->Id()   ,
                geofile        ,
                nodevector     ,
                dis            ,
                proc0map       );
          else
          {
            // standard case with direct support
            const int numnp = actele->NumNode();
            for (int inode=0; inode<numnp; ++inode)
            {
              if (myrank_==0) // proc0 can write its elements immediately
                Write(geofile, proc0map->LID(nodes[inode]->Id())+1);
              else // elements on other procs have to store their global node ids
                nodevector.push_back(nodes[inode]->Id());
            }
          }
          break;
        }
        case DRT::Element::nurbs9:
        {
          if(!writecp_)
            WriteNurbsCell(
                actele->Shape(),
                actele->Id()   ,
                geofile        ,
                nodevector     ,
                dis            ,
                proc0map       );
          else
          {
            // write subelements
            for (int isubele=0; isubele<GetNumSubEle(DRT::Element::quad9); ++isubele)
              for (int isubnode=0; isubnode<4; ++isubnode)
                if (myrank_==0) // proc0 can write its elements immidiately
                  Write(geofile, proc0map->LID(nodes[subquadmap[isubele][isubnode]]->Id())
                        +1);
                else // elements on other procs have to store their global node ids
                  nodevector.push_back(nodes[subquadmap[isubele][isubnode]]->Id());
          }
          break;
        }
        case DRT::Element::nurbs27:
        {
          if(!writecp_)
            WriteNurbsCell(
                actele->Shape(),
                actele->Id()   ,
                geofile        ,
                nodevector     ,
                dis            ,
                proc0map       );
          else
          {
            // write subelements
            for (int isubele=0; isubele<GetNumSubEle(DRT::Element::hex27); ++isubele)
              for (int isubnode=0; isubnode<8; ++isubnode)
                if (myrank_==0) // proc0 can write its elements immidiately
                  Write(geofile, proc0map->LID(nodes[subhexmap[isubele][isubnode]]->Id())
                        +1);
                else // elements on other procs have to store their global node ids
                  nodevector.push_back(nodes[subhexmap[isubele][isubnode]]->Id());
          }
          break;
        }
        break;
        default:
          dserror("don't know, how to write this element type as a cell");
        }
      }
    }

    // now do some communicative work for the parallel case:
    // proc 1 to proc n have to send their stored node connectivity to proc0
    // which does the writing

    WriteNodeConnectivityPar(geofile, dis, nodevector, proc0map);
  }
  return;
}


/*!
 * \brief communicate and write node connectivity in parallel case

 \author gjb
 \date 12/07
*/
void EnsightWriter::WriteNodeConnectivityPar(
  std::ofstream& geofile,
  const Teuchos::RCP<DRT::Discretization> dis,
  const std::vector<int>& nodevector,
  const Teuchos::RCP<Epetra_Map> proc0map) const
{
  // no we have communicate the connectivity infos from proc 1...proc n to proc 0

  std::vector<char> sblock; // sending block
  std::vector<char> rblock; // recieving block

  // create an exporter for communication
  DRT::Exporter exporter(dis->Comm());

  // pack my node ids into sendbuffer
  sblock.clear();

  DRT::PackBuffer data;
  DRT::ParObject::AddtoPack(data,nodevector);
  data.StartPacking();
  DRT::ParObject::AddtoPack(data,nodevector);
  swap( sblock, data() );

  // now we start the communication
  for (unsigned int pid=0;pid<static_cast<unsigned int>(dis->Comm().NumProc());++pid)
  {
    MPI_Request request;
    int         tag    =0;
    int         frompid=pid;
    int         topid  =0;
    int         length=sblock.size();

    //--------------------------------------------------
    // proc pid sends its values to proc 0
    if (myrank_==pid)
    {
      exporter.ISend(frompid,topid,
                     &(sblock[0]),sblock.size(),
                     tag,request);
    }

    //--------------------------------------------------
    // proc 0 receives from proc pid
    rblock.clear();
    if (myrank_ == 0)
    {
      exporter.ReceiveAny(frompid,tag,rblock,length);
      if(tag!=0)
      {
        dserror("Proc 0 received wrong message (ReceiveAny)");
      }
      exporter.Wait(request);
    }

    // for safety
    exporter.Comm().Barrier();

    //--------------------------------------------------
    // Unpack received block and write the data
    if (myrank_==0)
    {
      std::vector<char>::size_type index = 0;
      std::vector<int> nodeids;
      // extract data from recieved package
      while (index < rblock.size())
      {
        DRT::ParObject::ExtractfromPack(index,rblock,nodeids);
      }
      // compute node lid based on proc0map and write it to file
      for(int i=0;i<(int) nodeids.size();++i)
      {
        // using the same map as for the writing the node coordinates
        int id = (proc0map->LID(nodeids[i]))+1;
        Write(geofile, id);
      }
      nodeids.clear();
    } // end unpack

    // for safety
    exporter.Comm().Barrier();

  }// for pid

  return;
}


/*!
 * \brief parse all elements and get the global(!) number of elements for each distype
 * \author gjb
 * \date 01/08
 */
EnsightWriter::NumElePerDisType
EnsightWriter::GetNumElePerDisType(
  const Teuchos::RCP<DRT::Discretization> dis
  ) const
{
  const Epetra_Map* elementmap = dis->ElementRowMap();

  NumElePerDisType numElePerDisType;
  for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
  {
    DRT::Element* actele = dis->gElement(elementmap->GID(iele));
    const DRT::Element::DiscretizationType distype = actele->Shape();
    // update counter for current distype
    numElePerDisType[distype]++;
  }

  // in parallel case we have to sum up the local element distype numbers

  // determine maximum number of possible element discretization types
  DRT::Element::DiscretizationType numeledistypes = DRT::Element::max_distype;

  // write the final local numbers into a vector
  std::vector<int> myNumElePerDisType(numeledistypes);
  NumElePerDisType::const_iterator iter;
  for (iter=numElePerDisType.begin(); iter != numElePerDisType.end(); ++iter)
  {
    const DRT::Element::DiscretizationType distypeiter = iter->first;
    const int ne = iter->second;
    myNumElePerDisType[distypeiter]+=ne;
  }

  // wait for all procs before communication is started
  (dis->Comm()).Barrier();

  // form the global sum
  std::vector<int> globalnumeleperdistype(numeledistypes);
  (dis->Comm()).SumAll(&(myNumElePerDisType[0]),&(globalnumeleperdistype[0]),numeledistypes);

  // create return argument containing the global element numbers per distype
  NumElePerDisType globalNumElePerDisType;
  for(int i =0; i<numeledistypes ;++i)
  {
    if (globalnumeleperdistype[i]>0) // no entry when we have no element of this type
      globalNumElePerDisType[DRT::Element::DiscretizationType(i)]=globalnumeleperdistype[i];
  }

  return globalNumElePerDisType;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int EnsightWriter::GetNumEleOutput(
  const DRT::Element::DiscretizationType distype,
  const int numele) const
{
  return GetNumSubEle(distype)*numele;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int EnsightWriter::GetNumSubEle(
  const DRT::Element::DiscretizationType distype) const
{
  switch (distype)
  {
  case DRT::Element::hex18:
    return 4;
    break;
  case DRT::Element::hex27:
    return 8;
    break;
  case DRT::Element::nurbs27:
    return 8;
    break;
  case DRT::Element::quad9:
    return 4;
    break;
  case DRT::Element::nurbs9:
    return 4;
    break;
  case DRT::Element::line3:
    return 2;
    break;
  default:
    return 1; // no element splitting necessary
  }
}

/*!
 * \brief parse all elements and get the global ids of the elements for each distype
 */
EnsightWriter::EleGidPerDisType
EnsightWriter::GetEleGidPerDisType(
  const Teuchos::RCP<DRT::Discretization> dis,
  NumElePerDisType numeleperdistype
  ) const
{
  const Epetra_Map* elementmap = dis->ElementRowMap();

  EleGidPerDisType eleGidPerDisType;

  //allocate memory
  NumElePerDisType::const_iterator iter;
  for (iter=numElePerDisType_.begin(); iter != numElePerDisType_.end(); ++iter)
  {
    eleGidPerDisType[iter->first].reserve(iter->second);
  }

  for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
  {
    const int gid = elementmap->GID(iele);
    DRT::Element* actele = dis->gElement(gid);
    const DRT::Element::DiscretizationType distype = actele->Shape();
    // update counter for current distype
    eleGidPerDisType[distype].push_back(gid);
  }

  // in parallel case we have to provide the ele gids located on other procs as well
  EleGidPerDisType globaleleGidPerDisType;
  NumElePerDisType::const_iterator iterator;

  for (iterator=numElePerDisType_.begin(); iterator != numElePerDisType_.end(); ++iterator)
  {
    // wait for all procs before communication is started
    (dis->Comm()).Barrier();

    // no we have to communicate everything from proc 1...proc n to proc 0
    std::vector<char> sblock; // sending block
    std::vector<char> rblock; // recieving block

    // create an exporter for communication
    DRT::Exporter exporter(dis->Comm());

    // pack my element gids of this discretization type into sendbuffer
    sblock.clear();

    DRT::PackBuffer data;
    DRT::ParObject::AddtoPack(data,eleGidPerDisType[iterator->first]);
    data.StartPacking();
    DRT::ParObject::AddtoPack(data,eleGidPerDisType[iterator->first]);
    swap( sblock, data() );

    // now we start the communication
    for (unsigned int pid=0;pid<static_cast<unsigned int>(dis->Comm().NumProc());++pid)
    {
      MPI_Request request;
      int         tag    =0;
      int         frompid=pid;
      int         topid  =0;
      int         length=sblock.size();

      //--------------------------------------------------
      // proc pid sends its values to proc 0
      if (myrank_==pid)
      {
        exporter.ISend(frompid,topid,
                       &(sblock[0]),sblock.size(),
                       tag,request);
      }

      //--------------------------------------------------
      // proc 0 receives from proc pid
      rblock.clear();
      if (myrank_ == 0)
      {
        exporter.ReceiveAny(frompid,tag,rblock,length);
        if(tag!=0)
        {
          dserror("Proc 0 received wrong message (ReceiveAny)");
        }
        exporter.Wait(request);
      }

      // for safety
      exporter.Comm().Barrier();

      //--------------------------------------------------
      // Unpack received block and write the data
      if (myrank_==0)
      {
        std::vector<char>::size_type index = 0;
        std::vector<int> elegids;
        // extract data from recieved package
        while (index < rblock.size())
        {
          DRT::ParObject::ExtractfromPack(index,rblock,elegids);
        }
        for(int i=0;i<(int) elegids.size();++i)
        {
          globaleleGidPerDisType[iterator->first].push_back(elegids[i]);
        }
        elegids.clear();
      } // end unpack

    }// for pid
  } // for iter over type

    // note: this map is only filled on proc 0 !!!!
  return globaleleGidPerDisType;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::string EnsightWriter::GetEnsightString(
  const DRT::Element::DiscretizationType distype) const
{
  std::map<DRT::Element::DiscretizationType, std::string>::const_iterator entry;
  switch (distype)
  {
  case DRT::Element::hex18:
    entry = distype2ensightstring_.find(DRT::Element::hex8);
  case DRT::Element::hex27:
    entry = distype2ensightstring_.find(DRT::Element::hex8);
    break;
  case DRT::Element::quad9:
    entry = distype2ensightstring_.find(DRT::Element::quad4);
    break;
  case DRT::Element::tet10:
    entry = distype2ensightstring_.find(DRT::Element::tet10);
    break;
  default:
    entry = distype2ensightstring_.find(distype);
    break;
  }
  if (entry == distype2ensightstring_.end())
    dserror("no entry in distype2ensightstring_ found");
  return entry->second;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteResult(const std::string groupname,
                                const std::string name,
                                const ResultType restype,
                                const int numdf,
                                const int from/*=0*/,
                                const bool fillzeros/*=false*/)
{
  PostResult result(field_);
  bool foundit = false;
  while (result.next_result(groupname))
  {
    if (map_has_map(result.group(), groupname.c_str()))
    {
      foundit = true;
      break;
    }
  }
  if (!foundit) return;

  // new for file continuation
  bool multiple_files = false;

  // For NURBS control point output filename is extened by "_cp"
  std::string aux = "";
  if(writecp_)
  {
    aux = aux + "_cp";
  }

  // open file
  const std::string filename = filename_ + "_"+ field_->name() + aux + "."+ name;
  std::ofstream file;
  int startfilepos = 0;
  if (myrank_==0)
  {
    file.open(filename.c_str());
    startfilepos = file.tellp(); // file position should be zero, but we stay flexible
  }

  std::map<std::string, std::vector<std::ofstream::pos_type> > resultfilepos;
  int stepsize = 0;

  // distinguish between node- and element-based results
  switch(restype)
  {
  case dofbased:
  {
    if (myrank_==0)
      std::cout<<"writing dof-based field "<<name<<std::endl;
    // store information for later case file creation
    variableresulttypemap_[name] = "node";

    //const Epetra_Map* nodemap = field_->discretization()->NodeRowMap();
    //const int numnp = nodemap->NumGlobalElements();
    //int effnumdf = numdf;
    //if (numdf==2) effnumdf=3; // in 2D we still have to write a 3D vector with zero z-components!!!
    // get the number of bits to be written each time step
    //const int stepsize = 5*80+sizeof(int)+effnumdf*numnp*sizeof(float);

    WriteDofResultStep(file, result, resultfilepos, groupname, name, numdf, from, fillzeros);
    // how many bits are necessary per time step (we assume a fixed size)?
    if (myrank_==0)
    {
      stepsize = ((int) file.tellp())-startfilepos;
      if (stepsize <= 0) dserror("found invalid step size for result file");
    }
    else
      stepsize = 1; //use dummy value on other procs

    while (result.next_result(groupname))
    {
      if (map_has_map(result.group(), groupname.c_str()))
      {
        const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
        if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
        {
          FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
        }
        WriteDofResultStep(file, result, resultfilepos, groupname, name, numdf, from, fillzeros);
      }
    }
  }
  break;

  case nodebased:
  {
    if (myrank_==0)
      std::cout<<"writing node-based field "<<name<<std::endl;
    // store information for later case file creation
    variableresulttypemap_[name] = "node";

    WriteNodalResultStep(file, result, resultfilepos, groupname, name, numdf);
    // how many bits are necessary per time step (we assume a fixed size)?
    if (myrank_==0)
    {
      stepsize = ((int) file.tellp())-startfilepos;
      if (stepsize <= 0) dserror("found invalid step size for result file");
    }
    else
      stepsize = 1; //use dummy value on other procs

    while (result.next_result(groupname))
    {
      const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
      if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
      {
        FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
      }
      WriteNodalResultStep(file, result, resultfilepos, groupname, name, numdf);
    }
  }
  break;

  case elementdof:
  {
    if (myrank_==0)
      std::cout<<"writing element based field "<<name<<std::endl;
    // store information for later case file creation
    variableresulttypemap_[name] = "element";

    WriteElementDOFResultStep(file, result, resultfilepos, groupname, name, numdf, from);
    // how many bits are necessary per time step (we assume a fixed size)?
    if (myrank_==0)
    {
      stepsize = ((int) file.tellp())-startfilepos;
      if (stepsize <= 0) dserror("found invalid step size for result file");
    }
    else
      stepsize = 1; //use dummy value on other procs

    while (result.next_result(groupname))
    {
      const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
      if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
      {
        FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
      }
      WriteElementDOFResultStep(file, result, resultfilepos, groupname, name, numdf, from);
    }
  }
  break;

  case elementbased:
  {
    if (myrank_==0)
      std::cout<<"writing element-based field "<<name<<std::endl;
    // store information for later case file creation
    variableresulttypemap_[name] = "element";

    WriteElementResultStep(file, result, resultfilepos, groupname, name, numdf, from);
    // how many bits are necessary per time step (we assume a fixed size)?
    if (myrank_==0)
    {
      stepsize = ((int) file.tellp())-startfilepos;
      if (stepsize <= 0) dserror("found invalid step size for result file");
    }
    else
      stepsize = 1; //use dummy value on other procs

    while (result.next_result(groupname))
    {
      const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
      if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
      {
        FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
      }
      WriteElementResultStep(file, result, resultfilepos, groupname, name, numdf, from);
    }
  }
  break;

  case no_restype:
  case max_restype:
    dserror("found invalid result type");
    break;
  default:
    dserror("Invalid output type in WriteResult");
  } // end of switch(restype)

  // store information for later case file creation
  filesetmap_[name].push_back(file.tellp()/stepsize);// has to be done BEFORE writing the index table
  variablenumdfmap_[name] = numdf;
  variablefilenamemap_[name] = filename;
  // store solution times vector for later case file creation
  {
    PostResult res = PostResult(field_); // this is needed!
    std::vector<double> restimes = res.get_result_times(field_->name(),groupname);
    timesetmap_[name] = restimes;
  }

  // append index table
  WriteIndexTable(file, resultfilepos[name]);
  resultfilepos[name].clear();

  // close result file
  if (file.is_open())
    file.close();

  return;
}

void EnsightWriter::WriteResultOneTimeStep(PostResult& result,
                                const std::string groupname,
                                const std::string name,
                                const ResultType restype,
                                const int numdf,
                                bool firststep,
                                bool laststep,
                                const int from)
{

  if(not map_has_map(result.group(), groupname.c_str()))
    dserror("expected result: %s in step %i. Probably a return is missing here. But check!", groupname.c_str(), result.step());

  // FILE CONTINUATION NOT SUPPORTED DUE TO ITS COMPLEXITY FOR VARIABLE GEOMETRY ghamm 03.03.2014
  // guessing whether a new file is necessary should be possible with the current method because the file size is no hard limit
  // more tricky are the things that happen in FileSwitcher --> need adaptation for variable geometry
  bool multiple_files = false;

  // For NURBS control point output filename is extened by "_cp"
  std::string aux = "";
  if(writecp_)
  {
    aux = aux + "_cp";
  }

  // open file
  const std::string filename = filename_ + "_"+ field_->name() + aux + "."+ name;
  std::ofstream file;
  int startfilepos = 0;
  if (myrank_==0)
  {
    if( firststep )
    {
      file.open(filename.c_str());
      startfilepos = file.tellp(); // file position should be zero, but we stay flexible
    }
    else
    {
      file.open(filename.c_str(), std::ofstream::app | std::ofstream::binary);
      startfilepos = file.tellp(); // file position should be zero, but we stay flexible
    }
  }

  int stepsize = 0;

  // distinguish between node- and element-based results
  switch(restype)
  {
  case dofbased:
  {
    if (myrank_==0)
      std::cout<<"writing node-based field "<<name<<std::endl;
    // store information for later case file creation
    variableresulttypemap_[name] = "node";

//    const Epetra_Map* nodemap = field_->discretization()->NodeRowMap();
//    const int numnp = nodemap->NumGlobalElements();
//    int effnumdf = numdf;
//    if (numdf==2) effnumdf=3; // in 2D we still have to write a 3D vector with zero z-components!!!
//    // get the number of bits to be written each time step
//    const int stepsize = 5*80+sizeof(int)+effnumdf*numnp*sizeof(float);

    if (map_has_map(result.group(), groupname.c_str()))
    {
      WriteDofResultStep(file, result, resultfilepos_, groupname, name, numdf, from, false);

      // how many bits are necessary per time step (we assume a fixed size)?
      if (myrank_==0)
      {
        stepsize = ((int) file.tellp())-startfilepos;
        if (stepsize <= 0) dserror("found invalid step size for result file");
      }
      else
        stepsize = 1; //use dummy value on other procs

      const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
      if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
      {
        dserror("file continuation not supported for variable geometries");
        FileSwitcher(file, multiple_files, filesetmap_, resultfilepos_, stepsize, name, filename);
      }
    }
  }
  break;

  case nodebased:
  {
    if (myrank_==0)
      std::cout<<"writing node-based field "<<name<<std::endl;
    // store information for later case file creation
    variableresulttypemap_[name] = "node";

    WriteNodalResultStep(file, result, resultfilepos_, groupname, name, numdf);

    // how many bits are necessary per time step (we assume a fixed size)?
    if (myrank_==0)
    {
      stepsize = ((int) file.tellp())-startfilepos;
      if (stepsize <= 0) dserror("found invalid step size for result file");
    }
    else
      stepsize = 1; //use dummy value on other procs

    const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
    if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
    {
      dserror("file continuation not supported for variable geometries");
      FileSwitcher(file, multiple_files, filesetmap_, resultfilepos_, stepsize, name, filename);
    }
  }
  break;

  case elementdof:
  {
    if (myrank_==0)
      std::cout<<"writing element based field "<<name<<std::endl;
    // store information for later case file creation
    variableresulttypemap_[name] = "element";

    WriteElementDOFResultStep(file, result, resultfilepos_, groupname, name, numdf, from);

    // how many bits are necessary per time step (we assume a fixed size)?
    if (myrank_==0)
    {
      stepsize = ((int) file.tellp())-startfilepos;
      if (stepsize <= 0) dserror("found invalid step size for result file");
    }
    else
      stepsize = 1; //use dummy value on other procs

    const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
    if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
    {
      dserror("file continuation not supported for variable geometries");
      FileSwitcher(file, multiple_files, filesetmap_, resultfilepos_, stepsize, name, filename);
    }
  }
  break;

  case elementbased:
  {
    if (myrank_==0)
      std::cout<<"writing element-based field "<<name<<std::endl;
    // store information for later case file creation
    variableresulttypemap_[name] = "element";

    WriteElementResultStep(file, result, resultfilepos_, groupname, name, numdf, from);

    // how many bits are necessary per time step (we assume a fixed size)?
    if (myrank_==0)
    {
      stepsize = ((int) file.tellp())-startfilepos;
      if (stepsize <= 0) dserror("found invalid step size for result file");
    }
    else
      stepsize = 1; //use dummy value on other procs

    const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
    if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
    {
      dserror("file continuation not supported for variable geometries");
      FileSwitcher(file, multiple_files, filesetmap_, resultfilepos_, stepsize, name, filename);
    }
  }
  break;

  case no_restype:
  case max_restype:
    dserror("found invalid result type");
  break;
  default:
    dserror("Invalid output type in WriteResult");
  break;
  } // end of switch(restype)

  // store information for later case file creation
  timesetmap_[name].push_back(result.time());
  if(laststep)
  {
    filesetmap_[name].push_back((int)(timesetmap_[name].size()));
    variablenumdfmap_[name] = numdf;
    variablefilenamemap_[name] = filename;

    // append index table
    WriteIndexTable(file, resultfilepos_[name]);
    resultfilepos_[name].clear();
  }

  // close result file
  if (file.is_open())
    file.close();

  return;
}



void
EnsightWriter::WriteSpecialField (
    SpecialFieldInterface &special,
    PostResult& result,
    const ResultType  restype,
    const std::string &groupname,
    const std::vector<std::string> &fieldnames,
    const std::string &outinfo
    )
{
  const int numfiles = fieldnames.size();

  // new for file continuation
  std::vector<bool> multiple_files(numfiles);
  for (int i=0;i<numfiles;++i)
  {
    multiple_files[i] = false;
  }

  // open file
  std::vector<std::string> filenames(numfiles);
  for (int i=0;i<numfiles;++i)
  {
    filenames[i] = filename_ + "_"+ field_->name() + "."+ fieldnames[i];
  }

  std::vector<Teuchos::RCP<std::ofstream> > files(numfiles);
  std::vector<int> startfilepos(numfiles);
  for (int i=0;i<numfiles;++i)
    startfilepos[i] = 0;
  for (int i=0;i<numfiles;++i)
  {
    files[i] = Teuchos::rcp(new std::ofstream);

    if (myrank_==0)
    {
      files[i]->open(filenames[i].c_str());
      startfilepos[i] = files[i]->tellp(); // file position should be zero, but we stay flexible
    }
  }

  std::map<std::string, std::vector<std::ofstream::pos_type> > resultfilepos;
  std::vector<int> stepsize(numfiles);
  for (int i=0;i<numfiles;++i)
  {
    stepsize[i]=0;
  }

  if (myrank_==0)
  {
    if (restype == nodebased)
      std::cout << "writing node-based ";
    else if (restype == elementbased)
      std::cout << "writing element-based ";
    else
      std::cout << "writing [unknown type] ";
    std::cout << outinfo << std::endl;
  }

  // store information for later case file creation
  for (int i=0;i<numfiles;++i)
  {
    variableresulttypemap_[fieldnames[i]] = (restype == nodebased ? "node" : "element");
  }

  special(files,result,resultfilepos,groupname,fieldnames);

  // how many bits are necessary per time step (we assume a fixed size)?
  if (myrank_==0)
  {
    for (int i=0;i<numfiles;++i)
    {
      stepsize[i] = ((int) files[i]->tellp())-startfilepos[i];
      if (stepsize[i] <= 0) dserror("found invalid step size for result file");
    }
  }
  else
  {
    for (int i=0;i<numfiles;++i)
    {
      stepsize[i] = 1; //use dummy value on other procs
    }
  }

  while (result.next_result())
  {
    for (int i=0;i<numfiles;++i)
    {
      const int indexsize = 80+2*sizeof(int)+(files[i]->tellp()/stepsize[i]+2)*sizeof(long);
      if (static_cast<long unsigned int>(files[i]->tellp())+stepsize[i]+indexsize>= FILE_SIZE_LIMIT_)
      {
        bool mf = multiple_files[i];
        FileSwitcher(*(files[i]),mf,filesetmap_,resultfilepos,stepsize[i],fieldnames[i],filenames[i]);
      }
    }

    special(files,result,resultfilepos,groupname,fieldnames);
  }
  // store information for later case file creation

  const std::vector<int> numdfmap = special.NumDfMap();
  dsassert(static_cast<int>(numdfmap.size()) == numfiles, "Wrong number of components in NumDfMap.");
  for (int i=0;i<numfiles;++i)
  {
    filesetmap_[fieldnames[i]].push_back(files[i]->tellp()/stepsize[i]);// has to be done BEFORE writing the index table
    variablenumdfmap_[fieldnames[i]] = numdfmap[i];
    variablefilenamemap_[fieldnames[i]] = filenames[i];
  }

  // store solution times vector for later case file creation
  for (int i=0;i<numfiles;++i)
  {
    PostResult res = PostResult(field_); // this is needed!
    std::vector<double> restimes = res.get_result_times(field_->name(),groupname);
    timesetmap_[fieldnames[i]] = restimes;
  }

  //append index table
  for (int i=0;i<numfiles;++i)
  {
    WriteIndexTable(*(files[i]), resultfilepos[fieldnames[i]]);
    resultfilepos[fieldnames[i]].clear();
    if (files[i]->is_open()) files[i]->close();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::FileSwitcher(
  std::ofstream& file,
  bool& multiple_files,
  std::map<std::string, std::vector<int> >& filesetmap,
  std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
  const int stepsize,
  const std::string name,
  const std::string filename
  ) const
{
  if (myrank_==0)
  {
    std::ostringstream newfilename;

    if (multiple_files == false)
    {
      multiple_files = true;

      std::vector<int> numsteps;
      numsteps.push_back(file.tellp()/stepsize);
      filesetmap[name] = numsteps;

      // append index table
      WriteIndexTable(file, resultfilepos[name]);
      resultfilepos[name].clear();
      file.close();
      rename(filename.c_str(), (filename+"001").c_str());

      newfilename << filename << "002";
    }
    else
    {
      filesetmap[name].push_back(file.tellp()/stepsize);

      // append index table
      WriteIndexTable(file, resultfilepos[name]);
      resultfilepos[name].clear();
      file.close();

      newfilename << filename;
      newfilename.width(3);
      newfilename.fill('0');
      newfilename << filesetmap[name].size()+1;
    }
    file.open(newfilename.str().c_str());
  } // if (myrank_==0)
  return;
}


/*!
  \brief Write nodal values for one timestep for dof-based vectors
  Each node has to have the same number of dofs.
*/
void EnsightWriter::WriteDofResultStep(std::ofstream& file,
                                       PostResult& result,
                                       std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
                                       const std::string& groupname,
                                       const std::string& name,
                                       const int numdf,
                                       const int frompid,
                                       const bool fillzeros) const
{
  //-------------------------------------------
  // write some key words and read result data
  //-------------------------------------------

  std::vector<std::ofstream::pos_type>& filepos = resultfilepos[name];
  Write(file, "BEGIN TIME STEP");
  filepos.push_back(file.tellp());
  Write(file, "description");
  Write(file, "part");
  Write(file, field_->field_pos()+1);
  Write(file, "coordinates");

  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();
  const Epetra_Map* nodemap = dis->NodeRowMap(); //local node row map
  const int numnp = nodemap->NumGlobalElements();

  const Teuchos::RCP<Epetra_Vector> data = result.read_result(groupname);
  const Epetra_BlockMap& datamap = data->Map();

  // do stupid conversion into Epetra map
  Teuchos::RCP<Epetra_Map> epetradatamap;
  epetradatamap = Teuchos::rcp(new Epetra_Map(datamap.NumGlobalElements(),
                                     datamap.NumMyElements(),
                                     datamap.MyGlobalElements(),
                                     0,
                                     datamap.Comm()));

  // determine offset of dofs in case of multiple discretizations in
  // separate files (e.g. multi-scale problems). during calculation,
  // dofs are numbered consecutively for all discretizations. in the
  // post-processing phase, when only one discretization is called,
  // numbering always starts with 0, so a potential offset needs to be
  // taken into account.
  // NOTE 1: for the pressure result vector of FLUID calculations,
  //         this offset is 2 or 3, depending on the number of space dimensions.
  // NOTE 2: this command is only valid, if you use NOT MORE processors
  //         for filtering than for computation. Otherwise we have empty procs
  //         owning empty maps, and therefore epetradatamap->MinAllGID()
  //         will always return zero, resulting in a wrong offset value.
  //         This is the only reason, why not to use more (and empty) procs.
  //         All other code parts of post_drt_ensight can handle that.
  // NOTE 3: This check works for polynomial as well as NURBS discretizations!
  if (epetradatamap->NumMyElements()<1)
    dserror("Proc %d is empty. Do not use more procs for postprocessing than for calculation.",myrank_);
  int offset = epetradatamap->MinAllGID() - dis->DofRowMap()->MinAllGID();

  //switch between nurbs an others
  if(field_->problem()->SpatialApproximation()=="Nurbs" && !writecp_)
  {
    WriteDofResultStepForNurbs(
      file ,
      numdf,
      data,
      name,
      offset
      );
  }
  else if(field_->problem()->SpatialApproximation()=="Polynomial" or
          field_->problem()->SpatialApproximation()=="Meshfree" or
          field_->problem()->SpatialApproximation()=="HDG" or
          (field_->problem()->SpatialApproximation()=="Nurbs" && writecp_))
  {
    //------------------------------------------------------
    // each processor provides its result values for proc 0
    //------------------------------------------------------

    Teuchos::RCP<Epetra_Map> proc0datamap;
    proc0datamap = LINALG::AllreduceEMap(*epetradatamap,0);

    // contract result values on proc0 (proc0 gets everything, other procs empty)
    Epetra_Import proc0dataimporter(*proc0datamap,*epetradatamap);
    Teuchos::RCP<Epetra_Vector> proc0data = Teuchos::rcp(new Epetra_Vector(*proc0datamap));
    int err = proc0data->Import(*data,proc0dataimporter,Insert);
    if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

    const Epetra_BlockMap& finaldatamap = proc0data->Map();

    //------------------------------------------------------------------
    // each processor provides its dof global id information for proc 0
    //------------------------------------------------------------------

    // would be nice to have an Epetra_IntMultiVector, instead of casting to doubles
    Teuchos::RCP<Epetra_MultiVector> dofgidpernodelid = Teuchos::rcp(new Epetra_MultiVector(*nodemap,numdf));
    dofgidpernodelid->PutScalar(-1.0);

    const int mynumnp = nodemap->NumMyElements();
    for (int idf=0; idf<numdf; ++idf)
    {
      for (int inode=0; inode<mynumnp; inode++)
      {
        DRT::Node* n = dis->lRowNode(inode);

        const double dofgid = (double) dis->Dof(n, frompid + idf) + offset;
        if (dofgid > -1.0)
        {
          dofgidpernodelid->ReplaceMyValue(inode, idf, dofgid);
        }
        else
        {
          dserror("Error while creating Epetra_MultiVector dofgidperlocalnodeid");
        }
      }
    }

    // contract Epetra_MultiVector on proc0 (proc0 gets everything, other procs empty)
    Teuchos::RCP<Epetra_MultiVector> dofgidpernodelid_proc0 = Teuchos::rcp(new Epetra_MultiVector(*proc0map_,numdf));
    Epetra_Import proc0dofimporter(*proc0map_,*nodemap);
    err = dofgidpernodelid_proc0->Import(*dofgidpernodelid,proc0dofimporter,Insert);
    if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);


    //---------------
    // write results
    //---------------

    const int finalnumnode = proc0map_->NumGlobalElements();
    if (myrank_==0) // ensures pointer dofgids is valid
    {
      // care for rotationally symmetric periodic boundary conditions
      // note: only vector fields with numdf > 1 require checking!
      std::map<int,double> pbcslavenodemap;
      std::map<int,double>::iterator iter;
      if(numdf > 1)
        FLD::GetRelevantSlaveNodesOfRotSymPBC(pbcslavenodemap, dis);

      double* dofgids = (dofgidpernodelid_proc0->Values()); // columnwise data storage
      for (int idf=0; idf<numdf; ++idf)
      {
        for (int inode=0; inode<finalnumnode; inode++) // inode == lid of node because we use proc0map_
        {
          // local storage position of desired dof gid
          const int doflid = inode + (idf*numnp);
          // get the dof global id
          const int actdofgid = (int) (dofgids[doflid]);
          dsassert(actdofgid>= 0, "error while getting dof global id");
          // get the dof local id w.r.t. the finaldatamap
          int lid = finaldatamap.LID(actdofgid);
          if (lid > -1)
          {
            // is the current node a slave of a rot. symm. periodic boundary condition?
            const int nodegid = proc0map_->GID(inode);
            iter = pbcslavenodemap.find(nodegid);
            if (iter != pbcslavenodemap.end())
            {
              // this is the desired component of the rotated vector field result
              double value = FLD::GetComponentOfRotatedVectorField(idf,proc0data,lid,iter->second);
              Write(file, static_cast<float>(value));
            }
            else
            {
              // the standard case
              Write(file, static_cast<float>((*proc0data)[lid]));
            }
          }
          else
          {
            if(fillzeros)
              Write<float>(file, 0.);
            else
              dserror("received illegal dof local id: %d", lid);
          }
        }
      }// for idf

      // 2 component vectors in a 3d problem require a row of zeros.
      // do we really need this?
      if (numdf==2)
      {
        for (int inode=0; inode<numnp; inode++)
        {
          Write<float>(file, 0.);
        }
      }
    } // if (myrank_==0)
  }
  else
  {
    dserror("spatial approximation neither Nurbs nor Polynomial\n");
  }

  Write(file, "END TIME STEP");
  return;
}



/*!
  \brief Write nodal values for one timestep for node-based vectors
  Each node has to have the same number of dofs.
*/
void EnsightWriter::WriteNodalResultStep(std::ofstream& file,
                                         PostResult& result,
                                         std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
                                         const std::string& groupname,
                                         const std::string& name,
                                         const int numdf)
{
  const Teuchos::RCP<Epetra_MultiVector> data = result.read_multi_result(groupname);
  WriteNodalResultStep(file,data,resultfilepos,groupname,name,numdf);
}



/*!
  \brief Write nodal values for one timestep for node-based vectors
  Each node has to have the same number of dofs.
*/
void EnsightWriter::WriteNodalResultStep(std::ofstream& file,
                                         const Teuchos::RCP<Epetra_MultiVector>& data,
                                         std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
                                         const std::string& groupname,
                                         const std::string& name,
                                         const int numdf)
{
  //-------------------------------------------
  // write some key words and read result data
  //-------------------------------------------

  std::vector<std::ofstream::pos_type>& filepos = resultfilepos[name];
  Write(file, "BEGIN TIME STEP");
  filepos.push_back(file.tellp());
  Write(file, "description");
  Write(file, "part");
  Write(file, field_->field_pos()+1);
  Write(file, "coordinates");

  const Epetra_BlockMap& datamap = data->Map();

  //switch between nurbs an others
  if(field_->problem()->SpatialApproximation()=="Nurbs" && !writecp_)
  {
    WriteNodalResultStepForNurbs(
      file ,
      numdf,
      data,
      name,
      0
      );
  }
  else if(field_->problem()->SpatialApproximation()=="Polynomial" or
          field_->problem()->SpatialApproximation()=="Meshfree" or
          field_->problem()->SpatialApproximation()=="HDG" or
          (field_->problem()->SpatialApproximation()=="Nurbs" && writecp_))
  {
    // contract Epetra_MultiVector on proc0 (proc0 gets everything, other procs empty)
      Teuchos::RCP<Epetra_MultiVector> data_proc0 = Teuchos::rcp(new Epetra_MultiVector(*proc0map_,numdf));
    Epetra_Import proc0dofimporter(*proc0map_,datamap);
    int err = data_proc0->Import(*data,proc0dofimporter,Insert);
    if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

    //---------------
    // write results
    //---------------

    const int finalnumnode = proc0map_->NumGlobalElements();

    if (myrank_==0)
    {
      for (int idf=0; idf<numdf; ++idf)
      {
        Epetra_Vector* column = (*data_proc0)(idf);
        for (int inode=0; inode<finalnumnode; inode++) // inode == lid of node because we use proc0map_
        {
          Write(file, static_cast<float>((*column)[inode]));
        }
      }
    } // if (myrank_==0)

    // 2 component vectors in a 3d problem require a row of zeros.
    if (numdf==2)
    {
      for (int inode=0; inode<finalnumnode; inode++)
      {
        Write<float>(file, 0.);
      }
    }
  } // polynomial || meshfree
  else
  {
    dserror("spatial approximation neither Nurbs nor Polynomial\n");
  }

  Write(file, "END TIME STEP");
  return;
}


/*!
  \brief Write element dof values for one timestep

  Each element has to have the same number of dofs.
*/
void EnsightWriter::WriteElementDOFResultStep(
  std::ofstream& file,
  PostResult& result,
  std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
  const std::string& groupname,
  const std::string& name,
  const int numdof,
  const int from
  ) const
{
  //-------------------------------------------
  // write some key words and read result data
  //-------------------------------------------

  std::vector<std::ofstream::pos_type>& filepos = resultfilepos[name];
  Write(file, "BEGIN TIME STEP");
  filepos.push_back(file.tellp());
  Write(file, "description");
  Write(file, "part");
  Write(file, field_->field_pos()+1);

  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();
  const Epetra_Map* elementmap = dis->ElementRowMap(); //local node row map

  const Teuchos::RCP<Epetra_Vector> data = result.read_result(groupname);
  const Epetra_BlockMap& datamap = data->Map();

  // do stupid conversion into Epetra map
  Teuchos::RCP<Epetra_Map> epetradatamap;
  epetradatamap = Teuchos::rcp(new Epetra_Map(datamap.NumGlobalElements(),
                                     datamap.NumMyElements(),
                                     datamap.MyGlobalElements(),
                                     0,
                                     datamap.Comm()));
#if 0
  if (datamap->PointSameAs(*proc0map_))
    std::cout<<"INFO: proc0map and datamap are identical."<<std::endl;
  // check if the data is distributed over several processors
  bool isdistributed = (data->DistributedGlobal());
#endif

  //------------------------------------------------------
  // each processor provides its result values for proc 0
  //------------------------------------------------------

  Teuchos::RCP<Epetra_Map> proc0datamap;
  proc0datamap = LINALG::AllreduceEMap(*epetradatamap,0);

  // contract result values on proc0 (proc0 gets everything, other procs empty)
  Epetra_Import proc0dataimporter(*proc0datamap,*epetradatamap);
  Teuchos::RCP<Epetra_Vector> proc0data = Teuchos::rcp(new Epetra_Vector(*proc0datamap));
  int err = proc0data->Import(*data,proc0dataimporter,Insert);
  if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

  const Epetra_BlockMap& finaldatamap = proc0data->Map();

  //------------------------------------------------------------------
  // each processor provides its dof global id information for proc 0
  //------------------------------------------------------------------

  Teuchos::RCP<Epetra_MultiVector> dofgidperelementlid = Teuchos::rcp(new Epetra_MultiVector(*elementmap,numdof));
  dofgidperelementlid->PutScalar(-1.0);

  const int nummyelem = elementmap->NumMyElements();
  for (int idof=0; idof<numdof; ++idof)
  {
    for (int ielem=0; ielem<nummyelem; ielem++)
    {
      DRT::Element* n = dis->lRowElement(ielem);
      const double dofgid = (double) dis->Dof(n, from + idof);
      if (dofgid > -1.0)
      {
        dofgidperelementlid->ReplaceMyValue(ielem, idof, dofgid);
      }
      else
      {
        dserror("Error while creating Epetra_MultiVector dofgidperlocalnodeid");
      }
    }
  }

  // contract Epetra_MultiVector on proc0 (proc0 gets everything, other procs empty)
  Teuchos::RCP<Epetra_MultiVector> dofgidperelementlid_proc0 = Teuchos::rcp(new Epetra_MultiVector(*proc0map_,numdof));
  Epetra_Import proc0dofimporter(*proc0map_,*elementmap);
  err = dofgidperelementlid_proc0->Import(*dofgidperelementlid,proc0dofimporter,Insert);
  if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

  const int numglobelem = elementmap->NumGlobalElements();

  //-------------------------
  // specify the element type
  //-------------------------
  // loop over the different element types present
  EleGidPerDisType::const_iterator iter;
  for (iter=eleGidPerDisType_.begin(); iter != eleGidPerDisType_.end(); ++iter)
  {
    const std::string ensighteleString = GetEnsightString(iter->first);
    const int numelepertype = (iter->second).size();
    std::vector<int> actelegids(numelepertype);
    actelegids = iter->second;
    // write element type
    Write(file, ensighteleString);

    //---------------
    // write results
    //---------------
    if (myrank_==0)
    {
      if (eleGidPerDisType_.empty()==true) dserror("no element types available");
    }

    if (myrank_==0) // ensures pointer dofgids is valid
    {
      double* dofgids = (dofgidperelementlid_proc0->Values()); // columnwise data storage
      for (int idof=0; idof<numdof; ++idof)
      {
        for (int ielem=0; ielem<numelepertype; ielem++) // inode == lid of node because we use proc0map_
        {
          // local storage position of desired dof gid
          const int doflid = ielem + (idof*numglobelem);
          // get the dof global id
          const int actdofgid = (int) (dofgids[doflid]);
          dsassert(actdofgid>= 0, "error while getting dof global id");
          // get the dof local id w.r.t. the finaldatamap
          int lid = finaldatamap.LID(actdofgid);
          if (lid > -1)
          {
            Write(file, static_cast<float>((*proc0data)[lid]));
          }
          else
            dserror("received illegal dof local id: %d", lid);
        }
      }
    }// for idf

    // 2 component vectors in a 3d problem require a row of zeros.
    // do we really need this?
    if (numdof==2)
    {
      for (int ielem=0; ielem<numelepertype; ielem++)
      {
        Write<float>(file, 0.);
      }
    }

  } // eledistype


  Write(file, "END TIME STEP");
  return;
}



/*----------------------------------------------------------------------*/
/*!
  \brief Write element values for one timestep

  Each element has to have the same number of dofs.
  \author gjb
  \date 01/08
*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteElementResultStep(
  std::ofstream& file,
  PostResult& result,
  std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
  const std::string& groupname,
  const std::string& name,
  const int numdf,
  const int from
  )
{
  const Teuchos::RCP<Epetra_MultiVector> data = result.read_multi_result(groupname);
  WriteElementResultStep(file, data, resultfilepos, groupname, name, numdf, from);
}



/*----------------------------------------------------------------------*/
/*!
  \brief Write element values for one timestep

  Each element has to have the same number of dofs.
  \author gjb
  \date 01/08
*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteElementResultStep(
  std::ofstream& file,
  const Teuchos::RCP<Epetra_MultiVector>& data,
  std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
  const std::string &groupname,
  const std::string &name,
  const int numdf,
  const int from
  )
{

  //-------------------------------------------
  // write some key words and read result data
  //-------------------------------------------
  std::vector<std::ofstream::pos_type>& filepos = resultfilepos[name];
  Write(file, "BEGIN TIME STEP");
  filepos.push_back(file.tellp());
  Write(file, "description");
  Write(file, "part");
  Write(file, field_->field_pos()+1);

  const Epetra_BlockMap& datamap = data->Map();
  const int numcol = data->NumVectors();

  // do stupid conversion into Epetra map
  Teuchos::RCP<Epetra_Map> epetradatamap;
  epetradatamap = Teuchos::rcp(new Epetra_Map(datamap.NumGlobalElements(),
                                     datamap.NumMyElements(),
                                     datamap.MyGlobalElements(),
                                     0,
                                     datamap.Comm()));

  //------------------------------------------------------
  // each processor provides its result values for proc 0
  //------------------------------------------------------

  Teuchos::RCP<Epetra_Map> proc0datamap;
  proc0datamap = LINALG::AllreduceEMap(*epetradatamap,0);

  // contract result values on proc0 (proc0 gets everything, other procs empty)
  Epetra_Import proc0dataimporter(*proc0datamap,*epetradatamap);
  Teuchos::RCP<Epetra_MultiVector> proc0data = Teuchos::rcp(new Epetra_MultiVector(*proc0datamap,numcol));
  int err = proc0data->Import(*data,proc0dataimporter,Insert);
  if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

  const Epetra_BlockMap& finaldatamap = proc0data->Map();

  //-------------------------
  // specify the element type
  //-------------------------
  if (myrank_==0)
  {
    if (eleGidPerDisType_.empty()==true) dserror("no element types available");
  }
  // loop over the different element types present
  EleGidPerDisType::const_iterator iter;
  for (iter=eleGidPerDisType_.begin(); iter != eleGidPerDisType_.end(); ++iter)
  {
    const std::string ensighteleString = GetEnsightString(iter->first);

    int numsubele = GetNumSubEle(iter->first);

    const int numelepertype = (iter->second).size();
    std::vector<int> actelegids(numelepertype);
    actelegids = iter->second;
    // write element type
    Write(file, ensighteleString);

    //---------------
    // write results
    //---------------
    if (myrank_==0)
    {
      if (numdf+from > numcol) dserror("violated column range of Epetra_MultiVector: %d",numcol);
      for (int col=0; col<numdf; ++col)
      {
        //extract actual column
        Epetra_Vector* datacolumn = (*proc0data)(col+from);

        for (int iele=0; iele<numelepertype; iele++)
        {
          // extract element global id
          const int gid = actelegids[iele];
          // get the dof local id w.r.t. the finaldatamap
          //int lid = datamap.LID(gid);
          int lid = finaldatamap.LID(gid);
          if (lid > -1)
          {
            for(int i=0;i<numsubele;++i)
              Write(file, static_cast<float>((*datacolumn)[lid]));
          }
          else
            dserror("received illegal dof local id: %d", lid);
        }
      }
    } // if (myrank_==0)

    // 2 component vectors in a 3d problem require a row of zeros.
    if (numdf==2)
    {
      for (int iele=0; iele<numelepertype; iele++)
      {
        Write<float>(file, 0.);
      }
    }

  } // end iteration over eleGidPerDisType_;

    // finish writing the current time step
  Write(file, "END TIME STEP");
  return;

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteIndexTable(
  std::ofstream& file,
  const std::vector<std::ofstream::pos_type>& filepos
  ) const
{
  std::ofstream::pos_type lastpos = file.tellp();
  const unsigned steps = filepos.size();
  Write(file, steps);
  for (unsigned i=0; i<steps; ++i)
  {
    Write<long>(file, filepos[i]);
  }
  Write(file, 0);
  Write<long>(file, lastpos);
  Write(file, "FILE_INDEX");
  return;
}


/*----------------------------------------------------------------------*/
/*!
  \brief write strings of exactly 80 chars
*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteString(
  std::ofstream& stream,
  const std::string str) const
{
  // we need to write 80 bytes per string
  std::vector<char> s(str.begin(), str.end());
  while (s.size()<80)
  {
    s.push_back('\0');
  }
  stream.write(&s[0], 80);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::string EnsightWriter::GetVariableSection(
  std::map<std::string,std::vector<int> > filesetmap,
  std::map<std::string,int>               variablenumdfmap,
  std::map<std::string,std::string>       variablefilenamemap
) const
{
  std::stringstream str;

  std::map<std::string,int>::const_iterator variable;

  for (variable = variablenumdfmap.begin(); variable != variablenumdfmap.end(); ++variable)
  {
    const std::string key = variable->first;
    const int numdf = variable->second;
    const std::string filename = variablefilenamemap[key];

    // Get rid of path
    const size_t found_path = filename.find_last_of("/\\");
    const std::string filename_nopath = filename.substr(found_path+1);

    std::map<std::string,int>::const_iterator timeentry = timesetnumbermap_.find(key);
    if (timeentry == timesetnumbermap_.end())
      dserror("key not found!");
    const int timesetnumber = timeentry->second;

    std::map<std::string,int>::const_iterator entry1 = filesetnumbermap_.find(key);
    if (entry1 == filesetnumbermap_.end())
      dserror("key not found!");
    const int setnumber = entry1->second;

    std::map<std::string,std::vector<int> >::const_iterator entry2 = filesetmap.find(key);
    if (entry2 == filesetmap.end())
      dserror("filesetmap not defined for '%s'", key.c_str());

    const int numsubfilesteps = entry2->second.size();
    std::string filename_for_casefile;
    if (numsubfilesteps > 1)
    {
      filename_for_casefile = filename_nopath + "***";
    }
    else
    {
      filename_for_casefile = filename_nopath;
    }
    str << GetVariableEntryForCaseFile(numdf, setnumber, key, filename_for_casefile,timesetnumber);
  }

  return str.str();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::string EnsightWriter::GetVariableEntryForCaseFile(
  const int numdf,
  const unsigned int fileset,
  const std::string name,
  const std::string filename,
  const int timeset
  ) const
{
  std::stringstream str;

  // determine the type of this result variable (node-/element-based)
  std::map<std::string,std::string>::const_iterator entry = variableresulttypemap_.find(name);
  if (entry == variableresulttypemap_.end())
    dserror("key not found!");
  const std::string restypestring = entry->second;

  // create variable entry in the case-file
  switch (numdf)
  {
  case 9:
    str << "tensor asymm per "<<restypestring<<":\t"<< timeset <<"\t"<< fileset << "\t"<< name << "\t"<< filename << "\n";
    break;
  case 6:
    str << "tensor symm per "<<restypestring<<":\t"<< timeset <<"\t"<< fileset << "\t"<< name << "\t"<< filename << "\n";
    break;
  case 3:
  case 2:
    str << "vector per "<<restypestring<<":\t"<< timeset <<"\t"<< fileset << "\t"<< name << "\t"<< filename << "\n";
    break;
  case 1:
    str << "scalar per "<<restypestring<<":\t"<< timeset <<"\t"<< fileset << "\t"<< name << "\t"<< filename << "\n";
    break;
  default:
    dserror("unknown number of dof per node");
  };
  return str.str();
}


/*----------------------------------------------------------------------*/
/*!
  \brief create string for one TIME set in the case file
*/
/*----------------------------------------------------------------------*/
std::string EnsightWriter::GetTimeSectionString(
  const int timeset,
  const std::vector<double>& times
  ) const
{
  std::stringstream s;
  s << "time set:\t\t" << timeset << "\n"<< "number of steps:\t"<< times.size() << "\ntime values: ";
  for (unsigned i=0; i<times.size(); ++i)
  {
    s << std::setprecision(16) << times[i]<< " ";
    if (i%8==0&& i!=0)
    {
      s << "\n";
    }
  }
  s << "\n";
  return s.str();
}


/*----------------------------------------------------------------------*/
/*!
  \brief create string for the TIME section in the case file
*/
/*----------------------------------------------------------------------*/
std::string EnsightWriter::GetTimeSectionStringFromTimesets(
  const std::map<std::string,std::vector<double> >& timesetmap
  ) const
{
  std::stringstream s;
  std::map<std::string,std::vector<double> >::const_iterator timeset;
  std::set<int> donetimesets;

  for (timeset = timesetmap.begin(); timeset != timesetmap.end(); ++timeset)
  {
    std::string key = timeset->first;
    std::map<std::string,int>::const_iterator entry = timesetnumbermap_.find(key);
    if (entry == timesetnumbermap_.end())
      dserror("key not found!");
    const int timesetnumber = entry->second;
    const std::vector<double> soltimes = timeset->second;
    if (donetimesets.find(timesetnumber)==donetimesets.end()) // do not write redundant time sets
    {
      donetimesets.insert(timesetnumber);
      std::string outstring = GetTimeSectionString(timesetnumber, soltimes);
      s<< outstring<<std::endl;
    }
  }
  return s.str();
}

/*----------------------------------------------------------------------*/
/*!
  \brief create string for the FILE section in the case file
*/
/*----------------------------------------------------------------------*/
std::string EnsightWriter::GetFileSectionStringFromFilesets(
  const std::map<std::string,std::vector<int> >& filesetmap
  ) const
{
  std::stringstream s;
  std::map<std::string,std::vector<int> >::const_iterator fileset;

  // print filesets in increasing numbering, starting with "geo"

  std::map<std::string,int>::const_iterator entry = filesetnumbermap_.find("geo");
  if (entry == filesetnumbermap_.end())
    dserror("key 'geo' not found!");
  const int setnumber = entry->second;
  if (setnumber != 1) dserror ("geometry file must have file set number 1");
  fileset= filesetmap.find("geo");
  std::vector<int> stepsperfile = fileset->second;
  s << "file set:\t\t"<< setnumber << "\n";
  if (stepsperfile.size() == 1)
  {
    s << "number of steps:\t"<< stepsperfile[0] << "\n\n";
  }
  else
  {
    for (unsigned int j = 0; j < stepsperfile.size(); ++j)
    {
      s << "filename index:\t"<< 1+j << "\n";
      s << "number of steps:\t"<< stepsperfile[j] << "\n";
    }
    s << "\n";
  }

  for (fileset = filesetmap.begin(); fileset != filesetmap.end(); ++fileset)
  {
    std::string key = fileset->first;
    // skip geometry file since it was already considered above!
    if (key=="geo") continue;

    std::map<std::string,int>::const_iterator entry = filesetnumbermap_.find(key);
    if (entry == filesetnumbermap_.end())
      dserror("key not found!");
    const int setnumber = entry->second;
    std::vector<int> stepsperfile = fileset->second;
    s << "file set:\t\t"<< setnumber << "\n";
    if (stepsperfile.size() == 1)
    {
      s << "number of steps:\t"<< stepsperfile[0] << "\n\n";
    }
    else
    {
      for (unsigned int j = 0; j < stepsperfile.size(); ++j)
      {
        s << "filename index:\t"<< 1+j << "\n";
        s << "number of steps:\t"<< stepsperfile[j] << "\n";
      }
      s << "\n";
    }
  }
  return s.str();
}

/*----------------------------------------------------------------------*/
/*
    Write the coordinates for a Polynomial discretization
    The ccordinates of the vizualisation points (i.e. the corner
    nodes of elements displayed in paraview) are just the node
    coordinates of the nodes in the discretization.
*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteCoordinatesForPolynomialShapefunctions
(
  std::ofstream&                      geofile ,
  const Teuchos::RCP<DRT::Discretization> dis     ,
  Teuchos::RCP<Epetra_Map>&               proc0map
  )
{

  // refcountpointer to vector of all coordinates
  // distributed among all procs
  Teuchos::RCP<Epetra_MultiVector> nodecoords;

  const int NSD = 3; // number of space dimensions

  const Epetra_Map* nodemap = dis->NodeRowMap();
  const int numnp = nodemap->NumMyElements();
  nodecoords = Teuchos::rcp(new Epetra_MultiVector(*nodemap,3));

  // loop over the nodes on this proc and store the coordinate information
  for (int inode=0; inode<numnp; inode++)
  {
    int gid = nodemap->GID(inode);
    const DRT::Node* actnode = dis->gNode(gid);
    for (int isd=0; isd<NSD; ++isd)
    {
      double val = ((actnode->X())[isd]);
      nodecoords->ReplaceMyValue(inode, isd, val);
    }
  }

  // put all coordinate information on proc 0
  proc0map = LINALG::AllreduceEMap(*nodemap,0);

  // import my new values (proc0 gets everything, other procs empty)
  Epetra_Import proc0importer(*proc0map,*nodemap);
  Teuchos::RCP<Epetra_MultiVector> allnodecoords = Teuchos::rcp(new Epetra_MultiVector(*proc0map,3));
  int err = allnodecoords->Import(*nodecoords,proc0importer,Insert);
  if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

  // write the node coordinates (only proc 0)
  // ensight format requires x_1 .. x_n, y_1 .. y_n, z_1 ... z_n
  // this is fulfilled automatically due to Epetra_MultiVector usage (columnwise writing data)
  if (myrank_==0)
  {
    double* coords = allnodecoords->Values();
    int numentries = (3*(allnodecoords->GlobalLength()));
    dsassert(numentries == (3*nodemap->NumGlobalElements()),"proc 0 has not all of the node coordinates");
    if (nodeidgiven_)
    {
      // first write node global ids (default)
      for (int inode=0; inode<proc0map->NumGlobalElements(); ++inode)
      {
        Write(geofile,proc0map->GID(inode)+1);
        // gid+1 delivers the node numbering of the *.dat file starting with 1
      }
    }
    // now write the coordinate information
    for (int i=0; i<numentries; ++i)
    {
      Write(geofile, static_cast<float>(coords[i]));
    }
  }

  return;
}
