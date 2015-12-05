/*!----------------------------------------------------------------------
\file statmech_bilayer_gmsh.cpp
\brief gmsh output methods for statistical mechanics

<pre>
Maintainer: Dhrubajyoti Mukherjee
            mukherjee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15270
</pre>

*----------------------------------------------------------------------*/

#include "statmech_manager_bilayer.H"
#include "../drt_inpar/inpar_statmech.H"

#include <Teuchos_Time.hpp>
#include <iostream>
#include <string>
#include <sstream>

#include "../linalg/linalg_utils.H"
#include "../drt_discsh3/discsh3.H"
//#include "../drt_lib/drt_discret_faces.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 | writing Gmsh data for current step            (public)mukherjee 09/15|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::GmshOutput(const Epetra_Vector&          disrow,
                                           const std::ostringstream&            filename,
                                           const int&                           step,
                                           const double&                        time)
{
  /*the following method writes output data for Gmsh into file with name "filename"; all line elements are written;
   * the nodal displacements are handed over in the variable "dis"; note: in case of parallel computing only
   * processor 0 writes; it is assumed to have a fully overlapping column map and hence all the information about
   * all the nodal position; parallel output is now possible with the restriction that the nodes(processors) in question
   * are of the same machine*/


  //we need displacements also of ghost nodes and hence export displacement vector to column map format
  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(disrow, discol);

  // disrow.Print(std::cout);
  //  discol.Print(std::cout);

  // do output to file in c-style
  FILE* fp = NULL;

  // first processor starts by opening the file and writing the header, other processors have to wait
  if (!discret_->Comm().MyPID())
  {
    //open file to write output data into
    fp = fopen(filename.str().c_str(), "w");

    // write output to temporary std::stringstream;
    std::stringstream gmshfileheader;
    //      the beginning of the stream is "View \"" to indicate Gmsh that the following data is in order to create an image and
    //       * this command is followed by the name of that view displayed during it's shown in the video; in the following example
    //       * this name is for the 100th time step: Step00100; then the data to be presented within this view is written within { ... };
    //       * in the following example this data consists of scalar lines defined by the coordinates of their end points
    gmshfileheader <<"General.BackgroundGradient = 0;\n";
    gmshfileheader <<"View.LineType = 1;\n";
    gmshfileheader <<"View.LineWidth = 1.4;\n";
    gmshfileheader <<"View.PointType = 1;\n";
    gmshfileheader <<"View.PointSize = 3;\n";
    gmshfileheader <<"General.Antialiasing = 1;\n";
    gmshfileheader <<"General.ColorScheme = 1;\n";
    gmshfileheader <<"General.Color.Background = {255,255,255};\n";
    gmshfileheader <<"General.Color.Foreground = {255,255,255};\n";
    //gmshfileheader <<"General.Color.BackgroundGradient = {128,147,255};\n";
    gmshfileheader <<"General.Color.Foreground = {85,85,85};\n";
    gmshfileheader <<"General.Color.Text = {0,0,0};\n";
    gmshfileheader <<"General.Color.Axes = {0,0,0};\n";
    gmshfileheader <<"General.Color.SmallAxes = {0,0,0};\n";
    gmshfileheader <<"General.Color.AmbientLight = {25,25,25};\n";
    gmshfileheader <<"General.Color.DiffuseLight = {255,255,255};\n";
    gmshfileheader <<"General.Color.SpecularLight = {255,255,255};\n";
    gmshfileheader <<"View.ColormapAlpha = 1;\n";
    gmshfileheader <<"View.ColormapAlphaPower = 0;\n";
    gmshfileheader <<"View.ColormapBeta = 0;\n";
    gmshfileheader <<"View.ColormapBias = 0;\n";
    gmshfileheader <<"View.ColormapCurvature = 0;\n";
    gmshfileheader <<"View.ColormapInvert = 0;\n";
    gmshfileheader <<"View.ColormapNumber = 2;\n";
    gmshfileheader <<"View.ColormapRotation = 0;\n";
    gmshfileheader <<"View.ColormapSwap = 0;\n";
    gmshfileheader <<"View.ColorTable = {Black,Yellow,Blue,Orange,Red,Cyan,Purple,Brown,Green};\n";
    gmshfileheader << "View \" Step " << step << " \" {" << std::endl;

    //write content into file and close it
    fputs(gmshfileheader.str().c_str(), fp);
    fclose(fp);
  }
  //
  // wait for all processors to arrive at this point
  discret_->Comm().Barrier();


  //! \brief vector holding the NodeIDs of new faces
  std::vector<std::vector<int> >  new_faces_NodeID;
  new_faces_NodeID.clear();
  if(DRT::INPUT::IntegralValue<int>(statmechBilayerparams_,"BONDFLIP"))
  {
    std::map<int,std::vector<int> >::iterator it;
    for(it = new_faces_NodeID_.begin(); it != new_faces_NodeID_.end(); it++ )
    {
      std::vector<int> NewNodes; NewNodes.clear();
      NewNodes.push_back(it->second.at(0)); NewNodes.push_back(it->second.at(1));
      new_faces_NodeID.push_back(NewNodes);
    }
  }


  // loop over the participating processors each of which appends its part of the output to one output file
  for (int proc = 0; proc < discret_->Comm().NumProc(); proc++)
  {
//  facediscret_->PrintFaces(std::cout);
//  if (!discret_->Comm().MyPID())
    if (discret_->Comm().MyPID() == proc)
    {
      //open file again to append ("a") output data into
      fp = fopen(filename.str().c_str(), "a");
      // write output to temporary std::stringstream;
      std::stringstream gmshfilecontent;

      //looping through all elements on the processor
      for (int i=0; i<discret_->NumMyRowElements(); ++i)
      {
        // coordinates of nodes or binding spots (depending on element)
        LINALG::SerialDenseMatrix coord(3, 3, true); // 3 Nodes
        //getting pointer to current element
        DRT::Element* element = discret_->lRowElement(i);
        const DRT::ElementType & eot = element->ElementType();

        if(eot == DRT::ELEMENTS::DiscSh3Type::Instance())    // discrete shell element
        {
          for (int jd=0; jd<element->NumNode(); jd++)
          {
            DRT::Node * Node1 = NULL;
            DRT::Node * Node2 = NULL;
            Node1 = (element->Nodes())[jd];
            Node2 = (element->Nodes())[jd+1];
            if (jd==2)
            {
              Node1 = (element->Nodes())[2];
              Node2 = (element->Nodes())[0];
            }

            LINALG::Matrix <1,3> Node1Curr = GetSpatialPosition(Node1,discol);
            LINALG::Matrix <1,3> Node2Curr = GetSpatialPosition(Node2,discol);

            std::vector<int> PossibleLine1; PossibleLine1.clear();
            PossibleLine1.push_back(Node1->Id());   // connectivity Node1-Node2
            PossibleLine1.push_back(Node2->Id());

            std::vector<int> PossibleLine2; PossibleLine2.clear();
            PossibleLine2.push_back(Node2->Id());   // connectivity Node2-Node1
            PossibleLine2.push_back(Node1->Id());

            std::vector<std::vector<int> >::iterator it;
            it=find(new_faces_NodeID.begin(),new_faces_NodeID.end(),PossibleLine1);
            if(it!=new_faces_NodeID.end())
            {
              const double color_newbond=0.15; // yellow
              PrintGmshNewBondsToStream(gmshfilecontent,Node1Curr,Node2Curr,color_newbond);
            }
            else
            {
              std::vector<std::vector<int> >::iterator it_secondary;
              it_secondary=find(new_faces_NodeID.begin(),new_faces_NodeID.end(),PossibleLine2);
              if(it_secondary!=new_faces_NodeID.end())
              {
                const double color_newbond=0.15; // yellow
                PrintGmshNewBondsToStream(gmshfilecontent,Node1Curr,Node2Curr,color_newbond);
              }
              else
              {
                const double color_newbond=0.75; // purple
                PrintGmshNewBondsToStream(gmshfilecontent,Node1Curr,Node2Curr,color_newbond);
              }
            }
          }
        }

      } // discret_->NumMyColElements()


      //write content into file and close it (this way we make sure that the output is written serially)
      fputs(gmshfilecontent.str().c_str(), fp);
      fclose(fp);
    }
    discret_->Comm().Barrier();
  } // proc < discret_->Comm().NumProc()


  discret_->Comm().Barrier();

  // finish data section of this view by closing curly brackets
  if (discret_->Comm().MyPID() == 0)
  {
    fp = fopen(filename.str().c_str(), "a");
    std::stringstream gmshfileend;

    // add black dot for correct coloring...
//    if(periodlength_->at(0)==0.0)
      gmshfileend << "SP("<<periodlength_->at(0)<<",0,"<<periodlength_->at(2)<<"){0,0};"<<std::endl;
    // add green dot for correct coloring...
      gmshfileend << "SP("<<periodlength_->at(0)<<",0,"<<periodlength_->at(2)<<"){1,1};"<<std::endl;
    gmshfileend << "};" << std::endl;
    fputs(gmshfileend.str().c_str(), fp);

    fclose(fp);
  }
  // return simultaneously (not sure if really needed)
  discret_->Comm().Barrier();

  return;
} // STATMECH::StatMechManager::GmshOutput()


/*----------------------------------------------------------------------*
 |  print Gmsh Triangle to stringstream (private)        mukherjee 09/15|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::PrintGmshTriangleToStream(std::stringstream& gmshfilecontent,
                                                                 LINALG::SerialDenseMatrix coord)
{
  // "ST" is scalar triangle, followed by 3x coordinates (x,y,z) of vertices and color
  double color= 0.6;
  gmshfilecontent << "ST("<< std::scientific;
  gmshfilecontent << coord(0,0) << "," << coord(1,0) << "," << coord(2,0) << ",";
  gmshfilecontent << coord(0,1) << "," << coord(1,1) << "," << coord(2,1) << ",";
  gmshfilecontent << coord(0,2) << "," << coord(1,2) << "," << coord(2,2);
  gmshfilecontent << "){" << std::scientific;
  gmshfilecontent << color << "," << color << "," << color  << "};" << std::endl << std::endl;

  return;
}


/*----------------------------------------------------------------------*
 |  print Gmsh Line to stringstream (private)            mukherjee 09/15|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManagerBilayer::PrintGmshNewBondsToStream(std::stringstream& gmshfilecontent,
                                                             LINALG::Matrix <1,3> & MasterNodeSpatial,
                                                             LINALG::Matrix <1,3> & SlaveNodeSpatial,
                                                             const double color)
{
  // "ST" is scalar triangle, followed by 3x coordinates (x,y,z) of vertices and color

  // 1st line of Triangle
  gmshfilecontent << "SL(" << std::scientific;
  gmshfilecontent << MasterNodeSpatial(0) << "," << MasterNodeSpatial(1) << "," << MasterNodeSpatial(2) << "," << SlaveNodeSpatial(0) << "," << SlaveNodeSpatial(1) << "," << SlaveNodeSpatial(2);
  gmshfilecontent << ")" << "{" << std::scientific << color << ","<< color << "};" << std::endl;

  return;
}

