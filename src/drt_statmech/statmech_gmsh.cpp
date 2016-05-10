/*!----------------------------------------------------------------------
\file statmech_gmsh.cpp
\brief gmsh output methods for statistical mechanics

\maintainer Kei Müller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276

*----------------------------------------------------------------------*/

#include "statmech_manager.H"

#include <Teuchos_Time.hpp>

#include "../linalg/linalg_utils.H"
#include "../drt_beamcontact/beam3contact_manager.H"
#include "../drt_beam3/beam3.H"
#include "../drt_beam3/beam3r.H"
#include "../drt_beam3/beam3eb.H"
#include "../drt_truss3/truss3.H"
#include "../drt_spring3/spring3.H"
#include "../drt_truss3cl/truss3cl.H"
#include "../drt_beam3/beam3cl.H"
#include "../drt_torsion3/torsion3.H"
#include "../drt_rigidsphere/rigidsphere.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 | writing Gmsh data for current step                 public)cyron 01/09|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GmshOutput(const Epetra_Vector&                 disrow,
                                           const std::ostringstream&            filename,
                                           const int&                           step,
                                           const double&                        time,
                                           Teuchos::RCP<CONTACT::Beam3cmanager> beamcmanager)
{
  /*the following method writes output data for Gmsh into file with name "filename"; all line elements are written;
   * the nodal displacements are handed over in the variable "dis"; note: in case of parallel computing only
   * processor 0 writes; it is assumed to have a fully overlapping column map and hence all the information about
   * all the nodal position; parallel output is now possible with the restriction that the nodes(processors) in question
   * are of the same machine*/
  GmshPrepareVisualization(disrow);

  //we need displacements also of ghost nodes and hence export displacement vector to column map format
  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(disrow, discol);

  // do output to file in c-style
  FILE* fp = NULL;

  //number of solid elements by which a round line is depicted
  const int nline = 16;

  // first processor starts by opening the file and writing the header, other processors have to wait
  if (!discret_->Comm().MyPID())
  {
    //open file to write output data into
    fp = fopen(filename.str().c_str(), "w");

    // write output to temporary std::stringstream;
    std::stringstream gmshfileheader;
    /*the beginning of the stream is "View \"" to indicate Gmsh that the following data is in order to create an image and
     * this command is followed by the name of that view displayed during it's shown in the video; in the following example
     * this name is for the 100th time step: Step00100; then the data to be presented within this view is written within { ... };
     * in the following example this data consists of scalar lines defined by the coordinates of their end points*/
    gmshfileheader <<"General.BackgroundGradient = 0;\n";
    gmshfileheader <<"View.LineType = 1;\n";
    gmshfileheader <<"View.LineWidth = 1.4;\n";
    gmshfileheader <<"View.PointType = 1;\n";
    gmshfileheader <<"View.PointSize = 3;\n";
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

  // wait for all processors to arrive at this point
  discret_->Comm().Barrier();

  // loop over the participating processors each of which appends its part of the output to one output file
  for (int proc = 0; proc < discret_->Comm().NumProc(); proc++)
  {
    if (discret_->Comm().MyPID() == proc)
    {
      //open file again to append ("a") output data into
      fp = fopen(filename.str().c_str(), "a");
      // write output to temporary std::stringstream;
      std::stringstream gmshfilecontent;

      //looping through all elements on the processor
      for (int i=0; i<discret_->NumMyColElements(); ++i)
      {
        // coordinates of nodes or binding spots (depending on element)
        LINALG::SerialDenseMatrix coord(3, 2, true);
        //getting pointer to current element
        DRT::Element* element = discret_->lColElement(i);
        const DRT::ElementType & eot = element->ElementType();
//        bool PolarityConditions=true;
        bool IsBarbed=false; // Condition to indicate Barbed ends
//        bool IsPointed=false; // Condition to indicate Pointed ends

        // interpolated beam element
        if(eot==DRT::ELEMENTS::BeamCLType::Instance())
        {
          DRT::ELEMENTS::BeamCL* currele = NULL;
          currele = dynamic_cast<DRT::ELEMENTS::BeamCL*> (discret_->gElement(discret_->ElementColMap()->GID(i)));
          // Alternative way to get nodal positions
          DRT::Element* filelement = discret_->gElement(discret_->ElementColMap()->GID(i));
          DRT::Node* node0 = NULL;
          DRT::Node* node1 = NULL;
          std::vector<double> position(6);
          std::vector<int> dofnode0;
          std::vector<int> dofnode1;

          std::vector<LINALG::Matrix<1,2> > Ibp(2);
          for(int filament=0; filament<2; filament++)
             DRT::UTILS::shape_function_1D(Ibp[filament],currele->MyBindingPosition()[filament],currele->Shape());

          // determine positions of the interpolated nodes
          for(int ifil=0; ifil<coord.N(); ifil++)
          {
            node0 = discret_->gNode(filelement->NodeIds()[0+2*ifil]);
            node1 = discret_->gNode(filelement->NodeIds()[1+2*ifil]);

            dofnode0 = discret_->Dof(node0);
            dofnode1 = discret_->Dof(node1);
            position[0] = node0->X()[0] + discol[discret_->DofColMap()->LID(dofnode0[0])];
            position[1] = node0->X()[1] + discol[discret_->DofColMap()->LID(dofnode0[1])];
            position[2] = node0->X()[2] + discol[discret_->DofColMap()->LID(dofnode0[2])];
            position[3] = node1->X()[0] + discol[discret_->DofColMap()->LID(dofnode1[0])];
            position[4] = node1->X()[1] + discol[discret_->DofColMap()->LID(dofnode1[1])];
            position[5] = node1->X()[2] + discol[discret_->DofColMap()->LID(dofnode1[2])];

            //shift nodes back in case they have been shifted due to PBC
            UnshiftPositions(position);
            for(int id=0;id<3;id++)
             coord(id, ifil) = Ibp[ifil](0)*position[id]+Ibp[ifil](1)*position[id+3];
            //shift Crosslinker intperpolated nodal positions into Periodic Boundary Domain
            for(int fil=0;fil<2;fil++)
              for(int j=0;j<3;j++)
              {  if(coord(j,fil) > periodlength_->at(j))
                  coord(j,fil) -= periodlength_->at(j)*floor(coord(j,fil)/periodlength_->at(j));
                 if(coord(j,fil) < 0.0)
                   coord(j,fil) -= periodlength_->at(j)*floor(coord(j,fil)/periodlength_->at(j));
              }
          }
        }
        // interpolated truss-beam element
        else if(eot==DRT::ELEMENTS::Truss3CLType::Instance())
        {
          DRT::ELEMENTS::Truss3CL* currele = NULL;
          currele = dynamic_cast<DRT::ELEMENTS::Truss3CL*> (discret_->gElement(discret_->ElementColMap()->GID(i)));
          // Alternative way to get nodal positions
          DRT::Element* filelement = discret_->gElement(discret_->ElementColMap()->GID(i));
          // length of first element and second element required for hermite shape function interpolation
          double LengthofFilamentatRef;
          LINALG::Matrix<1,4>  Ibp;         //Hermite shape functions for beam3EB element
          std::vector<double> Tcurr_(6); // Vector containing tangent of Nodal Positions.
          DRT::Node* node0 = NULL;
          DRT::Node* node1 = NULL;
          std::vector<double> position(6);
          std::vector<double> InitialPosition(6);
          std::vector<int> dofnode0;
          std::vector<int> dofnode1;

          // determine positions of the interpolated nodes
          for(int ifil=0; ifil<coord.N(); ifil++)
          {
            node0 = discret_->gNode(filelement->NodeIds()[0+2*ifil]);
            node1 = discret_->gNode(filelement->NodeIds()[1+2*ifil]);

            dofnode0 = discret_->Dof(node0);
            dofnode1 = discret_->Dof(node1);

            InitialPosition[0] = node0->X()[0];
            InitialPosition[1] = node0->X()[1];
            InitialPosition[2] = node0->X()[2];
            InitialPosition[3] = node1->X()[0];
            InitialPosition[4] = node1->X()[1];
            InitialPosition[5] = node1->X()[2];

            // NodeShift Initial positions
            UnshiftPositions(InitialPosition);

            position[0] = node0->X()[0] + discol[discret_->DofColMap()->LID(dofnode0[0])];
            position[1] = node0->X()[1] + discol[discret_->DofColMap()->LID(dofnode0[1])];
            position[2] = node0->X()[2] + discol[discret_->DofColMap()->LID(dofnode0[2])];
            position[3] = node1->X()[0] + discol[discret_->DofColMap()->LID(dofnode1[0])];
            position[4] = node1->X()[1] + discol[discret_->DofColMap()->LID(dofnode1[1])];
            position[5] = node1->X()[2] + discol[discret_->DofColMap()->LID(dofnode1[2])];

            //Tangents at reference configuration
            std::vector<LINALG::Matrix<3,1> > Tref_(2);
            double norm2=0;
            for(int node = 0; node<2 ; node++)
            {
              for(int dof = 0; dof< 3 ; dof++ )
              {
                Tref_[node](dof) = InitialPosition[dof+3]- InitialPosition[dof];
              }
              norm2 = Tref_[node].Norm2();
              LengthofFilamentatRef=norm2;
              Tref_[node].Scale(1/norm2);
            }

            Tcurr_[0] = Tref_[0](0) + discol[discret_->DofColMap()->LID(dofnode0[3])];
            Tcurr_[1] = Tref_[0](1) + discol[discret_->DofColMap()->LID(dofnode0[4])];
            Tcurr_[2] = Tref_[0](2) + discol[discret_->DofColMap()->LID(dofnode0[5])];
            Tcurr_[3] = Tref_[1](0) + discol[discret_->DofColMap()->LID(dofnode1[3])];
            Tcurr_[4] = Tref_[1](1) + discol[discret_->DofColMap()->LID(dofnode1[4])];
            Tcurr_[5] = Tref_[1](2) + discol[discret_->DofColMap()->LID(dofnode1[5])];

            //shift nodes back in case they have been shifted due to PBC
            UnshiftPositions(position);

            DRT::UTILS::shape_function_hermite_1D(Ibp,currele->MyBindingPosition()[ifil],LengthofFilamentatRef,currele->Shape());
            for(int id=0;id<3;id++)
             coord(id, ifil) = Ibp(0)*position[id]+Ibp(1)*Tcurr_[id]+Ibp(2)*position[id+3]+Ibp(3)*Tcurr_[id+3];
            //shift Crosslinker intperpolated nodal positions into Periodic Boundary Domain
            for(int fil=0;fil<2;fil++)
              for(int j=0;j<3;j++)
              {  if(coord(j,fil) > periodlength_->at(j))
                  coord(j,fil) -= periodlength_->at(j)*floor(coord(j,fil)/periodlength_->at(j));
                 if(coord(j,fil) < 0.0)
                   coord(j,fil) -= periodlength_->at(j)*floor(coord(j,fil)/periodlength_->at(j));
              }
          }
        }
        else if(eot==DRT::ELEMENTS::Beam3ebType::Instance())  // Kirchhoff Type of beam elements
        {
          // this cast is necessary in order to use the method ->Tref()
          const DRT::ELEMENTS::Beam3eb* ele = dynamic_cast<const DRT::ELEMENTS::Beam3eb*>(element);
          //number of interpolated points in axial direction
          const int n_axial= statmechparams_.get<int>("GMSHNINTPT",10.0);
          // prepare storage for nodal coordinates
          int nnodes = element->NumNode();
          LINALG::SerialDenseMatrix nodalcoords(3,nnodes);
          LINALG::SerialDenseMatrix nodaltangents(3,nnodes);
          // Change the size of coord for interpolated coordinates
          coord.Reshape(3,n_axial);
          // Additional plot in case the filament ends are open to (de-)polymerization
          // i.e Polarity of filaments is taken into account.
          // In this case, plotting of Barbed ends and pointed ends.
//          if(PolarityConditions)
//          {
//            for (int i=0; i<(int)BarbedEnds_->size(); i++)
//            {
//              if (element->NodeIds()[0] == (int)(*BarbedEnds_)[i] || element->NodeIds()[1] == (int)(*BarbedEnds_)[i])
//                IsBarbed=true;
//            }
//          }

          // compute current nodal positions
          for (int i=0;i<3;++i)
          {
            for (int j=0;j<element->NumNode();++j)
            {
              double referenceposition = ((element->Nodes())[j])->X()[i];
              std::vector<int> dofnode = discret_->Dof((element->Nodes())[j]);
              double displacement = discol[discret_->DofColMap()->LID(dofnode[i])];
              nodalcoords(i,j) =  referenceposition + displacement;
              nodaltangents(i,j) =  ((ele->Tref())[j])(i) + discol[discret_->DofColMap()->LID(dofnode[3+i])];
            }
          }
          // Position vector for shifting broken elements to their original location
          std::vector<double> position(6);
          for(int i=0; i<3; i++)
          {
            position[i]=nodalcoords(i,0);
            position[i+3]=nodalcoords(i,1);
          }
          // Shift the broken filaments to their unbroken state or original state
          UnshiftPositions(position);

          // Position vector to break them again i.e. bring them in the periodic volume
          std::vector<double> ShiftPosition(n_axial*3);

          // Vector to interpolate the positions
          if (nnodes ==2)
          {
            LINALG::Matrix<12,1> disp_totlag(true);
            for (int i=0;i<3;i++)
            {
              disp_totlag(i)=position[i];
              disp_totlag(i+6)=position[i+3];
              disp_totlag(i+3)=nodaltangents(i,0);
              disp_totlag(i+9)=nodaltangents(i,1);
            }

            //Calculate axial positions within the element by using the Hermite interpolation of Kirchhoff beams
            for (int i=0;i<n_axial;i++)
            {
              double xi=-1.0 + i*2.0/(n_axial -1); // parameter coordinate of position vector on beam centerline
              LINALG::Matrix<3,1> r = ele->GetPos(xi, disp_totlag); //position vector on beam centerline

              for (int j=0;j<3;j++)
                ShiftPosition[3*i+j]=r(j);
            }

            // Shift the interpolated element back where they were, i.e. inside periodic volume
            ShiftPositions(ShiftPosition,n_axial);

            // Store the values in coord matrix for further plotting purposes
            for (int i=0;i<n_axial;i++)
              for (int j=0;j<3;j++)
                coord(j,i)=ShiftPosition[3*i+j];

          }
            else
              dserror("Only 2-noded Kirchhoff elements possible so far!");
        }

        else  // standard beam or truss element
        {
          for (int id=0; id<3; id++)
            for (int jd=0; jd<element->NumNode(); jd++)
            {
              double referenceposition = ((element->Nodes())[jd])->X()[id];
              std::vector<int> dofnode = discret_->Dof((element->Nodes())[jd]);
              double displacement = discol[discret_->DofColMap()->LID(dofnode[id])];
              coord(id, jd) = referenceposition + displacement;
            }
        }

        if (eot != DRT::ELEMENTS::RigidsphereType::Instance())  // beam element
        {
          //declaring variable for color of elements
          double color = 1.0;
          if (element->Id() >= basisnodes_)
          {
            if(eot != DRT::ELEMENTS::BeamCLType::Instance() &&
               eot != DRT::ELEMENTS::Torsion3Type::Instance() &&
               eot != DRT::ELEMENTS::Truss3CLType::Instance() &&
               element->NumNode()!=2)
                dserror("Crosslinker element has more than 2 nodes! No visualization has been implemented for such linkers!");
            //apply different colors for different crosslinkers
            if(crosslinkertype_!=Teuchos::null)
            {
              int crossLID = crosslinkermap_->LID((int)(*element2crosslink_)[discret_->ElementColMap()->LID(element->Id())]);
              if((*crosslinkertype_)[crossLID]==1.0) // active linker
              {
                color = 0.75;
                //std::cout<<"Proc "<<discret_->Comm().MyPID()<<" - SL active "<<crosslinkermap_->GID(crossLID)<<", element "<<element->Id()<<", Type "<<(*crosslinkertype_)[crossLID]<<std::endl;
              }
              else // standard single species crosslinker
              {
                color = 0.5;
                //std::cout<<"  Proc "<<discret_->Comm().MyPID()<<" - SL passive "<<crosslinkermap_->GID(crossLID)<<", element "<<element->Id()<<", Type "<<(*crosslinkertype_)[crossLID]<<std::endl;
              }
            }
            else // standard linker color (red) without any active linkers
              color = 0.5;
          }
          // different color for substrate filaments of motility assay
          if(networktype_==statmech_network_motassay && (*filamentnumber_)[discret_->NodeColMap()->LID(element->Nodes()[0]->Id())]<statmechparams_.get<int>("NUMSUBSTRATEFIL",0))
            color = 0.375;

          // highlight contacting elements
          if(beamcmanager!=Teuchos::null)
          {
            for(int j=0; j<(int)(beamcmanager->Pairs()).size(); j++)
              if(beamcmanager->Pairs()[j]->GetContactFlag() && (element->Id()==(beamcmanager->Pairs())[j]->Element1()->Id() || element->Id()==(beamcmanager->Pairs())[j]->Element2()->Id()))
                color = 1.0; //0.375;

            // loop over BTS pairs as well (if any), here we only need to check Element1 (beam ele)
            for (int j=0;j<(int)(beamcmanager->BTSPHPairs()).size(); j++)
            {
              if(beamcmanager->BTSPHPairs()[j]->GetContactFlag() && (element->Id()==(beamcmanager->BTSPHPairs())[j]->Element1()->Id() ))
                color = 1.0;
            }
          }

          //if no periodic boundary conditions are to be applied, we just plot the current element
          if (periodlength_->at(0) == 0.0)
          {
            // check whether the kinked visualization is to be applied
            bool kinked = CheckForKinkedVisual(element->Id());
            if (eot == DRT::ELEMENTS::Beam3Type::Instance() ||
                eot==DRT::ELEMENTS::Beam3rType::Instance() ||
                eot==DRT::ELEMENTS::BeamCLType::Instance() ||
                eot==DRT::ELEMENTS::Beam3ebType::Instance()||
                eot == DRT::ELEMENTS::Truss3Type::Instance() ||
                eot == DRT::ELEMENTS::Spring3Type::Instance() ||
                eot==DRT::ELEMENTS::Truss3CLType::Instance())
            {
              if (!kinked)
              {
                int numnode = element->NumNode();
                if(eot==DRT::ELEMENTS::BeamCLType::Instance() || eot==DRT::ELEMENTS::Truss3CLType::Instance())
                  numnode = 2;
                if(eot==DRT::ELEMENTS::Beam3ebType::Instance())
                {
                  const int n_axial= statmechparams_.get<int>("GMSHNINTPT",10.0);
                  for (int j=0; j<n_axial-1; j++)
                  {
                    //define output coordinates
                    LINALG::SerialDenseMatrix coordout(3,2);
                    for(int m=0; m<coordout.M(); m++)
                      for(int n=0; n<coordout.N(); n++)
                        coordout(m,n)=coord(m,j+n);


                    GmshWedge(nline,coordout,element,gmshfilecontent,color,IsBarbed);
                  }

                }
                else
                {
                  for (int j=0; j<numnode - 1; j++)
                  {
                    //define output coordinates
                    LINALG::SerialDenseMatrix coordout(3,2);
                    for(int m=0; m<coordout.M(); m++)
                      for(int n=0; n<coordout.N(); n++)
                        coordout(m,n)=coord(m,j+n);

                    GmshWedge(nline,coordout,element,gmshfilecontent,color);
                  }
                }
              }
              else
                GmshKinkedVisual(coord, 0.875, element->Id(), gmshfilecontent);
            }
            else if (eot == DRT::ELEMENTS::Torsion3Type::Instance())
            {
              double beadcolor = 0.75;
              for (int j=0; j<element->NumNode(); j++)
              {
                gmshfilecontent << "SP(" << std::scientific;
                gmshfilecontent << coord(0, j) << "," << coord(1, j) << ","<< coord(2, j);
                gmshfilecontent << ")" << "{" << std::scientific << beadcolor << ","<< beadcolor << "};" << std::endl;
              }
            }
            else
            {
              //nothing!
            }
          }
          //in case of periodic boundary conditions we have to take care to plot correctly an element broken at some boundary plane
          else
          {
            if(eot==DRT::ELEMENTS::Beam3ebType::Instance())
            {
              const int n_axial= statmechparams_.get<int>("GMSHNINTPT",10.0);
              for (int j=0; j<n_axial - 1; j++)
              {
                //define output coordinates
                LINALG::SerialDenseMatrix coordout(3,2);
                for(int m=0; m<coordout.M(); m++)
                  for(int n=0; n<coordout.N(); n++)
                    coordout(m,n)=coord(m,j+n);

                GmshOutputPeriodicBoundary(coordout, color, gmshfilecontent,element->Id(),false, IsBarbed);
              }
            }
              else
                GmshOutputPeriodicBoundary(coord, color, gmshfilecontent,element->Id(),false);
          }

        }//end if(eot != rigid sphere)

        else  // compute and stream gmsh output for rigid sphere
        {
          double color=0.5;

          DRT::ELEMENTS::Rigidsphere* currele = NULL;
          currele = dynamic_cast<DRT::ELEMENTS::Rigidsphere*> (discret_->gElement(discret_->ElementColMap()->GID(i)));

          double eleradius = currele->Radius();

          if(beamcmanager!=Teuchos::null)   // highlighting active contact and potential interaction
          {
            // loop over BTS pairs (if any), here we only need to check Element2 (Rigidsphere ele)
            for (int i=0;i<(int)(beamcmanager->BTSPHPairs()).size();++i)
            {
              // abbreviations
              int id2 =  (beamcmanager->BTSPHPairs())[i]->Element2()->Id();
              bool active = (beamcmanager->BTSPHPairs())[i]->GetContactFlag();

              // if element is member of an active contact pair, choose different color
              if ( currele->Id()==id2 && active) color = 1.0;
            }
          }

          double plotfactorthick = statmechparams_.get<double>("PlotFactorThick", 1.0);

          if (plotfactorthick==0.0)
          {
            // ********************** Visualization as a point ***********************************************

            // syntax for scalar point:  SP( coordinates x,y,z ){value at point (determines the color)}
            gmshfilecontent <<"SP(" << std::scientific << coord(0,0) << "," << coord(1,0) << "," << coord(2,0)
                                        << "){" << std::scientific  << color << "};"<< std::endl;

            // ***********************************************************************************************
          }
          else
          {
            double plotradius = eleradius*plotfactorthick;
            // ********************** Visualization as an icosphere ****************************************

            // for details see http://en.wikipedia.org/wiki/Icosahedron
            // and http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
            //
            // sphere is visualized as icosphere
            // the basic icosahedron consists of 20 equilateral triangles (12 vertices)
            // further refinement by subdividing the triangles

            // list storing the (x,y,z) coordinates of all vertices
            std::vector< std::vector<double> > vertexlist(12, std::vector<double>(3,0));

            // list storing the indices of the three vertices that define a triangular face
            std::vector< std::vector<int> > facelist(20, std::vector<int>(3,0));

            double normfac = sqrt( 1.0 + 0.25* pow(1+sqrt(5),2) );
            double c = 0.5*(1.0+sqrt(5))/normfac*plotradius;
            double d = 1/normfac*plotradius;

            // compute the final coordinates of the initial 12 vertices
            vertexlist[0][0]+=-d; vertexlist[0][1]+=c; vertexlist[0][2]+=0;
            vertexlist[1][0]+=d; vertexlist[1][1]+=c; vertexlist[1][2]+=0;
            vertexlist[2][0]+=-d; vertexlist[2][1]+=-c; vertexlist[2][2]+=0;
            vertexlist[3][0]+=d; vertexlist[3][1]+=-c; vertexlist[3][2]+=0;

            vertexlist[4][0]+=0; vertexlist[4][1]+=-d; vertexlist[4][2]+=c;
            vertexlist[5][0]+=0; vertexlist[5][1]+=d; vertexlist[5][2]+=c;
            vertexlist[6][0]+=0; vertexlist[6][1]+=-d; vertexlist[6][2]+=-c;
            vertexlist[7][0]+=0; vertexlist[7][1]+=d; vertexlist[7][2]+=-c;

            vertexlist[8][0]+=c; vertexlist[8][1]+=0; vertexlist[8][2]+=-d;
            vertexlist[9][0]+=c; vertexlist[9][1]+=0; vertexlist[9][2]+=d;
            vertexlist[10][0]+=-c; vertexlist[10][1]+=0; vertexlist[10][2]+=-d;
            vertexlist[11][0]+=-c; vertexlist[11][1]+=0; vertexlist[11][2]+=d;


            // fill initial facelist
            facelist[0][0]=0;   facelist[0][1]=11;    facelist[0][2]=5;
            facelist[1][0]=0;   facelist[1][1]=5;     facelist[1][2]=1;
            facelist[2][0]=0;   facelist[2][1]=1;     facelist[2][2]=7;
            facelist[3][0]=0;   facelist[3][1]=7;     facelist[3][2]=10;
            facelist[4][0]=0;   facelist[4][1]=10;    facelist[4][2]=11;

            facelist[5][0]=1;   facelist[5][1]=5;     facelist[5][2]=9;
            facelist[6][0]=5;   facelist[6][1]=11;    facelist[6][2]=4;
            facelist[7][0]=11;  facelist[7][1]=10;    facelist[7][2]=2;
            facelist[8][0]=10;  facelist[8][1]=7;     facelist[8][2]=6;
            facelist[9][0]=7;   facelist[9][1]=1;     facelist[9][2]=8;

            facelist[10][0]=3;  facelist[10][1]=9;    facelist[10][2]=4;
            facelist[11][0]=3;  facelist[11][1]=4;    facelist[11][2]=2;
            facelist[12][0]=3;  facelist[12][1]=2;    facelist[12][2]=6;
            facelist[13][0]=3;  facelist[13][1]=6;    facelist[13][2]=8;
            facelist[14][0]=3;  facelist[14][1]=8;    facelist[14][2]=9;

            facelist[15][0]=4;  facelist[15][1]=9;    facelist[15][2]=5;
            facelist[16][0]=2;  facelist[16][1]=4;    facelist[16][2]=11;
            facelist[17][0]=6;  facelist[17][1]=2;    facelist[17][2]=10;
            facelist[18][0]=8;  facelist[18][1]=6;    facelist[18][2]=7;
            facelist[19][0]=9;  facelist[19][1]=8;    facelist[19][2]=1;

            // level of refinement, num_faces = 20 * 4^(ref_level)
            int ref_level= 3;
            // refine the icosphere by calling GmshRefineIcosphere
            for (int p=0; p<ref_level; ++p)
            {
          //    std::cout << "\np= " << p << ", num_vertices= " << vertexlist.size() << ", num_faces= " << facelist.size() << std::endl;
              GmshRefineIcosphere(vertexlist,facelist,plotradius);
            }

            const double centercoord[] = {coord(0,0), coord(1,0), coord(2,0)};
            for (unsigned int i=0; i<facelist.size(); ++i)
              PrintGmshTriangleToStream(gmshfilecontent,vertexlist,facelist[i][0],facelist[i][1],facelist[i][2],color,centercoord);

            // ********************* end: visualization as an icosphere ***************************************
          }
        } // end: compute and stream gmsh output for rigid sphere

      } // discret_->NumMyColElements()
      //write content into file and close it (this way we make sure that the output is written serially)
      fputs(gmshfilecontent.str().c_str(), fp);
      fclose(fp);
    }
    discret_->Comm().Barrier();
  }

  // plot the periodic boundary box
  LINALG::Matrix<3,1> center;
  for(int i=0; i<(int)center.M(); i++)
    center(i) = periodlength_->at(i)/2.0;
  GmshOutputBox(0.0, &center, *periodlength_, &filename);
  // plot crosslink molecule diffusion and (partial) bonding
  double color = 0.125;
  GmshOutputCrosslinkDiffusion(color, &filename, disrow);
  // plot Neumann nodes
  GmshOutputPointNeumann(disrow, &filename, time, color);

  discret_->Comm().Barrier();

  // finish data section of this view by closing curly brackets
  if (discret_->Comm().MyPID() == 0)
  {
    fp = fopen(filename.str().c_str(), "a");
    std::stringstream gmshfileend;

    if(DRT::INPUT::IntegralValue<int>(statmechparams_,"GMSHNETSTRUCT") &&
       DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_structanaly)
    {
      // plot the rotated material triad (with axis length 0.5)
      for(int i=0; i<trafo_->M(); i++)
      {
        LINALG::SerialDenseMatrix coord(3,2);
        for(int j=0;j<coord.M(); j++)
        {
          coord(j,0) = cog_(j);
          coord(j,1) = cog_(j)+0.5*(*trafo_)(i,j);
        }
        GmshWedge(1,coord,discret_->lRowElement(0),gmshfileend,0.0,true,true);
      }
      // plot the cog
      std::vector<double> dimension(3,0.05);
      GmshOutputBox(0.75, &cog_, dimension, &filename, false);
      gmshfileend << "SP(" << std::scientific;
      gmshfileend << cog_(0)<<","<<cog_(1)<<","<<cog_(2)<<"){" << std::scientific << 0.75 << ","<< 0.75 <<"};"<<std::endl;

      // gmsh output of detected network structure volume
      int nline = 16;
      if(step>0)
        GmshNetworkStructVolume(nline, gmshfileend, 0.875);
    }
    // add black dot for correct coloring...
    if(periodlength_->at(0)==0.0)
      gmshfileend << "SP("<<periodlength_->at(0)<<",0,"<<periodlength_->at(2)<<"){0,0};"<<std::endl;
    gmshfileend << "};" << std::endl;
    fputs(gmshfileend.str().c_str(), fp);
    fclose(fp);
  }

  // return simultaneously (not sure if really needed)
  discret_->Comm().Barrier();

  return;
} // STATMECH::StatMechManager::GmshOutput()


/*----------------------------------------------------------------------*
 | gmsh output data in case of periodic boundary conditions             |
 |                                                    public)cyron 02/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GmshOutputPeriodicBoundary(const LINALG::SerialDenseMatrix& coord, const double& color, std::stringstream& gmshfilecontent, int eleid, bool ignoreeleid, bool IsBarbed)
{
  //number of solid elements by which a round line is depicted
  const int nline = 16;

  //number of spatial dimensions
  const int ndim = 3;
  // get Element Type of the first Element to determine the graphics output
  DRT::Element* element = discret_->gElement(eleid);

  bool dotline = false;
  bool kinked = false;

  if (ignoreeleid)
    dotline = true;
  else
  {
    // draw colored lines between two nodes of a beam or a truss element (meant for filaments/crosslinks/springs)
    const DRT::ElementType & eot = element->ElementType();
    if(element->ElementType().Name()=="Beam3rType")
      dotline = eot==DRT::ELEMENTS::Beam3rType::Instance();
    if (element->ElementType().Name() == "Beam3Type")
      dotline = eot == DRT::ELEMENTS::Beam3Type::Instance();
    else if(element->ElementType().Name()=="BeamCLType")
      dotline = eot==DRT::ELEMENTS::BeamCLType::Instance();
    else if(element->ElementType().Name()=="Beam3ebType")
      dotline = eot==DRT::ELEMENTS::Beam3ebType::Instance();
    else if(element->ElementType().Name()=="Truss3CLType")
      dotline = eot==DRT::ELEMENTS::Truss3CLType::Instance();
    else if (element->ElementType().Name() == "Truss3Type")
      dotline = eot == DRT::ELEMENTS::Truss3Type::Instance();
    else if (element->ElementType().Name() == "Spring3Type")
      dotline = eot == DRT::ELEMENTS::Spring3Type::Instance();
    // draw spheres at node positions ("beads" of the bead spring model)
    else if (eot == DRT::ELEMENTS::Torsion3Type::Instance())
    {
      double beadcolor = 0.75;
      for (int i=0; i<element->NumNode(); i++)
      {
        //writing element by nodal coordinates as a sphere
        gmshfilecontent << "SP(" << std::scientific;
        gmshfilecontent << coord(0,i) << "," << coord(1, i) << "," << coord(2,i);
        gmshfilecontent << ")" << "{" << std::scientific << beadcolor << ","<< beadcolor << "};" << std::endl;
      }
    }
    /* determine whether crosslink connects two filaments or occupies two binding spots on the same filament;
     * variable triggers different visualizations.*/
    kinked = CheckForKinkedVisual(element->Id());
  }

  if (dotline)
  {
    /*detect and save in vector "cut", at which boundaries the element is broken due to periodic boundary conditions;
     * the entries of cut have the following meaning: 0: element not broken in respective coordinate direction, 1:
     * element broken in respective coordinate direction (node 0 close to zero boundary and node 1 close to boundary
     * at periodlength_);  2: element broken in respective coordinate direction (node 1 close to zero boundary and node
     * 0 close to boundary at periodlength_);*/
    LINALG::SerialDenseMatrix cut= LINALG::SerialDenseMatrix(3, coord.N()-1, true);

    /* "coord" currently holds the shifted set of coordinates.
     * In order to determine the correct vector "dir" of the visualization at the boundaries,
     * a copy of "coord" with adjustments in the proper places is introduced*/
    for (int i=0; i<cut.N(); i++)
    {
      LINALG::SerialDenseMatrix unshift(3,2);

      int numshifts = 0;
      int shiftdof = -1;
      for (int dof=0; dof<ndim; dof++)
      {
        // initialize unshift with coord values
        unshift(dof,0) = coord(dof,i);
        unshift(dof,1) = coord(dof,i+1);
        if (fabs(coord(dof,i+1)-periodlength_->at(dof)-coord(dof,i)) < fabs(coord(dof,i+1) - coord(dof,i)))
        {
          cut(dof, i) = 1.0;
          shiftdof = dof;
          unshift(dof,1) -= periodlength_->at(dof);
          numshifts++;
        }
        if (fabs(coord(dof,1)+periodlength_->at(dof) - coord(dof,i)) < fabs(coord(dof,i+1)-coord(dof,i)))
        {
          cut(dof,i) = 2.0;
          shiftdof = dof;
          unshift(dof,1) += periodlength_->at(dof);
          numshifts++;
        }
      }

      // write special output for broken elements
      if (numshifts > 0)
      {
//        // Uncomment if check of oscillating node positions is needed
//        if(cut(2,0)==1.0)
//        {
//          gmshfilecontent << "SP(" << std::scientific;
//          gmshfilecontent << coord(0,i+1) << "," << coord(1, i+1) << "," << coord(2,i+1);
//          gmshfilecontent << ")" << "{" << std::scientific << 0.75 << ","<< 0.75 << "};" << std::endl;
//        }
//        else if(cut(2,0)==2.0)
//        {
//          gmshfilecontent << "SP(" << std::scientific;
//          gmshfilecontent << coord(0,i) << "," << coord(1, i) << "," << coord(2,i);
//          gmshfilecontent << ")" << "{" << std::scientific << 0.75 << ","<< 0.75 << "};" << std::endl;
//        }
        // directional vector
        LINALG::Matrix<3, 1> dir;
        for (int dof = 0; dof < ndim; dof++)
          dir(dof) = unshift(dof,1) - unshift(dof,0);
        dir.Scale(1.0/dir.Norm2());

        /* determine the intersection points of the line through unshift(:,0) and direction dir with the faces of the boundary cube
         * and sort them by distance. Thus, we obtain an order by which we have to shift the element back into the cube so that all
         * segments that arise by multiple shifts remain within the volume (see my notes on 6/12/2011).*/
        LINALG::Matrix<3,2> lambdaorder;
        lambdaorder.PutScalar(1e6);
        // collect lambdas
        for(int dof=0; dof<(int)lambdaorder.M(); dof++)
        {
          switch((int)cut(dof,i))
          {
            case 1:
            {
              lambdaorder(dof,0) = -coord(dof, i) / dir(dof);
              lambdaorder(dof,1) = dof;
            }
            break;
            case 2:
            {
              lambdaorder(dof,0) = (periodlength_->at(dof) - coord(dof,i)) / dir(dof);
              lambdaorder(dof,1) = dof;
            }
            break;
            default:
            {
              lambdaorder(dof,1) = dof;
            }
            break;
          }
        }
        // sort the lambdas (ascending values) and indices accordingly
        // in case of multiple shifts
        if(numshifts>1)
        {
          for(int j=0; j<(int)lambdaorder.M()-1; j++)
            for(int k=j+1; k<(int)lambdaorder.M(); k++)
              if(lambdaorder(k,0)<lambdaorder(j,0))
              {
                double temp = lambdaorder(j,0);
                int tempindex = (int)lambdaorder(j,1);
                lambdaorder(j,0) = lambdaorder(k,0);
                lambdaorder(j,1) = lambdaorder(k,1);
                lambdaorder(k,0) = temp;
                lambdaorder(k,1) = tempindex;
              }
        }
        else  // for a single shift (the majority of broken elements), just put the index and the lambda of the broken dof in front
          for(int n=0; n<(int)lambdaorder.N(); n++)
          {
            double tmp = lambdaorder(shiftdof,n);
            lambdaorder(0,n) = tmp;
          }

        // calculate segment lambdas
        for(int dof=numshifts-1; dof>0; dof--)
          lambdaorder(dof,0) -= lambdaorder(dof-1,0);

        // the idea is to gradually shift the matrix "unshift" back into the volume and, while doing so,
        // calculate the segments except for the last one
        // determine closest boundary component-wise
        for(int shift=0; shift<numshifts; shift++)
        {
          //second point
          for(int j=0 ;j<unshift.M(); j++)
            unshift(j,1) = unshift(j,0) + lambdaorder(shift,0)*dir(j);
          //GmshOutput
          if(!ignoreeleid)
            GmshWedge(nline,unshift,element,gmshfilecontent,color);
          else
            GmshWedge(nline,unshift,element,gmshfilecontent,color,true);

          int currshift = (int)lambdaorder(shift,1);
          // shift the coordinates of the second point
          if(cut(currshift,i)==1.0)
            unshift(currshift,1) += periodlength_->at(currshift);
          else if(cut(currshift,i)==2.0)
            unshift(currshift,1) -= periodlength_->at(currshift);
          // make second point the first and calculate new second point in the next iteration!
          for(int j=0; j<unshift.M(); j++)
            unshift(j,0) = unshift(j,1);
        }
        //last segment
        // the last segment
        for(int dof=0; dof<unshift.M(); dof++)
          unshift(dof,1) = coord(dof,i+1);
        if(!ignoreeleid)
          GmshWedge(nline,unshift,element,gmshfilecontent,color);
        else
          GmshWedge(nline,unshift,element,gmshfilecontent,color,true);
      }
      else // output for continuous elements
      {
        if (!kinked)
        {
          if(!ignoreeleid)
            GmshWedge(nline,coord,element,gmshfilecontent,color);
          else
            GmshWedge(nline,coord,element,gmshfilecontent,color,true);
        }
        else
          GmshKinkedVisual(coord, 0.875, element->Id(), gmshfilecontent);
      }
    }
  }
  return;
} // STATMECH::StatMechManager::GmshOutputPeriodicBoundary()

/*----------------------------------------------------------------------*
 | plot the periodic boundary box                  (public) mueller 7/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GmshOutputBox(double boundarycolor, LINALG::Matrix<3,1>* boxcenter,
                                              std::vector<double>& dimension,
                                              const std::ostringstream *filename,
                                              bool barrier,
                                              int pid)
{
  // plot the periodic box in case of periodic boundary conditions (first processor)
  if (periodlength_->at(0) > 0.0 && discret_->Comm().MyPID() == pid)
  {
    FILE *fp = fopen(filename->str().c_str(), "a");
    std::stringstream gmshfilefooter;
    // get current period length

    double xmin = (*boxcenter)(0)-dimension[0]/2.0;
    double xmax = (*boxcenter)(0)+dimension[0]/2.0;
    double ymin = (*boxcenter)(1)-dimension[1]/2.0;
    double ymax = (*boxcenter)(1)+dimension[1]/2.0;
    double zmin = (*boxcenter)(2)-dimension[2]/2.0;
    double zmax = (*boxcenter)(2)+dimension[2]/2.0;

    // define boundary lines
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmin << "," << ymin << "," << zmin << "," << xmax << ","<< ymin << "," << zmin;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << std::endl;
    // line 2
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmax << "," << ymin << "," << zmin << "," << xmax << "," << ymax << "," << zmin;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << std::endl;
    // line 3
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmax << "," << ymax << "," << zmin << "," << xmax << "," << ymax << "," << zmax;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << std::endl;
    // line 4
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmax << "," << ymax << "," << zmax << "," << xmin << "," << ymax << "," << zmax;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << std::endl;
    // line 5
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmin << "," << ymax << "," << zmax << "," << xmin << "," << ymin << "," << zmax;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << std::endl;
    // line 6
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmin << "," << ymin << "," << zmax << "," << xmin << ","<< ymin << "," << zmin;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << std::endl;
    // line 7
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmin << "," << ymin << "," << zmin << "," << xmin << ","<< ymax << "," << zmin;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << std::endl;
    // line 8
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmin << "," << ymax << "," << zmin << "," << xmax << "," << ymax << "," << zmin;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << std::endl;
    // line 9
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmin << "," << ymax << "," << zmin << "," << xmin << "," << ymax << "," << zmax;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << std::endl;
    // line 10
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmax << "," << ymin << "," << zmin << "," << xmax << "," << ymin << "," << zmax;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << std::endl;
    // line 11
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmax << "," << ymin << "," << zmax << "," << xmax << "," << ymax<< "," << zmax;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << std::endl;
    // line 12
    gmshfilefooter << "SL(" << std::scientific;
    gmshfilefooter << xmax << "," << ymin << "," << zmax << "," << xmin << "," << ymin<< "," << zmax;
    gmshfilefooter << ")" << "{" << std::scientific << boundarycolor << ","<< boundarycolor << "};" << std::endl;

    fputs(gmshfilefooter.str().c_str(), fp);
    fclose(fp);
  }
  // wait for Proc 0 to catch up to the others
  if(barrier)
    discret_->Comm().Barrier();
}// STATMECH::StatMechManager::GmshOutputBoundaryBox

/*----------------------------------------------------------------------*
 | Gmsh output for crosslink molecule diffusion    (public) mueller 7/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GmshOutputCrosslinkDiffusion(double color, const std::ostringstream *filename, const Epetra_Vector& disrow)
{
  // export row displacement to column map format
  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(disrow, discol);

  //In case of the 4-noded Crosslinker Element we need the following vectors
  Teuchos::RCP<Epetra_MultiVector> bspotpositions = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));
  Teuchos::RCP<Epetra_MultiVector> bspotrotations = Teuchos::null;
  if(statmechparams_.get<double>("ILINK",0.0)>0.0)
    bspotrotations = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));

  GetBindingSpotPositions(discol, bspotpositions, bspotrotations);

  if (!discret_->Comm().MyPID())
  {
    FILE *fp = fopen(filename->str().c_str(), "a");
    //special visualization for crosslink molecules with one/two bond(s); going through the Procs

    //fp = fopen(filename->str().c_str(), "a");
    std::stringstream gmshfilebonds;

    // Print Free monomers
    /*
    for (int i=0; i<MonomerBindingStatus_->MyLength(); i++)
    {
      double beadcolor2 = 3*color;
      //writing element by nodal coordinates as a sphere
      gmshfilebonds << "SP(" << std::scientific;
      gmshfilebonds<< (*VisualizeMonomerPositions_)[0][i]<< "," << (*VisualizeMonomerPositions_)[1][i] << "," << (*VisualizeMonomerPositions_)[2][i];
      gmshfilebonds << ")" << "{" << std::scientific << beadcolor2 << "," << beadcolor2 << "};" << std::endl;
    }*/

    // first, just update positions: redundant information on all procs
    for (int i=0; i<numbond_->MyLength(); i++)
    {
      switch ((int) (*numbond_)[i])
      {
        // free linkers
        case 0:
        {
          /*
          double beadcolor = 5*color;
          //writing element by nodal coordinates as a sphere
          gmshfilebonds << "SP(" << std::scientific;
          gmshfilebonds<< (*visualizepositions_)[0][i]<< "," << (*visualizepositions_)[1][i] << "," << (*visualizepositions_)[2][i];
          gmshfilebonds << ")" << "{" << std::scientific << beadcolor << "," << beadcolor << "};" << std::endl;
          */

        }
        break;
        // crosslink molecule with one bond
        case 1:
        {
          // determine position of bpspotLID entry
          int bspotLID = bspotcolmap_->LID((int)(*crosslinkerbond_)[0][i]);
          if((*crosslinkerbond_)[0][i]<-0.9)
            bspotLID = bspotcolmap_->LID((int)(*crosslinkerbond_)[1][i]);

          LINALG::SerialDenseMatrix coord(3, 2, true);
          double length = 0.0;
          for (int j=0; j<coord.M(); j++)
          {
            coord(j, 0) = (*bspotpositions)[j][bspotLID];
            coord(j, 1) = (*visualizepositions_)[j][i];
            length += (coord(j,1)-coord(j,0))*(coord(j,1)-coord(j,0));
          }
          double beadcolor = 2*color; //blue
          // in case of periodic boundary conditions
          if (periodlength_->at(0) > 0.0)
          {
            // get arbitrary element (we just need it to properly visualize)
            DRT::Element* tmpelement=discret_->lRowElement(0);
            GmshOutputPeriodicBoundary(coord, 2*color, gmshfilebonds, tmpelement->Id(), true);
            // visualization of "real" crosslink molecule positions (attention: shifted by one step since StatMechUpdate is called before the Newton scheme)
            //beadcolor = 0.0; //black
            //gmshfilebonds << "SP(" << scientific;
            //gmshfilebonds << (*crosslinkerpositions_)[0][i] << ","<< (*crosslinkerpositions_)[1][i] << ","<< (*crosslinkerpositions_)[2][i];
            //gmshfilebonds << ")" << "{" << scientific << beadcolor << ","<< beadcolor << "};" << std::endl;
          }
          else
          {
            gmshfilebonds << "SL(" << std::scientific;
            gmshfilebonds << coord(0, 0) << "," << coord(1, 0) << ","<< coord(2, 0) << "," << coord(0, 1) << "," << coord(1, 1)<< "," << coord(2, 1);
            gmshfilebonds << ")" << "{" << std::scientific << 2*color << ","<< 2*color << "};" << std::endl;
            gmshfilebonds << "SP(" << std::scientific;
            gmshfilebonds << coord(0, 1) << "," << coord(1, 1) << ","<< coord(2, 1);
            gmshfilebonds << ")" << "{" << std::scientific << beadcolor << ","<< beadcolor << "};" << std::endl;
          }
        }
        break;
        // crosslinker element: crosslink molecule (representation) position (Proc 0 only)
        case 2:
        {
          // actual crosslinker element connecting two filaments (self-binding kinked crosslinkers are visualized in GmshKinkedVisual())
          //TODO reactivate if passive linkers are needed
//          if((*searchforneighbours_)[i] > 0.9)
//          {

          bool crosslinked = true;
          // linker attached to binding sites on ONE filament -> different visualization, hence crosslinked = false!
          if(linkermodel_ == statmech_linker_stdintpol || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)
          {
            if((*filamentnumber_)[discret_->NodeColMap()->LID((int)(*bspot2nodes_)[0][bspotcolmap_->GID((*crosslinkerbond_)[0][i])])] == (*filamentnumber_)[discret_->NodeColMap()->LID((int)(*bspot2nodes_)[0][bspotcolmap_->GID((*crosslinkerbond_)[1][i])])] && statmechparams_.get<double>("K_ON_SELF",0.0)>0.0)
              crosslinked = false;
          }
          else
          {
            if((*filamentnumber_)[bspotcolmap_->LID((int)(*crosslinkerbond_)[0][i])]==(*filamentnumber_)[bspotcolmap_->LID((int)(*crosslinkerbond_)[1][i])] && statmechparams_.get<double>("K_ON_SELF",0.0)>0.0)
              crosslinked = false;
          }

          if(crosslinked)
          {
            double beadcolor = 4* color ; //red
            if(crosslinkertype_!=Teuchos::null)
            {
              if((*crosslinkertype_)[i]==1.0)
              {
                beadcolor = 6*color;
                //std::cout<<"Proc "<<discret_->Comm().MyPID()<<" - SP active "<<i<<", element "<<(*crosslink2element_)[i]<<", Type "<<(*crosslinkertype_)[i]<<std::endl;
              }
              //else
              //  std::cout<<"  Proc "<<discret_->Comm().MyPID()<<" - SP passive "<<i<<", element "<<(*crosslink2element_)[i]<<", Type "<<(*crosslinkertype_)[i]<<std::endl;
            }
            gmshfilebonds << "SP(" << std::scientific;
            gmshfilebonds << (*visualizepositions_)[0][i] << ","<< (*visualizepositions_)[1][i] << ","<< (*visualizepositions_)[2][i];
            gmshfilebonds << ")" << "{" << std::scientific << beadcolor << ","<< beadcolor << "};" << std::endl;
          }
//          }
          //TODO reactivate if passive linkers are needed
//          else  // passive crosslink molecule
//          {
//            if((linkermodel_ == statmech_linker_stdintpol || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol) &&
//                statmechparams_.get<double>("K_ON_SELF",0.0)>0.0)
//              dserror("Passified links not implemented for BeamCL elements!");
//            // determine position of nodeGID entry
//            int bspotLID = bspotcolmap_->LID((int)(*crosslinkerbond_)[0][i]);
//            if((*crosslinkerbond_)[0][i]<-0.9)
//              bspotLID = bspotcolmap_->LID((int)(*crosslinkerbond_)[1][i]);
//
//            DRT::Node *node = discret_->lColNode(bspotLID);
//            LINALG::SerialDenseMatrix coord(3, 2, true);
//            for (int j=0; j<coord.M(); j++)
//            {
//              int dofgid = discret_->Dof(node)[j];
//              coord(j, 0) = node->X()[j] + discol[dofgid];
//              coord(j, 1) = (*visualizepositions_)[j][i];
//            }
//
//            double beadcolor = 3*color;
//            // in case of periodic boundary conditions
//            if (periodlength_->at(0) > 0.0)
//            {
//              // get arbitrary element (we just need it to properly visualize)
//              DRT::Element* tmpelement=discret_->lRowElement(0);
//              GmshOutputPeriodicBoundary(coord, 3*color, gmshfilebonds, tmpelement->Id(), true);
//              // visualization of "real" crosslink molecule positions (attention: shifted by one step since StatMechUpdate is called before the Newton scheme)
//              //beadcolor = 0.0; //black
//              //gmshfilebonds << "SP(" << scientific;
//              //gmshfilebonds << (*crosslinkerpositions_)[0][i] << ","<< (*crosslinkerpositions_)[1][i] << ","<< (*crosslinkerpositions_)[2][i];
//              //gmshfilebonds << ")" << "{" << scientific << beadcolor << ","<< beadcolor << "};" << std::endl;
//            }
//            else
//            {
//              gmshfilebonds << "SL(" << std::scientific;
//              gmshfilebonds << coord(0, 0) << "," << coord(1, 0) << ","<< coord(2, 0) << "," << coord(0, 1) << "," << coord(1, 1)<< "," << coord(2, 1);
//              gmshfilebonds << ")" << "{" << std::scientific << 3*color << ","<< 3*color << "};" << std::endl;
//              gmshfilebonds << "SP(" << std::scientific;
//              gmshfilebonds << coord(0, 1) << "," << coord(1, 1) << ","<< coord(2, 1);
//              gmshfilebonds << ")" << "{" << std::scientific << beadcolor << ","<< beadcolor << "};" << std::endl;
//            }
//          }
        }
          break;
        default:
          continue;
      }
    }
    fputs(gmshfilebonds.str().c_str(), fp);
    fclose(fp);
  }

  discret_->Comm().Barrier();
}// GmshOutputCrosslinkDiffusion

/*----------------------------------------------------------------------*
 | Special Gmsh output for crosslinkers occupying two binding spots on  |
 | the same filament                              (public) mueller 07/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GmshKinkedVisual(const LINALG::SerialDenseMatrix& coord, const double& color, int eleid, std::stringstream& gmshfilecontent)
{
  /* We need a third point in order to visualize the crosslinker.
   * It marks the location of the kink.
   */
  std::vector<double> thirdpoint(3,0.0);
  // get the element
  DRT::Element* element = discret_->gElement(eleid);

  // calculate unit tangent
  LINALG::Matrix<3,1> t(true);
  for (int j=0; j<coord.M(); j++)
    t(j) = coord(j, coord.N() - 1) - coord(j, 0);
  t.Scale(1.0/t.Norm2());

  // calculate normal via cross product: [0 0 1]x[tx ty tz]
  LINALG::Matrix<3,1> n(true);
  n(0) = -t(1);
  n(1) = t(0);
  // norm it since the cross product does not keep the length
  n.Scale(1.0/n.Norm2());

  // by modulo operation involving the node IDs
  double alpha = fmod((double) (element->Nodes()[element->NumNode() - 1]->Id() + element->Nodes()[0]->Id()), 2*M_PI);

  // rotate the normal by alpha
  RotationAroundFixedAxis(t,n,alpha);

  // calculation of the third point lying in the direction of the rotated normal
  // height of the third point above the filament
  double h = statmechparams_.get<double> ("R_LINK", 0.0)/3.0;
  for (int j=0; j<3; j++)
    thirdpoint.at(j) = (coord(j, 0) + coord(j, 1)) / 2.0 + h * n(j);

  gmshfilecontent << "SL(" << std::scientific << coord(0,0) << ","<< coord(1,0) << "," << coord(2,0) << ","
                  << thirdpoint.at(0) << ","<< thirdpoint.at(1) << "," << thirdpoint.at(2) << ")"
                  << "{" << std::scientific<< color << "," << color << "};" << std::endl;
  gmshfilecontent << "SL(" << std::scientific << thirdpoint.at(0) << "," << thirdpoint.at(1)<< "," << thirdpoint.at(2) << ","
                  << coord(0, 1) << "," << coord(1,1) << "," << coord(2,1) << ")"
                  << "{" << std::scientific<< color << "," << color << "};" << std::endl;
  gmshfilecontent << "SP(" << std::scientific
                  << thirdpoint.at(0) << "," << thirdpoint.at(1) << ","<< thirdpoint.at(2)
                  << ")" << "{" << std::scientific << color << ","<< color << "};" << std::endl;
  return;
}// STATMECH::StatMechManager::GmshKinkedVisual

/*----------------------------------------------------------------------*
 | prepare visualization vector for Gmsh Output  (private) mueller 08/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GmshPrepareVisualization(const Epetra_Vector& dis)
{
  double ronebond = statmechparams_.get<double> ("R_LINK", 0.0) / 2.0;

  // column map displacement and displacement increment

  Epetra_Vector discol(*(discret_->DofColMap()), true);
  LINALG::Export(dis, discol);

  //In case of the 4-noded Crosslinker Element we need the following vectors
  Teuchos::RCP<Epetra_MultiVector> bspotpositions = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));
//  Teuchos::RCP<Epetra_MultiVector> bspotrotations = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));
  Teuchos::RCP<Epetra_MultiVector> bspotrotations = Teuchos::null;
  if(statmechparams_.get<double>("ILINK",0.0)>0.0)
    bspotrotations = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,3,true));
  GetBindingSpotPositions(discol, bspotpositions, bspotrotations);

//  // get binding spot triads
  Teuchos::RCP<Epetra_MultiVector> bspottriadscol = Teuchos::rcp(new Epetra_MultiVector(*bspotcolmap_,4,true));
  GetBindingSpotTriads(bspotrotations, bspottriadscol);

  if (discret_->Comm().MyPID() == 0)
  {
//    // Capture motion of Free monomers at different time steps
//    for (int i=0; i<MonomerBindingStatus_->MyLength(); i++)
//    {
//      for (int j=0; j<VisualizeMonomerPositions_->NumVectors(); j++)
//        (*VisualizeMonomerPositions_)[j][i] = (*MonomerPositions_)[j][i];
//    }
    for (int i=0; i<numbond_->MyLength(); i++)
    {
      switch((int)(*numbond_)[i])
      {
        // diffusion
        case 0:
        {
          for (int j=0; j<visualizepositions_->NumVectors(); j++)
            (*visualizepositions_)[j][i] = (*crosslinkerpositions_)[j][i];
        }
        break;
        // one bonded crosslink molecule
        case 1:
        {
          // determine position of bspotLID entry
          int bspotLID = bspotcolmap_->LID((int)(*crosslinkerbond_)[0][i]);
          if((*crosslinkerbond_)[0][i]<-0.9)
            bspotLID = bspotcolmap_->LID((int)(*crosslinkerbond_)[1][i]);

          //calculate unit tangent
          LINALG::Matrix<3,1> bspotpos0;
          for (int j=0; j<3; j++)
            bspotpos0(j,0) = (*bspotpositions)[j][bspotLID];

          // first and second vector of the nodal triad
          LINALG::Matrix<3, 1> tangent;
          LINALG::Matrix<3, 1> normal;
          // rotation angle
          double alpha = 0.0;

          LINALG::Matrix<3,3> bspottriad;
          // auxiliary variable for storing a triad in quaternion form
          LINALG::Matrix<4, 1> qnode;
          // triad of node on first filament which is affected by the new crosslinker
          for (int j=0; j<4; j++)
            qnode(j) = (*bspottriadscol)[j][bspotLID];
          LARGEROTATIONS::quaterniontotriad(qnode, bspottriad);

          for(int j=0; j<(int)tangent.M(); j++)
          {
            tangent(j) = bspottriad(j,0);
            normal(j) = bspottriad(j,1);
          }
          // rotation angle
          if(filamentmodel_ == statmech_filament_helical)
            alpha = (*bspotorientations_)[bspotLID];
          else if(networktype_ == statmech_network_motassay)
            alpha = M_PI/2.0;
          else
            alpha = fmod((double) i, 2.0*M_PI);

          // rotate the normal by alpha and store the new direction into the same Matrix ("normal")
          if(alpha!=0.0)
            RotationAroundFixedAxis(tangent,normal,alpha);

          // calculation of the visualized point lying in the direction of the rotated normal
          // For Beam3eb elements this value is not calculated properly.
          // But it's not important since it's only for visualization purposes.
          // Besides singly bound linkers do not have mechanical contribution.
          for (int j=0; j<visualizepositions_->NumVectors(); j++)
            (*visualizepositions_)[j][i] = bspotpos0(j,0) + ronebond*normal(j,0);

        }

        break;
        case 2:
        {
          // actual crosslinker element (not kinked)
          //TODO reactivate if passive linkers are needed
//          if((*searchforneighbours_)[i]>0.9)
//          {
            // loop over filament node components
            for (int j=0; j<visualizepositions_->NumVectors(); j++)
            {
              std::vector<double> dofbspotpositions(crosslinkerbond_->NumVectors(), 0.0);
              // loop over filament node GIDs
              for (int k=0; k<crosslinkerbond_->NumVectors(); k++)
              {
                int bspotGID = (int) (*crosslinkerbond_)[k][i];
                dofbspotpositions.at(k) = (*bspotpositions)[j][bspotGID];
              }
              /* Check if the crosslinker element is broken/discontinuous; if so, reposition the second nodal value.
               * It does not matter which value is shifted as long the shift occurs in a consistent way.*/
              if (periodlength_->at(j) > 0.0)
              {
                (*visualizepositions_)[j][i] = dofbspotpositions.at(0);
                for (int k=0; k<1; k++)
                {
                  // shift position if it is found to be outside the boundary box
                  if (fabs(dofbspotpositions.at(k+1) - periodlength_->at(j) - dofbspotpositions.at(k))< fabs(dofbspotpositions.at(k+1) - dofbspotpositions.at(k)))
                    dofbspotpositions.at(k+1) -= periodlength_->at(j);
                  if (fabs(dofbspotpositions.at(k+1) + periodlength_->at(j) - dofbspotpositions.at(k))< fabs(dofbspotpositions.at(k+1) - dofbspotpositions.at(k)))
                    dofbspotpositions.at(k+1) += periodlength_->at(j);

                  (*visualizepositions_)[j][i] += dofbspotpositions.at(k+1);
                }
                // new crosslink molecule position
                (*visualizepositions_)[j][i] /= 2.0;
              }
              else
                (*visualizepositions_)[j][i] /= 2.0;
            }
//          }
        //TODO reactivate if passive linkers are needed
//          else  // passive crosslink molecule
//          {
//            if((linkermodel_ == statmech_linker_stdintpol  || linkermodel_ == statmech_linker_activeintpol || linkermodel_ == statmech_linker_bellseqintpol)&& statmechparams_.get<double>("K_ON_SELF",0.0)>0.0)
//              dserror("Passified links not implemented for BeamCL elements!");
//            // determine position of bspotLID entry
//            int bspotLID = bspotcolmap_->LID((int)(*crosslinkerbond_)[0][i]);
//            if ((*crosslinkerbond_)[0][i] < -0.9)
//              bspotLID = bspotcolmap_->LID((int)(*crosslinkerbond_)[1][i]);
//
//            const DRT::Node *node0 = discret_->lColNode(bspotLID);
//            //calculate unit tangent
//            LINALG::Matrix<3,1> nodepos0;
//
//            for (int j=0; j<3; j++)
//            {
//              int dofgid0 = discret_->Dof(node0)[j];
//              nodepos0(j,0) = node0->X()[j] + discol[discret_->DofColMap()->LID(dofgid0)];
//            }
//
//            // first and second vector of the nodal triad
//            LINALG::Matrix<3, 1> tangent;
//            LINALG::Matrix<3, 1> normal;
//            // rotation angle
//            double alpha = 0.0;
//
//            // determine tangent and normal direction when helical binding spot geometry is applied
//            if(filamentmodel_ == statmech_filament_helical)
//            {
//              LINALG::Matrix<3,3> bspottriad;
//              // auxiliary variable for storing a triad in quaternion form
//              LINALG::Matrix<4, 1> qnode;
//              // triad of node on first filament which is affected by the new crosslinker
//              for (int j=0; j<4; j++)
//                qnode(j) = (*bspottriadscol)[j][bspotLID];
//              LARGEROTATIONS::quaterniontotriad(qnode, bspottriad);
//
//              for(int j=0; j<(int)tangent.M(); j++)
//              {
//                tangent(j) = bspottriad(j,0);
//                normal(j) = bspottriad(j,1);
//              }
//              // rotation angle
//              alpha = (*bspotorientations_)[bspotLID];
//            }
//            else  // conventional case
//            {
//              // choose a second (neighbour) node
//              const DRT::Node *node1 = NULL;
//              int currfilament = (int)(*filamentnumber_)[bspotLID];
//              if(bspotLID < basisnodes_-1)
//              {
//                if((*filamentnumber_)[bspotLID+1]==currfilament)
//                  node1 = discret_->lColNode(bspotLID+1);
//                else
//                  node1 = discret_->lColNode(bspotLID-1);
//              }
//              if(bspotLID == basisnodes_-1)
//                if((*filamentnumber_)[bspotLID-1]==currfilament)
//                  node1 = discret_->lColNode(bspotLID-1);
//
//              //calculate unit tangent
//              for (int j=0; j<3; j++)
//              {
//                int dofgid1 = discret_->Dof(node1)[j];
//                double nodeposj1 = node1->X()[j] + discol[discret_->DofColMap()->LID(dofgid1)];
//                tangent(j) = nodeposj1 - nodepos0(j,0);
//              }
//              tangent.Scale(1 / tangent.Norm2());
//              // calculate normal via cross product: [0 0 1]x[tx ty tz]
//              normal.Clear();
//              normal(0) = -tangent(1);
//              normal(1) = tangent(0);
//              // norm it since the cross product does not keep the length
//              normal.Scale(1 / normal.Norm2());
//              // random angle
//              // by modulo operation involving the crosslink molecule number
//              alpha = fmod((double) i, 2*M_PI);
//            }
//
//            // rotate the normal by alpha and store the new direction into the same Matrix ("normal")
//            RotationAroundFixedAxis(tangent,normal,alpha);
//
//            // calculation of the visualized point lying in the direction of the rotated normal
//            for (int j=0; j<visualizepositions_->NumVectors(); j++)
//              (*visualizepositions_)[j][i] = nodepos0(j,0) + ronebond*normal(j,0);
//          }
        }
        break;
      }
    }
    if (periodlength_->at(0) > 0.0)
      CrosslinkerPeriodicBoundaryShift(visualizepositions_);
  }
  // synchronize results
  Teuchos::RCP<Epetra_MultiVector> visualizepositionstrans = Teuchos::rcp(new Epetra_MultiVector(*transfermap_, 3, true));
  CommunicateMultiVector(visualizepositionstrans, visualizepositions_);
  // debug std::couts
  //std::cout<<*visualizepositions_<<std::endl;
}//GmshPrepareVisualization

/*----------------------------------------------------------------------*
 | wedge output for two-noded beams                        cyron   11/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GmshWedge(const int& n,
                                          const Epetra_SerialDenseMatrix& coord,
                                          DRT::Element* thisele,
                                          std::stringstream& gmshfilecontent,
                                          const double color,
                                          bool ignoreeleid,
                                          bool drawsphere,
                                          bool IsBarbed)
{
  //if this element is a line element capable of providing its radius get that radius
  double radius = 0.0;
  if(!ignoreeleid)
  {
    const DRT::ElementType & eot = thisele->ElementType();
    if(eot == DRT::ELEMENTS::Beam3Type::Instance())
      radius = sqrt(sqrt(4 * ((dynamic_cast<DRT::ELEMENTS::Beam3*>(thisele))->Izz()) / M_PI));
    else if(eot == DRT::ELEMENTS::Beam3rType::Instance())
      radius = sqrt(sqrt(4 * ((dynamic_cast<DRT::ELEMENTS::Beam3r*>(thisele))->Izz()) / M_PI));
    else if(eot == DRT::ELEMENTS::BeamCLType::Instance())
      radius = sqrt(sqrt(4 * ((dynamic_cast<DRT::ELEMENTS::BeamCL*>(thisele))->Izz()) / M_PI));
    else if(eot == DRT::ELEMENTS::Beam3ebType::Instance())
      radius = sqrt(sqrt(4 * ((dynamic_cast<DRT::ELEMENTS::Beam3eb*>(thisele))->Izz()) / M_PI));
    else if(eot == DRT::ELEMENTS::Truss3CLType::Instance())
      radius = sqrt((dynamic_cast<DRT::ELEMENTS::Truss3CL*>(thisele))->CSec() / M_PI);
    else if(eot == DRT::ELEMENTS::Truss3Type::Instance())
      radius = sqrt((dynamic_cast<DRT::ELEMENTS::Truss3*>(thisele))->CSec() / M_PI);
    else if(eot == DRT::ELEMENTS::Spring3Type::Instance())
      radius = sqrt(4.75166e-06 / M_PI); // Hardcode radius if for Spring3 Element
    else
      dserror("thisele is not a line element providing its radius. Check your input file and your defines flags!");
    // case: crosslinker
    if(thisele->Id()>basisnodes_)
      radius = sqrt(statmechparams_.get<double>("ALINK",4.75166e-06) / M_PI); //default value according to diss. Tharmann
  }
  //line elements are plotted by a factor PlotFactorThick thicker than they are actually to allow for better visibility in gmsh pictures
  radius *= statmechparams_.get<double>("PlotFactorThick", 1.0);

  //solid lines plotted if PlotFactorThick != 0; otherwise output utilizes gmsh line visualization (lines width with fixed number of pixels, not a physical diameter)
  if(radius > 0.0)
  {

    // some local variables
    LINALG::Matrix<3,6> prism;
    LINALG::Matrix<3,1> axis;
    LINALG::Matrix<3,1> radiusvec1;
    LINALG::Matrix<3,1> radiusvec2;
    LINALG::Matrix<3,1> auxvec;
    LINALG::Matrix<3,1> theta;
    LINALG::Matrix<3,3> R;


    // compute three dimensional angle theta
    for (int j=0;j<3;++j)
      axis(j) = coord(j,1) - coord(j,0);
    double norm_axis = axis.Norm2();
    for (int j=0;j<3;++j)
      theta(j) = axis(j) / norm_axis * 2 * M_PI / n;

    // Compute rotation matirx R from rotation angle theta
    LARGEROTATIONS::angletotriad(theta,R);

    // Now the first prism will be computed via two radiusvectors, that point from each of
    // the nodes to two points on the beam surface. Further prisms will be computed via a
    // for-loop, where the second node of the previous prism is used as the first node of the
    // next prism, whereas the central points (=nodes) stay  identic for each prism. The
    // second node will be computed by a rotation matrix and a radiusvector.

    // compute radius vector for first surface node of first prims
    for (int j=0;j<3;++j) auxvec(j) = coord(j,0) + norm_axis;

    // radiusvector for point on surface
    radiusvec1(0) = auxvec(1)*axis(2) - auxvec(2)*axis(1);
    radiusvec1(1) = auxvec(2)*axis(0) - auxvec(0)*axis(2);
    radiusvec1(2) = auxvec(0)*axis(1) - auxvec(1)*axis(0);

    // initialize all prism points to nodes
    for (int j=0;j<3;++j)
    {
      prism(j,0) = coord(j,0);
      prism(j,1) = coord(j,0);
      prism(j,2) = coord(j,0);
      prism(j,3) = coord(j,1);
      prism(j,4) = coord(j,1);
      prism(j,5) = coord(j,1);
    }

    // get first point on surface for node1 and node2
    for (int j=0;j<3;++j)
    {
      prism(j,1) += radiusvec1(j) / radiusvec1.Norm2() * radius;
      prism(j,4) += radiusvec1(j) / radiusvec1.Norm2() * radius;
    }

    // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
    radiusvec2.Multiply(R,radiusvec1);

    // get second point on surface for node1 and node2
    for(int j=0;j<3;j++)
    {
      prism(j,2) += radiusvec2(j) / radiusvec2.Norm2() * radius;
      prism(j,5) += radiusvec2(j) / radiusvec2.Norm2() * radius;
    }

    // now first prism is built -> put coordinates into filecontent-stream
    // Syntax for gmsh input file
    // SI(x,y--,z,  x+.5,y,z,    x,y+.5,z,   x,y,z+.5, x+.5,y,z+.5, x,y+.5,z+.5){1,2,3,3,2,1};
    // SI( coordinates of the six corners ){colors}
    gmshfilecontent << "SI("<< std::scientific;
    gmshfilecontent << prism(0,0) << "," << prism(1,0) << "," << prism(2,0) << ",";
    gmshfilecontent << prism(0,1) << "," << prism(1,1) << "," << prism(2,1) << ",";
    gmshfilecontent << prism(0,2) << "," << prism(1,2) << "," << prism(2,2) << ",";
    gmshfilecontent << prism(0,3) << "," << prism(1,3) << "," << prism(2,3) << ",";
    gmshfilecontent << prism(0,4) << "," << prism(1,4) << "," << prism(2,4) << ",";
    gmshfilecontent << prism(0,5) << "," << prism(1,5) << "," << prism(2,5);
    gmshfilecontent << "){" << std::scientific;
    gmshfilecontent << color << "," << color << "," << color << "," << color << "," << color << "," << color << "};" << std::endl;

    // now the other prisms will be computed
    for (int sector=0;sector<n-1;++sector)
    {
      // initialize for next prism
      // some nodes of last prims can be taken also for the new next prism
      for (int j=0;j<3;++j)
      {
        prism(j,1)=prism(j,2);
        prism(j,4)=prism(j,5);
        prism(j,2)=prism(j,0);
        prism(j,5)=prism(j,3);
      }

      // old radiusvec2 is now radiusvec1; radiusvec2 is set to zero
      for (int j=0;j<3;++j)
      {
        radiusvec1(j) = radiusvec2(j);
        radiusvec2(j) = 0.0;
      }

      // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
      radiusvec2.Multiply(R,radiusvec1);


      // get second point on surface for node1 and node2
      for (int j=0;j<3;++j)
      {
        prism(j,2) += radiusvec2(j) / radiusvec2.Norm2() * radius;
        prism(j,5) += radiusvec2(j) / radiusvec2.Norm2() * radius;
      }

      // put coordinates into filecontent-stream
      // Syntax for gmsh input file
      // SI(x,y--,z,  x+.5,y,z,    x,y+.5,z,   x,y,z+.5, x+.5,y,z+.5, x,y+.5,z+.5){1,2,3,3,2,1};
      // SI( coordinates of the six corners ){colors}
      gmshfilecontent << "SI("<< std::scientific;
      gmshfilecontent << prism(0,0) << "," << prism(1,0) << "," << prism(2,0) << ",";
      gmshfilecontent << prism(0,1) << "," << prism(1,1) << "," << prism(2,1) << ",";
      gmshfilecontent << prism(0,2) << "," << prism(1,2) << "," << prism(2,2) << ",";
      gmshfilecontent << prism(0,3) << "," << prism(1,3) << "," << prism(2,3) << ",";
      gmshfilecontent << prism(0,4) << "," << prism(1,4) << "," << prism(2,4) << ",";
      gmshfilecontent << prism(0,5) << "," << prism(1,5) << "," << prism(2,5);
      gmshfilecontent << "){" << std::scientific;
      gmshfilecontent << color << "," << color << "," << color << "," << color << "," << color << "," << color << "};" << std::endl;
    }
  }
  //no thickness >0 specified; plot elements as line segements without physical volume
  else
  {
    if(IsBarbed)
    {
      double BarbedColor= 0.2;
      gmshfilecontent << "SL(" << std::scientific;
      gmshfilecontent << coord(0,0) << "," << coord(1,0) << "," << coord(2,0) << ","
          << coord(0,1) << "," << coord(1,1) << "," << coord(2,1);
      gmshfilecontent << ")" << "{" << std::scientific << BarbedColor << ","<< BarbedColor << "};" << std::endl;
    }

    {
      gmshfilecontent << "SL(" << std::scientific;
      gmshfilecontent << coord(0,0) << "," << coord(1,0) << "," << coord(2,0) << ","
          << coord(0,1) << "," << coord(1,1) << "," << coord(2,1);
      gmshfilecontent << ")" << "{" << std::scientific << color << ","<< color << "};" << std::endl;
    }
  }
  // crosslink molecules are marked with an additional small ball if they are plotted as volumeless lines
  if(ignoreeleid && drawsphere)
  {
    double beadcolor = color;
    gmshfilecontent << "SP(" << std::scientific;
    gmshfilecontent << coord(0,1) << "," << coord(1,1) << ","<< coord(2,1);
    gmshfilecontent << ")" << "{" << std::scientific << beadcolor << ","<< beadcolor << "};" << std::endl;
  }
  return;
}//GmshWedge

/*----------------------------------------------------------------------*
 |  print Gmsh Triangle to stringstream (private)            grill 03/14|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::PrintGmshTriangleToStream(std::stringstream& gmshfilecontent,
                                                      const std::vector< std::vector< double > > &vertexlist,
                                                      int i, int j, int k, double color,
                                                      const double centercoord[])
{
  // "ST" is scalar triangle, followed by 3x coordinates (x,y,z) of vertices and color
  gmshfilecontent << "ST("<< std::scientific;
  gmshfilecontent << centercoord[0] + vertexlist[i][0] << "," << centercoord[1] + vertexlist[i][1] << "," << centercoord[2] + vertexlist[i][2] << ",";
  gmshfilecontent << centercoord[0] + vertexlist[j][0] << "," << centercoord[1] + vertexlist[j][1] << "," << centercoord[2] + vertexlist[j][2] << ",";
  gmshfilecontent << centercoord[0] + vertexlist[k][0] << "," << centercoord[1] + vertexlist[k][1] << "," << centercoord[2] + vertexlist[k][2];
  gmshfilecontent << "){" << std::scientific;
  gmshfilecontent << color << "," << color << "," << color  << "};" << std::endl << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 |  Refine Icosphere                     (private)            grill 03/14|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GmshRefineIcosphere(std::vector< std::vector<double> > &vertexlist,
                                                  std::vector< std::vector<int> > &facelist,
                                                  double radius)
{
  int num_faces_old = facelist.size();
  std::vector<double> newvertex(3,0.0);
  double scalefac = 0.0;
  std::vector<int> newface(3,0);

  // subdivide each face into four new triangular faces:
  /*               /_\
   *              /_V_\
   */
  for (int i=0; i<num_faces_old; ++i)
  {
    int oldvertices[] = {facelist[i][0],facelist[i][1],facelist[i][2]};

    // compute, normalize and store new vertices in vertexlist
    // subdivide all three edges (all connections of three old vertices)
    for (int j=0; j<3; ++j)
    {
      for (int k=j+1; k<3; ++k)
      {
        newvertex[0]=0.5* (vertexlist[oldvertices[j]][0] + vertexlist[oldvertices[k]][0]);
        newvertex[1]=0.5* (vertexlist[oldvertices[j]][1] + vertexlist[oldvertices[k]][1]);
        newvertex[2]=0.5* (vertexlist[oldvertices[j]][2] + vertexlist[oldvertices[k]][2]);

        // scale new vertex to lie on sphere with given radius
        scalefac = radius / sqrt( pow(newvertex[0],2) + pow(newvertex[1],2) + pow(newvertex[2],2) );
        for (int q=0; q<3; ++q) newvertex[q] *= scalefac;

        vertexlist.push_back(newvertex);
      }
    }

    int len_vertexlist = (int) vertexlist.size();
    // add four new triangles to facelist
    newface[0]=oldvertices[0];  newface[1]=len_vertexlist-3; newface[2]=len_vertexlist-2;
    facelist.push_back(newface);
    newface[0]=oldvertices[1];  newface[1]=len_vertexlist-3; newface[2]=len_vertexlist-1;
    facelist.push_back(newface);
    newface[0]=oldvertices[2];  newface[1]=len_vertexlist-2; newface[2]=len_vertexlist-1;
    facelist.push_back(newface);
    newface[0]=len_vertexlist-3;  newface[1]=len_vertexlist-2; newface[2]=len_vertexlist-1;
    facelist.push_back(newface);
  }

  // erase the old faces
  facelist.erase(facelist.begin(), facelist.begin()+num_faces_old);

  return;
}

/*----------------------------------------------------------------------*
 | Gmsh Output of detected network structure volume        mueller 12/10|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GmshNetworkStructVolume(const int& n, std::stringstream& gmshfilecontent, const double color)
{
  if(DRT::INPUT::IntegralValue<INPAR::STATMECH::StatOutput>(statmechparams_, "SPECIAL_OUTPUT")==INPAR::STATMECH::statout_structanaly)
  {
    std::cout<<"Visualizing test volume: ";
    switch(structuretype_)
    {
      // either cluster or homogeneous network
      case 0:
      {
        std::cout<<"Cluster"<<std::endl;
        // draw three octagons/hexadecagon lying in the base planes
        for(int i=0; i<3; i++)//spatial comp i
          for(int j=0; j<3; j++)//spatial comp j
            if(j<i)
              for(int k=0; k<3; k++)//spatial comp k perp. to the others
                if(k!=i && k!=j)
                {
                  double radius = characlength_/2.0;
                  // some local variables
                  LINALG::SerialDenseMatrix edge(3,2,true);
                  LINALG::Matrix<3,1> radiusvec1;
                  LINALG::Matrix<3,1> radiusvec2;
                  LINALG::Matrix<3,1> auxvec;
                  LINALG::Matrix<3,1> theta;
                  LINALG::Matrix<3,3> R;

                  // compute three dimensional angle theta
                  // unit vector perpendicular to jk-plane
                  LINALG::Matrix<3,1> axis;
                  for (int l=0;l<3;++l)
                    if(l==k)
                      axis(l) = 1.0;
                    else
                      axis(l) = 0.0;
                  for (int l=0;l<3;++l)
                    theta(l) = axis(l) * 2 * M_PI / n;

                  // Compute rotation matrix R from rotation angle theta
                  LARGEROTATIONS::angletotriad(theta,R);

                  // compute radius vector for first surface node of first edges
                  auxvec.Clear();
                  for (int l=0;l<3;++l)
                    if(l==j)
                      auxvec(l) = 1.0;

                  // radiusvector for point on surface
                  radiusvec1(0) = auxvec(1)*axis(2) - auxvec(2)*axis(1);
                  radiusvec1(1) = auxvec(2)*axis(0) - auxvec(0)*axis(2);
                  radiusvec1(2) = auxvec(0)*axis(1) - auxvec(1)*axis(0);

                  // initialize all edge points to nodes
                  for (int l=0;l<3;++l)
                  {
                    edge(l,0) = cog_(l);
                    edge(l,1) = cog_(l);
                  }

                  // get first point on surface for node1 and node2
                  for (int l=0;l<3;++l)
                    edge(l,0) += radiusvec1(l) / radiusvec1.Norm2() * radius;


                  // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
                  radiusvec2.Multiply(R,radiusvec1);

                  // get second point on surface for node1 and node2
                  for(int l=0;l<3;l++)
                    edge(l,1) += radiusvec2(l) / radiusvec2.Norm2() * radius;

                  // shift the coordinates according to periodic boundary conditions
                  LINALG::SerialDenseMatrix edgeshift = edge;
                  for(int l=0; l<(int)edge.M(); l++)
                    for(int m=0; m<(int)edge.N(); m++)
                    {
                      if(edge(l,m)>periodlength_->at(l))
                        edgeshift(l,m) -= periodlength_->at(l);
                      else if(edge(l,m)<0.0)
                        edgeshift(l,m) += periodlength_->at(l);
                    }
                  // write only the edge of the triangle (i.e. the line connecting two corners of the octagon/hexadecagon)
                  GmshOutputPeriodicBoundary(edgeshift, color, gmshfilecontent, discret_->lRowElement(0)->Id(),true, false);

                  // now the other edges will be computed
                  for (int sector=0;sector<n-1;++sector)
                  {
                    // initialize for next edge
                    for (int l=0;l<3;++l)
                    {
                      edge(l,0)=edge(l,1);
                      edge(l,1)=cog_(l);
                    }

                    // old radiusvec2 is now radiusvec1; radiusvec2 is set to zero
                    radiusvec1 = radiusvec2;
                    radiusvec2.Clear();

                    // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
                    radiusvec2.Multiply(R,radiusvec1);

                    // get second point on surface for node1 and node2
                    for (int l=0;l<3;++l)
                      edge(l,1) += radiusvec2(l) / radiusvec2.Norm2() * radius;
                    edgeshift = edge;
                    for(int l=0; l<(int)edge.M(); l++)
                      for(int m=0; m<(int)edge.N(); m++)
                      {
                        if(edge(l,m)>periodlength_->at(l))
                          edgeshift(l,m) -= periodlength_->at(l);
                        else if(edge(l,m)<0.0)
                          edgeshift(l,m) += periodlength_->at(l);
                      }

                    GmshOutputPeriodicBoundary(edgeshift, color, gmshfilecontent, discret_->lRowElement(0)->Id(),true);
                  }
                }
      }
      break;
      // bundle network
      case 1:
      {
        std::cout<<"Bundle"<<std::endl;

        double radius = characlength_;
        LINALG::Matrix<3,2> coord;
        for(int i=0; i<(int)coord.M(); i++)
          for(int j=0; j<(int)coord.N(); j++)
          coord(i,j) = testvolumepos_[j](i);

        // some local variables
        LINALG::Matrix<3,4> edges;
        LINALG::Matrix<3,1> axis;
        LINALG::Matrix<3,1> radiusvec1;
        LINALG::Matrix<3,1> radiusvec2;
        LINALG::Matrix<3,1> auxvec;
        LINALG::Matrix<3,1> theta;
        LINALG::Matrix<3,3> R;

        // compute three dimensional angle theta
        for (int j=0;j<3;++j)
          axis(j) = coord(j,1) - coord(j,0);
        double norm_axis = axis.Norm2();
        for (int j=0;j<3;++j)
          theta(j) = axis(j) / norm_axis * 2 * M_PI / n;

        // Compute rotation matrix R from rotation angle theta
        LARGEROTATIONS::angletotriad(theta,R);

        // compute radius vector for first surface node of first edges
        for (int j=0;j<3;++j) auxvec(j) = coord(j,0) + norm_axis;

        // radiusvector for point on surface
        radiusvec1(0) = auxvec(1)*axis(2) - auxvec(2)*axis(1);
        radiusvec1(1) = auxvec(2)*axis(0) - auxvec(0)*axis(2);
        radiusvec1(2) = auxvec(0)*axis(1) - auxvec(1)*axis(0);

        // initialize all edge points to nodes
        for (int j=0;j<3;++j)
        {
          edges(j,0) = coord(j,0);
          edges(j,1) = coord(j,1);
          edges(j,2) = coord(j,1);
          edges(j,3) = coord(j,0);
        }

        // get first point on surface for node1 and node2
        for (int j=0;j<3;++j)
        {
          edges(j,0) += radiusvec1(j) / radiusvec1.Norm2() * radius;
          edges(j,1) += radiusvec1(j) / radiusvec1.Norm2() * radius;
        }

        // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
        radiusvec2.Multiply(R,radiusvec1);

        // get second point on surface for node1 and node2
        for(int j=0;j<3;j++)
        {
          edges(j,2) += radiusvec2(j) / radiusvec2.Norm2() * radius;
          edges(j,3) += radiusvec2(j) / radiusvec2.Norm2() * radius;
        }
        // write only the edge of the triangle (i.e. the line connecting two corners of the octagon/hexadecagon)
        // and the connecting lines between the prism bases -> rectangle
        for(int j=1; j<(int)edges.N()+1; j++)
        {
          int jlow = (j-1)%(int)edges.N();
          int jhigh= j%(int)edges.N();
          int numsections = 10;
          // original (unshifted) edges
          LINALG::SerialDenseMatrix curredge(3,2);
          for(int k=0; k<curredge.M(); k++)
          {
            curredge(k,0) = edges(k, jlow);
            curredge(k,1) = edges(k, jhigh);
          }

          GmshNetworkStructVolumePeriodic(curredge, numsections, gmshfilecontent,color-0.5);
          /*gmshfilecontent << "SL(" << std::scientific;
          gmshfilecontent << edges(0,ilow) << "," << edges(1,ilow) << "," << edges(2,ilow) << ",";
          gmshfilecontent << edges(0,ihigh) << "," << edges(1,ihigh) << "," << edges(2,ihigh);
          gmshfilecontent << ")" << "{" << std::scientific << color-0.125 << ","<< color-0.125 << "};" << std::endl;*/
        }

        // now the other edges will be computed
        for (int sector=0;sector<n-1;++sector)
        {
          // initialize for next edge
          for (int j=0;j<3;++j)
          {
            edges(j,0)=edges(j,3);
            edges(j,1)=edges(j,2);
            edges(j,2)=coord(j,1);
            edges(j,3)=coord(j,0);
          }

          // old radiusvec2 is now radiusvec1; radiusvec2 is set to zero
          radiusvec1 = radiusvec2;
          radiusvec2.Clear();

          // compute radiusvec2 by rotating radiusvec1 with rotation matrix R
          radiusvec2.Multiply(R,radiusvec1);

          // get second point on surface for node1 and node2
          for (int j=0;j<3;++j)
          {
            edges(j,2) += radiusvec2(j) / radiusvec2.Norm2() * radius;
            edges(j,3) += radiusvec2(j) / radiusvec2.Norm2() * radius;
          }

          // put coordinates into filecontent-stream
          for(int j=1; j<(int)edges.N()+1; j++)
          {
            int jlow = (j-1)%(int)edges.N();
            int jhigh= j%(int)edges.N();
            int numsections = 10;
            LINALG::SerialDenseMatrix curredge(3,2);
            for(int k=0; k<curredge.M(); k++)
            {
              curredge(k,0) = edges(k,jlow);
              curredge(k,1) = edges(k,jhigh);
            }

            GmshNetworkStructVolumePeriodic(curredge, numsections, gmshfilecontent,color-0.5);
            /*gmshfilecontent << "SL(" << std::scientific;
            gmshfilecontent << edges(0,ilow) << "," << edges(1,ilow) << "," << edges(2,ilow) << ",";
            gmshfilecontent << edges(0,ihigh) << "," << edges(1,ihigh) << "," << edges(2,ihigh);
            gmshfilecontent << ")" << "{" << std::scientific << color-0.125 << ","<< color-0.125 << "};" << std::endl;*/
          }
        }
      }
      break;
      // layer
      case 2:
      {
        std::cout<<"Layer ("<<testvolumepos_.size()<<"-noded)"<<std::endl;
        double halfthick = characlength_/2.0;
        // compute normal
        // cross product v_1 x v_2, plane normal
        LINALG::Matrix<3,1> normal;
        LINALG::Matrix<3,1> firstdir = testvolumepos_[1];
        LINALG::Matrix<3,1> secdir = testvolumepos_[1];
        firstdir -= testvolumepos_[0];
        secdir -= testvolumepos_[2];
        normal(0) = firstdir(1)*secdir(2) - firstdir(2)*secdir(1);
        normal(1) = firstdir(2)*secdir(0) - firstdir(0)*secdir(2);
        normal(2) = firstdir(0)*secdir(1) - firstdir(1)*secdir(0);
        normal.Scale(1.0/normal.Norm2());
        // upper and lower bound
        // componentwise comparison between two points to determine whether or not their connection is an edge
        // of the volume to visualize
        for(int i=0; i<(int)testvolumepos_.size(); i++)
          for(int j=0; j<(int)testvolumepos_.size(); j++)
            if(j>i) // consider entries above diagonal only
              for(int k=0; k<3; k++)
                if(fabs((testvolumepos_[i])(k)-(testvolumepos_[j])(k))<1e-7)
                {
                  int numsections = 10;

                  LINALG::SerialDenseMatrix loweredge(3,2);
                  LINALG::SerialDenseMatrix upperedge(3,2);
                  LINALG::SerialDenseMatrix connection0(3,2);
                  LINALG::SerialDenseMatrix connection1(3,2);
                  for(int l=0; l<loweredge.M(); l++)
                  {
                    loweredge(l,0) = testvolumepos_[i](l)-halfthick*normal(l);
                    loweredge(l,1) = testvolumepos_[j](l)-halfthick*normal(l);
                    upperedge(l,0) = testvolumepos_[i](l)+halfthick*normal(l);
                    upperedge(l,1) = testvolumepos_[j](l)+halfthick*normal(l);
                    connection0(l,0) = loweredge(l,0);
                    connection0(l,1) = upperedge(l,0);
                    connection1(l,0) = loweredge(l,1);
                    connection1(l,1) = upperedge(l,1);
                  }
                  // (long) edges
                  GmshNetworkStructVolumePeriodic(loweredge,numsections,gmshfilecontent,color-0.5);
                  GmshNetworkStructVolumePeriodic(upperedge,numsections,gmshfilecontent,color-0.5);
                  // connecting edges
                  GmshNetworkStructVolumePeriodic(connection0,4,gmshfilecontent,color-0.5);
                  GmshNetworkStructVolumePeriodic(connection1,4,gmshfilecontent,color-0.5);
                  // draw continous box also (for now)
                  // upper edge
                  gmshfilecontent << "SL(" << std::scientific;
                  gmshfilecontent << testvolumepos_[i](0)+halfthick*normal(0) << "," << testvolumepos_[i](1)+halfthick*normal(1) << "," << testvolumepos_[i](2)+halfthick*normal(2) << ",";
                  gmshfilecontent << testvolumepos_[j](0)+halfthick*normal(0) << "," << testvolumepos_[j](1)+halfthick*normal(1) << "," << testvolumepos_[j](2)+halfthick*normal(2);
                  gmshfilecontent << ")" << "{" << std::scientific << color-0.25 << ","<< color-0.25 << "};" << std::endl;
                  // lower edge
                  gmshfilecontent << "SL(" << std::scientific;
                  gmshfilecontent << testvolumepos_[i](0)-halfthick*normal(0) << "," << testvolumepos_[i](1)-halfthick*normal(1) << "," << testvolumepos_[i](2)-halfthick*normal(2) << ",";
                  gmshfilecontent << testvolumepos_[j](0)-halfthick*normal(0) << "," << testvolumepos_[j](1)-halfthick*normal(1) << "," << testvolumepos_[j](2)-halfthick*normal(2);
                  gmshfilecontent << ")" << "{" << std::scientific << color-0.25 << ","<< color-0.25 << "};" << std::endl;
                  // connections
                  gmshfilecontent << "SL(" << std::scientific;
                  gmshfilecontent << testvolumepos_[i](0)+halfthick*normal(0) << "," << testvolumepos_[i](1)+halfthick*normal(1) << "," << testvolumepos_[i](2)+halfthick*normal(2) << ",";
                  gmshfilecontent << testvolumepos_[i](0)-halfthick*normal(0) << "," << testvolumepos_[i](1)-halfthick*normal(1) << "," << testvolumepos_[i](2)-halfthick*normal(2);
                  gmshfilecontent << ")" << "{" << std::scientific << color-0.25 << ","<< color-0.25 << "};" << std::endl;
                  gmshfilecontent << "SL(" << std::scientific;
                  gmshfilecontent << testvolumepos_[j](0)+halfthick*normal(0) << "," << testvolumepos_[j](1)+halfthick*normal(1) << "," << testvolumepos_[j](2)+halfthick*normal(2) << ",";
                  gmshfilecontent << testvolumepos_[j](0)-halfthick*normal(0) << "," << testvolumepos_[j](1)-halfthick*normal(1) << "," << testvolumepos_[j](2)-halfthick*normal(2);
                  gmshfilecontent << ")" << "{" << std::scientific << color-0.25 << ","<< color-0.25 << "};" << std::endl;
                }
      }
      break;
      case 3:
        std::cout<<"Homogeneous network"<<std::endl;
      break;
    }
  }
}//GmshNetworkStructVolume()

/*----------------------------------------------------------------------*
 | Gmsh periodic bound. output for test volumes  (protected)mueller 1/11|
 *----------------------------------------------------------------------*/
void STATMECH::StatMechManager::GmshNetworkStructVolumePeriodic(const Epetra_SerialDenseMatrix& coord, const int numsections, std::stringstream& gmshfilecontent,const double color)
{
  // direction of edge

  LINALG::Matrix<3, 1> dir;
  for(int k=0; k<(int)dir.M(); k++)
    dir(k) = coord(k,1)-coord(k,0);
  dir.Scale(1.0/dir.Norm2());
  // divide given edge into intervals for periodic boundary shift and visualization

  for(int k=0; k<numsections; k++)
  {
    // k-th intersection point
    LINALG::SerialDenseMatrix linepart(3,2);
    for(int l=0; l<(int)linepart.M(); l++)
    {
      linepart(l,0) = coord(l,0)+(coord(l,1)-coord(l,0))/(double)numsections*(double)k;
      linepart(l,1) = coord(l,0)+(coord(l,1)-coord(l,0))/(double)numsections*(double)(k+1);
    }

    // determine how the line is broken
    // Periodic boundary shift
    for(int l=0; l<linepart.M(); l++)
      for(int m=0; m<linepart.N(); m++)
      {
        if(linepart(l,m)>periodlength_->at(l))
          while(linepart(l,m)>periodlength_->at(l))
            linepart(l,m) -= periodlength_->at(l);
        if(linepart(l,m)<0.0)
          while(linepart(l,m)<0.0)
            linepart(l,m) += periodlength_->at(l);
      }

    LINALG::Matrix<3,1> cut;
    cut.Clear();
    for (int dof=0; dof<linepart.M(); dof++)
    {
      // point 0 close to zero boundary
      if (fabs(linepart(dof,1)-periodlength_->at(dof)-linepart(dof,0)) < fabs(linepart(dof,1) - linepart(dof,0)))
        cut(dof) = 1.0;
      // point 1 close to zero boundary
      if (fabs(linepart(dof,1)+periodlength_->at(dof) - linepart(dof,0)) < fabs(linepart(dof,1)-linepart(dof,0)))
        cut(dof) = 2.0;
    }

    // write special output for broken elements
    if (cut(0) + cut(1) + cut(2) > 0.0)
    {
      //from node 0 to nearest boundary where element is broken you get by vector X + lambda0*dir
      double lambda0 = 1e4;
      for (int dof = 0; dof < linepart.M(); dof++)
      {
        if (cut(dof) == 1.0)
        {
          if (fabs(-linepart(dof, 0) / dir(dof)) < fabs(lambda0))
            lambda0 = -linepart(dof, 0) / dir(dof);
        }
        else if (cut(dof) == 2.0)
        {
          if (fabs((periodlength_->at(dof) - linepart(dof, 0)) / dir(dof)) < fabs(lambda0))
            lambda0 = (periodlength_->at(dof) - linepart(dof, 0)) / dir(dof);
        }
      }
      //from node 1 to nearest boundary where element is broken you get by vector X + lambda1*dir
      double lambda1 = 1e4;
      for (int dof = 0; dof < linepart.M(); dof++)
      {
        if (cut(dof) == 2.0)
        {
          if (fabs(-linepart(dof, 1) / dir(dof)) < fabs(lambda1))
            lambda1 = -linepart(dof, 1) / dir(dof);
        }
        else if (cut(dof) == 1.0)
        {
          if (fabs((periodlength_->at(dof) - linepart(dof, 1)) / dir(dof)) < fabs(lambda1))
            lambda1 = (periodlength_->at(dof) - linepart(dof, 1)) / dir(dof);
        }
      }
      //define output coordinates for broken elements, first segment
      LINALG::SerialDenseMatrix coordout=linepart;
      for(int l=0 ;l<coordout.M(); l++)
        coordout(l,1) = linepart(l,0) + lambda0*dir(l);
      GmshWedge(1,coordout,discret_->lRowElement(0),gmshfilecontent,color,true);

      //define output coordinates for broken elements, second segment
      for(int l=0; l<coordout.M(); l++)
      {
        coordout(l,0) = linepart(l,1);
        coordout(l,1) = linepart(l,1)+lambda1*dir(l);
      }
      GmshWedge(1,coordout,discret_->lRowElement(0),gmshfilecontent,color,true);
    }
    else // output for continuous elements
      GmshWedge(1,linepart,discret_->lRowElement(0),gmshfilecontent,color,true,false);
  }
}//GmshNetworkStructVolumePeriodic()

/*----------------------------------------------------------------------*
 | Output Point Neumann force nodes              (private) mueller 09/14|
 *----------------------------------------------------------------------*/
void  STATMECH::StatMechManager::GmshOutputPointNeumann(const Epetra_Vector&      disrow,
                                                        const std::ostringstream* filename,
                                                        const double&             time,
                                                        const double&             color)
{
  // do output to file in c-style
  if(!nbcnodesets_.empty())
  {
    // export row displacement to column map format
    Epetra_Vector discol(*(discret_->DofColMap()), true);
    LINALG::Export(disrow, discol);

    //In case of the 4-noded Crosslinker Element we need the following vectors
    Teuchos::RCP<Epetra_MultiVector> nodepositions = Teuchos::rcp(new Epetra_MultiVector(*discret_->NodeColMap(),3,true));
    GetNodalBindingSpotPositionsFromDisVec(discol, nodepositions, Teuchos::null);

    FILE* fp = NULL;
    fp = fopen(filename->str().c_str(), "a");
    std::stringstream gmshfilecontent;

    double beadcolor = 7.0*color;
    int nodesetindex = timeintervalstep_-bctimeindex_;
    if((int)nbcnodesets_.size()==1)
      nodesetindex = 0;
    else
    {
      if(timeintervalstep_>1)
        nodesetindex--;
    }

    INPAR::STATMECH::NBCType nbctype = DRT::INPUT::IntegralValue<INPAR::STATMECH::NBCType>(statmechparams_, "NBCTYPE");
    switch(nbctype)
    {
      case INPAR::STATMECH::nbctype_constcreep:
      break;
      case INPAR::STATMECH::nbctype_randompointforce:
      {
        //writing element by nodal coordinates as a sphere
        if(discret_->NodeRowMap()->MyGID(nbcnodesets_[nodesetindex][0]))
        {
          LINALG::Matrix<3,1> position(true);
          for(int i=0; i<(int)position.M(); i++)
            position(i) = (*nodepositions)[i][bspotcolmap_->LID(nbcnodesets_[nodesetindex][0])];
          std::vector<double> boxsize(3, 1.5*statmechparams_.get<double>("R_LINK",0.1));
          GmshOutputBox(beadcolor,&position,boxsize,filename, false, discret_->Comm().MyPID());
          gmshfilecontent << "SP(" << std::scientific;
          gmshfilecontent<< position(0)<< ","
                         << position(1) << ","
                         << position(2) << ")" << "{" << std::scientific << beadcolor << "," << beadcolor << "};" << std::endl;
          // direction and amplitude of the force
          int oscdir = statmechparams_.get<int>("DBCDISPDIR",0)-1;
          int curvenum = statmechparams_.get<int>("NBCCURVENUMBER", 0)-1;
          double nbcamp = statmechparams_.get<double>("NBCFORCEAMP",0.0);
          LINALG::Matrix<3,1> tip(position);
          tip(oscdir) += 0.25*nbcamp*(DRT::Problem::Instance()->Curve(curvenum).f(time));
          gmshfilecontent << "SL(" << std::scientific;
          gmshfilecontent << position(0) << "," << position(1) << "," << position(2) << "," << tip(0) << ","<< tip(1) << "," << tip(2);
          gmshfilecontent << ")" << "{" << std::scientific << beadcolor << ","<< beadcolor << "};" << std::endl;
          std::vector<double> tipsize(3,0.025);
          GmshOutputBox(beadcolor, &tip, tipsize, filename, false, discret_->Comm().MyPID());
        }
      }
      break;
      default: break;
    }
    fputs(gmshfilecontent.str().c_str(), fp);
    fclose(fp);
  }
  return;
}

/*----------------------------------------------------------------------*
 | check the visualization mode of a crosslinker depending on the       |
 | filament(s) it is bound to.                    (public) mueller 07/10|
 *----------------------------------------------------------------------*/
bool STATMECH::StatMechManager::CheckForKinkedVisual(int eleid)
{
  // if element is a crosslinker
  if(eleid>basisnodes_)
  {
    DRT::Element *element = discret_->lColElement(discret_->ElementColMap()->LID(eleid));
    // note: choose first and last node since these two nodes definitely tell how the crosslinker element is bound
    int nodelid0 = element->NodeIds()[0];
    int nodelid1 = element->NodeIds()[element->NumNode()-1];
    if((*filamentnumber_)[nodelid0]==(*filamentnumber_)[nodelid1])
      return true;
  }
  return false;
}//STATMECH::StatMechManager::CheckForKinkedVisual
