/*----------------------------------------------------------------------*/
/*! \file

\brief Create parametrized structured grid.

\level 0


*/
/*----------------------------------------------------------------------*/


#include "lib_gridgenerator.H"
#include "lib_elementdefinition.H"
#include "rebalance.H"
#include "lib_utils_factory.H"
#include "lib_utils_parallel.H"
#include "lib_discret.H"
#include "io_pstream.H"

#include <algorithm>

namespace DRT
{
  namespace GRIDGENERATOR
  {
    // forward declarations
    Teuchos::RCP<DRT::Element> CreateHexElement(int eleid, int nodeoffset, int myrank,
        DRT::INPUT::LineDefinition* linedef, std::array<int, 3> interval, std::string elementtype,
        std::string distype);

    Teuchos::RCP<DRT::Element> CreateWedgeElement(int eleid, int nodeoffset, int myrank,
        DRT::INPUT::LineDefinition* linedef, std::array<int, 3> interval, std::string elementtype,
        std::string distype);

    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    void CreateRectangularCuboidDiscretization(DRT::Discretization& dis,
        const DRT::GRIDGENERATOR::RectangularCuboidInputs& inputData, bool outputFlag)
    {
      const Epetra_Comm& comm = dis.Comm();
      const int myrank = comm.MyPID();
      const int numproc = comm.NumProc();

      DRT::INPUT::ElementDefinition ed;
      ed.SetupValidElementLines();

      // safety checks
      for (int i = 0; i < 3; ++i)
      {
        if (inputData.bottom_corner_point_[i] >= inputData.top_corner_point_[i])
          dserror("lower bound in domain reader must be smaller than upper bound");

        if (inputData.interval_[i] <= 0)
          dserror("intervals in domain reader must be greater than zero");
      }

      Teuchos::RCP<Epetra_Map> nodeRowMap;
      Teuchos::RCP<Epetra_Map> nodeColMap;
      Teuchos::RCP<Epetra_Map> elementRowMap;
      Teuchos::RCP<Epetra_Map> elementColMap;

      // Create initial (or final) map of row elements
      DRT::Element::DiscretizationType distype_enum = DRT::StringToDistype(inputData.distype_);
      int numnewele = inputData.interval_[0] * inputData.interval_[1] * inputData.interval_[2];
      if (inputData.autopartition_)  // linear map
      {
        int scale = 1;
        if (distype_enum == DRT::Element::wedge6 or distype_enum == DRT::Element::wedge15)
        {
          scale = 2;
        }
        elementRowMap = Teuchos::rcp(new Epetra_Map(scale * numnewele, 0, comm));
      }
      else  // fancy final box map
      {
        // Error for invalid element types!!!
        if (distype_enum != DRT::Element::hex8 and distype_enum != DRT::Element::hex20 and
            distype_enum != DRT::Element::hex27)
        {
          dserror("This map-partition is only available for HEX-elements!");
        }

        std::vector<int> factors;
        int nproc = numproc;
        for (int fac = 2; fac < nproc + 1;)
        {
          if (nproc % fac == 0)
          {
            factors.push_back(fac);
            nproc /= fac;
          }
          else
          {
            fac++;
          }
        }
        if (nproc != 1) dserror("Could not split numproc.");

        unsigned int subdivisions[] = {1, 1, 1};
        const double dinterval[] = {static_cast<double>(inputData.interval_[0]),
            static_cast<double>(inputData.interval_[1]),
            static_cast<double>(inputData.interval_[2])};
        for (std::vector<int>::const_reverse_iterator fac = factors.rbegin(); fac != factors.rend();
             ++fac)
        {
          const double ratios[] = {dinterval[0] / subdivisions[0], dinterval[1] / subdivisions[1],
              dinterval[2] / subdivisions[2]};
          if (ratios[0] >= ratios[1] && ratios[0] >= ratios[2])
            subdivisions[0] *= *fac;
          else if (ratios[1] >= ratios[0] && ratios[1] >= ratios[2])
            subdivisions[1] *= *fac;
          else if (ratios[2] >= ratios[0] && ratios[2] >= ratios[1])
            subdivisions[2] *= *fac;
        }

        if (myrank == 0 && outputFlag)
          IO::cout << "Determined domain subdivision to: " << subdivisions[0] << "x"
                   << subdivisions[1] << "x" << subdivisions[2] << "\n";

        unsigned int xranges[subdivisions[0] + 1ul];
        for (size_t i = 0; i < subdivisions[0] + 1ul; ++i)
          xranges[i] =
              std::max(0, std::min(inputData.interval_[0],
                              static_cast<int>(round(i * dinterval[0] / subdivisions[0]))));

        unsigned int yranges[subdivisions[1] + 1ul];
        for (size_t i = 0; i < subdivisions[1] + 1ul; ++i)
          yranges[i] =
              std::max(0, std::min(inputData.interval_[1],
                              static_cast<int>(round(i * dinterval[1] / subdivisions[1]))));

        unsigned int zranges[subdivisions[2] + 1ul];
        for (size_t i = 0; i < subdivisions[2] + 1ul; ++i)
          zranges[i] =
              std::max(0, std::min(inputData.interval_[2],
                              static_cast<int>(round(i * dinterval[2] / subdivisions[2]))));

        const unsigned int mysection[] = {myrank % subdivisions[0],
            (myrank / subdivisions[0]) % subdivisions[1],
            myrank / (subdivisions[0] * subdivisions[1])};
        const int nummynewele = (xranges[mysection[0] + 1] - xranges[mysection[0]]) *
                                (yranges[mysection[1] + 1] - yranges[mysection[1]]) *
                                (zranges[mysection[2] + 1] - zranges[mysection[2]]);
        int mynewele[nummynewele];

        size_t idx = 0;
        for (size_t iz = zranges[mysection[2]]; iz < zranges[mysection[2] + 1]; ++iz)
          for (size_t iy = yranges[mysection[1]]; iy < yranges[mysection[1] + 1]; ++iy)
            for (size_t ix = xranges[mysection[0]]; ix < xranges[mysection[0] + 1]; ++ix)
              mynewele[idx++] = (iz * inputData.interval_[1] + iy) * inputData.interval_[0] + ix;

        elementRowMap = Teuchos::rcp(new Epetra_Map(-1, nummynewele, mynewele, 0, comm));
      }

      // Create the actual elements according to the row map
      for (int lid = 0; lid < elementRowMap->NumMyElements(); ++lid)
      {
        int eleid = elementRowMap->GID(lid);
        dsassert(eleid >= 0, "Missing gid");

        // For the time being we support old and new input facilities. To
        // smooth transition.
        DRT::INPUT::LineDefinition* linedef =
            ed.ElementLines(inputData.elementtype_, inputData.distype_);
        if (linedef == NULL)
          dserror("a matching line definition is needed for %s %s", inputData.elementtype_.c_str(),
              inputData.distype_.c_str());

        std::istringstream eleargstream(inputData.elearguments_);
        if (not linedef->Read(eleargstream, &inputData.distype_))
        {
          IO::cout << "\n"
                   << eleid << " " << inputData.elementtype_ << " " << inputData.distype_ << " ";
          linedef->Print(IO::cout.cout_replacement());
          IO::cout << "\n";
          dserror("failed to read element %d %s %s", eleid, inputData.elementtype_.c_str(),
              inputData.distype_.c_str());
        }

        // Create specified elemnts
        switch (distype_enum)
        {
          case DRT::Element::hex8:
          case DRT::Element::hex20:
          case DRT::Element::hex27:
          {
            Teuchos::RCP<DRT::Element> ele =
                CreateHexElement(eleid, inputData.node_gid_of_first_new_node_, myrank, linedef,
                    inputData.interval_, inputData.elementtype_, inputData.distype_);
            // add element to discretization
            dis.AddElement(ele);
            break;
          }
          case DRT::Element::wedge6:
          case DRT::Element::wedge15:
          {
            Teuchos::RCP<DRT::Element> ele = DRT::GRIDGENERATOR::CreateWedgeElement(eleid,
                inputData.node_gid_of_first_new_node_, myrank, linedef, inputData.interval_,
                inputData.elementtype_, inputData.distype_);
            dis.AddElement(ele);
            break;
          }
          default:
            dserror(
                "The discretization type %s, is not implemented. Currently only HEX(8,20,27) and "
                "WEDGE(6,15) are implemented for the box geometry generation.",
                inputData.distype_.c_str());
        }
      }

      // redistribute the elements
      if (inputData.autopartition_)
      {
        std::tie(nodeRowMap, nodeColMap) =
            REBALANCE::RebalanceNodeMaps(Teuchos::rcp(&dis, false), elementRowMap, comm.NumProc());
      }
      else  // do not destroy our manual partitioning
      {
        Teuchos::RCP<const Epetra_CrsGraph> graph =
            REBALANCE::BuildGraph(Teuchos::rcp(&dis, false), elementRowMap);
        nodeRowMap = Teuchos::rcp(new Epetra_Map(
            -1, graph->RowMap().NumMyElements(), graph->RowMap().MyGlobalElements(), 0, comm));
        nodeColMap = Teuchos::rcp(new Epetra_Map(
            -1, graph->ColMap().NumMyElements(), graph->ColMap().MyGlobalElements(), 0, comm));
      }


      // now we have all elements in a linear map roweles
      // build reasonable maps for elements from the
      // already valid and final node maps
      // note that nothing is actually redistributed in here
      dis.BuildElementRowColumn(*nodeRowMap, *nodeColMap, elementRowMap, elementColMap);

      // we can now export elements to resonable row element distribution
      dis.ExportRowElements(*elementRowMap);

      // export to the column map / create ghosting of elements
      dis.ExportColumnElements(*elementColMap);

      // Create the nodes according to their elements
      // number of nodes per direction
      const size_t nx = 2 * inputData.interval_[0] + 1;
      const size_t ny = 2 * inputData.interval_[1] + 1;
      int maxgid = -1;

      // as we are using the redistributed row node map, the nodes are directly created on the
      // correct processors

      // Compute midpoint for rotations of the box geometry
      std::vector<double> coordm(3, 0.0);
      if (inputData.rotation_angle_[0] != 0.0 || inputData.rotation_angle_[1] != 0.0 ||
          inputData.rotation_angle_[2] != 0.0)
      {
        coordm[0] = (inputData.top_corner_point_[0] + inputData.bottom_corner_point_[0]) / 2.;
        coordm[1] = (inputData.top_corner_point_[1] + inputData.bottom_corner_point_[1]) / 2.;
        coordm[2] = (inputData.top_corner_point_[2] + inputData.bottom_corner_point_[2]) / 2.;
      }

      for (int lid = 0; lid < nodeRowMap->NumMyElements(); ++lid)
      {
        const int gid = nodeRowMap->GID(lid);
        maxgid = std::max(gid, maxgid);

        const int posid = gid - inputData.node_gid_of_first_new_node_;
        dsassert(posid >= 0, "Tried to access a node gid that was not on this proc");
        size_t i = posid % nx;
        size_t j = (posid / nx) % ny;
        size_t k = posid / (nx * ny);

        double coords[3];
        coords[0] = static_cast<double>(i) / (2 * inputData.interval_[0]) *
                        (inputData.top_corner_point_[0] - inputData.bottom_corner_point_[0]) +
                    inputData.bottom_corner_point_[0];
        coords[1] = static_cast<double>(j) / (2 * inputData.interval_[1]) *
                        (inputData.top_corner_point_[1] - inputData.bottom_corner_point_[1]) +
                    inputData.bottom_corner_point_[1];
        coords[2] = static_cast<double>(k) / (2 * inputData.interval_[2]) *
                        (inputData.top_corner_point_[2] - inputData.bottom_corner_point_[2]) +
                    inputData.bottom_corner_point_[2];

        // If set perform rotations, applied in the order, x,y,z-axis
        for (int rotaxis = 0; rotaxis < 3; ++rotaxis)
        {
          if (inputData.rotation_angle_[rotaxis] != 0.0)
          {
            // add rotation around mitpoint here.
            double dx[3];
            dx[0] = coords[0] - coordm[0];
            dx[1] = coords[1] - coordm[1];
            dx[2] = coords[2] - coordm[2];

            double calpha = cos(inputData.rotation_angle_[rotaxis] * M_PI / 180);
            double salpha = sin(inputData.rotation_angle_[rotaxis] * M_PI / 180);

            coords[0] = coordm[0];  //+ calpha*dx[0] + salpha*dx[1];
            coords[1] = coordm[1];  //+ -salpha*dx[0] + calpha*dx[1];
            coords[2] = coordm[2];

            coords[(rotaxis + 1) % 3] +=
                calpha * dx[(rotaxis + 1) % 3] + salpha * dx[(rotaxis + 2) % 3];
            coords[(rotaxis + 2) % 3] +=
                calpha * dx[(rotaxis + 2) % 3] - salpha * dx[(rotaxis + 1) % 3];
            coords[rotaxis] += dx[rotaxis];
          }
        }

        Teuchos::RCP<DRT::Node> node = Teuchos::rcp(new DRT::Node(gid, coords, myrank));
        dis.AddNode(node);
      }
      dis.ExportColumnNodes(*nodeColMap);
    }

    /*----------------------------------------------------------------------*
     | create HEX type elements for the partition                           |
     *----------------------------------------------------------------------*/
    Teuchos::RCP<DRT::Element> CreateHexElement(int eleid, int nodeOffset, int myrank,
        DRT::INPUT::LineDefinition* linedef, std::array<int, 3> interval, std::string elementtype,
        std::string distype)
    {
      // Reserve nodeids for this element type
      std::vector<int> nodeids(DRT::UTILS::getNumberOfElementNodes(DRT::StringToDistype(distype)));

      // current element position
      const size_t ex = 2 * (eleid % interval[0]);
      const size_t ey = 2 * ((eleid / interval[0]) % interval[1]);
      const size_t ez = 2 * (eleid / (interval[0] * interval[1]));

      // number of nodes per direction
      const size_t nx = 2 * interval[0] + 1;
      const size_t ny = 2 * interval[1] + 1;

      switch (nodeids.size())
      {
        case 27:
          nodeids[20] = nodeOffset + (ez * ny + ey + 1) * nx + ex + 1;
          nodeids[21] = nodeOffset + ((ez + 1) * ny + ey) * nx + ex + 1;
          nodeids[22] = nodeOffset + ((ez + 1) * ny + ey + 1) * nx + ex + 2;
          nodeids[23] = nodeOffset + ((ez + 1) * ny + ey + 2) * nx + ex + 1;
          nodeids[24] = nodeOffset + ((ez + 1) * ny + ey + 1) * nx + ex;
          nodeids[25] = nodeOffset + ((ez + 2) * ny + ey + 1) * nx + ex + 1;
          nodeids[26] = nodeOffset + ((ez + 1) * ny + ey + 1) * nx + ex + 1;
          [[fallthrough]];
        case 20:
          nodeids[8] = nodeOffset + (ez * ny + ey) * nx + ex + 1;
          nodeids[9] = nodeOffset + (ez * ny + ey + 1) * nx + ex + 2;
          nodeids[10] = nodeOffset + (ez * ny + ey + 2) * nx + ex + 1;
          nodeids[11] = nodeOffset + (ez * ny + ey + 1) * nx + ex;
          nodeids[12] = nodeOffset + ((ez + 1) * ny + ey) * nx + ex;
          nodeids[13] = nodeOffset + ((ez + 1) * ny + ey) * nx + ex + 2;
          nodeids[14] = nodeOffset + ((ez + 1) * ny + ey + 2) * nx + ex + 2;
          nodeids[15] = nodeOffset + ((ez + 1) * ny + ey + 2) * nx + ex;
          nodeids[16] = nodeOffset + ((ez + 2) * ny + ey) * nx + ex + 1;
          nodeids[17] = nodeOffset + ((ez + 2) * ny + ey + 1) * nx + ex + 2;
          nodeids[18] = nodeOffset + ((ez + 2) * ny + ey + 2) * nx + ex + 1;
          nodeids[19] = nodeOffset + ((ez + 2) * ny + ey + 1) * nx + ex;
          [[fallthrough]];
        case 8:
          nodeids[0] = nodeOffset + (ez * ny + ey) * nx + ex;
          nodeids[1] = nodeOffset + (ez * ny + ey) * nx + ex + 2;
          nodeids[2] = nodeOffset + (ez * ny + ey + 2) * nx + ex + 2;
          nodeids[3] = nodeOffset + (ez * ny + ey + 2) * nx + ex;
          nodeids[4] = nodeOffset + ((ez + 2) * ny + ey) * nx + ex;
          nodeids[5] = nodeOffset + ((ez + 2) * ny + ey) * nx + ex + 2;
          nodeids[6] = nodeOffset + ((ez + 2) * ny + ey + 2) * nx + ex + 2;
          nodeids[7] = nodeOffset + ((ez + 2) * ny + ey + 2) * nx + ex;
          break;
        default:
          dserror("The number of nodeids: %d, does not correspond to a supported HEX-element.",
              nodeids.size());
          break;
      }
      // let the factory create a matching empty element
      Teuchos::RCP<DRT::Element> ele = DRT::UTILS::Factory(elementtype, distype, eleid, myrank);
      ele->SetNodeIds(nodeids.size(), &(nodeids[0]));
      ele->ReadElement(elementtype, distype, linedef);
      return ele;
    }

    /*----------------------------------------------------------------------*
     | Create WEDGE type elements for the partition.                        |
     | For even eleids -> create 1st part of HEX equivalent, odd -> 2nd     |
     | part of HEX equivalent.                                              |
     | Wedges aligned in z-direction.                                       |
     *----------------------------------------------------------------------*/
    Teuchos::RCP<DRT::Element> CreateWedgeElement(int eleid, int nodeoffset, int myrank,
        DRT::INPUT::LineDefinition* linedef, std::array<int, 3> interval, std::string elementtype,
        std::string distype)
    {
      // Reserve nodeids for this element type
      std::vector<int> nodeids(DRT::UTILS::getNumberOfElementNodes(DRT::StringToDistype(distype)));

      // HEX-equivalent element
      int hex_equiv_eleid = int(eleid / 2);

      // current element position
      const size_t ex = 2 * (hex_equiv_eleid % interval[0]);
      const size_t ey = 2 * ((hex_equiv_eleid / interval[0]) % interval[1]);
      const size_t ez = 2 * (hex_equiv_eleid / (interval[0] * interval[1]));

      // number of nodes per direction
      const size_t nx = 2 * interval[0] + 1;
      const size_t ny = 2 * interval[1] + 1;

      // Create 2 elements for every hex element. Even - Odd pairs.
      if (eleid % 2 == 0)  // Even - elements
      {
        switch (nodeids.size())
        {
          case 15:
            nodeids[6] = nodeoffset + (ez * ny + ey) * nx + ex + 1;             // HEX-eqvi: 8
            nodeids[8] = nodeoffset + (ez * ny + ey + 1) * nx + ex;             // HEX-eqvi: 11
            nodeids[9] = nodeoffset + ((ez + 1) * ny + ey) * nx + ex;           // HEX-eqvi: 12
            nodeids[10] = nodeoffset + ((ez + 1) * ny + ey) * nx + ex + 2;      // HEX-eqvi: 13
            nodeids[11] = nodeoffset + ((ez + 1) * ny + ey + 2) * nx + ex;      // HEX-eqvi: 15
            nodeids[12] = nodeoffset + ((ez + 2) * ny + ey) * nx + ex + 1;      // HEX-eqvi: 16
            nodeids[14] = nodeoffset + ((ez + 2) * ny + ey + 1) * nx + ex;      // HEX-eqvi: 19
            nodeids[7] = nodeoffset + (ez * ny + ey + 1) * nx + ex + 1;         // HEX-eqvi: 20
            nodeids[13] = nodeoffset + ((ez + 2) * ny + ey + 1) * nx + ex + 1;  // HEX-eqvi: 25
            [[fallthrough]];
          case 6:
            nodeids[0] = nodeoffset + (ez * ny + ey) * nx + ex;            // HEX-eqvi: 0
            nodeids[1] = nodeoffset + (ez * ny + ey) * nx + ex + 2;        // HEX-eqvi: 1
            nodeids[2] = nodeoffset + (ez * ny + ey + 2) * nx + ex;        // HEX-eqvi: 3
            nodeids[3] = nodeoffset + ((ez + 2) * ny + ey) * nx + ex;      // HEX-eqvi: 4
            nodeids[4] = nodeoffset + ((ez + 2) * ny + ey) * nx + ex + 2;  // HEX-eqvi: 5
            nodeids[5] = nodeoffset + ((ez + 2) * ny + ey + 2) * nx + ex;  // HEX-eqvi: 7
            break;
            //---------------------
          default:
            dserror("The number of nodeids: %d, does not correspond to a supported WEDGE-element.",
                nodeids.size());
            break;
        }
      }
      else  // Odd - elements
      {
        switch (nodeids.size())
        {
          case 15:
            nodeids[6] = nodeoffset + (ez * ny + ey + 1) * nx + ex + 2;         // HEX-eqvi: 9
            nodeids[7] = nodeoffset + (ez * ny + ey + 2) * nx + ex + 1;         // HEX-eqvi: 10
            nodeids[9] = nodeoffset + ((ez + 1) * ny + ey) * nx + ex + 2;       // HEX-eqvi: 13
            nodeids[10] = nodeoffset + ((ez + 1) * ny + ey + 2) * nx + ex + 2;  // HEX-eqvi: 14
            nodeids[11] = nodeoffset + ((ez + 1) * ny + ey + 2) * nx + ex;      // HEX-eqvi: 15
            nodeids[12] = nodeoffset + ((ez + 2) * ny + ey + 1) * nx + ex + 2;  // HEX-eqvi: 17
            nodeids[13] = nodeoffset + ((ez + 2) * ny + ey + 2) * nx + ex + 1;  // HEX-eqvi: 18
            nodeids[8] = nodeoffset + (ez * ny + ey + 1) * nx + ex + 1;         // HEX-eqvi: 20
            nodeids[14] = nodeoffset + ((ez + 2) * ny + ey + 1) * nx + ex + 1;  // HEX-eqvi: 25
          case 6:
            nodeids[0] = nodeoffset + (ez * ny + ey) * nx + ex + 2;            // HEX-eqvi: 1
            nodeids[1] = nodeoffset + (ez * ny + ey + 2) * nx + ex + 2;        // HEX-eqvi: 2
            nodeids[2] = nodeoffset + (ez * ny + ey + 2) * nx + ex;            // HEX-eqvi: 3
            nodeids[3] = nodeoffset + ((ez + 2) * ny + ey) * nx + ex + 2;      // HEX-eqvi: 5
            nodeids[4] = nodeoffset + ((ez + 2) * ny + ey + 2) * nx + ex + 2;  // HEX-eqvi: 6
            nodeids[5] = nodeoffset + ((ez + 2) * ny + ey + 2) * nx + ex;      // HEX-eqvi: 7
            break;
            //---------------------
          default:
            dserror("The number of nodeids: %d, does not correspond to a supported WEDGE-element.",
                nodeids.size());
            break;
        }
      }

      // let the factory create a matching empty element
      Teuchos::RCP<DRT::Element> ele = DRT::UTILS::Factory(elementtype, distype, eleid, myrank);
      ele->SetNodeIds(nodeids.size(), &(nodeids[0]));
      ele->ReadElement(elementtype, distype, linedef);
      return ele;
    }

  }  // namespace GRIDGENERATOR
}  // namespace DRT
