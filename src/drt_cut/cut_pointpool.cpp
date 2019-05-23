/*---------------------------------------------------------------------*/
/*!
\file cut_pointpool.cpp

\brief PointPool, stores a points in the cut and decides if points are merged or new points are
created

\level 3

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249

*----------------------------------------------------------------------*/

#include "cut_pointpool.H"
#include "cut_output.H"
#include "../drt_lib/drt_globalproblem.H"


/*-----------------------------------------------------------------------------------------*
 * If a point with the coordinates "x" does not exists, it creates a new point correspondingly
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::Point* GEO::CUT::OctTreeNode::NewPoint(const double* x, Edge* cut_edge, Side* cut_side,
    double tolerance, Pointpool_MergeStrategy merge_strategy)
{
  // check if the point already exists
#if CUT_CREATION_INFO
  bool new_point = false;
#endif
  LINALG::Matrix<3, 1> px(x);

  Point* p = GetPoint(x, cut_edge, cut_side, tolerance, merge_strategy);

  if (p == NULL)
  {
    p = &*CreatePoint(points_.size(), x, cut_edge, cut_side, tolerance);  // create the point
#if CUT_CREATION_INFO
    new_point = true;
#endif

#if 1
    if (points_.size() % 1000 == 0)  // split the node starting from level 0
    {
      Split(0);
    }
#endif
  }


#if CUT_CREATION_INFO
  // if it was merged
  if (not new_point)
  {
    std::stringstream info;
    info << "// Another point was merged in this one with\n";
    info << "// Initial coordinates" << std::setprecision(15) << px << std::endl;
    // merged_to point coordinates
    LINALG::Matrix<3, 1> nx;
    p->Coordinates(nx.A());
    px.Update(-1, nx, 1);
    info << "// Merged with tolerance of " << std::setprecision(15) << px.Norm2() << std::endl;
    p->AddAdditionalCreationInfo(info.str());
  }

  else
  {
    p->AddAdditionalCreationInfo(
        "// This point was really created with the exact coordinates given");
  }
#endif

  return p;
}

/*-----------------------------------------------------------------------------------------*
 * Get the point with the specified coordinates "x" from the pointpool
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::Point* GEO::CUT::OctTreeNode::GetPoint(const double* x, Edge* cut_edge, Side* cut_side,
    double tolerance, Pointpool_MergeStrategy merge_strategy)
{
  // try to find the point in all of the 8 children nodes
  if (not IsLeaf())
  {
    // stop finding the point when point is not included in the current bounding box
    if (!bb_->Within(1.0, x)) return NULL;

    for (int i = 0; i < 8; ++i)
    {
      Point* p = nodes_[i]->GetPoint(x, cut_edge, cut_side, tolerance, merge_strategy);
      if (p != NULL)
      {
        return p;
      }
    }
  }
  else
  {
    LINALG::Matrix<3, 1> px(x);
    LINALG::Matrix<3, 1> nx;

    double tol = TOPOLOGICAL_TOLERANCE * norm_;


    switch (merge_strategy)
    {
      case Pointpool_MergeStrategy::SelfCutLoad:
      {
        tol = TOPOLOGICAL_TOLERANCE * norm_ * NODAL_POINT_TOLERANCE_SELFCUT_SCALE;
        break;
      }
      case Pointpool_MergeStrategy::InitialLoad:
      {
        // when we are loading the geometry into the cut  we want the distance to be large than
        // topological tolerance. Otherwise we might experience huge problems
        tol = TOPOLOGICAL_TOLERANCE * norm_ * NODAL_POINT_TOLERANCE_SCALE;
        // safety check
#ifdef NODAL_POINT_TOLERANCE_SCALE
        // this should happen with both cut_side and cut_edge equal to NULL
        // note: be carefull with other cases, apart from mesh loading
        if (cut_side and cut_edge)
        {
          dserror("Scaling is %lf for non-NULL cut_side and cut_edge. This should not be possible!",
              NODAL_POINT_TOLERANCE_SCALE);
        }
#endif
        break;
      }
      case Pointpool_MergeStrategy::NormalCutLoad:
      {
        // no additional scale
        break;
      }
      default:
        dserror("Unknown merge strategy is equal to %d", merge_strategy);
    }

    // linear search for the node in the current leaf
    std::vector<Point*> merge_candidates;  // canditates for merge

    for (RCPPointSet::iterator i = points_.begin(); i != points_.end(); ++i)
    {
      Point* p = &**i;

      p->Coordinates(nx.A());
      nx.Update(-1, px, 1);
      if (nx.Norm2() <= tol)
      {
        merge_candidates.push_back(p);
      }
    }


    // if there are merge candidates
    if (merge_candidates.size() > 0)
    {
      // first consider topologically connnected points
      // and try to find any candidate there
      std::vector<Point*> topological_candidates;
      Point* mymerge = NULL;
      int match_number = 0;
      // base tolerance, all the candidates should fit into that one
      for (unsigned int i = 0; i < merge_candidates.size(); ++i)
      {
        // if this point was cut but same side and edge.
        // can happen if this edge is shared between multiple sides, we want the intersection point
        // to be consistent and not to merged somewhere else, as then it will cause problem
        if (merge_candidates[i]->IsCut(cut_side, cut_edge))
        {
          topological_candidates.push_back(merge_candidates[i]);
          mymerge = merge_candidates[i];
          match_number++;
        }
      }

      if (match_number >= 2)
      {
        std::stringstream err_msg;
        err_msg << "More then one merge candidate in the pointpool, that both are intersected by "
                   "the same side"
                   "and edge and fit into tolerance. This should not happen. Some points are "
                   "probably created twice or not merged";
        // do gmsh side and edge dump for analysis
        std::ofstream file("pointpool_conflict_topologically_connected.pos");
        for (std::vector<Point*>::iterator it = topological_candidates.begin();
             it != topological_candidates.end(); ++it)
        {
#if CUT_CREATION_INFO
          std::cout << "Id" << (*it)->Id()
                    << ((*it)->IsReallyCut(std::pair<Side*, Edge*>(cut_side, cut_edge))
                               ? " is really cut"
                               : "connection is added later")
                    << std::endl;
#endif
          (*it)->DumpConnectivityInfo();
        }

        if (cut_side)
        {
          GEO::CUT::OUTPUT::GmshNewSection(file, "CurrentIntersectionSide");
          GEO::CUT::OUTPUT::GmshSideDump(file, cut_side, false, NULL);
          GEO::CUT::OUTPUT::GmshEndSection(file, false);
        }
        if (cut_edge)
        {
          GEO::CUT::OUTPUT::GmshNewSection(file, "CurrentIntersectionEdge");
          GEO::CUT::OUTPUT::GmshEdgeDump(file, cut_edge, false, NULL);
          GEO::CUT::OUTPUT::GmshEndSection(file, false);
        }

        file.close();
        dserror(err_msg.str());
      }

      // if there are no points intersected by same side and edge but there  are still other points
      // that fit into the tolerance
      if (mymerge == NULL)
      {
        // there only unknown non-topologically connected case
        if (merge_candidates.empty()) return NULL;

        // if nothing was intersected before return first matching point. else do nothing
        mymerge = merge_candidates[0];

#if EXTENDED_CUT_DEBUG_OUTPUT
        // If there are more matching candidates we basically don't really know where to merge.
        // For now we just merge into the closest one
        // Otherwise One possible idea is is would be to
        // merge into the point that is "the most close". To do this we just need to find minimum
        // merge_candidates[i] ->Coordinates( nx.A() )
        // min (nx  - px).Norm2()  and merge into that one.
        if (merge_candidates.size() > 1)
        {
          std::cout << "NOTE: More then one merge candidate in the pointpool that fit into the "
                       "tolerance. There are possibly"
                    << merge_candidates.size() << "merging candidates. Mergin it into "
                    << mymerge->Id() << std::endl;

          // do gmsh side and edge dump for analysis
          std::ofstream file("pointpool_conflict_side_no_shared_cut.pos", std::ios_base::app);

          if (cut_side)
          {
            GEO::CUT::OUTPUT::GmshNewSection(file, "InterSides");
            GEO::CUT::OUTPUT::GmshSideDump(file, cut_side, false, NULL);
            GEO::CUT::OUTPUT::GmshEndSection(file, false);
          }
          if (cut_edge)
          {
            GEO::CUT::OUTPUT::GmshNewSection(file, "InterEdges");
            GEO::CUT::OUTPUT::GmshEdgeDump(file, cut_edge, false, NULL);
            GEO::CUT::OUTPUT::GmshEndSection(file, false);
          }

          file.close();
          std::cout << "Actual point coordinates, (the real one, before merging) ar "
                    << std::setprecision(15) << px << std::endl;
        }
#endif
      }

      // after this we found proper megin candidates
      if (mymerge != NULL)
      {  // just a safety check

        mymerge->Coordinates(nx.A());
        nx.Update(-1, px, 1);
        if (not mymerge->IsCut(cut_side, cut_edge))
        {
          if (cut_edge != NULL)
          {
            // Here we need to set this point local coordinates on the edge, based on the unmerged
            // point, that it why "t" is called explicitely before AddEdge
            mymerge->t(cut_edge, px);
            mymerge->AddEdge(cut_edge);
          }
          if (cut_side != NULL)
          {
            mymerge->AddSide(cut_side);
          }
          if ((cut_side != NULL) && (cut_edge != NULL))
          {
            mymerge->AddPair(cut_side, cut_edge);
#if CUT_CREATION_INFO
            std::pair<Side*, Edge*> int_pair = std::make_pair(cut_side, cut_edge);
            mymerge->AddMergedPair(int_pair, &px);
#endif
          }
        }

        // This is done because we want to merge cut_mesh into the normal mesh first
        if (merge_strategy == Pointpool_MergeStrategy::InitialLoad)
        {
          mymerge->MovePoint(x);
        }
        return mymerge;
      }
    }
  }
  return NULL;
}


/*-----------------------------------------------------------------------------------------*
 * Get the point with the specified coordinates "x" from the pointpool
 *-----------------------------------------------------------------------------------------*/
Teuchos::RCP<GEO::CUT::Point> GEO::CUT::OctTreeNode::CreatePoint(
    unsigned newid, const double* x, Edge* cut_edge, Side* cut_side, double tolerance)
{
  if (not IsLeaf())
  {
    // call recursively CreatePoint for the child where the Point shall lie in
    Teuchos::RCP<Point> p = Leaf(x)->CreatePoint(newid, x, cut_edge, cut_side, tolerance);
    // add the pointer not only in the leaf but also on the current level
    AddPoint(x, p);
    return p;
  }
  else
  {
    // create a new point and add the point at the lowest level
    Teuchos::RCP<Point> p = GEO::CUT::CreatePoint(newid, x, cut_edge, cut_side, tolerance);
    AddPoint(x, p);
    return p;
  }
}


/*-----------------------------------------------------------------------------------------*
 * Simply insert p into the pointpool and correspondingly modify the boundingbox size
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::OctTreeNode::AddPoint(const double* x, Teuchos::RCP<Point> p)
{
  points_.insert(p);  // insert the point in the pointpool
  bb_->AddPoint(x);   // modify the boundingbox size
}


/*-----------------------------------------------------------------------------------------*
 * get the leaf where the point with the given coordinates lies in
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::OctTreeNode* GEO::CUT::OctTreeNode::Leaf(const double* x)
{
  // navigate to the right one of the 8 children nodes
  //
  //    z <0            1  |  3         z > 0        5  |  7
  //        ____ y     ____|____                    ____|____
  //       |               |                            |
  //       |            2  |  4                      6  |  8
  //       x
  //

  int idx = 0;
  if (x[0] > splitpoint_(0))  // add an index of one to move in x direction
  {
    idx += 1;
  }
  if (x[1] > splitpoint_(1))  // add an index of two to move in y direction
  {
    idx += 2;
  }
  if (x[2] > splitpoint_(2))  // add an index of four to move in z direction
  {
    idx += 4;
  }
  return &*nodes_[idx];
}


/*-----------------------------------------------------------------------------------------*
 * split the current boounding box (tree-node)
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::OctTreeNode::Split(int level)
{
  // We must not end up with a OctTreeNode that holds just nodes from the
  // cutter mesh. However, there is no real way to test this right now.

  if (points_.size() > 125)  /// 125 = 1/8 *1000 -> see NewPoint
  {
    LINALG::Matrix<3, 1> x;
    bool first = true;

    for (RCPPointSet::iterator i = points_.begin(); i != points_.end(); ++i)
    {
      Point* p = &**i;
      if (first)
      {
        first = false;
        p->Coordinates(splitpoint_.A());
      }
      else
      {
        p->Coordinates(x.A());
        splitpoint_.Update(1, x, 1);
      }
    }

    splitpoint_.Scale(1. / points_.size());

    for (int i = 0; i < 8; ++i)
    {
      nodes_[i] = Teuchos::rcp(new OctTreeNode(norm_));
    }

    // avoid empty room (room not covered by boundary boxes)
    for (int i = 0; i < 8; ++i)
    {
      // always have the split point in all boxes
      nodes_[i]->bb_->AddPoint(splitpoint_);

      // always have the outmost point in each box
      double x[3];
      bb_->CornerPoint(i, x);
      Leaf(x)->bb_->AddPoint(x);
    }

    for (RCPPointSet::iterator i = points_.begin(); i != points_.end(); ++i)
    {
      Teuchos::RCP<Point> p = *i;
      double x[3];
      p->Coordinates(x);
      Leaf(x)->AddPoint(x, p);
    }

    for (int i = 0; i < 8; ++i)
    {
      nodes_[i]->Split(level + 1);
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 * collect all edges
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::OctTreeNode::CollectEdges(const BoundingBox& edgebox, plain_edge_set& edges)
{
  if (not IsLeaf())
  {
    if (edgebox.Within(norm_, *bb_))
    {
      for (int i = 0; i < 8; ++i)
      {
        nodes_[i]->CollectEdges(edgebox, edges);
      }
    }
  }
  else
  {
    Teuchos::RCP<BoundingBox> sbox = Teuchos::rcp(BoundingBox::Create());
    for (RCPPointSet::iterator i = points_.begin(); i != points_.end(); ++i)
    {
      Point* p = &**i;
      const plain_edge_set& sds = p->CutEdges();
      for (plain_edge_set::const_iterator i = sds.begin(); i != sds.end(); ++i)
      {
        Edge* s = *i;
        if (edges.count(s) == 0)
        {
          sbox->Assign(*s);
          if (sbox->Within(norm_, edgebox))
          {
            edges.insert(s);
          }
        }
      }
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 * collect all sides
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::OctTreeNode::CollectSides(const BoundingBox& sidebox, plain_side_set& sides)
{
  if (not IsLeaf())
  {
    if (sidebox.Within(norm_, *bb_))
    {
      for (int i = 0; i < 8; ++i)
      {
        nodes_[i]->CollectSides(sidebox, sides);
      }
    }
  }
  else
  {
    Teuchos::RCP<BoundingBox> sbox = Teuchos::rcp(BoundingBox::Create());
    for (RCPPointSet::iterator i = points_.begin(); i != points_.end(); ++i)
    {
      Point* p = &**i;
      const plain_side_set& sds = p->CutSides();
      for (plain_side_set::const_iterator i = sds.begin(); i != sds.end(); ++i)
      {
        Side* s = *i;
        if (sides.count(s) == 0)
        {
          sbox->Assign(*s);
          if (sbox->Within(norm_, sidebox))
          {
            sides.insert(s);
          }
        }
      }
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 * collect all elements near the sidebox
 * (unused, does not work properly when there is no point adjacent to elements in a tree's leaf,
 * e.g. when side lies within an element)
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::OctTreeNode::CollectElements(const BoundingBox& sidebox, plain_element_set& elements)
{
  // see REMARK in cut_mesh.cpp
  dserror(
      "collecting elements via the OctTreeNode does not find all possible element-side "
      "intersections");

  if (not IsLeaf())
  {
    if (sidebox.Within(
            norm_, *bb_))  // within check is a check of overlap between the 2 bounding boxes
    {
      for (int i = 0; i < 8; ++i)
      {
        nodes_[i]->CollectElements(sidebox, elements);
      }
    }
  }
  else
  {
    Teuchos::RCP<BoundingBox> elementbox = Teuchos::rcp(BoundingBox::Create());
    for (RCPPointSet::iterator i = points_.begin(); i != points_.end(); ++i)
    {
      Point* p = &**i;
      const plain_element_set& els = p->Elements();

      // add all elements adjacent to the current point
      // REMARK: this does not find all elements that have an overlap with the sidebox!!!
      for (plain_element_set::const_iterator i = els.begin(); i != els.end(); ++i)
      {
        Element* e = *i;
        if (elements.count(e) == 0)
        {
          elementbox->Assign(*e);
          if (elementbox->Within(norm_, sidebox))
          {
            elements.insert(e);
          }
        }
      }
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 * reset the Point::Position of outside points
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::OctTreeNode::ResetOutsidePoints()
{
  for (RCPPointSet::iterator i = points_.begin(); i != points_.end(); ++i)
  {
    Point* p = &**i;
    if (p->Position() == Point::outside)
    {
      p->Position(Point::undecided);
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 * print the tree at a given level
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::OctTreeNode::Print(int level, std::ostream& stream)
{
  if (not IsLeaf())
  {
    for (int i = 0; i < 8; ++i)
    {
      nodes_[i]->Print(level + 1, stream);
    }
  }
  else
  {
    for (RCPPointSet::iterator i = points_.begin(); i != points_.end(); ++i)
    {
      Point* p = &**i;
      p->Plot(stream);
    }
    stream << "\n";
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::PointPool::PointPool(double norm)
    : tree_(norm), probdim_(DRT::Problem::Instance()->NDim())
{
}


// find node of the octree where current point resides, search start from "this" node and goes
// downwards
GEO::CUT::OctTreeNode* GEO::CUT::OctTreeNode::FindNode(const double* coord, Point* p)
{
  if (p)
  {
    if (not IsLeaf())
    {
      // stop finding the point when point is not included in the current bounding box
      if (!bb_->Within(1.0, coord)) return NULL;
      for (int i = 0; i < 8; ++i)
      {
        GEO::CUT::OctTreeNode* node = (nodes_[i]->FindNode(coord, p));
        if (node != NULL)
        {
          return node;
        }
      }
    }
    // if it is leaf
    else
    {
      RCPPointSet::iterator i = points_.begin();
      for (; i != points_.end(); ++i)
      {
        Point* b = &**i;
        if (b == p)
        {
          // return current node
          return this;
        }
      }
      if (i == points_.end())
      {
        return NULL;
      }
    }
  }
  else
  {
    dserror("Invalid point");
    return NULL;
  }
  return NULL;
}
