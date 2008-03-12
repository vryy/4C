#ifdef CCADISCRET

#include <iostream>

#include "drt_node.H"
#include "drt_element.H"
#include "drt_elementgroup.H"
#include "drt_discret.H"
#include "drt_utils_local_connectivity_matrices.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::EGROUP::SystemEvaluatorBase::SystemEvaluatorBase(Teuchos::RCP<DRT::Discretization> dis,
                                                      const Teuchos::ParameterList& params,
                                                      Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
                                                      Teuchos::RCP<Epetra_Vector> systemvector)
  : dis_(dis),
    params_(params),
    systemmatrix_(systemmatrix),
    systemvector_(systemvector)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::EGROUP::ElementGroupManager::ElementGroupManager(const Discretization& dis)
  : dis_(dis)
{
  GroupElements();

  std::vector<int> mysize;
  mysize.push_back(aligned_.Size());
  mysize.push_back(others_.Size());

  std::vector<int> globalsize(mysize.size());

  dis_.Comm().SumAll(&mysize[0],&globalsize[0],mysize.size());

  if (dis.Comm().MyPID()==0)
  {
    std::cout << "aligned hex8 elements: " << globalsize[0] << "\n"
              << "other elements       : " << globalsize[1] << "\n\n";
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::EGROUP::ElementGroupManager::Evaluate(Evaluator& evaluator)
{
  evaluator.Evaluate(others_);
  evaluator.Evaluate(aligned_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::EGROUP::ElementGroupManager::GroupElements()
{
  int numele = dis_.NumMyColElements();
  for (int i=0; i<numele; ++i)
  {
    DRT::Element* actele = dis_.lColElement(i);
    DRT::Element::DiscretizationType dt = actele->Shape();

    switch (dt)
    {
    case DRT::Element::hex8:
      FindHex8Align(actele);
      break;
    default:
      others_.Add(actele);
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool DRT::EGROUP::ElementGroupManager::FindHex8Align(DRT::Element* ele)
{
  blitz::Array<double,2> xyze(3,8,blitz::ColumnMajorArray<2>());

  DRT::Node** nodes = ele->Nodes();
  for (int inode=0; inode<8; ++inode)
  {
    const double* x = nodes[inode]->X();
    xyze(0,inode) = x[0];
    xyze(1,inode) = x[1];
    xyze(2,inode) = x[2];
  }

  blitz::Range _ = blitz::Range::all();

  blitz::Array<double,1> edge0(3);
  blitz::Array<double,1> edge1(3);
  blitz::Array<double,1> edge2(3);
  blitz::Array<double,1> edge3(3);

  blitz::Array<double,1> edge4(3);
  blitz::Array<double,1> edge5(3);
  blitz::Array<double,1> edge6(3);
  blitz::Array<double,1> edge7(3);

  blitz::Array<double,1> edge8(3);
  blitz::Array<double,1> edge9(3);
  blitz::Array<double,1> edge10(3);
  blitz::Array<double,1> edge11(3);

  edge0  = xyze(_,0) - xyze(_,1);
  edge1  = xyze(_,1) - xyze(_,2);
  edge2  = xyze(_,2) - xyze(_,3);
  edge3  = xyze(_,0) - xyze(_,3);

  edge4  = xyze(_,4) - xyze(_,5);
  edge5  = xyze(_,5) - xyze(_,6);
  edge6  = xyze(_,6) - xyze(_,7);
  edge7  = xyze(_,4) - xyze(_,7);

  edge8  = xyze(_,0) - xyze(_,4);
  edge9  = xyze(_,1) - xyze(_,5);
  edge10 = xyze(_,2) - xyze(_,6);
  edge11 = xyze(_,3) - xyze(_,7);

  if (Xdir(edge0))
  {
    if (Ydir(edge1))
    {
      if (Xdir(edge2) and Ydir(edge3) and
          Xdir(edge4) and Ydir(edge5) and Xdir(edge6) and Ydir(edge7) and
          Zdir(edge8) and Zdir(edge9) and Zdir(edge10) and Zdir(edge11))
      {
        // XYZ
//         std::cout << "element: " << ele->Id()
//                   << " scale (" << edge0(0) << "," << edge1(1) << "," << edge8(2) << ")\n";

        aligned_.Add(ele);
        return true;
      }
    }
    else if (Zdir(edge1))
    {
      if (Xdir(edge2) and Zdir(edge3) and
          Xdir(edge4) and Zdir(edge5) and Xdir(edge6) and Zdir(edge7) and
          Ydir(edge8) and Ydir(edge9) and Ydir(edge10) and Ydir(edge11))
      {
        // XZY
//         std::cout << "element: " << ele->Id()
//                   << " scale (" << edge0(0) << "," << edge8(1) << "," << edge1(2) << ")\n";

        DRT::Node** nodes = ele->Nodes();
        int* nodeids = const_cast<int*>(ele->NodeIds());

        Rotate(nodeids[0],nodeids[3],nodeids[7],nodeids[4]);
        Rotate(nodeids[1],nodeids[2],nodeids[6],nodeids[5]);

        Rotate(nodes[0],nodes[3],nodes[7],nodes[4]);
        Rotate(nodes[1],nodes[2],nodes[6],nodes[5]);

        aligned_.Add(ele);
        return true;
      }
    }
  }
  else if (Ydir(edge0))
  {
    if (Xdir(edge1))
    {
      if (Ydir(edge2) and Xdir(edge3) and
          Ydir(edge4) and Xdir(edge5) and Ydir(edge6) and Xdir(edge7) and
          Zdir(edge8) and Zdir(edge9) and Zdir(edge10) and Zdir(edge11))
      {
        // YXZ
//         std::cout << "element: " << ele->Id()
//                   << " scale (" << edge1(0) << "," << edge0(1) << "," << edge8(2) << ")\n";

        DRT::Node** nodes = ele->Nodes();
        int* nodeids = const_cast<int*>(ele->NodeIds());

        Rotate(nodeids[0],nodeids[1],nodeids[2],nodeids[3]);
        Rotate(nodeids[4],nodeids[5],nodeids[6],nodeids[7]);

        Rotate(nodes[0],nodes[1],nodes[2],nodes[3]);
        Rotate(nodes[4],nodes[5],nodes[6],nodes[7]);

        aligned_.Add(ele);
        return true;
      }
    }
    else if (Zdir(edge1))
    {
      if (Ydir(edge2) and Zdir(edge3) and
          Ydir(edge4) and Zdir(edge5) and Ydir(edge6) and Zdir(edge7) and
          Xdir(edge8) and Xdir(edge9) and Xdir(edge10) and Xdir(edge11))
      {
        // YZX
//         std::cout << "element: " << ele->Id()
//                   << " scale (" << edge8(0) << "," << edge0(1) << "," << edge1(2) << ")\n";

        DRT::Node** nodes = ele->Nodes();
        int* nodeids = const_cast<int*>(ele->NodeIds());

        Rotate(nodeids[4],nodeids[3],nodeids[1]);
        Rotate(nodeids[5],nodeids[7],nodeids[2]);

        Rotate(nodes[4],nodes[3],nodes[1]);
        Rotate(nodes[5],nodes[7],nodes[2]);

        aligned_.Add(ele);
        return true;
      }
    }
  }
  else if (Zdir(edge0))
  {
    if (Xdir(edge1))
    {
      if (Zdir(edge2) and Xdir(edge3) and
          Zdir(edge4) and Xdir(edge5) and Zdir(edge6) and Xdir(edge7) and
          Ydir(edge8) and Ydir(edge9) and Ydir(edge10) and Ydir(edge11))
      {
        // ZXY
//         std::cout << "element: " << ele->Id()
//                   << " scale (" << edge1(0) << "," << edge8(1) << "," << edge0(2) << ")\n";

        DRT::Node** nodes = ele->Nodes();
        int* nodeids = const_cast<int*>(ele->NodeIds());

        Rotate(nodeids[1],nodeids[3],nodeids[4]);
        Rotate(nodeids[2],nodeids[7],nodeids[5]);

        Rotate(nodes[1],nodes[3],nodes[4]);
        Rotate(nodes[2],nodes[7],nodes[5]);

        aligned_.Add(ele);
        return true;
      }
    }
    else if (Ydir(edge1))
    {
      if (Zdir(edge2) and Ydir(edge3) and
          Zdir(edge4) and Ydir(edge5) and Zdir(edge6) and Ydir(edge7) and
          Xdir(edge8) and Xdir(edge9) and Xdir(edge10) and Xdir(edge11))
      {
        // ZYX
//         std::cout << "element: " << ele->Id()
//                   << " scale (" << edge8(0) << "," << edge1(1) << "," << edge0(2) << ")\n";

        DRT::Node** nodes = ele->Nodes();
        int* nodeids = const_cast<int*>(ele->NodeIds());

        Rotate(nodeids[0],nodeids[1],nodeids[5],nodeids[4]);
        Rotate(nodeids[3],nodeids[2],nodeids[6],nodeids[7]);

        Rotate(nodes[0],nodes[1],nodes[5],nodes[4]);
        Rotate(nodes[3],nodes[2],nodes[6],nodes[7]);

        aligned_.Add(ele);
        return true;
      }
    }
  }

  others_.Add(ele);
  return false;
}


#endif
