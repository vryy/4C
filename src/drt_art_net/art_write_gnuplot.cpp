/*!----------------------------------------------------------------------
\file art_write_gnuplot.cpp
\brief Method to print the arteries in a way that could be displayed by
\gnuplot

\level 3

\maintainer Lena Yoshihara

*----------------------------------------------------------------------*/


#include "art_write_gnuplot.H"
#include "../drt_lib/drt_utils.H"
#include <sstream>

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                    ismail 08/09|
 |                                                                      |
 |                                                                      |
 |       ------> (direction of the flow)                                |
 |       1                 2                 3                 4        |
 |       +-----------------o-----------------o-----------------+        |
 |       ^        ^                 ^                 ^        ^        |
 |    ___|____    |                 |                 |     ___|____    |
 |   [DPOINT 1]   |                 |                 |    [DPOINT 2]   |
 |             ___|___           ___|___           ___|___              |
 |            [DLINE 1]         [DLINE 1]         [DLINE 1]             |
 |                                                                      |
 | ...................................................................  |
 |                                                                      |
 | The gnuplot format exporter will export the results (DOFs) of each   |
 | artery in a different file.                                          |
 | Each artery is defined as a set of elements that belong to a similar |
 | design line (DLINE)                                                  |
 |                                                                      |
 | Therefore, ArtWriteGnuplotWrapper will check how many arteries are   |
 | there to export and generate the associated condition which will     |
 | export it.                                                           |
 |                                                                      |
 | For now we will consider that each artery must have a ascending      |
 | node numbering in the direction of the flow. This could be improved  |
 | later! ;)                                                            |
 |                                                                      |
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

ART::UTILS::ArtWriteGnuplotWrapper::ArtWriteGnuplotWrapper(
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::ParameterList& params)
    : discret_(actdis)
{
  // -------------------------------------------------------------------
  // Get all gnuplot export conditions
  // -------------------------------------------------------------------
  std::vector<DRT::Condition*> myConditions;
  discret_->GetCondition("ArtWriteGnuplotCond", myConditions);
  int numofcond = myConditions.size();

  // -------------------------------------------------------------------
  // if gnuplot export conditions exist then create the classes
  // which will export the files
  // -------------------------------------------------------------------
  if (numofcond > 0 && discret_->Comm().MyPID() == 0)
  {
    // Start by creating a map of classes that will export the wanted arteries
    for (unsigned int i = 0; i < myConditions.size(); i++)
    {
      // ---------------------------------------------------------------
      // Read in the artery number and the nodes assosiated with the
      // condition
      // ---------------------------------------------------------------
      int Artery_Number;
      Artery_Number = myConditions[i]->GetInt("ArteryNumber");
      const std::vector<int>* nodes = myConditions[i]->Nodes();

      // ---------------------------------------------------------------
      // Sort all nodes so such that inlet node is the first and outlet
      // node is the last
      // ---------------------------------------------------------------

      // step (1) find both inlet and outlet nodes
      DRT::Node* ndi = NULL;  // ith node
      DRT::Node* ndl = NULL;  // last node

      for (unsigned int n = 0; n < nodes->size(); n++)
      {
        DRT::Node* nd = actdis->gNode((*nodes)[n]);
        if (nd->GetCondition("ArtInOutCond"))
        {
          std::string TerminalType =
              *(nd->GetCondition("ArtInOutCond")->Get<std::string>("terminaltype"));
          if (TerminalType == "inlet")
            ndi = nd;
          else
            ndl = nd;
        }
      }

      if (ndl == NULL) dserror("artery %d has no outlet node!", Artery_Number);
      if (ndi == NULL) dserror("artery %d has no inlet node!", Artery_Number);


      // loop over all nodes
      std::vector<int>* sorted_nodes = new std::vector<int>;
      DRT::Element** Elements = ndi->Elements();

      DRT::Element* Elem_i;
      if (ndi->NumElement() != 1)
        dserror("artery %d must have one element connected to the inlet node!", Artery_Number);

      Elem_i = &Elements[0][0];

      sorted_nodes->push_back(ndi->Id());

      for (unsigned int n = 0; n < nodes->size() - 2; n++)
      {
        // find the next node!
        if (Elem_i->Nodes()[0]->Id() != ndi->Id())
          ndi = Elem_i->Nodes()[0];
        else
          ndi = Elem_i->Nodes()[1];
        if (ndi->NumElement() != 2)
          dserror(
              "artery %d must have two elements connected to any internal node!", Artery_Number);

        // find the next element
        Elements = ndi->Elements();

        if (Elements[0][0].Id() != Elem_i->Id())
          Elem_i = Elements[0];
        else
          Elem_i = Elements[1];
        sorted_nodes->push_back(ndi->Id());
      }

      sorted_nodes->push_back(ndl->Id());

      // ---------------------------------------------------------------
      // Allocate the gnuplot export condition
      // ---------------------------------------------------------------
      Teuchos::RCP<ArtWriteGnuplot> artgnu_c = Teuchos::rcp(new ArtWriteGnuplot(Artery_Number));


      // ---------------------------------------------------------------
      // Sort the export ondition in a map and check whether the
      // condition exists more than once, which shouldn't be allowed
      // ---------------------------------------------------------------
      bool inserted = agmap_.insert(std::make_pair(Artery_Number, artgnu_c)).second;
      bool inserted2 = agnode_map_.insert(std::make_pair(Artery_Number, sorted_nodes)).second;

      if (!inserted || !inserted2)
        dserror("Each artery must have a unique artery number, please correct your input file\n");

      std::cout << "----------------------------------------------------------" << std::endl;
      std::cout << "Artery[" << Artery_Number << "] has the following sorted nodes" << std::endl;
      for (unsigned int n = 0; n < sorted_nodes->size(); n++)
      {
        std::cout << (*sorted_nodes)[n] << "\t";
      }
      std::cout << std::endl;
      std::cout << "----------------------------------------------------------" << std::endl;
    }
  }
  // throw;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Write (public)                                          ismail 08/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//

void ART::UTILS::ArtWriteGnuplotWrapper::Write(Teuchos::ParameterList& params)
{
  //----------------------------------------------------------------------
  // Exit if the function accessed by a non-master processor
  //----------------------------------------------------------------------
  if (discret_->Comm().MyPID() == 0)
  {
    // -------------------------------------------------------------------
    // loop over all conditions and export the arteries values
    // -------------------------------------------------------------------
    std::map<const int, Teuchos::RCP<class ArtWriteGnuplot>>::iterator mapiter;

    // defining a constant that will have the artery number
    int art_num;
    for (mapiter = agmap_.begin(); mapiter != agmap_.end(); mapiter++)
    {
      art_num = mapiter->first;
      mapiter->second->ArtWriteGnuplot::Write(discret_, params, agnode_map_[art_num]);
    }
  }
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                    ismail 08/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
ART::UTILS::ArtWriteGnuplot::ArtWriteGnuplot(int ArteryNum) : ArteryNum_(ArteryNum)
{
  // -------------------------------------------------------------------
  // Create the file with the following name
  // artery[ArteryNum]_.art
  // -------------------------------------------------------------------
  std::stringstream out;
  std::string str, Numb_str;
  char* cstr;
  out << ArteryNum;
  Numb_str = out.str();
  str.clear();
  str = "xxx";
  str += Numb_str;
  str += "_";
  str += ".art";
  cstr = new char[str.size() + 1];
  strcpy(cstr, str.c_str());
  fout_ = Teuchos::rcp(new std::ofstream(cstr));
  delete[] cstr;

  // Avoid warning on unused variable
  (void)ArteryNum_;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                    ismail 08/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
ART::UTILS::ArtWriteGnuplot::ArtWriteGnuplot() {}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
/*----------------------------------------------------------------------*
 |  Constructor (public)                                    ismail 08/09|
 *----------------------------------------------------------------------*/
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>//
void ART::UTILS::ArtWriteGnuplot::Write(Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::ParameterList& params, const std::vector<int>* nodes)
{
  // defining the Length
  double L = 0.0;
  double dL, time;
  int ElemNum;

  for (unsigned int i = 0; i < nodes->size() - 1; i++)
  {
    // get the elements connected to the node
    if (!discret->HaveGlobalNode((*nodes)[i]))
    {
      int proc = discret->Comm().MyPID();
      dserror("Global Node (%d) doesn't exist on processor (%d)\n", (*nodes)[i], proc);
      exit(1);
    }

    //    DRT::Node * nd = discret->lColNode((*nodes)[i]);
    DRT::Node* nd = discret->gNode((*nodes)[i]);
    DRT::Element** ele = nd->Elements();

    // get element location vector, dirichlet flags and ownerships
    std::vector<int> lm;
    std::vector<int> lmstride;
    Teuchos::RCP<std::vector<int>> lmowner = Teuchos::rcp(new std::vector<int>);
    const int* ele_nodes = ele[0][0].NodeIds();

    if (ele_nodes[0] == (*nodes)[i])
      ElemNum = 0;
    else
      ElemNum = 1;

    ele[ElemNum][0].LocationVector(*discret, lm, *lmowner, lmstride);

    // get node coordinates and number of elements per node
    LINALG::Matrix<3, 2> xyze;
    for (int inode = 0; inode < 2; inode++)
    {
      //      const double* x = discret->lColNode((*nodes)[i+inode])->X();
      const double* x = discret->gNode((*nodes)[i + inode])->X();
      xyze(0, inode) = x[0];
      xyze(1, inode) = x[1];
      xyze(2, inode) = x[2];
    }
    // calculate Length of the element
    dL = sqrt(pow(xyze(0, 0) - xyze(0, 1), 2) + pow(xyze(1, 0) - xyze(1, 1), 2) +
              pow(xyze(2, 0) - xyze(2, 1), 2));

    // get the degrees of freedom
    Teuchos::RCP<const Epetra_Vector> qanp = discret->GetState("qanp");
    std::vector<double> myqanp(lm.size());
    DRT::UTILS::ExtractMyValues(*qanp, myqanp, lm);

    // get the current simulation time
    time = params.get<double>("total time");

    // export the degrees of freedom
    (*fout_) << time << "\t" << L << "\t";
    for (unsigned int j = 0; j < lm.size() / 2; j++)
    {
      (*fout_) << myqanp[j] << "\t";
    }
    (*fout_) << nd->Id() << std::endl;
    // Update L
    L += dL;
    // export the dof of the final node
    if (i == nodes->size() - 2)
    {
      (*fout_) << time << "\t" << L << "\t";
      for (unsigned int j = lm.size() / 2; j < lm.size(); j++)
      {
        (*fout_) << myqanp[j] << "\t";
      }
      (*fout_) << nd->Id() << std::endl;
    }
  }
  (*fout_) << std::endl;
}
