/*
 * post_drt_q9neconvert.cpp
 *
 *  Created on: Jun 4, 2009
 *      Author: wiesner
 */

#ifdef CCADISCRET
#ifdef D_FLUID2

#include "../post_drt_common/post_drt_common.H"
#include "../drt_f2/fluid2.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_node.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include <Teuchos_RCP.hpp>
#include <blitz/array.h>

// classes

class Converter
{
public:
	Converter(PostField* field, PostProblem& problem);	// constructor
	virtual ~Converter(){};

	virtual void ConvertResults();

	virtual void write_vel_vector_result(string result_name, PostResult* result);
	virtual void write_pre_vector_result(string result_name, PostResult* result);

protected:
	PostField* field_;
	RCP<DRT::Discretization> sdis_;		// source discretization (fluid)
	RCP<DRT::Discretization> tvdis_;	// target discretization for velocity
	RCP<DRT::Discretization> tpdis_;	// target discretization for pressure

	RCP<IO::OutputControl> velcontrol_;
	RCP<IO::OutputControl> precontrol_;
	RCP<IO::DiscretizationWriter> velwriter_;
	RCP<IO::DiscretizationWriter> prewriter_;

	std::map<int, int > nodeids_;		// map (global id -> local id)
};

class FSIConverter : public Converter
{
public:
	FSIConverter(PostField* field, PostField* strufield, PostProblem& problem);

	virtual void ConvertResults();

	virtual void write_vel_vector_result(string result_name, PostField* field, PostResult* result);
	virtual void write_pre_vector_result(string result_name, PostField* field, PostResult* result);

protected:
	RCP<DRT::Discretization> sdisstru_;	// source discretization (structure)
	RCP<IO::OutputControl> strucontrol_;
	RCP<IO::DiscretizationWriter> struwriter_;
};

/////////////////////////////////////////////
// implementation


Converter::Converter(PostField* field, PostProblem& problem)
	: field_(field)
{
	// Membervariablen belegen
	sdis_ = field->discretization();

	tvdis_ = Teuchos::rcp(new DRT::Discretization("none", Teuchos::rcp(new Epetra_SerialComm())));
	tpdis_ = Teuchos::rcp(new DRT::Discretization("none", Teuchos::rcp(new Epetra_SerialComm())));

	// Anzahl der Knoten
	int numnodes = sdis_->NumGlobalNodes();

	for (int i=0; i<numnodes; ++i)
	{
		DRT::Node* actnode = sdis_->lRowNode(i);
		nodeids_[actnode->Id()] = i;		// map: globale NodeID -> lokale Id

		/*coords( i, 0) = actnode->X()[0];
		coords( i, 1) = actnode->X()[1];*/

		/////////////////////////
		// Koordinaten auslesen
		double x[2];
		x[0] = actnode->X()[0];
		x[1] = actnode->X()[1];

		int gid = actnode->Id();		// global Node ID

		/////////////////////////
		// Knoten erzeugen
		{	// target velocity discretization
			RCP<DRT::Node> newnode = rcp(new DRT::Node(gid,x,0));	// globale NodeID wird übernommen!
			tvdis_->AddNode(newnode);
		}
		if(sdis_->NumDof(actnode)==3)
		{
			RCP<DRT::Node> newnode = rcp(new DRT::Node(gid,x,0));
			tpdis_->AddNode(newnode);
		}
	}

	// create elements
	int numelements = sdis_->NumGlobalElements();
	for (int i=0; i<numelements; ++i)
	{
		// aktuelles Element aus Source Discretization bestimmen
		DRT::Element* actele = sdis_->lRowElement(i);
		DRT::ELEMENTS::Fluid2* actfluid = dynamic_cast<DRT::ELEMENTS::Fluid2*>(actele);

		if(actfluid==NULL) dserror("error: no fluid2 element?");
		if(actfluid->Shape()!=DRT::Element::quad9) dserror("no quad9 fluid2 elment?");
		if(actfluid->DisMode()!=DRT::ELEMENTS::Fluid2::dismod_taylorhood) dserror("no Taylor Hood element?");

		// 9 Knoten für aktuelles Element bestimmen
		std::vector<int> sourcenodeids;
		copy(actele->NodeIds(), actele->NodeIds()+9, back_inserter(sourcenodeids));

		// Elemente für Geschwindigkeitsdiskretisierung
		{
			Teuchos::RCP<DRT::ELEMENTS::Fluid2> f2 = rcp(new DRT::ELEMENTS::Fluid2(actfluid->Id(),0));
			f2->SetNodeIds(9,&sourcenodeids[0]);
			tvdis_->AddElement(f2);
		}
		// Elemente für Druckdiskretisierung
		{
			Teuchos::RCP<DRT::ELEMENTS::Fluid2> f2 = rcp(new DRT::ELEMENTS::Fluid2(actfluid->Id(),0));
			f2->SetNodeIds(4,&sourcenodeids[0]);
			tpdis_->AddElement(f2);
		}
	}

	tvdis_->FillComplete();
	tpdis_->FillComplete();

	string vfilename = problem.basename();
	vfilename.append(".fluid.vel");
	string pfilename = problem.basename();
	pfilename.append(".fluid.pre");

	velcontrol_ = rcp(new IO::OutputControl(tvdis_->Comm(), "none", "Polynomial", "generated", vfilename, 2, 0, 1000));
	velwriter_ = rcp(new IO::DiscretizationWriter(tvdis_, velcontrol_));
	velwriter_->WriteMesh(0, 0);

	precontrol_ = rcp(new IO::OutputControl(tpdis_->Comm(), "none", "Polynomial", "generated", pfilename, 2, 0, 1000));
	prewriter_ = rcp(new IO::DiscretizationWriter(tpdis_,precontrol_));
	prewriter_->WriteMesh(0, 0);

}

void Converter::ConvertResults()
{
	PostResult result = PostResult(field_);
	while( result.next_result())
	{
		velwriter_->NewStep(result.step(), result.time());
		prewriter_->NewStep(result.step(), result.time());
		cout << "next step... " << endl;

		if (map_has_map(result.group(), "velnp"))
		{
			cout << "Export velnp..." << endl;
			write_vel_vector_result("velnp", &result);
		}
		if (map_has_map(result.group(), "pressure"))
		{
			cout << "Export pressure..." << endl;
			write_pre_vector_result("pressure", &result);
		}

		MAP_ITERATOR iter;
		init_map_iterator(&iter,result.group());

		while (next_map_node(&iter))
		{
			// We do not support multiple definitions of the same name here. We just
			// use the map node to obtain the key string. Afterward we can use normal
			// map functions to find out if this key names an element vector group.
			MAP_NODE* node = iterator_get_node(&iter);
			char* key = node->key;
			if (map_has_map(result.group(),key))
			{
				string skey = string(key);
				if(skey.find("vel__")==0)
				{
					cout << "write vel vector " << skey.substr(5,skey.length()-5) << endl;
					write_vel_vector_result(skey.c_str(),  &result);
				}
				if(skey.find("pre__")==0)
				{
					cout << "write pre vector " << skey.substr(5,skey.length()-5) << endl;
					write_pre_vector_result(skey.c_str(), &result);
				}
			}
		}
	}
}

void Converter::write_vel_vector_result(string result_name, PostResult* result)
{
	// Daten aus source auslesen
	RCP<Epetra_Vector> srcdata = result->read_result(result_name);
	const Epetra_BlockMap& srcmap = srcdata->Map();

	RCP<Epetra_Vector> dstdata = rcp(new Epetra_Vector(*tvdis_->DofRowMap()));	// Zielvektor
	const Epetra_BlockMap& dstmap = dstdata->Map();

	// alle Knoten in tvdis durchiterieren, jeweils der DOFs auslesen und die entsprechenden Werte aus srcdata eintragen
	int numnodes = tvdis_->NumGlobalNodes();
	for (int i=0; i<numnodes; ++i)
	{
		// Vorgehen: Iteration über lokale IDs
		// 1) Quellknoten bestimmen + zugehörige DOFs (global)
		// 2) Zielknoten bestimmen (passen hier automatisch zusammen!) + zugehörige DOFs (global)
		// 3) Werte in Vektoren übertragen:
		//    hierbei globale DOFs in lokale DOFs umwandeln und an entsprechende Position in Vektor schreiben
		//    Alternative: verwende ReplaceMyValues...

		// Source Knoten
		DRT::Node* srcnode = sdis_->lRowNode(i);
		std::vector<int> srcdof = sdis_->Dof(srcnode);

		// Zielknoten
		DRT::Node* dstnode = tvdis_->lRowNode(i);
		std::vector<int> dstdof = tvdis_->Dof(dstnode); // enthält globale DOF Ids

		(*dstdata)[dstmap.LID(dstdof[0])] = (*srcdata)[srcmap.LID(srcdof[0])];
		(*dstdata)[dstmap.LID(dstdof[1])] = (*srcdata)[srcmap.LID(srcdof[1])];
	}

	velwriter_->WriteVector(result_name,dstdata);
}

void Converter::write_pre_vector_result(string result_name, PostResult* result)
{

	RCP<Epetra_Vector> srcdata = result->read_result(result_name);
	const Epetra_BlockMap& srcmap = srcdata->Map();
	//cout << srcmap << endl;

	RCP<Epetra_Vector> dstdata = rcp(new Epetra_Vector(*tpdis_->DofRowMap()));
	const Epetra_BlockMap& dstmap = dstdata->Map();
	//cout << dstmap << endl;

	// global: Schleife über alle Knoten
	int numnodes = sdis_->NumGlobalNodes();
	for (int i=0; i<numnodes; ++i)
	{
		// Source Knoten
		DRT::Node* srcnode = sdis_->lRowNode(i);

		if(sdis_->NumDof(srcnode)<3) continue;	// Druckknoten erkennen

		int lid = nodeids_[srcnode->Id()];

		std::vector<int> srcdof = sdis_->Dof(srcnode);

		// Zielknoten
		DRT::Node* dstnode = tpdis_->gNode(i);	// globale Node ID bei srcnodes und dstnodes gleich!
		std::vector<int> dstdof = tpdis_->Dof(dstnode);

		// debug
		cout << "lid: " << lid << " i: " << i << " SRCID: " << srcnode->Id() << " DOFS: " << srcdof[0] << "," << srcdof[1] << "," << srcdof[2] << " DSTID: " << dstnode->Id() << " DSTDOFS: " << dstdof[0] << "," << dstdof[1] << "," << dstdof[2] << endl;

		(*dstdata)[dstmap.LID(dstdof[0])] = (*srcdata)[srcmap.LID(srcdof[2])];
	}
	prewriter_->WriteVector(result_name, dstdata);

}

/////////////////////////////////////////////
// FSI Converter
FSIConverter::FSIConverter(PostField* field, PostField* strufield, PostProblem& problem)
 : Converter(field, problem)
{
	cout << "This is FSI-Converter" << endl;
}

void FSIConverter::ConvertResults()
{
	cout << "ConvertResults" << endl;
	PostResult result = PostResult(field_);
	while( result.next_result())
	{
		velwriter_->NewStep(result.step(), result.time());
		prewriter_->NewStep(result.step(), result.time());
		cout << "next step... " << endl;

		if (map_has_map(result.group(), "velnp"))
		{
			cout << "Export velnp..." << endl;
			write_vel_vector_result("velnp", field_, &result);
		}
		if (map_has_map(result.group(), "dispnp"))
		{
			cout << "Export dispnp..." << endl;
			write_vel_vector_result("dispnp", field_, &result);
		}
		if (map_has_map(result.group(), "pressure"))
		{
			cout << "Export pressure..." << endl;
			write_pre_vector_result("pressure", field_, &result);
		}

		MAP_ITERATOR iter;
		init_map_iterator(&iter,result.group());

		while (next_map_node(&iter))
		{
			// We do not support multiple definitions of the same name here. We just
			// use the map node to obtain the key string. Afterward we can use normal
			// map functions to find out if this key names an element vector group.
			MAP_NODE* node = iterator_get_node(&iter);
			char* key = node->key;
			if (map_has_map(result.group(),key))
			{
				string skey = string(key);
				if(skey.find("vel__")==0)
				{
					cout << "write vel vector " << skey.substr(5,skey.length()-5) << endl;
					write_vel_vector_result(skey.c_str(), field_, &result);
				}
				if(skey.find("pre__")==0)
				{
					cout << "write pre vector " << skey.substr(5,skey.length()-5) << endl;
					write_pre_vector_result(skey.c_str(), field_, &result);
				}
			}
		}

	}
}

// Erzeuge Vektor in Geschwindigkeitsdiskretisierung FSI
void FSIConverter::write_vel_vector_result(string result_name, PostField* field, PostResult* result)
{
	RCP<Epetra_Vector> srcdata = result->read_result(result_name);
	//const Epetra_BlockMap& srcmap = srcdata->Map();
	const Epetra_BlockMap& srcmap = *(sdis_->DofRowMap());

	RCP<Epetra_Vector> dstdata = rcp(new Epetra_Vector(*tvdis_->DofRowMap()));
	const Epetra_BlockMap& dstmap = dstdata->Map();

	// global: Schleife über alle Knoten
	int numnodes = sdis_->NumGlobalNodes();
	for (int i=0; i<numnodes; ++i)
	{
		// Source Knoten
		DRT::Node* srcnode = sdis_->lRowNode(i);
		int lid = nodeids_[srcnode->Id()];
		std::vector<int> srcdof = sdis_->Dof(srcnode);
		//cout << *(dis_->DofRowMap()) << endl;
		/*const double* pos1 = srcnode->X();
		cout << "Pos: " << pos1[0] << ", " << pos1[1] << ", " << pos1[2] << endl;*/

		// Zielknoten
		DRT::Node* dstnode = tvdis_->lRowNode(lid);
		std::vector<int> dstdof = tvdis_->Dof(dstnode);
		//cout << *(veldis_->DofRowMap()) << endl;
		/*const double* pos2 = srcnode->X();
				cout << "Pos: " << pos2[0] << ", " << pos2[1] << ", " << pos2[2] << endl;*/

		/*for(int j=0; j<dstdof.size(); j++)
			cout << "SRCDOF: " << srcdof[j] << "(" << srcmap.LID(srcdof[j]) << ") -> DSTDof: " << dstdof[j] << "(" << dstmap.LID(dstdof[j]) << ")  |   ";
		cout << endl;*/

		// nur Geschwindigkeitsdaten/Displacements... kopieren
		// TODO das passt noch gar nicht!
		(*dstdata)[dstmap.LID(dstdof[0])] = (*srcdata)[srcmap.LID(srcdof[0])];
		(*dstdata)[dstmap.LID(dstdof[1])] = (*srcdata)[srcmap.LID(srcdof[1])];
	}
	//cout << "DSTDATA: " << *dstdata << endl;
	velwriter_->WriteVector(result_name, dstdata);

}

// spezielle Funktion: extrahiert nur skalare Druckdaten für FSI Fall
void FSIConverter::write_pre_vector_result(string result_name, PostField* field, PostResult* result)
{
	RCP<Epetra_Vector> srcdata = result->read_result(result_name);
	const Epetra_BlockMap& srcmap = srcdata->Map();
	//const Epetra_BlockMap& srcmap = *(dis_->DofRowMap());
	//cout << srcmap << endl;

	RCP<Epetra_Vector> dstdata = rcp(new Epetra_Vector(*tpdis_->DofRowMap()));
	const Epetra_BlockMap& dstmap = dstdata->Map();
	//cout << dstmap << endl;

	// global: Schleife über alle Knoten
	int numnodes = sdis_->NumGlobalNodes();
	int npnode = 0;	// lokaler Zähler für Druckknoten
	for (int i=0; i<numnodes; ++i)
	{
		// Source Knoten
		DRT::Node* srcnode = sdis_->lRowNode(i);
		const double* pos = srcnode->X();

		if(sdis_->NumDof(srcnode)<3) continue;

		std::vector<int> srcdof = sdis_->Dof(srcnode);

		// Zielknoten
		DRT::Node* dstnode = tpdis_->gNode(srcnode->Id());
		const double* pos2 = dstnode->X();
		std::vector<int> dstdof = tpdis_->Dof(dstnode);

		//if(pos[0]!=pos2[0] || pos[1]!=pos2[1] || pos[2]!=pos2[2])
		//	dserror("Knoten haben nicht die gleichen Koordinaten?");

		// globaler SourceDof: START AUSKOMMENTIERT
		/*int globalsrcdof = srcdof[2];
		int globaldstdof = dstdof[0];


		double srcvalue = (*srcdata)[srcmap.LID(globalsrcdof)];
		dstdata->ReplaceGlobalValue(globalsrcdof,0,srcvalue);*/ //END AUSKOMMENTIERT

		/*for(int j=0; j<dstdof.size(); j++)
		{
			cout << "SRCDOF: " << srcdof[j] << " (" << srcmap.LID(srcdof[j])  << ") -> DSTDOF: " << dstdof[j] << " ("<< dstmap.LID(dstdof[0]) << ")   |   ";
		}
		cout << endl;*/

		//(*dstdata)[dstmap.LID(dstdof[0])] = (*srcdata)[srcmap.LID(srcdof[2])];	// geändert
		//cout << srcmap.GID(npnode) << " -> " << dstmap.LID(dstdof[0]) << endl;
		//(*dstdata)[dstmap.LID(dstdof[0])] = (*srcdata)[srcmap.GID(npnode)]; // Druck als X Komponente sonst +2 in dstdata
		(*dstdata)[dstmap.LID(dstdof[0])] = (*srcdata)[npnode]; // die Zeile scheint zu funktionieren! dafür geht fsidebugwriter nicht mehr :-(

		npnode++;
	}
	//cout << "DSTDATA: " << *dstdata << endl;
	prewriter_->WriteVector(result_name, dstdata);
}

/////////////////////////////////////////////
// Hauptroutine
int main(int argc, char** argv)
{
	Teuchos::CommandLineProcessor My_CLP;
	My_CLP.setDocString("Post DRT quad9 (non equal) conversion program\n");

	PostProblem problem(My_CLP, argc, argv);

	switch (problem.Problemtype())
	{
	case prb_fluid:
	{
		// einfaches Fluidproblem
		PostField* field = problem.get_discretization(0);
		Converter conv(field,problem);
		conv.ConvertResults();
		break;
	}
	case prb_fsi:
	{
		// FSI Problem

		PostField* fluidfield = NULL;
		PostField* structurefield = NULL;

		cout << "FSI: number of discretizations: " << problem.num_discr() << endl;
		for(int i=0; i< problem.num_discr(); i++)
		{
			PostField* field = problem.get_discretization(i);
			cout << "Feld: " << i << " " << field->name() << " Nodes: " << field->num_nodes() << " Elements: " << field->num_elements() << endl;
			if(field->name()=="fluid")
				fluidfield = field;
			if(field->name()=="structure")
				structurefield = field;
		}

		if(fluidfield == NULL) dserror ("no fluid field?");

		FSIConverter conv(fluidfield,structurefield,problem);
		conv.ConvertResults();
		break;
	}
	default:
		dserror("problem type %d not yet supported",problem.Problemtype());
	}

	return 0;
}


#endif
#endif
