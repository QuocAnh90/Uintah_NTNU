/*
 *  Salmon.cc:  The Geometry Viewer
 *
 *  Written by:
 *   Steven G. Parker
 *   Department of Computer Science
 *   University of Utah
 *   April 1994
 *
 *  Copyright (C) 1994 SCI Group
 */

#include <Modules/Salmon/Salmon.h>
#include <Modules/Salmon/Renderer.h>
#include <Modules/Salmon/Roe.h>
#include <Classlib/HashTable.h>
#include <Classlib/NotFinished.h>
#include <Comm/MessageTypes.h>
#include <Dataflow/Connection.h>
#include <Dataflow/HelpUI.h>
#include <Dataflow/ModuleHelper.h>
#include <Dataflow/ModuleList.h>
#include <Datatypes/GeometryComm.h>
#include <Devices/DBContext.h>
#include <Geom/Geom.h>
#include <Geom/HeadLight.h>
#include <iostream.h>

static Module* make_Salmon(const clString& id)
{
    return new Salmon(id);
}

static RegisterModule db1("Geometry", "Salmon", make_Salmon);
static RegisterModule db2("Dave", "Salmon", make_Salmon);

Salmon::Salmon(const clString& id)
: Module("Salmon", id, Sink), max_portno(0)
{
    // Add a headlight
    lighting.lights.add(new HeadLight(Color(1,1,1)));
    // Create the input port
    add_iport(new GeometryIPort(this, "Geometry", GeometryIPort::Atomic));
    default_matl=new Material(Color(.1,.1,.1), Color(.6,0,0),
			      Color(.7,.7,.7), 10);
    busy_bit=1;
    have_own_dispatch=1;
}

Salmon::~Salmon()
{
}

Module* Salmon::clone(int deep)
{
    return new Salmon(*this, deep);
}

void Salmon::do_execute()
{
    while(1){
	if(mailbox.nitems() == 0){
	    // See if anything needs to be redrawn...
	    for(int i=0;i<roe.size();i++)
		roe[i]->redraw_if_needed();
	}
	busy_bit=0;
	MessageBase* msg=mailbox.receive();
	busy_bit=1;
	GeometryComm* gmsg=(GeometryComm*)msg;
	switch(msg->type){
#ifdef OLDUI
	case MessageTypes::DoDBCallback:
	    {
		DBCallback_Message* cmsg=(DBCallback_Message*)msg;
		cmsg->mcb->perform(cmsg->context, cmsg->which,
				   cmsg->value, cmsg->delta, cmsg->cbdata);
	    }
	    break;
#endif
	case MessageTypes::ExecuteModule:
	    // We ignore these messages...
	    break;
	case MessageTypes::RoeRedraw:
	    {
		RedrawMessage* rmsg=(RedrawMessage*)msg;
		for(int i=0;i<roe.size();i++){
		    Roe* r=roe[i];
		    if(r->id == rmsg->rid){
			r->redraw();
			break;
		    }
		}
	    }
	    break;
	case MessageTypes::RoeMouse:
	    {
		RoeMouseMessage* rmsg=(RoeMouseMessage*)msg;
		for(int i=0;i<roe.size();i++){
		    Roe* r=roe[i];
		    if(r->id == rmsg->rid){
			(r->*(rmsg->handler))(rmsg->action, rmsg->x, rmsg->y);
			break;
		    }
		}
	    }
	    break;
	case MessageTypes::GeometryInit:
	    initPort(gmsg->reply);
	    break;	
	case MessageTypes::GeometryAddObj:
	    addObj(gmsg->portno, gmsg->serial, gmsg->obj, gmsg->name);
	    break;
	case MessageTypes::GeometryDelObj:
	    delObj(gmsg->portno, gmsg->serial);
	    break;
	case MessageTypes::GeometryDelAll:
	    delAll(gmsg->portno);
	    break;
	case MessageTypes::GeometryFlush:
	    flushViews();
	    break;
	default:
	    cerr << "Salomon: Illegal Message type: " << msg->type << endl;
	    break;
	}
	delete msg;
    }
}


int Salmon::should_execute()
{
    // See if there is new data upstream...
    int changed=0;
    for(int i=0;i<iports.size();i++){
	IPort* port=iports[i];
	for(int c=0;c<port->nconnections();c++){
	    Module* mod=port->connection(c)->iport->get_module();
	    if(mod->sched_state == SchedNewData){
		sched_state=SchedNewData;
		changed=1;
		break;
	    }
	}
    }
    return changed;
}

void Salmon::initPort(Mailbox<GeomReply>* reply)
{
    reply->send(GeomReply(max_portno++, &busy_bit));
}

void Salmon::flushViews()
{
    for (int i=0; i<topRoe.size(); i++) {
	topRoe[i]->redraw();
    }
}

void Salmon::addObj(int portno, int serial, GeomObj *obj,
		    const clString& name)
{
    HashTable<int, SceneItem*>* serHash;
    if (!portHash.lookup(portno, serHash)) {
	// need to make this table
	serHash = new HashTable<int, SceneItem*>;
	portHash.insert(portno, serHash);
    }
    clString pname(name+" ("+to_string(portno)+")");
    SceneItem* si=new SceneItem(obj, pname);
    serHash->insert(serial, si);
    for (int i=0; i<roe.size(); i++)
	roe[i]->itemAdded(si);
}

void Salmon::delObj(int portno, int serial)
{
    HashTable<int, SceneItem*>* serHash;
    if (portHash.lookup(portno, serHash)) {
	SceneItem *si;
	if(serHash->lookup(serial, si)){
	    serHash->remove(serial);
	    for (int i=0; i<roe.size(); i++)
		roe[i]->itemDeleted(si);
	    delete si->obj;
	    delete si;
	} else {
	    cerr << "Error deleting object, object not in database...(serial=" << serial << ")" << endl;
	}
    } else {
	cerr << "Error deleting object, port not in database...(port=" << portno << ")" << endl;
    }
}

void Salmon::delAll(int portno)
{

    HashTable<int, SceneItem*>* serHash;
    if (portHash.lookup(portno, serHash)) {
	HashTableIter<int, SceneItem*> iter(serHash);
	for (iter.first(); iter.ok();) {
	    SceneItem* si=iter.get_data();
	    int serial=iter.get_key();
	    ++iter; // We have to increment before we nuke the
	            // current element...
	    serHash->remove(serial);
	    for (int i=0; i<roe.size(); i++) {
		roe[i]->itemDeleted(si);
	    }
	    delete si->obj;
	    delete si;
	}
    }
}

void Salmon::addTopRoe(Roe *r)
{
    topRoe.add(r);
}

void Salmon::delTopRoe(Roe *r)
{
    for (int i=0; i<topRoe.size(); i++) {
	if (r==topRoe[i]) topRoe.remove(i);
    }
} 

#ifdef OLDUI
void Salmon::spawnIndCB(CallbackData*, void*)
{
  topRoe.add(new Roe(this));
  topRoe[topRoe.size()-1]->SetTop();
  GeomItem *item;
  for (int i=0; i<topRoe[0]->geomItemA.size(); i++) {
      item=topRoe[0]->geomItemA[i];
      topRoe[topRoe.size()-1]->itemAdded(item->geom, item->name);
  }
//  printFamilyTree();
}
#endif

Salmon::Salmon(const Salmon& copy, int deep)
: Module(copy, deep)
{
    NOT_FINISHED("Salmon::Salmon");
}

void Salmon::connection(ConnectionMode mode, int which_port, int)
{
    if(mode==Disconnected){
	remove_iport(which_port);
    } else {
	add_iport(new GeometryIPort(this, "Geometry", GeometryIPort::Atomic));
    }
}

void Salmon::tcl_command(TCLArgs& args, void* userdata)
{
    if(args.count() < 2){
	args.error("Salmon needs a minor command");
	return;
    }
    if(args[1] == "addroe"){
	if(args.count() != 3){
	    args.error("addroe must have a RID");
	    return;
	}
	Roe* r=new Roe(this, args[2]);
	roe.add(r);
    } else if(args[1] == "listrenderers"){
	Array1<clString> list;
	AVLTreeIter<clString, make_Renderer> iter(Renderer::get_db());
	for(iter.first();iter.ok();++iter){
	    list.add(iter.get_key());
	}
	args.result(args.make_list(list));
    } else {
	Module::tcl_command(args, userdata);
    }
}

void Salmon::execute()
{
    // Never gets called...
}

RedrawMessage::RedrawMessage(const clString& rid)
: MessageBase(MessageTypes::RoeRedraw), rid(rid)
{
}

RedrawMessage::~RedrawMessage()
{
}

SceneItem::SceneItem(GeomObj* obj, const clString& name)
: obj(obj), name(name)
{
}

SceneItem::~SceneItem()
{
}
