/*
 *  DenseMatrixSQLReader.cc:
 *
 *  Written by:
 *   Jason V. Morgan
 *   December 21, 2000
 *
 */

#include <Dataflow/Network/Module.h>
#include <Core/Malloc/Allocator.h>
#include <Dataflow/Ports/MatrixPort.h>
#include <Core/Datatypes/DenseMatrix.h>
#include <Core/GuiInterface/GuiVar.h>

#include <Packages/Morgan/Core/Dbi/Dbi.h>
#include <Packages/Morgan/Core/Dbi/Dbd.h>
#include <Packages/Morgan/share/share.h>
#include <stdlib.h>
#include <memory>

using Morgan::Dbi::connect;
using namespace Morgan::Dbi;
using namespace SCIRun; 
using namespace std;

namespace Morgan {
namespace Modules {

class MorganSHARE DenseMatrixSQLReader : public Module {
public:
  DenseMatrixSQLReader(const string& id);

  virtual ~DenseMatrixSQLReader();

  virtual void execute();

  virtual void tcl_command(TCLArgs&, void*);

private:

  // data from the user interface
  GuiString database_ui;
  GuiString hostname_ui;
  GuiString port_ui;
  GuiString username_ui;
  GuiString password_ui;
  GuiString sql_ui;

  // the actual output port
  MatrixOPort* omatrix;
};

extern "C" MorganSHARE Module* make_DenseMatrixSQLReader(const string& id) {
  return new DenseMatrixSQLReader(id);
}

#define init_var(name) name##_ui(#name, id, this)
DenseMatrixSQLReader::DenseMatrixSQLReader(const string& id)
  : Module("DenseMatrixSQLReader", id, Source, "Readers", "Morgan"),
    init_var(database),
    init_var(hostname),
    init_var(port),
    init_var(username),
    init_var(password),
    init_var(sql)
{
}

DenseMatrixSQLReader::~DenseMatrixSQLReader(){
}

void DenseMatrixSQLReader::execute() {
  omatrix = (MatrixOPort *)get_oport("Field");

  string dbase=database_ui.get();
  string hname = hostname_ui.get();
  string pnum = port_ui.get();
  string uname = username_ui.get();
  string pword = password_ui.get();

 cerr << "dbase is: " << dbase << endl;
 cerr << "hname is: " << hname << endl;
 cerr << "pnum is: " << pnum << endl;
 cerr << "uname is: " << uname << endl;
 cerr << "pword is: " << pword << endl;

  auto_ptr<Dbd> dbd(connect(dbase.c_str(),
                              hname.c_str(),
                              atoi(pnum.c_str()), 
                              uname.c_str(), 
                              pword.c_str()));

        
/*     system("/home/sci/butson/test.pl"); */
    if(!dbd.get()) {
        fprintf(stderr, "Unable to connect to database\n");
        return; // could not connect
    }

    int rows = 0;


    // count the number of rows
    if(!dbd->execute(sql_ui.get().c_str())) {
        fprintf(stderr, "Can't execute SQL query %s\n", sql_ui.get().c_str());
    }

    while(!dbd->at_end()) {
        ++rows;
        dbd->next_row();
    }

    if(!dbd->execute(sql_ui.get().c_str())) {
        fprintf(stderr, "Can't execute SQL query %s\n", sql_ui.get().c_str());
    }

    int cols = dbd->cols();
    DenseMatrix* dm = scinew DenseMatrix(rows, cols);

    for(int row = 0 ; row < rows ; ++row) {
        for(int col = 0 ; col < cols ; ++col) {
            string strval;
            dbd->fetch(strval, col);
            dm->put(row, col, atof(strval.c_str()));
        }
        dbd->next_row();
    }

/*
    for(int row = 0 ; row < dm->nrows() ; ++row) {
        for(int col = 0 ; col < dm->ncols() ; ++col) {
            cerr << dm->get(row, col) << "\t";
        }
        cerr << endl;
    }
*/

    MatrixHandle omat = dm;
    omatrix->send(omat);
}

void DenseMatrixSQLReader::tcl_command(TCLArgs& args, void* userdata) {
  Module::tcl_command(args, userdata);
}

} // End namespace Modules
} // End namespace Morgan


