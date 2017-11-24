#include <QtGui>

//#include "renderarea.h"
#include "controlpad.h"
#include "window.h"

Window::Window(PatchDBServer<PolyMesh>* ds, DB_Test<PolyMesh>* dt)
{
    //
    dbserver = ds;
    dbtest = dt;
    
    controlpad = new ControlPad(ds, dt);
    
    QHBoxLayout * layout = new QHBoxLayout;
    layout->addWidget(controlpad);
    
    setLayout(layout);
    
    setWindowTitle(tr("Pattern Viewer"));

}

