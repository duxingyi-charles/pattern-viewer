#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include "db_test.h"


class ControlPad;


class Window : public QWidget
{
    Q_OBJECT
    
public:
    Window(PatchDBServer<PolyMesh>*, DB_Test<PolyMesh>*);
    
private slots:

    
private:
    ControlPad* controlpad;
    
    PatchDBServer<PolyMesh> * dbserver;
    DB_Test<PolyMesh> * dbtest;
    
    
};

#endif
