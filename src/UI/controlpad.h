#ifndef CONTROLPAD_H
#define CONTROLPAD_H

#include <QWidget>

#include "db_test.h"
//#include "my_collapse.h"

//#include "common_header.h"

QT_BEGIN_NAMESPACE
class QLabel;
class QSpinBox;
class QCheckBox;
QT_END_NAMESPACE
class RenderArea;
class ExpandVar;




class ControlPad : public QWidget
{
    Q_OBJECT
    
public:
    ControlPad(PatchDBServer<PolyMesh>*, DB_Test<PolyMesh>*, QWidget *parent = 0);
    
public slots:
    
    void sideChanged();
    void indexChanged();
    void collapseChanged();
    void expandChanged();
    
    void patchChanged();
    void varVecChanged();
    
    
    
private:
    int side;
    int index;
    bool isCollapse;
    bool isExpand;
    
    PatchDBServer<PolyMesh> * dbserver;
    DB_Test<PolyMesh> * dbtest;
    
    RenderArea * raw_renderarea;
    RenderArea * collapse_renderarea;
    RenderArea * expand_renderarea;
    
    QLabel   * sideLabel;
    QSpinBox * sideSpinBox;
    
    QLabel   * indexLabel;
    QSpinBox * indexSpinBox;

    QCheckBox * collapseCheckBox;
    QCheckBox * expandCheckBox;
        
    ExpandVar * expandVars;
    
    
};





#endif