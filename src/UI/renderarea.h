#ifndef RENDERAREA_H
#define RENDERAREA_H

#include <QWidget>

#include "db_test.h"
#include "my_collapse.h"



class RenderArea : public QWidget
{
    Q_OBJECT
    
public:
    RenderArea(DB_Test<PolyMesh>*, QWidget *parent = 0);
    
    QSize minimumSizeHint() const;   //?
    QSize sizeHint() const;          //?

public slots:
    void set_patch(Patch<PolyMesh>);
    
protected:
    void paintEvent(QPaintEvent *event);   //?
    
private:
    DB_Test<PolyMesh> * dbtest;
    
    Patch<PolyMesh> patch;

};



#endif