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


//Window::Window(PatchDBServer<PolyMesh>* ds, DB_Test<PolyMesh>* dt)
//{
//    //
//    dbserver = ds;
//    dbtest = dt;
//    
//    //renderArea
//    l_renderArea = new RenderArea(ds, dt);
//    r_renderArea = new RenderArea(ds, dt, true); //enable collapse
//    
//    //collapse checkbox
//    collapseCheckBox = new QCheckBox(tr("&Collaspe"));
//    collapseCheckBox->setCheckState(Qt::CheckState::Unchecked);
//    connect(collapseCheckBox, SIGNAL(stateChanged(int)), this, SLOT(collapseChanged()));
//    
//    //nSide
//    nSideSpinBox = new QSpinBox;
//    nSideSpinBox->setRange(3,20);
//    nSideSpinBox->setSpecialValueText(tr("3 (triangle)"));
//    
//    nSideLabel = new QLabel(tr("Polygon &Sides:"));
//    nSideLabel->setBuddy(nSideSpinBox);
//    
//    connect(nSideSpinBox, SIGNAL(valueChanged(int)), this, SLOT(nSideChanged()));
//    
//    //index
//    indexSpinBox = new QSpinBox;
//    indexSpinBox->setRange(0, 100);
//    indexSpinBox->setSpecialValueText(tr("0 (first pattern)"));
//    
//    indexLabel = new QLabel(tr("Pattern &Index:"));
//    indexLabel->setBuddy(indexSpinBox);
//    
//    connect(indexSpinBox, SIGNAL(valueChanged(int)), this, SLOT(indexChanged()));
//    
//    //
//    nSideChanged();
//    indexChanged();
//    collapseChanged();
//    
//    
//    //layout
//    QHBoxLayout *renderLayout = new QHBoxLayout;
//    renderLayout->addWidget(l_renderArea);
//    renderLayout->addWidget(r_renderArea);
//    
//    QVBoxLayout *layout = new QVBoxLayout;
//    layout->addLayout(renderLayout);
//    layout->addWidget(collapseCheckBox);
//    layout->addWidget(nSideLabel);
//    layout->addWidget(nSideSpinBox);
//    layout->addWidget(indexLabel);
//    layout->addWidget(indexSpinBox);
//    
//    setLayout(layout);
//    
//    setWindowTitle(tr("Pattern Viewer"));
//    
//}

//void Window::nSideChanged()
//{
//    int nside = nSideSpinBox->value();
//    dbserver->myload(nside);
//    indexSpinBox->setRange(0, dbserver->numberOfCachedPatches(nside)-1);
//    indexSpinBox->setValue(0);
//    l_renderArea->set_nSide(nside);
//    l_renderArea->set_index(0);
//    Qt::CheckState state = collapseCheckBox->checkState();
//    if (state == Qt::CheckState::Checked) {
//        r_renderArea->set_nSide(nside);
//        r_renderArea->set_index(0);
//    }
//}
//
//
//void Window::indexChanged()
//{
//    int idx = indexSpinBox->value();
//    l_renderArea->set_index(idx);
//    Qt::CheckState state = collapseCheckBox->checkState();
//    if (state == Qt::CheckState::Checked) {
//        r_renderArea->set_index(idx);
//    }
//}
//
//
//void Window::collapseChanged()
//{
//    Qt::CheckState state = collapseCheckBox->checkState();
//    if (state == Qt::CheckState::Checked) {
//        int nside = nSideSpinBox->value();
//        int idx = indexSpinBox->value();
//        r_renderArea->set_nSide(nside);
//        r_renderArea->set_index(idx);
//        r_renderArea->show();
//    }
//    else {
//        r_renderArea->hide();
//    }
//}


