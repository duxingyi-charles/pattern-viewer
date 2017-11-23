#include "controlpad.h"
#include "renderarea.h"
#include "expandvar.h"

#include "db_test.h"
#include "my_collapse.h"

#include <QtGui>



ControlPad::ControlPad(PatchDBServer<PolyMesh>* ds, DB_Test<PolyMesh>* dt, QWidget * parent)
: QWidget(parent)
{
    //
    dbserver = ds;
    dbtest   = dt;
    
    isCollapse = false;
    isExpand = false;
    
    //render
    raw_renderarea      = new RenderArea(dt);
    collapse_renderarea = new RenderArea(dt);
    expand_renderarea   = new RenderArea(dt);
    
    //side
    sideSpinBox = new QSpinBox;
    sideSpinBox->setRange(3, 20);
    sideSpinBox->setSpecialValueText(tr("3(triangle)"));
    
    sideLabel = new QLabel(tr("Pattern &Sides:"));
    sideLabel->setBuddy(sideSpinBox);
    
    connect(sideSpinBox, SIGNAL(valueChanged(int)), this, SLOT(sideChanged()));
    
    //index
    indexSpinBox = new QSpinBox;
    indexSpinBox->setRange(0, 100);
    indexSpinBox->setSpecialValueText(tr("0(first pattern)"));
    
    indexLabel = new QLabel(tr("Pattern &Index:"));
    indexLabel->setBuddy(indexSpinBox);
    
    connect(indexSpinBox, SIGNAL(valueChanged(int)), this, SLOT(indexChanged()));
    
    //collapse
    collapseCheckBox = new QCheckBox(tr("&Collapse"));
    collapseCheckBox->setCheckState(Qt::CheckState::Unchecked);
    connect(collapseCheckBox, SIGNAL(stateChanged(int)), this, SLOT(collapseChanged()));
    
    //expand
    expandCheckBox = new QCheckBox(tr("Expand"));
    expandCheckBox->setCheckState(Qt::CheckState::Unchecked);
    connect(expandCheckBox, SIGNAL(stateChanged(int)), this, SLOT(expandChanged()));
    
    expandVars = new ExpandVar;
    connect(expandVars, SIGNAL(valueChanged()), this, SLOT(varVecChanged()));
    
    
    //layout
    QHBoxLayout * topLayout = new QHBoxLayout;
    topLayout->addWidget(raw_renderarea);
    topLayout->addWidget(collapse_renderarea);
    topLayout->addWidget(expand_renderarea);
    
    QVBoxLayout *bottomLeftLayout = new QVBoxLayout;
    QHBoxLayout *sideLayout = new QHBoxLayout;
    sideLayout->addWidget(sideLabel);
    sideLayout->addWidget(sideSpinBox);
    QHBoxLayout *indexLayout = new QHBoxLayout;
    indexLayout->addWidget(indexLabel);
    indexLayout->addWidget(indexSpinBox);
    bottomLeftLayout->addLayout(sideLayout);
    bottomLeftLayout->addLayout(indexLayout);
    
    QVBoxLayout *bottomRightLayout = new QVBoxLayout;
    bottomRightLayout->addWidget(collapseCheckBox);
    bottomRightLayout->addWidget(expandCheckBox);
    bottomRightLayout->addWidget(expandVars);
    
    QHBoxLayout * bottomLayout = new QHBoxLayout;
    bottomLayout->addLayout(bottomLeftLayout);
    bottomLayout->addLayout(bottomRightLayout);
    
    QVBoxLayout * layout = new QVBoxLayout;
    layout->addLayout(topLayout);
    layout->addLayout(bottomLayout);
    
    setLayout(layout);
    
    //
    sideChanged();
    collapseChanged();
    expandChanged();
    
}


void ControlPad::sideChanged()
{
    side = sideSpinBox->value();
    
    dbserver->myload(side);
    indexSpinBox->setRange(0, dbserver->numberOfCachedPatches(side)-1);
    indexSpinBox->setValue(0);
    indexChanged();
}


void ControlPad::indexChanged()
{
    index = indexSpinBox->value();
    patchChanged();

}


void ControlPad::collapseChanged()
{
    Qt::CheckState state = collapseCheckBox->checkState();
    isCollapse = (state == Qt::CheckState::Checked ? true : false);
    if (isCollapse) {
        Patch<PolyMesh> tpatch = dbserver->at(side, index);
        vcg::tri::pl::PolyMesh & mesh = tpatch.getMesh();
        dbtest->initCornerBit(tpatch);
        
        vcg::tri::myPolychordCollapse<PolyMesh>::CollapseAllPolychords(mesh, dbtest->cornerBit, true);
        vcg::tri::Allocator<PolyMesh>::CompactFaceVector(mesh);
        vcg::tri::Allocator<PolyMesh>::CompactVertexVector(mesh);
        assert(mesh.vert.size() == mesh.VN());
        assert(mesh.face.size() == mesh.FN());
        
        //setPatch
        std::deque<vcg::face::Pos<PolyMesh::FaceType>> corners;
        vcg::tri::pl::n_corners_type nCorners;
        std::vector<std::pair<vcg::tri::pl::valence_type,vcg::tri::pl::num_type>> sigVec;
        
        dbtest->patchPrepare(mesh, corners, nCorners);
        tpatch.getSingularityVec(sigVec);
        
        Patch<PolyMesh> collapsed_patch;
        
        collapsed_patch.setPatch(mesh, corners, nCorners, mesh.FN(), mesh.VN(), sigVec);
//        collapsed_patch.reconstructPatch();

        collapse_renderarea->set_patch(collapsed_patch);
        collapse_renderarea->show();
    }
    else {
        collapse_renderarea->hide();
    }
}

void ControlPad::expandChanged()
{
    Qt::CheckState state = expandCheckBox->checkState();
    isExpand = (state == Qt::CheckState::Checked ? true : false);
    if (isExpand) {
        Patch<PolyMesh> patch = dbserver->at(side, index);
        int nVar = patch.getBoundaryVarVec().size();
        
        expandVars->reset(nVar);
        expand_renderarea->set_patch(patch);
        
        expandVars->show();
        expand_renderarea->show();
    }
    else {
        expandVars->hide();
        expand_renderarea->hide();
    }
}


void ControlPad::patchChanged()
{
    Patch<PolyMesh> raw_patch = dbserver->at(side, index);
    raw_renderarea->set_patch(raw_patch);
    
    if (isCollapse) {
        collapseChanged();
    }
    
    if (isExpand) {
        expandChanged();
    }
}


void ControlPad::varVecChanged()
{
    Patch<PolyMesh> patch = dbserver->at(side, index);
    std::vector<vcg::tri::pl::num_type> new_varVec;
    expandVars->getvarVec(new_varVec);
    
    patch.setBoundaryVarVec(new_varVec);
//    patch.reconstructPatch();
    
    expand_renderarea->set_patch(patch);
}

