#include <QtGui>
#include <cmath>
#include "renderarea.h"

#define PI 3.1415926

RenderArea::RenderArea(DB_Test<PolyMesh> *dt, QWidget *parent)
    : QWidget(parent)
{
    dbtest = dt;
    
    setBackgroundRole(QPalette::Base);  //?
    setAutoFillBackground(true);        //?
}



QSize RenderArea::minimumSizeHint() const
{
    return QSize(100, 100);
}


QSize RenderArea::sizeHint() const
{
    return QSize(400, 200);
}


//SLOT

void RenderArea::set_patch(Patch<PolyMesh> p)
{
    patch = p;
    update();
}



void RenderArea::paintEvent(QPaintEvent * /* event */) //??
{
    ///test
//    std::cout << std::endl << this << " paintEvent ..." << std::endl;
    ///
    patch.reconstructPatch();

    QPainter painter(this);
    painter.save();
    painter.translate(width()/2, height()/2);
    
    ///////
    unsigned int radius = 90;
    dbtest->geometrize(patch, radius);
    
    PolyMesh & mesh = patch.getMesh();
    dbtest->initCornerBit(patch);
    
    //mesh to vector<QPoint>
    std::vector<QPointF> vertices(mesh.VN());
    for (int i=0; i<mesh.vert.size(); ++i) {
        auto &p = mesh.vert[i].P();
        vertices[i].setX(-p[0]);
        vertices[i].setY(-p[1]);
    }
    
    //calc VVAdj
    std::vector<std::vector<int>> Adj;
    dbtest->calc_VVAdj(mesh, Adj);
    
    //draw edges
    for (int i=0; i<vertices.size(); ++i) {
        for (int j=i; j<vertices.size(); ++j) {
            if (Adj[i][j] == 1) {
                painter.drawLine(vertices[i], vertices[j]);
                ///test
//                std::cout << "v[" << i << "]" << "(" << vertices[i].x() << ", " << vertices[i].y() << ")"
//                << " v[" << j << "] " <<  "(" << vertices[j].x() << ", " << vertices[j].y() << ")"
//                << std::endl;
                ///
            }
        }
    }
    
    //highlight corners and vertices
    std::vector<int> vVec;
    dbtest->calc_valenceVec(mesh, vVec);
    QColor v_color;
    for (int i=0; i<vertices.size(); ++i) {
        //change color according to valence
        if (vVec[i] >=2 && vVec[i] <= 6) {
            v_color.setHsv(420-60*vVec[i], 255, 255);
        }
        else {
            v_color.setHsv(0, 255, 255);  //red
        }
        if (vVec[i] == 4) {
            v_color.setHsv(0, 0, 0);   //black
        }
        painter.setPen(QPen(v_color, 5));
        painter.drawPoint(vertices[i]);
        
        if (mesh.vert[i].IsUserBit(dbtest->cornerBit)) { //corners
            painter.setPen(Qt::red);
            painter.drawEllipse(vertices[i], 3, 3);
        }
    }
    //
    
    painter.restore();
    
    painter.drawRect(QRect(0, 0, width() - 1, height() - 1)); //?
    
    
    ///test
//    std::cout << this << " paintEvent end" << std::endl << std::endl;

}

