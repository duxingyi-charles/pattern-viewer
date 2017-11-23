#include "expandvar.h"

#include <QtGui>


ExpandVar::ExpandVar(QWidget *parent) :
    QWidget(parent)
{
    size = 0;
    
    //
    QVBoxLayout *layout = new QVBoxLayout;
    
    for (int i=0; i<MAX_N_VAR; ++i) {
        QHBoxLayout *hlayout = new QHBoxLayout;
        QLabel * label = new QLabel(tr(std::to_string(i).data()));
        label->hide();
        varLabelVec.push_back(label);
        hlayout->addWidget(label);
        
        QSpinBox * sbox = new QSpinBox;
        sbox->setRange(1, 10);
        sbox->setValue(1);
        connect(sbox, SIGNAL(valueChanged(int)), this, SLOT(set_varVec()));
        sbox->hide();
        varSpinBoxVec.push_back(sbox);
        varVec.push_back(1);
        hlayout->addWidget(sbox);
        
        layout->addLayout(hlayout);
    }
    
    setLayout(layout);
    
    
}




int ExpandVar::getsize()
{
    return size;
}

void ExpandVar::getvarVec(std::vector<vcg::tri::pl::num_type> &vec)
{
    vec.clear();
    for (int i=0; i<size; ++i) {
        vec.push_back(varVec[i]);
    }
}



//SLOT

void ExpandVar::reset(int s)
{
    //check input
    if (s > MAX_N_VAR) {
        return;
    }
    
    //
    size = s;
    
    for (int i=0; i<size; ++i) {
        varLabelVec[i]->show();
        varSpinBoxVec[i]->show();
        varSpinBoxVec[i]->setValue(1);
        varVec[i] = 1;
    }
    
    for (int i=size; i<MAX_N_VAR; ++i) {
        varLabelVec[i]->hide();
        varSpinBoxVec[i]->hide();
    }
}



void ExpandVar::set_varVec()
{
    for (int i=0; i<size; ++i) {
        varVec[i] = varSpinBoxVec[i]->value();
    }
    emit valueChanged();
}