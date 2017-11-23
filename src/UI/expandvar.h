#ifndef EXPAND_VAR_H
#define EXPAND_VAR_H

//#include "common_header.h"
#include "db_test.h"

#include <QWidget>

#define MAX_N_VAR 21   //max boudaryVarVec.size()

class QLabel;
class QSpinBox;


class ExpandVar : public QWidget
{
    Q_OBJECT
    
public:
    ExpandVar(QWidget * parent = 0);
    
    int getsize();
    void getvarVec(std::vector<vcg::tri::pl::num_type> &);
    
public slots:
    void reset(int);
    void set_varVec();
    
signals:
    void valueChanged();
    
    
private:
    int size;
    std::vector<vcg::tri::pl::num_type> varVec;
    std::vector<QLabel*>    varLabelVec;
    std::vector<QSpinBox*>  varSpinBoxVec;
    
    
};


#endif
