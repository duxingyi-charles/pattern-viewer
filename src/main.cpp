#include <QApplication>
#include "window.h"


#include "db_test.h"




int main(int argc, char* argv[]){
    
    
    //db init
    string dbFile = "/Users/dohoney/patches.db";
    
    PatchDBServer<PolyMesh> dbserver;
    
    dbserver.setFilename(dbFile);
    
    dbserver.initialize();
    
    
    //db_test
    DB_Test<PolyMesh> dbtest;

    
    //UI
    QApplication app(argc, argv);
    Window window(&dbserver, &dbtest);
    
    window.show();
    
    app.exec();
    
    
    //db close
    dbserver.close();
    
    return 1;
}