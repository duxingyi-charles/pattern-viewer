#include <QApplication>
#include <iostream>
#include "window.h"
#include "db_test.h"




int main(int argc, char* argv[]){
    
    
    //db init
    if (argc != 2) {
        std::cout << "Usage: ./PatternViewer [database_file]" << std::endl;
        return 0;
    }
    string dbFile = argv[1];
    
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