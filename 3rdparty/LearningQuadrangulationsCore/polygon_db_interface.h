#ifndef POLYGON_DB_INTERFACE
#define POLYGON_DB_INTERFACE

#include <iostream>
#include <string>
#include <sstream>
#include <list>
#include <utility>
#include <sqlite3.h>
#include "mesh_support.h"

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief The PolygonDBInterface class defines the base operations to access the database.
 *
 * Usage:
 * \code
 * std::string dbFilename;
 * ...
 * PolygonDBInterface polygonDB;
 * polygonDB.setFilename(dbFilename);
 * polygonDB.setReadWrite(true);
 * polygonDB.initialize();
 * ... // getters and/or setters here
 * polygonDB.close();
 * \endcode
 */
class PolygonDBInterface {
public:
  /**
   * @brief POLYGONDBFILENAME The name of the DB file.
   * @return
   */
  static const std::string & POLYGONDBFILENAME() {
    static const std::string polygonDBFilename("polygons.db");
    return polygonDBFilename;
  }

  /**
   * @brief POLYGON_TABLE_NAME Name of the table storing the polygons for the collapsed patches.
   * @return
   */
  static const std::string & POLYGON_TABLE_NAME() {
    static const std::string polygonTableName("Polygons");
    return polygonTableName;
  }

  /**
   * @brief INDEX_POSTFIX Postfixed name of the index used to improve queries performance.
   * @return
   */
  static const std::string & INDEX_POSTFIX() {
    static const std::string indexPostifix("_index");
    return indexPostifix;
  }

  /**
   * @brief COLUMN_TOPOLOGY Header of topology column.
   * @return
   */
  static const std::string & COLUMN_TOPOLOGY() {
    static const std::string columnTopology("topology");
    return columnTopology;
  }

  /**
   * @brief COLUMN_MESH_FILENAME Header of mesh column.
   * @return
   */
  static const std::string & COLUMN_MESH_FILENAME() {
    static const std::string columnMeshFilename("mesh");
    return columnMeshFilename;
  }

  /**
   * @brief COLUMN_POLYGON Header of polygon (POSes) column.
   * @return
   */
  static const std::string & COLUMN_POLYGON() {
    static const std::string columnPolygon("polygon");
    return columnPolygon;
  }

  /**
   * @brief COLUMN_N_FACES Header of number of faces column.
   * @return
   */
  static const std::string & COLUMN_N_FACES() {
    static const std::string columnNFaces("nFaces");
    return columnNFaces;
  }

  /**
   * @brief PolygonDBInterface Default constructor.
   */
  PolygonDBInterface() : _db(0), _queryType(EMPTY), _stmt(0), _readwrite(false) { }

  /**
   * @brief PolygonDBInterface Constructor setting the DB filename.
   * @param dbFilename
   */
  PolygonDBInterface(const std::string &dbFilename) : _dbFilename(dbFilename), _db(0), _queryType(EMPTY),
                                                      _stmt(0), _readwrite(false) { }

  /**
   * @brief PolygonDBInterface Safe copy constructor.
   * @param polygonDB
   */
  PolygonDBInterface(const PolygonDBInterface &polygonDB) : _dbFilename(polygonDB.getFilename()), _db(0),
                                                            _queryType(EMPTY), _stmt(0),
                                                            _readwrite(polygonDB._readwrite) { }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief PolygonDBInterface Move constructor.
   * @param polygonDB The PolygonDBInterface will be moved and reset.
   */
  PolygonDBInterface(PolygonDBInterface&& polygonDB) : _dbFilename(std::move(polygonDB._dbFilename)), _db(polygonDB._db),
                                                       _queryType(EMPTY), _stmt(polygonDB._stmt),
                                                       _readwrite(polygonDB._readwrite) {
    polygonDB._dbFilename.clear();
    polygonDB._db = 0;
    polygonDB._queryType = EMPTY;
    polygonDB._stmt = 0;
    polygonDB._readwrite = false;
  }
#endif

  /**
   * @brief ~PolygonDBInterface Destructor.
   */
  ~PolygonDBInterface() {
    close();
  }

  /**
   * @brief operator = Copy assignment operator.
   * @param polygonDB The PolygonDBInterface to copy.
   * @return A reference to this.
   */
  PolygonDBInterface & operator = (const PolygonDBInterface &polygonDB) {
    if (&polygonDB != this) {
      close();
      _dbFilename = polygonDB.getFilename();
      _db = 0;
      _queryType = EMPTY;
      _stmt = 0;
      _readwrite = polygonDB._readwrite;
    }
    return *this;
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief operator = Move assignment operator.
   * @param polygonDB The PolygonDBInterface to move.
   * @return A reference to this.
   */
  PolygonDBInterface & operator = (PolygonDBInterface&& polygonDB) {
    if (&polygonDB != this) {
      close();
      _dbFilename = std::move(polygonDB._dbFilename);
      _db = polygonDB._db;
      _queryType = polygonDB._queryType;
      _stmt = polygonDB._stmt;
      _readwrite = polygonDB._readwrite;
      polygonDB._dbFilename.clear();
      polygonDB._db = 0;
      polygonDB._queryType = EMPTY;
      polygonDB._stmt = 0;
      polygonDB._readwrite = false;
    }
    return *this;
  }
#endif

  /**
   * @brief getFilename return current DB filename.
   * @return
   */
  const std::string & getFilename() const {
    return _dbFilename;
  }

  /**
   * @brief setFilename sets the DB filename.
   * @param dbFilename
   */
  void setFilename(const std::string &dbFilename) {
    _dbFilename = dbFilename;
  }

  /**
   * @brief setReadWrite makes SQLite open the DB in read-write mode. Call this before initialize.
   */
  void setReadWrite(const bool readWrite) {
    _readwrite = readWrite;
  }

  /**
   * @brief initialize opens the DB connection.
   */
  void initialize() {
    if (_db != 0) {
      std::cout << "DB already initialized. Skipped." << std::endl;
      return;
    }
    if (_dbFilename.empty()) {
      std::cout << "No DB filename to open." << std::endl;
      return;
    }

    // open the database
    int rc = 0;
    int flags = SQLITE_OPEN_FULLMUTEX | SQLITE_OPEN_SHAREDCACHE;
    if (_readwrite)
      flags |= SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE;
    else
      flags |= SQLITE_OPEN_READONLY;
    rc = sqlite3_open_v2(_dbFilename.c_str(), &_db, flags, 0);
    if (rc != SQLITE_OK) {
      sqlite3_close(_db);
      _db = 0;
      std::cout << "Polygon db file '" << _dbFilename << "' not open." << std::endl;
      return;
    }
    if (_readwrite)
      createPolygonsTable();
  }

  /**
   * @brief isInitialized returns true if a DB connection has been already created.
   * @return
   */
  bool isInitialized() const {
    return _db != 0;
  }

  /**
   * @brief numberOfRows returns how many rows are extracted by a SELECT query for patches.
   * @param topology The topology string of the SELECT query.
   * @return The number of rows extracted by the query.
   */
  count_type numberOfRows(const std::string &topology) {
    // check if the database is open
    if (_db == 0) {
      std::cout << "Not initialized. Impossible to prepare the count query." << std::endl;
      return 0;
    }

    // reset the previous statement
    finalize();

    // query
    std::stringstream stream;
    stream << "SELECT count(*) FROM " << POLYGON_TABLE_NAME();
    stream << " INDEXED BY " << POLYGON_TABLE_NAME() << INDEX_POSTFIX();
    stream << " WHERE " << COLUMN_TOPOLOGY() << " = '" << topology << "'";
    stream << " ORDER BY ";
    stream << COLUMN_MESH_FILENAME() << ", ";
    stream << COLUMN_N_FACES() << ", ";
    stream << COLUMN_POLYGON();
    stream << ';';
    sqlite3_prepare_v2(_db, stream.str().c_str(), stream.str().size() + 1, &_stmt, 0);
    sqlite3_step(_stmt);

    // get the number of rows
    sqlite3_int64 number = sqlite3_column_int64(_stmt, 0);

    // terminate
    finalize();

    // return
    return number;
  }

  /**
   * @brief querySelect prepares the query SELECT to get the rows with COLUMN_TOPOLOGY() = topology.
   * @param topology The topology key string to search.
   */
  void querySelect(const std::string &topology) {
    // check if the database is open
    if (_db == 0) {
      std::cout << "Not initialized. Impossible to prepare the select query." << std::endl;
      return;
    }

    // reset the previous statement
    finalize();

    // query
    std::stringstream stream;
    stream << "SELECT ";
    stream << COLUMN_MESH_FILENAME();
    stream << ", " << COLUMN_POLYGON();
    stream << " FROM " << POLYGON_TABLE_NAME();
    stream << " INDEXED BY " << POLYGON_TABLE_NAME() << INDEX_POSTFIX();
    stream << " WHERE " << COLUMN_TOPOLOGY() << " = '" << topology << "'";
    stream << " ORDER BY ";
    stream << COLUMN_MESH_FILENAME() << ", ";
    stream << COLUMN_N_FACES() << ", ";
    stream << COLUMN_POLYGON();
    stream << ';';
    sqlite3_prepare_v2(_db, stream.str().c_str(), stream.str().size() + 1, &_stmt, 0);

    _queryType = SELECT;
  }

  /**
   * @brief getNextRow retrieves data from the next row of the current SELECT statement.
   * @param meshFilename The mesh filename.
   * @param polygonStr The polygon's corners' Pos in a string format.
   * @return true if there is a new row, false otherwise.
   */
  bool getNextRow(std::string &meshFilename, std::string &polygonStr) {
    // init parameters
    meshFilename.clear();
    polygonStr.clear();

    // check if the database is open
    if (_db == 0) {
      std::cout << "Not initialized. Impossible to retrieve data." << std::endl;
      return false;
    }

    // check query type
    if (_queryType != SELECT) {
      std::cout << "Wrong query. Impossible to retrive data with a non-SELECT statement." << std::endl;
      return false;
    }

    // check if there is any other row
    if (sqlite3_step(_stmt) != SQLITE_ROW)
      return false;

    // read mesh filename
    std::stringstream stream;
    stream << sqlite3_column_text(_stmt, 0);
    meshFilename = stream.str();

    // read polygon string
    stream.str(std::string());
    stream << sqlite3_column_text(_stmt, 1);
    polygonStr = stream.str();

    return true;
  }

  /**
   * @brief queryInsert prepares a INSERT query using named parameters. Use insertNewRow() to put a new row data into the db.
   */
  void queryInsert() {
    // check if the database is open
    if (_db == 0) {
      std::cout << "Not initialized. Impossible to prepare the insert query." << std::endl;
      return;
    }

    // reset the previous statement
    finalize();

    // check if writable
    if (!_readwrite) {
      std::cout << "Not writable. Impossible to prepare the insert query." << std::endl;
      return;
    }

    // begin transaction
    beginTransaction();

    // insert query
    std::stringstream stream;
    stream << "INSERT OR IGNORE INTO " << POLYGON_TABLE_NAME();
    stream << " VALUES (";
    stream <<   '$' << COLUMN_TOPOLOGY();
    stream << ", $" << COLUMN_MESH_FILENAME();
    stream << ", $" << COLUMN_POLYGON();
    stream << ", $" << COLUMN_N_FACES();
    stream << ");";
    sqlite3_prepare_v2(_db, stream.str().c_str(), stream.str().size() + 1, &_stmt, 0);

    _queryType = INSERT;
  }

  /**
   * @brief insertNewRow inserts a new row by binding named parameters to values.
   * @param topology The topology key string.
   * @param meshFilename The mesh filename.
   * @param polygonStr The polygon's corners' Pos in a string format.
   * @param nFaces The number of faces in the polygon used for sorting.
   */
  void insertNewRow(const std::string &topology, const std::string &meshFilename, const std::string &polygonStr, const num_type nFaces) {
    std::string str;
    int parameterIndex = 0;
    int rc = SQLITE_OK;

    // check if the database is open
    if (_db == 0) {
      std::cout << "Not initialized. Impossible to insert data." << std::endl;
      return;
    }

    // check query type
    if (_queryType != INSERT) {
      std::cout << "Wrong query. Impossible to insert data with a non-INSERT statement." << std::endl;
      return;
    }

    // reset the statement
    sqlite3_reset(_stmt);

    // bind parameters to values
    str = std::string("$") + COLUMN_TOPOLOGY();
    parameterIndex = sqlite3_bind_parameter_index(_stmt, str.c_str());
    rc = sqlite3_bind_text(_stmt, parameterIndex, topology.c_str(), -1, SQLITE_STATIC);
    assert(rc == SQLITE_OK);
    str = std::string("$") + COLUMN_MESH_FILENAME();
    parameterIndex = sqlite3_bind_parameter_index(_stmt, str.c_str());
    rc = sqlite3_bind_text(_stmt, parameterIndex, meshFilename.c_str(), -1, SQLITE_STATIC);
    assert(rc == SQLITE_OK);
    str = std::string("$") + COLUMN_POLYGON();
    parameterIndex = sqlite3_bind_parameter_index(_stmt, str.c_str());
    rc = sqlite3_bind_text(_stmt, parameterIndex, polygonStr.c_str(), -1, SQLITE_STATIC);
    assert(rc == SQLITE_OK);
    str = std::string("$") + COLUMN_N_FACES();
    parameterIndex = sqlite3_bind_parameter_index(_stmt, str.c_str());
    rc = sqlite3_bind_int64(_stmt, parameterIndex, nFaces);
    assert(rc == SQLITE_OK);

    // execute
    sqlite3_step(_stmt);
  }

  /**
   * @brief finalize finalizes any prepared statement.
   */
  void finalize() {
    // finalize stmt
    if (_stmt != 0) {
      sqlite3_finalize(_stmt);
      _stmt = 0;
      if (_queryType == INSERT)
        endTransaction();
    }
    _queryType = EMPTY;
  }

  /**
   * @brief close closes the DB connection.
   */
  void close() {
    finalize();
    // close the db
    if (_db != 0)
      sqlite3_close(_db);
    _db = 0;
  }

  /**
   * @brief PolygonCornersToString represents the corners' Pos of a polygon in a string format.
   * @param mesh The input mesh.
   * @param corners The corners's Pos of the polygon.
   * @param polygonStr The output string.
   */
  template < typename PolyMeshType >
  static void PolygonCornersToString(const PolyMeshType &mesh,
                                     const std::deque<typename MeshSupport<PolyMeshType>::PosType> &corners,
                                     std::string &polygonStr) {
    std::stringstream stream;
    for (size_t c = 0; c < corners.size(); ++c) {
      // print face index
      stream << POSFACELABEL << LABLINDEXSEP;
      stream << std::setw(2 * sizeof(num_type)) << std::setfill('0') << std::hex << std::noshowbase;
      stream << (num_type)vcg::tri::Index(mesh, corners[c].F());
      stream << INDSEPARATOR;
      // print edge index
      stream << POSEDGELABEL << LABLINDEXSEP;
      stream << std::setw(1) << std::hex << std::noshowbase;
      stream << corners[c].E();
      stream << INDSEPARATOR;
      // print vertex index
      stream << POSVERTLABEL << LABLINDEXSEP;
      stream << std::setw(2 * sizeof(num_type)) << std::setfill('0') << std::hex << std::noshowbase;
      stream << (num_type)vcg::tri::Index(mesh, corners[c].V());
      if (c < corners.size() - 1)
        stream << POSSEPARATOR;
    }
    polygonStr = stream.str();
  }

  /**
   * @brief PolygonCornersFromString converts the corners' Pos of a polygon from a string format.
   * @param mesh The input mesh.
   * @param corners The output corners's Pos of the polygon.
   * @param polygonStr The input string.
   */
  template < typename PolyMeshType >
  static void PolygonCornersFromString(PolyMeshType &mesh,
                                       std::deque<typename MeshSupport<PolyMeshType>::PosType> &corners,
                                       const std::string &polygonStr) {
    corners.clear();
    if (polygonStr.empty() || mesh.IsEmpty())
      return;

    typename MeshSupport<PolyMeshType>::PosType cornerPos;
    std::stringstream stream;
    size_t pos = 0, faceIndex = 0, vertexIndex = 0;
    int edgeIndex = 0;
    for (size_t i = 0; i < polygonStr.length(); i += 1/*POSFACELABEL*/ + 1/*LABINDEXSEP*/
                                                  + 2 * sizeof(num_type)/*face index*/
                                                  + 1/*INDSEPARATOR*/
                                                  + 1/*POSEDGELABEL*/ + 1/*LABINDEXSEP*/
                                                  + 1/*edge index*/
                                                  + 1/*INDSEPARATOR*/
                                                  + 1/*POSVERTLABEL*/ + 1/*LABINDEXSEP*/
                                                  + 2 * sizeof(num_type)/*vertex index*/
                                                  + 1/*POSSEPARATOR*/) {
      // read face index
      pos = i + 1/*POSFACELABEL*/ + 1/*LABINDEXSEP*/;
      stream.str(polygonStr.substr(pos, 2 * sizeof(num_type)) + ' ');
      stream >> std::hex >> faceIndex;
      // read edge index
      pos += 2 * sizeof(num_type)/*face index*/ + 1/*INDSEPARATOR*/ + 1/*POSEDGELABEL*/ + 1/*LABINDEXSEP*/;
      stream.str(polygonStr.substr(pos, 1) + ' ');
      stream >> edgeIndex;
      // read vertex index
      pos += 1/*edge index*/ + 1/*INDSEPARATOR*/ + 1/*POSVERTLABEL*/ + 1/*LABINDEXSEP*/;
      stream.str(polygonStr.substr(pos, 2 * sizeof(num_type)) + ' ');
      stream >> std::hex >> vertexIndex;
      // push the corner
      cornerPos.Set(&mesh.face[faceIndex], edgeIndex, &mesh.vert[vertexIndex]);
      corners.push_back(cornerPos);
    }
  }

private:
  /**
   * @brief The QueryType enum defines the type of the current query.
   */
  enum QueryType {
    EMPTY,
    SELECT,
    INSERT
  };

  /**
   * @brief beginTransaction starts a new INSERT transaction.
   */
  void beginTransaction() {
    if (_db == 0)
      return;
    if (_stmt != 0)
      sqlite3_finalize(_stmt);
    _stmt = 0;
    std::string str("BEGIN TRANSACTION;");
    sqlite3_prepare_v2(_db, str.c_str(), str.size() + 1, &_stmt, 0);
    sqlite3_step(_stmt);
    sqlite3_finalize(_stmt);
    _stmt = 0;
  }

  /**
   * @brief endTransaction ends an ongoing INSERT transaction.
   */
  void endTransaction() {
    if (_db == 0)
      return;
    if (_stmt != 0)
      sqlite3_finalize(_stmt);
    _stmt = 0;
    std::string str("COMMIT TRANSACTION;");
    sqlite3_prepare_v2(_db, str.c_str(), str.size() + 1, &_stmt, 0);
    sqlite3_step(_stmt);
    sqlite3_finalize(_stmt);
    _stmt = 0;
  }

  /**
   * @brief createPolygonsTable
   */
  void createPolygonsTable() {
    if (_db == 0)
      return;

    // finalize stmt
    if (_stmt != 0)
      sqlite3_finalize(_stmt);
    _stmt = 0;

    // create table (if it does not exist)
    std::stringstream stream;
    stream << "CREATE TABLE IF NOT EXISTS " << POLYGON_TABLE_NAME();
    stream << " (" << COLUMN_TOPOLOGY() << " TEXT NOT NULL";
    stream << ", " << COLUMN_MESH_FILENAME() << " TEXT NOT NULL";
    stream << ", " << COLUMN_POLYGON() << " TEXT NOT NULL";
    stream << ", " << COLUMN_N_FACES() << " INTEGER CHECK(" << COLUMN_N_FACES() << " > 0)";
    stream << ", PRIMARY KEY";
    stream << " (" << COLUMN_TOPOLOGY();
    stream << ", " << COLUMN_MESH_FILENAME();
    stream << ", " << COLUMN_POLYGON();
    stream << ')';
    stream << ')';
    stream << ';';
    sqlite3_prepare_v2(_db, stream.str().c_str(), stream.str().size() + 1, &_stmt, 0);
    sqlite3_step(_stmt);
    sqlite3_finalize(_stmt);
    _stmt = 0;
    // create index (if it does not exist)
    stream.str(std::string());
    stream << "CREATE INDEX IF NOT EXISTS " << POLYGON_TABLE_NAME() << INDEX_POSTFIX();
    stream << " ON " << POLYGON_TABLE_NAME();
    stream << " (";
    stream << COLUMN_TOPOLOGY();
    stream << ", " << COLUMN_MESH_FILENAME();
    stream << ", " << COLUMN_N_FACES();
    stream << ", " << COLUMN_POLYGON();
    stream << ");";
    sqlite3_prepare_v2(_db, stream.str().c_str(), stream.str().size() + 1, &_stmt, 0);
    sqlite3_step(_stmt);
    sqlite3_finalize(_stmt);
    _stmt = 0;
  }

  std::string     _dbFilename;    ///< Name of the DB file.
  sqlite3        *_db;            ///< DB accessed by sqlite API.
  QueryType       _queryType;     ///< Current query type.
  sqlite3_stmt   *_stmt;          ///< SQLite statement.
  bool            _readwrite;     ///< true if the DB must be opened in read-write mode, false if in read-only mode (default).
};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // POLYGON_DB_INTERFACE

