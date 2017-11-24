#ifndef POLYGON_REGISTER_H
#define POLYGON_REGISTER_H

#include <queue>
#include <string>
#include <sstream>
#include <iomanip>
#include <sqlite3.h>
#include "mesh_support.h"
#if __cplusplus >= 201103L || _MSC_VER >= 1800
  #include "string_support.h"
#endif

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief The PolygonRegister class is responsible of registering the visited polygons on a given mesh.
 * It is able to check if a given polygon has been already visited.
 */
template < typename PolyMeshType >
class PolygonRegister {
public:
  typedef typename MeshSupport<PolyMeshType>::PosType   PosType;

  /**
   * @brief REGISTER_TABLE_NAME is the name of the main table.
   * @return The name of the main table.
   */
  static const std::string & REGISTER_TABLE_NAME() {
    static const std::string registerTableName("Polygons");
    return registerTableName;
  }

  /**
   * @brief REGISTER_COLUMN_NAME is the name of the main column.
   * @return The name of the main column.
   */
  static const std::string & REGISTER_COLUMN_NAME() {
    static const std::string registerColumnName("polygon");
    return registerColumnName;
  }

  /**
   * @brief PolygonRegister Default constructor.
   */
  PolygonRegister() : _db(0), _stmt(0), _bufferSize(std::numeric_limits<count_type>::max()) { }

  /**
   * @brief PolygonRegister Constructor which set the database filename specified by filename.
   * @param filename
   */
  PolygonRegister(const std::string &filename) : _dbFilename(filename), _db(0), _stmt(0),
                                                 _bufferSize(std::numeric_limits<count_type>::max()) { }

  /**
   * @brief PolygonRegister Safe copy constructor.
   * @param polygonRegister
   */
  PolygonRegister(const PolygonRegister &polygonRegister) : _dbFilename(polygonRegister._dbFilename),
                                                            _db(0), _stmt(0),
                                                            _bufferSize(polygonRegister._bufferSize) { }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief PolygonRegister Move constructor.
   * @param polygonRegister The PolygonRegister to move.
   */
  PolygonRegister(PolygonRegister&& polygonRegister) : _dbFilename(std::move(polygonRegister._dbFilename)),
                                                       _db(polygonRegister._db),
                                                       _stmt(polygonRegister._stmt),
                                                       _keysToInsert(std::move(polygonRegister._keysToInsert)),
                                                       _bufferSize(polygonRegister._bufferSize) {
    polygonRegister._dbFilename.clear();
    polygonRegister._db = 0;
    polygonRegister._stmt = 0;
    polygonRegister._keysToInsert.clear();
    polygonRegister._bufferSize = std::numeric_limits<count_type>::max();
  }
#endif

  /**
   * @brief ~PolygonRegister destructor.
   */
  virtual ~PolygonRegister() {
    finalize();
  }

  /**
   * @brief operator = Copy assignment operator.
   * @param polygonRegister The PolygonRegister to copy.
   * @return A reference to this.
   */
  PolygonRegister & operator = (const PolygonRegister &polygonRegister) {
    if (&polygonRegister != this) {
      // close the db
      if (_stmt != 0)
        sqlite3_finalize(_stmt);
      _stmt = 0;
      if (_db != 0)
        sqlite3_close(_db);
      _db = 0;
      _dbFilename = polygonRegister._dbFilename;
      _keysToInsert.clear();
      _bufferSize = polygonRegister._bufferSize;
    }
    return *this;
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief operator = Move assignment operator.
   * @param polygonRegister The PolygonRegister to move.
   * @return A reference to this.
   */
  PolygonRegister & operator = (PolygonRegister&& polygonRegister) {
    if (&polygonRegister != this) {
      _dbFilename = std::move(polygonRegister._dbFilename);
      _db = polygonRegister._db;
      _stmt = polygonRegister._stmt;
      _keysToInsert = std::move(polygonRegister._keysToInsert);
      _bufferSize = polygonRegister._bufferSize;
      polygonRegister._dbFilename.clear();
      polygonRegister._db = 0;
      polygonRegister._stmt = 0;
      polygonRegister._keysToInsert.clear();
      polygonRegister._bufferSize = std::numeric_limits<count_type>::max();
    }
    return *this;
  }
#endif

  /**
   * @brief getFilename returns the current DB filename.
   * @return The DB filename.
   */
  virtual std::string getFilename() const {
    return _dbFilename;
  }

  /**
   * @brief setFilename
   * @note Call this before Initialize().
   * @param filename
   */
  virtual void setFilename(const std::string &filename) {
    _dbFilename = filename;
  }

  /**
   * @brief initialize just does all the work to make the underlying container ready to use.
   */
  virtual void initialize() {
    int rc = 0;
    std::stringstream stream;
    std::string str;

    if (_db != 0) {
      std::cout << "DB already initialized. Skipped." << std::endl;
      return;
    }
    if (_stmt != 0)
      sqlite3_finalize(_stmt);
    _stmt = 0;
    // open the database
    rc = sqlite3_open_v2(_dbFilename.c_str(), &_db,
                         SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_PRIVATECACHE | SQLITE_OPEN_NOMUTEX,
                         0);
    if (rc != SQLITE_OK) {
      std::cout << (_dbFilename.empty() || _dbFilename == ":memory:" ? "Temporary db" : _dbFilename) << " not open." << std::endl;
      sqlite3_close(_db);
      _db = 0;
    } else {
      // create the table
      createTable();
    }
  }

  /**
   * @brief isInitialized returns true if a DB connection has been already created.
   * @return
   */
  virtual bool isInitialized() const {
    return _db != 0;
  }

  /**
   * @brief resetDB drops all the tables into the database, if open.
   */
  virtual void resetDB() {
    std::queue<std::string> dropStringQueue;
    std::stringstream stream;
    std::string str;
    int rc = 0;

    if (_db == 0) {
      std::cout << "Resetting a NULL database. Skipped." << std::endl;
      return;
    }
    if (_stmt != 0)
      sqlite3_finalize(_stmt);
    _stmt = 0;

    // drop all
    str = "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name;";
    sqlite3_prepare_v2(_db, str.c_str(), str.size()+1, &_stmt, 0);
    do {
      rc = sqlite3_step(_stmt);
      if (rc == SQLITE_ROW) {
        stream.str(std::string());
        const unsigned char *name = sqlite3_column_text(_stmt, 0);
        stream << "DROP TABLE " << name << ";";
        dropStringQueue.push(stream.str());
      }
    } while (rc != SQLITE_DONE);
    sqlite3_finalize(_stmt);
    _stmt = 0;
    while (!dropStringQueue.empty()) {
      str = dropStringQueue.front();
      dropStringQueue.pop();
      sqlite3_prepare_v2(_db, str.c_str(), str.size()+1, &_stmt, 0);
      sqlite3_step(_stmt);
      sqlite3_finalize(_stmt);
      _stmt = 0;
    }
    // recreate table
    createTable();
  }

  /**
   * @brief setBufferSize sets the maximum size of the buffer.
   * @param bufferSize
   */
  virtual void setBufferSize(const count_type bufferSize = std::numeric_limits<count_type>::max()) {
    _bufferSize = bufferSize;
  }

  /**
   * @brief BuildKey construct a string key representing a polygon (boundary) to be used as unique value in the set
   * of the registered polygons.
   * @param mesh is used for obtaining the index of the vertices.
   * @param corners is a std::vector of vcg::face::Pos which are the corners of the boundary.
   * @param key is a std::string being the final key value.
   */
  static void BuildKey(const PolyMeshType &mesh, const std::deque<PosType> &corners, std::string &key) {
    PosType start;
    std::stringstream stream;

    // store the number of sides
    stream << std::setw(2 * sizeof(n_corners_type)) << std::setfill('0') << std::hex << std::noshowbase;
    stream << corners.size();

    // if it has at most 2 sides, use all the vertices
    if (corners.size() <= 2) {
      start = corners.front();
      do {
        stream << std::setw(2 * sizeof(num_type)) << std::setfill('0') << std::hex << std::noshowbase;
        stream << (num_type)vcg::tri::Index(mesh, start.V());
        start.FlipV();
        start.FlipE();
        while (!start.IsBorder() && vcg::tri::IsMarked(mesh, start.FFlip())) {
          start.FlipF();
          start.FlipE();
        }
      } while (start != corners.front());

    } else {  // else, use only the corners
      for (size_t j = 0; j < corners.size(); ++j) {
        stream << std::setw(2 * sizeof(num_type)) << std::setfill('0') << std::hex << std::noshowbase;
        stream << (num_type)vcg::tri::Index(mesh, corners[j].V());
      }
    }

    // copy the final string
    key = stream.str();
  }

  /**
   * @brief isRegistered checks if a given polygon (represented by its boundary) has been already registered.
   * If it isn't, then this method register the polygon.
   * @param key is the unique key obtained with BuildKey().
   * @return true if already registered, false otherwise.
   */
  virtual bool isRegistered(const std::string &key) {
    int parameterIndex = 0;
    int rc = 0;
    bool isAlreadyRegistered;

    // check if the database is open
    if (_db == 0) {
      std::cout << "NULL database. Not registered." << std::endl;
      return false;
    }

    // if found in the buffer, return
#if __cplusplus >= 201103L || _MSC_VER >= 1800
    Hasher::HashValue hashValue = Hasher::StringToHashValue(key);
    if (_keysToInsert.find(hashValue) != _keysToInsert.end())
#else
    if (_keysToInsert.find(key) != _keysToInsert.end())
#endif
      return true;

    // else search in the DB
    if (_stmt == 0)
      querySearch();
    sqlite3_reset(_stmt);
    std::string str(std::string("$") + REGISTER_COLUMN_NAME());
    parameterIndex = sqlite3_bind_parameter_index(_stmt, str.c_str());
#if __cplusplus >= 201103L || _MSC_VER >= 1800
    rc = sqlite3_bind_int64(_stmt, parameterIndex, hashValue);
#else
    rc = sqlite3_bind_text(_stmt, parameterIndex, key.c_str(), -1, SQLITE_STATIC);
#endif
    assert(rc == SQLITE_OK);
    rc = sqlite3_step(_stmt);
    isAlreadyRegistered = rc == SQLITE_ROW;

    if (!isAlreadyRegistered) {
      // register it in the buffer first
      try {
#if __cplusplus >= 201103L || _MSC_VER >= 1800
        _keysToInsert.insert(hashValue);
#else
        _keysToInsert.insert(key);
#endif
      } catch(std::bad_alloc &e) {
        std::cout << "bad_alloc caught: " << e.what() << std::endl;
        writeAllIntoTheDB();
        try {
#if __cplusplus >= 201103L || _MSC_VER >= 1800
          _keysToInsert.insert(hashValue);
#else
          _keysToInsert.insert(key);
#endif
        } catch (std::bad_alloc &ba) {
          std::cerr << "bad_alloc caught: " << ba.what() << ". Something's gone wrong. Exiting." << std::endl;
          sqlite3_close(_db);
          exit(1);
        }
      }

      // check the buffer size and flush
      if (_keysToInsert.size() >= _bufferSize)
        writeAllIntoTheDB();
    }

    return isAlreadyRegistered;
  }

  /**
   * @brief finalize just does all the work to make the underlying container consistent when it is not used anymore.
   */
  virtual void finalize() {
    _keysToInsert.clear();
    // close the db
    if (_stmt != 0)
      sqlite3_finalize(_stmt);
    _stmt = 0;
    if (_db != 0)
      sqlite3_close(_db);
    _db = 0;
  }



private:
  /**
   * @brief querySearch starts a query to select rows.
   */
  void querySearch() {
    if (_db == 0)
      return;
    if (_stmt != 0)
      sqlite3_finalize(_stmt);
    _stmt = 0;
    std::stringstream stream;
    stream << "SELECT " << REGISTER_COLUMN_NAME() << " FROM " << REGISTER_TABLE_NAME();
    stream << " WHERE " << REGISTER_COLUMN_NAME() << " = $" << REGISTER_COLUMN_NAME() << ";";
    sqlite3_prepare_v2(_db, stream.str().c_str(), stream.str().size()+1, &_stmt, 0);
  }

  /**
   * @brief writeAllIntoTheDB flushes the buffer writing all the registered keys into the DB.
   */
  void writeAllIntoTheDB() {
    if (_db == 0) {
      std::cout << "Database not open. Impossible to write keys." << std::endl;
      return;
    }
    if (_stmt != 0)
      sqlite3_finalize(_stmt);
    _stmt = 0;

    // write patch
    int parameterIndex = 0;
    int rc = SQLITE_OK;
    std::stringstream stream;
    std::string str("BEGIN TRANSACTION;");
    sqlite3_prepare_v2(_db, str.c_str(), str.size() + 1, &_stmt, 0);
    sqlite3_step(_stmt);
    sqlite3_finalize(_stmt);
    _stmt = 0;
    stream << "INSERT INTO " << REGISTER_TABLE_NAME();
    stream << " VALUES ($" << REGISTER_COLUMN_NAME() << ");";
    sqlite3_prepare_v2(_db, stream.str().c_str(), stream.str().size()+1, &_stmt, 0);
    str = std::string("$") + REGISTER_COLUMN_NAME();
    parameterIndex = sqlite3_bind_parameter_index(_stmt, str.c_str());
    while (!_keysToInsert.empty()) {
      // insert new row
      sqlite3_reset(_stmt);
#if __cplusplus >= 201103L || _MSC_VER >= 1800
      rc = sqlite3_bind_int64(_stmt, parameterIndex, *_keysToInsert.begin());
#else
      rc = sqlite3_bind_text(_stmt, parameterIndex, _keysToInsert.begin()->c_str(), -1, SQLITE_STATIC);
#endif
      assert(rc == SQLITE_OK);
      sqlite3_step(_stmt);
      _keysToInsert.erase(_keysToInsert.begin());
    }
    sqlite3_finalize(_stmt);
    _stmt = 0;
    str = "COMMIT TRANSACTION;";
    sqlite3_prepare_v2(_db, str.c_str(), str.size() + 1, &_stmt, 0);
    sqlite3_step(_stmt);
    sqlite3_finalize(_stmt);
    _stmt = 0;
  }

  /**
   * @brief createTable creates the tables in the DB.
   */
  void createTable() {
    if (_db == 0) {
      std::cout << "Database not open. Impossible to create table." << std::endl;
      return;
    }
    if (_stmt != 0)
      sqlite3_finalize(_stmt);
    _stmt = 0;

    // set cache size to 10MB
    std::string str = "PRAGMA cache_size = 10240;";
    sqlite3_prepare_v2(_db, str.c_str(), str.size()+1, &_stmt, 0);
    sqlite3_step(_stmt);
    sqlite3_finalize(_stmt);
    _stmt = 0;
    // create table
    std::stringstream stream;
    stream << "CREATE TABLE IF NOT EXISTS " << REGISTER_TABLE_NAME();
    stream << " (" << REGISTER_COLUMN_NAME();
#if __cplusplus >= 201103L || _MSC_VER >= 1800
    stream << " INTEGER ";
#else
    stream << " TEXT ";
#endif
    stream << "PRIMARY KEY)";
    if (SQLITE_VERSION_NUMBER >= 3008002)
      stream << " WITHOUT ROWID";
    stream << ';';
    str = stream.str();
    sqlite3_prepare_v2(_db, str.c_str(), str.size()+1, &_stmt, 0);
    sqlite3_step(_stmt);
    sqlite3_finalize(_stmt);
    _stmt = 0;
  }



  std::string       _dbFilename;        ///< The filename of the DB.
  sqlite3          *_db;                ///< The Db connection.
  sqlite3_stmt     *_stmt;              ///< The DB query statement.

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  typedef Hasher::HashValue             HashValue;
  typedef std::unordered_set<HashValue> RegisterBuffer;
#else
  typedef std::set<std::string>         RegisterBuffer;
#endif
  typedef RegisterBuffer::iterator      RegisterBufferIterator;

  RegisterBuffer    _keysToInsert;      ///< Buffer to keep keys in memory before writing them all in a block in the DB.
  count_type        _bufferSize;        ///< Maximum size of the map.
};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // POLYGON_REGISTER_H
