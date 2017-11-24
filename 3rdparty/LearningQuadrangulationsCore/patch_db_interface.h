#ifndef PATCH_DB_INTERFACE_H
#define PATCH_DB_INTERFACE_H

#include <iostream>
#include <string>
#include <sstream>
#include <list>
#include <utility>
#include <sqlite3.h>
#include "default_types.h"

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief The PatchDBDataPack struct wraps all DB data.
 */
struct PatchDBDataPack {
  std::string     topology;
  std::string     boundaryIDs;
  std::string     singularities;
  num_type        nSingularities;
  n_corners_type  nCorners;
  num_type        nFaces;
  num_type        nVertices;
  count_type      nOccurrences;
  n_corners_type  nConcaveCorners;
  PatchDBDataPack() : nSingularities(0)
                    , nCorners(0)
                    , nFaces(0)
                    , nVertices(0)
                    , nOccurrences(0)
                    , nConcaveCorners(0) { }
};

/**
 * @brief The PatchDBInterface class defines the base operations to access the database.
 *
 * Usage:
 * \code
 * std::string dbFilename;
 * ...
 * PatchDBInterface patchDB;
 * patchDB.setFilename(dbFilename);
 * patchDB.setReadWrite(true);
 * patchDB.initialize();
 * ... // getters and/or setters here
 * patchDB.close();
 * \endcode
 */
class PatchDBInterface {
public:
  /**
   * @brief PATCHDBFILENAME The name of the DB file.
   * @return
   */
  static const std::string & PATCHDBFILENAME() {
    static const std::string patchDBFilename("patches.db");
    return patchDBFilename;
  }

  /**
   * @brief PATCH_TABLE_NAME Name of the table storing the collapsed patches.
   * @return
   */
  static const std::string & PATCH_TABLE_NAME() {
    static const std::string patchTableName("Patches");
    return patchTableName;
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
   * @brief COLUMN_BOUNDARY_IDS Header of boundary IDs column.
   * @return
   */
  static const std::string & COLUMN_BOUNDARY_IDS() {
    static const std::string columnBoundaryIDs("boundaryIDs");
    return columnBoundaryIDs;
  }

  /**
   * @brief COLUMN_SINGULARITIES Header of singularities column.
   * @return
   */
  static const std::string & COLUMN_SINGULARITIES() {
    static const std::string columnSingularities("singularities");
    return columnSingularities;
  }

  /**
   * @brief COLUMN_N_SINGULARITIES Header of number of singularities column.
   * @return
   */
  static const std::string & COLUMN_N_SINGULARITIES() {
    static const std::string columnNSingularities("nSingularities");
    return columnNSingularities;
  }

  /**
   * @brief COLUMN_N_CORNERS Header of number of corners column.
   * @return
   */
  static const std::string & COLUMN_N_CORNERS() {
    static const std::string columnNCorners("nCorners");
    return columnNCorners;
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
   * @brief COLUMN_N_VERTICES Header of number of vertices column.
   * @return
   */
  static const std::string & COLUMN_N_VERTICES() {
    static const std::string columnNVertices("nVertices");
    return columnNVertices;
  }

  /**
   * @brief COLUMN_N_OCCURRENCES Header of number of occurrences column.
   * @return
   */
  static const std::string & COLUMN_N_OCCURRENCES() {
    static const std::string columnNOccurrences("nOccurrences");
    return columnNOccurrences;
  }

  /**
   * @brief COLUMN_N_CONCAVE_CORNERS Header of number of concave corners column.
   * @return
   */
  static const std::string & COLUMN_N_CONCAVE_CORNERS() {
    static const std::string columnNConcaveCorners("nConcaveCorners");
    return columnNConcaveCorners;
  }

  /**
   * @brief PatchDBInterface Default constructor.
   */
  PatchDBInterface() : _db(0), _queryType(EMPTY), _stmt(0), _readwrite(false) { }

  /**
   * @brief PatchDBInterface Constructor setting the DB filename.
   * @param dbFilename
   */
  PatchDBInterface(const std::string &dbFilename) : _dbFilename(dbFilename), _db(0), _queryType(EMPTY),
                                                    _stmt(0), _readwrite(false) { }

  /**
   * @brief PatchDBInterface Safe copy constructor.
   * @param patchDB
   */
  PatchDBInterface(const PatchDBInterface &patchDB) : _dbFilename(patchDB.getFilename()), _db(0),
                                                      _queryType(EMPTY), _stmt(0),
                                                      _readwrite(patchDB._readwrite) { }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief PatchDBInterface Move constructor.
   * @param patchDB The PatchDBInterface will be moved and reset.
   */
  PatchDBInterface(PatchDBInterface&& patchDB) : _dbFilename(std::move(patchDB._dbFilename)), _db(patchDB._db),
                                                 _queryType(EMPTY), _stmt(patchDB._stmt),
                                                 _readwrite(patchDB._readwrite) {
    patchDB._dbFilename.clear();
    patchDB._db = 0;
    patchDB._queryType = EMPTY;
    patchDB._stmt = 0;
    patchDB._readwrite = false;
  }
#endif

  /**
   * @brief ~PatchDBInterface Destructor.
   */
  ~PatchDBInterface() {
    close();
  }

  /**
   * @brief operator = Copy assignment operator.
   * @param patchDB The PatchDBInterface to copy.
   * @return A reference to this.
   */
  PatchDBInterface & operator = (const PatchDBInterface &patchDB) {
    if (&patchDB != this) {
      close();
      _dbFilename = patchDB._dbFilename;
      _db = 0;
      _queryType = EMPTY;
      _stmt = 0;
      _readwrite = patchDB._readwrite;
    }
    return *this;
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief operator = Move assignment operator.
   * @param patchDB The PatchDBInterface to move.
   * @return A reference to this.
   */
  PatchDBInterface & operator = (PatchDBInterface&& patchDB) {
    if (&patchDB != this) {
      close();
      _dbFilename = std::move(patchDB._dbFilename);
      _db = patchDB._db;
      _queryType = patchDB._queryType;
      _stmt = patchDB._stmt;
      _readwrite = patchDB._readwrite;
      patchDB._dbFilename.clear();
      patchDB._db = 0;
      patchDB._queryType = EMPTY;
      patchDB._stmt = 0;
      patchDB._readwrite = false;
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
      std::cout << "Patch db file '" << _dbFilename << "' not open." << std::endl;
      return;
    }
    if (_readwrite)
      createPatchesTable();
  }

  /**
   * @brief isInitialized returns true if a DB connection has been already created.
   * @return
   */
  bool isInitialized() const {
    return _db != 0;
  }

  /**
   * @brief getCornerNumbers gives the list of distinct corner numbers that patches have.
   * @param cornerNumbers A list of numbers of corners.
   */
  void getCornerNumbers(std::list<n_corners_type> &cornerNumbers) {
    cornerNumbers.clear();

    // check if the database is open
    if (_db == 0) {
      std::cout << "Not initialized. Impossible to prepare the corner query." << std::endl;
      return;
    }

    // reset the previous statement
    finalize();

    // query
    std::stringstream stream;
    stream << "SELECT DISTINCT " << COLUMN_N_CORNERS() << " FROM " << PATCH_TABLE_NAME();
    stream << " INDEXED BY " << PATCH_TABLE_NAME() << INDEX_POSTFIX();
    stream << ';';
    sqlite3_prepare_v2(_db, stream.str().c_str(), stream.str().size() + 1, &_stmt, 0);

    // get the number of corners
    while (sqlite3_step(_stmt) == SQLITE_ROW)
      cornerNumbers.push_back(static_cast<n_corners_type>(sqlite3_column_int64(_stmt, 0)));

    // terminate
    finalize();
  }

  /**
   * @brief numberOfRows returns how many rows are extracted by a query for patches.
   * @param nCorners The number of corners. If < 0, then all patches are counted.
   * @param limit Only return 'limit' DB rows.
   * @param offset Skip the first 'offset' DB rows.
   * @param fromEnd true to return the less frequent patches, false the most frequent.
   * @return The number of rows extracted by the query.
   */
  count_type numberOfRows(const num_type nCorners = -1,
                          const count_type limit = std::numeric_limits<count_type>::max(),
                          const count_type offset = 0,
                          const bool fromEnd = false) {
    // check if the database is open
    if (_db == 0) {
      std::cout << "Not initialized. Impossible to prepare the count query." << std::endl;
      return 0;
    }

    // reset the previous statement
    finalize();

    // query
    std::stringstream stream;
    stream << "SELECT count(*) FROM " << PATCH_TABLE_NAME();
    stream << " INDEXED BY " << PATCH_TABLE_NAME() << INDEX_POSTFIX();
    if (nCorners >= 0)
      stream << " WHERE " << COLUMN_N_CORNERS() << " = " << nCorners;
    stream << " ORDER BY ";
    if (nCorners < 0) {
      if (fromEnd)
        stream << COLUMN_N_CORNERS() << " DESC, ";
      else
        stream << COLUMN_N_CORNERS() << " ASC, ";
    }
    if (fromEnd)
      stream << COLUMN_N_SINGULARITIES() << " DESC, " << COLUMN_SINGULARITIES() << " DESC, " << COLUMN_N_OCCURRENCES() << " ASC, " << COLUMN_N_FACES() << " DESC";
    else
      stream << COLUMN_N_SINGULARITIES() << " ASC, " << COLUMN_SINGULARITIES() << " ASC, " << COLUMN_N_OCCURRENCES() << " DESC, " << COLUMN_N_FACES() << " ASC";
    stream << " LIMIT " << (limit > std::numeric_limits<sqlite3_int64>::max() ? std::numeric_limits<sqlite3_int64>::max() : limit)
           << " OFFSET " << (offset > std::numeric_limits<sqlite3_int64>::max() ? std::numeric_limits<sqlite3_int64>::max() : offset);
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
   * @brief queryFind looks up the record with a given topology key, if it exists. Use findRecord() to search a desired topology.
   */
  void queryFind() {
    // check if the database is open
    if (_db == 0) {
      std::cout << "Not initialized. Impossible to prepare the find query." << std::endl;
      return;
    }

    // reset the previous statement
    finalize();

    // query
    std::stringstream stream;
    stream << "SELECT " << COLUMN_N_OCCURRENCES() << " FROM " << PATCH_TABLE_NAME();
    stream << " WHERE " << COLUMN_TOPOLOGY() << " = $" << COLUMN_TOPOLOGY() << ";";
    sqlite3_prepare_v2(_db, stream.str().c_str(), stream.str().size() + 1, &_stmt, 0);

    _queryType = FIND;
  }

  /**
   * @brief findRecord looks up the record with the given topology.
   * @param topology The topology key to look up.
   * @return The number of occurrences, or 0 if it does not exist.
   */
  count_type findRecord(const std::string &topology) {
    // check if the database is open
    if (_db == 0) {
      std::cout << "Not initialized. Impossible to find data." << std::endl;
      return 0;
    }

    // check query type
    if (_queryType != FIND) {
      std::cout << "Wrong query. Impossible to find data with a non-FIND statement." << std::endl;
      return 0;
    }

    // reset the statement
    sqlite3_reset(_stmt);

    // bind parameter
    std::string str(std::string("$") + COLUMN_TOPOLOGY());
    int parameterIndex = sqlite3_bind_parameter_index(_stmt, str.c_str());
    int rc = sqlite3_bind_text(_stmt, parameterIndex, topology.c_str(), -1, SQLITE_STATIC);
    assert(rc == SQLITE_OK);

    // execute
    if (sqlite3_step(_stmt) != SQLITE_ROW)
      return 0;

    // return number of occurrences
    return static_cast<count_type>(sqlite3_column_int64(_stmt, 0));
  }

  /**
   * @brief querySelect selects the rows with COLUMN_N_CORNERS() == nCorners. It is possible to restrict the set of results
   * by using limit and offset. Use getNextRow() to retrieve the results one at a time.
   * @param nCorners The number of corners.
   * @param limit The maximum number of rows to return.
   * @param offset The number of first rows to skip.
   * @param fromEnd true to return the less frequent patches, false the most frequent.
   */
  void querySelect(const num_type nCorners = -1,
                   const count_type limit = std::numeric_limits<count_type>::max(),
                   const count_type offset = 0,
                   const bool fromEnd = false) {
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
    stream << COLUMN_TOPOLOGY();
    stream << ", " << COLUMN_BOUNDARY_IDS();
    stream << ", " << COLUMN_SINGULARITIES();
    stream << ", " << COLUMN_N_SINGULARITIES();
    stream << ", " << COLUMN_N_CORNERS();
    stream << ", " << COLUMN_N_FACES();
    stream << ", " << COLUMN_N_VERTICES();
    stream << ", " << COLUMN_N_OCCURRENCES();
    stream << ", " << COLUMN_N_CONCAVE_CORNERS();
    stream << " FROM " << PATCH_TABLE_NAME();
    stream << " INDEXED BY " << PATCH_TABLE_NAME() << INDEX_POSTFIX();
    if (nCorners >= 0)
      stream << " WHERE " << COLUMN_N_CORNERS() << " = " << nCorners;
    stream << " ORDER BY ";
    if (nCorners < 0) {
      if (fromEnd)
        stream << COLUMN_N_CORNERS() << " DESC, ";
      else
        stream << COLUMN_N_CORNERS() << " ASC, ";
    }
    if (fromEnd)
      stream << COLUMN_N_SINGULARITIES() << " DESC, " << COLUMN_SINGULARITIES() << " DESC, " << COLUMN_N_OCCURRENCES() << " ASC, " << COLUMN_N_FACES() << " DESC";
    else
      stream << COLUMN_N_SINGULARITIES() << " ASC, " << COLUMN_SINGULARITIES() << " ASC, " << COLUMN_N_OCCURRENCES() << " DESC, " << COLUMN_N_FACES() << " ASC";
    stream << " LIMIT " << (limit > std::numeric_limits<sqlite3_int64>::max() ? std::numeric_limits<sqlite3_int64>::max() : limit)
           << " OFFSET " << (offset > std::numeric_limits<sqlite3_int64>::max() ? std::numeric_limits<sqlite3_int64>::max() : offset);
    stream << ';';
    sqlite3_prepare_v2(_db, stream.str().c_str(), stream.str().size() + 1, &_stmt, 0);

    _queryType = SELECT;
  }

  /**
   * @brief getNextRow retrieves data from the next row of the current SELECT statement.
   * @param data
   * @return true if there is a new row, false otherwise.
   */
  bool getNextRow(PatchDBDataPack &data) {
    // init parameters
    data.topology.clear();
    data.boundaryIDs.clear();
    data.singularities.clear();
    data.nSingularities = 0;
    data.nFaces = 0;
    data.nVertices = 0;
    data.nOccurrences = 0;
    data.nConcaveCorners = 0;

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

    // read topology
    std::stringstream stream;
    stream << sqlite3_column_text(_stmt, 0);
    data.topology = stream.str();

    // read boundary equations
    stream.str(std::string());
    stream << sqlite3_column_text(_stmt, 1);
    data.boundaryIDs = stream.str();

    // read singularities
    stream.str(std::string());
    stream << sqlite3_column_text(_stmt, 2);
    data.singularities = stream.str();

    // read number of singularities
    data.nSingularities = static_cast<num_type>(sqlite3_column_int64(_stmt, 3));

    // read number of faces
    data.nCorners = static_cast<n_corners_type>(sqlite3_column_int64(_stmt, 4));

    // read number of faces
    data.nFaces = static_cast<num_type>(sqlite3_column_int64(_stmt, 5));

    // read number of vertexes
    data.nVertices = static_cast<num_type>(sqlite3_column_int64(_stmt, 6));

    // read number of occurrences
    data.nOccurrences = static_cast<count_type>(sqlite3_column_int64(_stmt, 7));

    // read number of concave corners
    data.nConcaveCorners = static_cast<n_corners_type>(sqlite3_column_int64(_stmt, 8));

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
    stream << "INSERT OR IGNORE INTO " << PATCH_TABLE_NAME();
    stream << " VALUES (";
    stream <<   '$' << COLUMN_TOPOLOGY();
    stream << ", $" << COLUMN_BOUNDARY_IDS();
    stream << ", $" << COLUMN_SINGULARITIES();
    stream << ", $" << COLUMN_N_SINGULARITIES();
    stream << ", $" << COLUMN_N_CORNERS();
    stream << ", $" << COLUMN_N_FACES();
    stream << ", $" << COLUMN_N_VERTICES();
    stream << ", $" << COLUMN_N_OCCURRENCES();
    stream << ", $" << COLUMN_N_CONCAVE_CORNERS();
    stream << ");";
    sqlite3_prepare_v2(_db, stream.str().c_str(), stream.str().size() + 1, &_stmt, 0);

    _queryType = INSERT;
  }

  /**
   * @brief insertNewRow inserts a new row by binding named parameters to values.
   * @param data
   */
  void insertNewRow(const PatchDBDataPack &data) {
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
    rc = sqlite3_bind_text(_stmt, parameterIndex, data.topology.c_str(), -1, SQLITE_STATIC);
    assert(rc == SQLITE_OK);
    str = std::string("$") + COLUMN_BOUNDARY_IDS();
    parameterIndex = sqlite3_bind_parameter_index(_stmt, str.c_str());
    rc = sqlite3_bind_text(_stmt, parameterIndex, data.boundaryIDs.c_str(), -1, SQLITE_STATIC);
    assert(rc == SQLITE_OK);
    str = std::string("$") + COLUMN_SINGULARITIES();
    parameterIndex = sqlite3_bind_parameter_index(_stmt, str.c_str());
    rc = sqlite3_bind_text(_stmt, parameterIndex, data.singularities.c_str(), -1, SQLITE_STATIC);
    assert(rc == SQLITE_OK);
    str = std::string("$") + COLUMN_N_SINGULARITIES();
    parameterIndex = sqlite3_bind_parameter_index(_stmt, str.c_str());
    rc = sqlite3_bind_int64(_stmt, parameterIndex, data.nSingularities);
    assert(rc == SQLITE_OK);
    str = std::string("$") + COLUMN_N_CORNERS();
    parameterIndex = sqlite3_bind_parameter_index(_stmt, str.c_str());
    rc = sqlite3_bind_int64(_stmt, parameterIndex, data.nCorners);
    assert(rc == SQLITE_OK);
    str = std::string("$") + COLUMN_N_FACES();
    parameterIndex = sqlite3_bind_parameter_index(_stmt, str.c_str());
    rc = sqlite3_bind_int64(_stmt, parameterIndex, data.nFaces);
    assert(rc == SQLITE_OK);
    str = std::string("$") + COLUMN_N_VERTICES();
    parameterIndex = sqlite3_bind_parameter_index(_stmt, str.c_str());
    rc = sqlite3_bind_int64(_stmt, parameterIndex, data.nVertices);
    assert(rc == SQLITE_OK);
    str = std::string("$") + COLUMN_N_OCCURRENCES();
    parameterIndex = sqlite3_bind_parameter_index(_stmt, str.c_str());
    rc = sqlite3_bind_int64(_stmt, parameterIndex, data.nOccurrences);
    assert(rc == SQLITE_OK);
    str = std::string("$") + COLUMN_N_CONCAVE_CORNERS();
    parameterIndex = sqlite3_bind_parameter_index(_stmt, str.c_str());
    rc = sqlite3_bind_int64(_stmt, parameterIndex, data.nConcaveCorners);
    assert(rc == SQLITE_OK);

    // execute
    sqlite3_step(_stmt);
  }

  /**
   * @brief queryUpdate prepares a UPDATE query using named parameters. Use updateRow() to update a given row data.
   */
  void queryUpdate() {
    // check if the database is open
    if (_db == 0) {
      std::cout << "Not initialized. Impossible to prepare the update query." << std::endl;
      return;
    }

    // reset the previous statement
    finalize();

    // check if writable
    if (!_readwrite) {
      std::cout << "Not writable. Impossible to prepare the update query." << std::endl;
      return;
    }

    // begin transaction
    beginTransaction();

    // insert query
    std::stringstream stream;
    stream << "UPDATE OR IGNORE " << PATCH_TABLE_NAME();
    stream << " SET " << COLUMN_N_OCCURRENCES() << " = " << COLUMN_N_OCCURRENCES() << " + $" << COLUMN_N_OCCURRENCES();
    stream << " WHERE " << COLUMN_TOPOLOGY() << " = $" << COLUMN_TOPOLOGY() << ';';
    sqlite3_prepare_v2(_db, stream.str().c_str(), stream.str().size() + 1, &_stmt, 0);

    _queryType = UPDATE;
  }

  /**
   * @brief updateRow updates the row identified by topology, if exists, with COLUMN_N_OCCURRENCES() + nOccurrencesToAdd,
   * by binding named parameters.
   * @param topology
   * @param nOccurrencesToAdd
   */
  void updateRow(const std::string &topology, const count_type nOccurrencesToAdd) {
    std::string str;
    int parameterIndex = 0;
    int rc = SQLITE_OK;

    // check if the database is open
    if (_db == 0) {
      std::cout << "Not initialized. Impossible to update data." << std::endl;
      return;
    }

    // check query type
    if (_queryType != UPDATE) {
      std::cout << "Wrong query. Impossible to update data with a non-UPDATE statement." << std::endl;
      return;
    }

    // reset the statement
    sqlite3_reset(_stmt);

    // bind parameters to values
    str = std::string("$") + COLUMN_TOPOLOGY();
    parameterIndex = sqlite3_bind_parameter_index(_stmt, str.c_str());
    rc = sqlite3_bind_text(_stmt, parameterIndex, topology.c_str(), -1, SQLITE_STATIC);
    assert(rc == SQLITE_OK);
    str = std::string("$") + COLUMN_N_OCCURRENCES();
    parameterIndex = sqlite3_bind_parameter_index(_stmt, str.c_str());
    rc = sqlite3_bind_int64(_stmt, parameterIndex, nOccurrencesToAdd);
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
      if (_queryType == INSERT || _queryType == UPDATE)
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

private:
  /**
   * @brief The QueryType enum defines the type of the current query.
   */
  enum QueryType {
    EMPTY,
    FIND,
    SELECT,
    INSERT,
    UPDATE
  };

  /**
   * @brief beginTransaction starts a new INSERT/UPDATE transaction.
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
   * @brief endTransaction ends an ongoing INSERT/UPDATE transaction.
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
   * @brief createPatchesTable
   */
  void createPatchesTable() {
    if (_db == 0)
      return;

    // finalize stmt
    if (_stmt != 0)
      sqlite3_finalize(_stmt);
    _stmt = 0;

    // create table (if it does not exist)
    std::stringstream stream;
    stream << "CREATE TABLE IF NOT EXISTS " << PATCH_TABLE_NAME();
    stream << " (" << COLUMN_TOPOLOGY() << " TEXT PRIMARY KEY";
    stream << ", " << COLUMN_BOUNDARY_IDS() << " TEXT NOT NULL";
    stream << ", " << COLUMN_SINGULARITIES() << " TEXT NOT NULL";
    stream << ", " << COLUMN_N_SINGULARITIES() << " INTEGER CHECK (" << COLUMN_N_SINGULARITIES() << " >= 0)";
    stream << ", " << COLUMN_N_CORNERS() << " INTEGER CHECK (" << COLUMN_N_CORNERS() << " >= 0)";
    stream << ", " << COLUMN_N_FACES() << " INTEGER CHECK (" << COLUMN_N_FACES() << " > 0)";
    stream << ", " << COLUMN_N_VERTICES() << " INTEGER CHECK (" << COLUMN_N_VERTICES() << " > 0)";
    stream << ", " << COLUMN_N_OCCURRENCES() << " INTEGER CHECK (" << COLUMN_N_OCCURRENCES() << " > 0)";
    stream << ", " << COLUMN_N_CONCAVE_CORNERS() << " INTEGER CHECK (" << COLUMN_N_CONCAVE_CORNERS() << " >= 0)";
    stream << ")";
    stream << ";";
    sqlite3_prepare_v2(_db, stream.str().c_str(), stream.str().size() + 1, &_stmt, 0);
    sqlite3_step(_stmt);
    sqlite3_finalize(_stmt);
    _stmt = 0;
    // create index (if it does not exist)
    stream.str(std::string());
    stream << "CREATE INDEX IF NOT EXISTS " << PATCH_TABLE_NAME() << INDEX_POSTFIX();
    stream << " ON " << PATCH_TABLE_NAME();
    stream << " (" << COLUMN_N_CORNERS() << " ASC";
    stream << ", " << COLUMN_N_SINGULARITIES() << " ASC";
    stream << ", " << COLUMN_SINGULARITIES() << " ASC";
    stream << ", " << COLUMN_N_OCCURRENCES() << " DESC";
    stream << ", " << COLUMN_N_FACES() << " ASC";
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

#endif // PATCH_DB_INTERFACE_H
