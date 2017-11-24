#ifndef PATCH_LEARNER_H
#define PATCH_LEARNER_H

#include "patch_db_interface.h"
#include "polygon_db_interface.h"
#include "patch.h"
#if __cplusplus >= 201103L || _MSC_VER >= 1800
  #include "string_support.h"
#endif

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief The PatchLearner class is used to acquire and store in the DB a number of patches.
 *
 * Usage:
 * \code
 * std::string dbFilename = ...;
 * Patch<PolyMeshType> patch;
 * ...  // set the patch
 * PatchLearner<PolyMeshType> patchLearner(dbFilename);
 * patchLearner.initialize();
 * patcheLearner.acquirePatch(patch);
 * patcheLearner.finalize();
 * \endcode
 */
template < typename PolyMeshType >
class PatchLearner {
public:
  typedef pl::Patch<PolyMeshType>     Patch;

  /**
   * @brief PatchLearner Default constructor.
   */
  PatchLearner() : _nNewPatches(0),
                   _insertBufferSize(std::numeric_limits<count_type>::max()),
                   _updateBufferSize(std::numeric_limits<count_type>::max()),
                   _counterSize(0) { }

  /**
   * @brief PatchLearner Constructor setting the DB filename.
   * @param dbFilename
   */
  PatchLearner(const std::string &dbFilename) : _patchDB(dbFilename),
                                                _nNewPatches(0),
                                                _insertBufferSize(std::numeric_limits<count_type>::max()),
                                                _updateBufferSize(std::numeric_limits<count_type>::max()),
                                                _counterSize(0) { }

  /**
   * @brief PatchLearner Copy constructor.
   * @param cpl The PatchLearner to copy.
   */
  PatchLearner(const PatchLearner &cpl) : _patchDB(cpl._patchDB),
                                          _patchCounter(cpl._patchCounter),
                                          _nNewPatches(cpl._nNewPatches),
                                          _patchDataToInsert(cpl._patchDataToInsert),
                                          _patchDataToUpdate(cpl._patchDataToUpdate),
                                          _insertBufferSize(cpl._insertBufferSize),
                                          _updateBufferSize(cpl._updateBufferSize),
                                          _polygonDB(cpl._polygonDB),
                                          _polygonDataToInsert(cpl._polygonDataToInsert),
                                          _meshFilename(cpl._meshFilename),
                                          _counterSize(cpl._counterSize) { }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief PatchLearner Move constructor.
   * @param cpl The PatchLearner to move.
   */
  PatchLearner(PatchLearner&& cpl) : _patchDB(std::move(cpl._patchDB)),
                                     _patchCounter(std::move(cpl._patchCounter)),
                                     _nNewPatches(cpl._nNewPatches),
                                     _patchDataToInsert(std::move(cpl._patchDataToInsert)),
                                     _patchDataToUpdate(std::move(cpl._patchDataToUpdate)),
                                     _insertBufferSize(cpl._insertBufferSize),
                                     _updateBufferSize(cpl._updateBufferSize),
                                     _polygonDB(std::move(cpl._polygonDB)),
                                     _polygonDataToInsert(std::move(cpl._polygonDataToInsert)),
                                     _meshFilename(std::move(cpl._meshFilename)),
                                     _counterSize(cpl._counterSize) {
    cpl._patchCounter.clear();
    cpl._nNewPatches = 0;
    cpl._patchDataToInsert.clear();
    cpl._patchDataToUpdate.clear();
    cpl._polygonDataToInsert.clear();
    _meshFilename.clear();
  }
#endif

  /**
   * @brief ~PatchLearner Destructor.
   */
  virtual ~PatchLearner() {
    finalize();
  }

  /**
   * @brief operator = Copy assignment operator.
   * @param cpl The PatchLearner to copy.
   * @return A reference to this.
   */
  PatchLearner & operator = (const PatchLearner &cpl) {
    if (&cpl != this) {
      flushAndClose();
      _patchDB = cpl._patchDB;
      _patchCounter = cpl._patchCounter;
      _nNewPatches = cpl._nNewPatches;
      _patchDataToInsert = cpl._patchDataToInsert;
      _patchDataToUpdate = cpl._patchDataToUpdate;
      _insertBufferSize = cpl._insertBufferSize;
      _updateBufferSize = cpl._updateBufferSize;
      _polygonDB = cpl._polygonDB;
      _polygonDataToInsert = cpl._polygonDataToInsert;
      _meshFilename = cpl._meshFilename;
      _counterSize = cpl._counterSize;
    }
    return *this;
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief operator = Move assignment operator.
   * @param cpl The PatchLearner to move.
   * @return A reference to this.
   */
  PatchLearner & operator = (PatchLearner&& cpl) {
    if (&cpl != this) {
      flushAndClose();
      _patchDB = std::move(cpl._patchDB);
      _patchCounter = std::move(cpl._patchCounter);
      _nNewPatches = cpl._nNewPatches;
      _patchDataToInsert = std::move(cpl._patchDataToInsert);
      _patchDataToUpdate = std::move(cpl._patchDataToUpdate);
      _insertBufferSize = cpl._insertBufferSize;
      _updateBufferSize = cpl._updateBufferSize;
      _polygonDB = std::move(cpl._polygonDB);
      _polygonDataToInsert = std::move(cpl._polygonDataToInsert);
      _meshFilename = std::move(cpl._meshFilename);
      _counterSize = cpl._counterSize;
      cpl._patchCounter.clear();
      cpl._nNewPatches = 0;
      cpl._patchDataToInsert.clear();
      cpl._patchDataToUpdate.clear();
      cpl._meshFilename.clear();
    }
    return *this;
  }
#endif

  /**
   * @brief getNumberOfPatches gives the number of template patches found.
   * @return The number of template patches.
   */
  virtual count_type getNumberOfPatches() const {
    return _patchCounter.size();
  }

  /**
   * @brief getNumberOfNewPatches gives the number of template patches added to the DB.
   * @return The number of new template patches.
   */
  virtual count_type getNumberOfNewPatches() const {
    return _nNewPatches;
  }

  /**
   * @brief setPatchDBFilename sets the DB filename.
   * @param dbFilename
   */
  virtual void setPatchDBFilename(const std::string &dbFilename) {
    flushAndClose();
    _patchDB.setFilename(dbFilename);
  }

  /**
   * @brief setPolygonDBFilename sets the Poltgon DB filename.
   * @param dbFilename
   */
  virtual void setPolygonDBFilename(const std::string &dbFilename) {
    flushAndClose();
    _polygonDB.setFilename(dbFilename);
  }

  /**
   * @brief setMeshFilename sets the current mesh filename.
   * @param meshFilename
   */
  virtual void setMeshFilename(const std::string &meshFilename) {
    _polygonDataToInsert.clear();
    _polygonDB.close();
    _meshFilename = meshFilename;
  }

  /**
   * @brief initialize opens the DB and resets temporary data structures.
   */
  virtual void initialize() {
    _patchDB.setReadWrite(true);
    _patchDB.initialize();
    _patchDB.queryFind();
    if (!_polygonDB.getFilename().empty() && !_meshFilename.empty() && _counterSize > 0) {
      _polygonDB.setReadWrite(true);
      _polygonDB.initialize();
    }
    _patchCounter.clear();
    _nNewPatches = 0;
  }

  /**
   * @brief isInitialized returns true if a DB connection has been already created.
   * @return
   */
  virtual bool isInitialized() const {
    return _patchDB.isInitialized();
  }

  /**
   * @brief setInsertBufferSize sets the maximum size of the INSERT map.
   * @param bufferSize
   */
  virtual void setInsertBufferSize(const count_type bufferSize = std::numeric_limits<count_type>::max()) {
    _insertBufferSize = bufferSize;
  }

  /**
   * @brief setUpdateBufferSize sets the maximum size of the UPDATE map.
   * @param bufferSize
   */
  virtual void setUpdateBufferSize(const count_type bufferSize = std::numeric_limits<count_type>::max()) {
    _updateBufferSize = bufferSize;
  }

  /**
   * @brief setPolygonCounterSize sets the maximum size of the counter for each topology key.
   * @param counterSize
   */
  virtual void setPolygonCounterSize(const count_type counterSize = 0) {
    _counterSize = counterSize;
  }

  /**
   * @brief acquirePatch stores a patch and updates its probability.
   * @param patch
   */
  virtual void acquirePatch(Patch &patch) {
    if (!_patchDB.isInitialized()) {
      std::cout << "DB not initialized. Skipped acquisition." << std::endl;
      return;
    }

    // check input
    if (patch.isNull()) {
      std::cout << "NULL patch to acquire. Skipped acquisition." << std::endl;
      return;
    }

    PatchDBUpdateDataMapIterator updIt = _patchDataToUpdate.find(patch.getTopologyStr());
    if (updIt != _patchDataToUpdate.end())
      ++updIt->second;
    else {
#if __cplusplus >= 201103L || _MSC_VER >= 1800
      HashValue hashValue = Hasher::StringToHashValue(patch.getTopologyStr());
      PatchDBInsertDataMapIterator insIt = _patchDataToInsert.find(hashValue);
#else
      PatchDBInsertDataMapIterator insIt = _patchDataToInsert.find(patch.getTopologyStr());
#endif
        if (insIt != _patchDataToInsert.end()) {
//        ++insIt->second.nOccurrences;
        ///dxy test
        insIt->second.nOccurrences += patch.getNumberOfOccurrences();
        ///
        }
      else {
#if __cplusplus >= 201103L || _MSC_VER >= 1800
        if (_patchCounter.find(hashValue) == _patchCounter.end())
          _patchCounter.insert(hashValue);
#else
        if (_patchCounter.find(patch.getTopologyStr()) == _patchCounter.end())
          _patchCounter.insert(patch.getTopologyStr());
#endif
        // wrap data
        PatchDBDataPack data;
        data.topology = patch.getTopologyStr();
        patch.getBoundaryIDsStr(data.boundaryIDs);
        data.singularities = patch.getSingularityStr();
        data.nSingularities = patch.getNumberOfSingularities();
        data.nCorners = patch.getNumberOfCorners();
        data.nFaces = patch.getNumberOfFaces();
        data.nVertices = patch.getNumberOfVertices();
//        data.nOccurrences = 1;
          ///dxy test
          data.nOccurrences = patch.getNumberOfOccurrences();
          ///
        data.nConcaveCorners = patch.getNumberOfConcaveCorners();
        // check if present in the DB
        if (_patchDB.findRecord(data.topology) > 0) {
          // insert into update-cache
          try {
            _patchDataToUpdate.insert(std::make_pair(data.topology, data.nOccurrences));
          } catch (std::bad_alloc &e) {
            std::cout << "bad_alloc caught: " << e.what() << std::endl;
            flushUpdateBuffer();
            _patchDB.queryFind();
            try {
              _patchDataToUpdate.insert(std::make_pair(data.topology, data.nOccurrences));
            } catch (std::bad_alloc &ba) {
              std::cerr << "bad_alloc caught: " << ba.what() << ". Something's gone wrong. Exiting." << std::endl;
              flushAndClose();
              exit(1);
            }
          }
        } else {
          ++_nNewPatches;
          // insert into insert-cache
          try {
#if __cplusplus >= 201103L || _MSC_VER >= 1800
            _patchDataToInsert.insert(std::make_pair(hashValue, data));
#else
            _patchDataToInsert.insert(std::make_pair(data.topology, data));
#endif
          } catch (std::bad_alloc &e) {
            std::cout << "bad_alloc caught: " << e.what() << std::endl;
            flushInsertBuffer();
            _patchDB.queryFind();
            try {
#if __cplusplus >= 201103L || _MSC_VER >= 1800
              _patchDataToInsert.insert(std::make_pair(hashValue, data));
#else
              _patchDataToInsert.insert(std::make_pair(data.topology, data));
#endif
            } catch (std::bad_alloc &ba) {
              std::cerr << "bad_alloc caught: " << ba.what() << ". Something's gone wrong. Exiting." << std::endl;
              flushAndClose();
              exit(1);
            }
          }
        }
        // check the buffer size and flush
        bool resetQueryFind = false;
        if (_patchDataToUpdate.size() > _updateBufferSize) {
          flushUpdateBuffer();
          resetQueryFind = true;
        }
        if (_patchDataToInsert.size() > _insertBufferSize) {
          flushInsertBuffer();
          resetQueryFind = true;
        }
        if (resetQueryFind)
          _patchDB.queryFind();
      }
    }
  }

  /**
   * @brief bindTopologyToPolygon stores a link between a patch's topology and a polygon where the patch has been found.
   * @param topology The patch's topology string.
   * @param polygonStr The polygon's corners' Pos in string format.
   * @param nFaces The number of faces in the polygon used for sorting.
   */
  virtual void bindTopologyToPolygon(const std::string &topology, const std::string &polygonStr, const num_type nFaces) {
    if (!_polygonDB.isInitialized() || _counterSize == 0 || _meshFilename.empty() || topology.empty() || polygonStr.empty() || nFaces <= 0)
      return;

    PolygonData pData;
    pData.polygonStr = polygonStr;
    pData.nFaces = nFaces;

    // insert in the buffer
    PolygonDBInsertDataMapIterator pIt = _polygonDataToInsert.find(topology);
    if (pIt != _polygonDataToInsert.end()) {
      PolygonDataContainerIterator dIt = std::upper_bound(pIt->second.begin(), pIt->second.end(), pData);
      pIt->second.insert(dIt, pData);
    } else {
      pIt = _polygonDataToInsert.insert(std::make_pair(topology, PolygonDataContainer(1, pData))).first;
    }

    // check size
    if (pIt->second.size() > _counterSize)
      pIt->second.pop_back();
  }

  /**
   * @brief finalize writes all the in-memory patches into the DB and closes the DB connection.
   */
  virtual void finalize() {
    flushAndClose();
  }

private:
  /**
   * @brief flushAndClose writes all the in-memory patches into the DB and closes the DB connection.
   */
  void flushAndClose() {
    if (_patchDB.isInitialized())
      writeAllIntoTheDB();
    _patchDB.close();
    _polygonDB.close();
    _polygonDataToInsert.clear();
    _meshFilename.clear();
  }

  /**
   * @brief flushUpdateBuffer writes all updated data from the UPDATE buffer.
   */
  void flushUpdateBuffer() {
    if (!_patchDB.isInitialized()) {
      std::cout << "Database not open. Impossible to flush out patches." << std::endl;
      return;
    }
    // update patches
    _patchDB.queryUpdate();
    while (!_patchDataToUpdate.empty()) {
      _patchDB.updateRow(_patchDataToUpdate.begin()->first, _patchDataToUpdate.begin()->second);
      _patchDataToUpdate.erase(_patchDataToUpdate.begin());
    }
  }

  /**
   * @brief flushInsertBuffer writes all new data from the INSERT buffer.
   */
  void flushInsertBuffer() {
    if (!_patchDB.isInitialized()) {
      std::cout << "Database not open. Impossible to flush out patches." << std::endl;
      return;
    }
    // insert patches
    _patchDB.queryInsert();
    while (!_patchDataToInsert.empty()) {
      _patchDB.insertNewRow(_patchDataToInsert.begin()->second);
      _patchDataToInsert.erase(_patchDataToInsert.begin());
    }
  }

  /**
   * @brief flushPolygonBuffer writes all polygon data into the DB.
   */
  void flushPolygonBuffer() {
    if (!_polygonDB.isInitialized())
      return;
    PolygonDBInsertDataMapIterator pIt;
    _polygonDB.queryInsert();
    while (!_polygonDataToInsert.empty()) {
      pIt = _polygonDataToInsert.begin();
      for (size_t p = 0; p < pIt->second.size(); ++p)
        _polygonDB.insertNewRow(pIt->first, _meshFilename, pIt->second[p].polygonStr, pIt->second[p].nFaces);
      _polygonDataToInsert.erase(pIt);
    }
  }

  /**
   * @brief writeAllIntoTheDB performs all the writings into the DB, if any.
   */
  void writeAllIntoTheDB() {
    // update patches
    flushUpdateBuffer();
    // insert patches
    flushInsertBuffer();
    // insert polygons
    flushPolygonBuffer();
  }

  /**
   * @brief The PolygonData struct stores polygon's data to store into the polygon DB.
   */
  struct PolygonData {
    std::string     polygonStr;   ///< The string reprensenting the boundary of a polygon.
    num_type        nFaces;       ///< The number of faces for the sort criterion.

    /**
     * @brief PolygonData Default constructor.
     */
    PolygonData() : nFaces(0) { }

    /**
     * @brief PolygonData Copy constructor.
     * @param p The PolygonData to copy.
     */
    PolygonData(const PolygonData &p) : polygonStr(p.polygonStr), nFaces(p.nFaces) { }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
    /**
     * @brief PolygonData Move constructor.
     * @param p The PolygonData to move.
     */
    PolygonData(PolygonData&& p) : polygonStr(std::move(p.polygonStr)), nFaces(p.nFaces) { }
#endif

    /**
     * @brief operator = Copy assignment.
     * @param p The PolygonData to copy.
     * @return A reference to this.
     */
    PolygonData & operator = (const PolygonData &p) {
      if (&p != this) {
        polygonStr = p.polygonStr;
        nFaces = p.nFaces;
      }
      return *this;
    }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
    /**
     * @brief operator = Move assignment.
     * @param p The PolygonData to move.
     * @return A reference to this.
     */
    PolygonData & operator = (PolygonData&& p) {
      if (&p != this) {
        polygonStr = std::move(p.polygonStr);
        nFaces = p.nFaces;
      }
      return *this;
    }
#endif

    /**
     * @brief operator < defines the order between two PolygonData.
     * @param p The right PolygonData.
     * @return true if this comes first than p, false otherwise.
     */
    bool operator < (const PolygonData &p) const {
      if (nFaces != p.nFaces)
        return nFaces < p.nFaces;
      return polygonStr < p.polygonStr;
    }
  };

  typedef std::deque<PolygonData>                                 PolygonDataContainer;
  typedef typename PolygonDataContainer::iterator                 PolygonDataContainerIterator;
#if __cplusplus >= 201103L || _MSC_VER >= 1800
  typedef Hasher::HashValue                                       HashValue;
  typedef std::unordered_set<HashValue>                           PatchCounter;
  typedef std::unordered_map<HashValue, PatchDBDataPack>          PatchDBInsertDataMap;
  typedef std::unordered_map<std::string, count_type>             PatchDBUpdateDataMap;
  typedef std::unordered_map<std::string, PolygonDataContainer>   PolygonDBInsertDataMap;
#else
  typedef std::set<std::string>                                   PatchCounter;
  typedef std::map<std::string, PatchDBDataPack>                  PatchDBInsertDataMap;
  typedef std::map<std::string, count_type>                       PatchDBUpdateDataMap;
  typedef std::map<std::string, PolygonDataContainer>             PolygonDBInsertDataMap;
#endif
  typedef PatchDBInsertDataMap::iterator                          PatchDBInsertDataMapIterator;
  typedef PatchDBUpdateDataMap::iterator                          PatchDBUpdateDataMapIterator;
  typedef typename PolygonDBInsertDataMap::iterator               PolygonDBInsertDataMapIterator;

  PatchDBInterface            _patchDB;             ///< Database for patches.
  PatchCounter                _patchCounter;        ///< Set to store template patches in order to count how many of them are inserted.
  count_type                  _nNewPatches;         ///< How many new patches are inserted in the DB.
  PatchDBInsertDataMap        _patchDataToInsert;   ///< Map to keep patches in memory before inserting them all in a block in the DB.
  PatchDBUpdateDataMap        _patchDataToUpdate;   ///< Map to keep patches in memory before updating them all in a block in the DB.
  count_type                  _insertBufferSize;    ///< Maximum size of the INSERT map.
  count_type                  _updateBufferSize;    ///< Maximum size of the UPDATE map.
  PolygonDBInterface          _polygonDB;           ///< Database for template patches' source.
  PolygonDBInsertDataMap      _polygonDataToInsert; ///< Map to keep polygons in memory before inserting them al in a block in the DB.
  std::string                 _meshFilename;        ///< Current mesh filename.
  count_type                  _counterSize;         ///< Count limit for polygon insertions (for each topology key).
};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // PATCH_LEARNER_H
