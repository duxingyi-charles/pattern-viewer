#ifndef PATCH_DB_SERVER_H
#define PATCH_DB_SERVER_H

#include <deque>
#include <vector>
#include "patch_db_interface.h"
#include "patch.h"

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief The PatchDBServer class extracts patches from the DB and cache them.
 */
template < typename PolyMeshType >
class PatchDBServer {
public:
  typedef pl::Patch<PolyMeshType>         Patch;

  /**
   * @brief PatchDBServer Default constructor.
   */
  PatchDBServer() { }

  /**
   * @brief PatchDBServer Constructor binding to a DB file.
   * @param dbFilename The DB filename.
   */
  PatchDBServer(const std::string &dbFilename) : _patchDB(dbFilename) { }

  /**
   * @brief PatchDBServer Copy contructor.
   * @param ps The PatchDBServer to copy.
   */
  PatchDBServer(const PatchDBServer &ps) : _patchDB(ps._patchDB),
                                           _cache(ps._cache),
                                           _cached(ps._cached) { }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief PatchDBServer Move constructor.
   * @param cps The PatchDBServer to move.
   */
  PatchDBServer(PatchDBServer&& ps) : _patchDB(std::move(ps._patchDB)),
                                      _cache(std::move(ps._cache)),
                                      _cached(std::move(ps._cached)) {
    ps._cache.clear();
    ps._cached.clear();
  }
#endif

  /**
   * @brief operator = Copy assignment.
   * @param cps The PatchDBServer to copy.
   * @return A reference to this.
   */
  PatchDBServer & operator= (const PatchDBServer &ps) {
    if (this != &ps) {
      _patchDB = ps._patchDB;
      _cache = ps._cache;
      _cached = ps._cached;
    }
    return *this;
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief operator = Move assignment.
   * @param cps The PatchDBServer to move.
   * @return A reference to this.
   */
  PatchDBServer & operator= (PatchDBServer&& ps) {
    if (this != &ps) {
      _patchDB = std::move(ps._patchDB);
      _cache = std::move(ps._cache);
      _cached = std::move(ps._cached);
    }
    return *this;
  }
#endif

  /**
   * @brief getFilename return current DB filename.
   * @return
   */
  const std::string & getFilename() const {
    return _patchDB.getFilename();
  }

  /**
   * @brief setFilename sets the DB filename.
   * @param dbFilename
   */
  void setFilename(const std::string &dbFilename) {
    _patchDB.setFilename(dbFilename);
  }

  /**
   * @brief initialize sets up the DB connection.
   */
  void initialize() {
    _patchDB.setReadWrite(false);
    _patchDB.initialize();
  }

  /**
   * @brief isInitialized
   * @return true if it is initialized, false otherwise.
   */
  bool isInitialized() const {
    return _patchDB.isInitialized();
  }

  /**
   * @brief close clears everything and closes the DB connection.
   */
  void close() {
    _cache.clear();
    _cached.clear();
    _patchDB.close();
  }

  /**
   * @brief load caches all patches with n number of corners from the DB.
   * @param nCorners The number of corners.
   */
  void load(const n_corners_type nCorners) const {
    // check DB
    if (!_patchDB.isInitialized()) {
      std::cout << "DB not initialized. Impossible to cache." << std::endl;
      return;
    }

    // resize
    if (nCorners >= _cache.size()) {
      _cache.resize(nCorners + 1);
      _cached.resize(nCorners + 1, false);
    } else {
      _cache[nCorners].clear();
      _cached[nCorners] = false;
    }

    // cache the patches
    Patch patch;
    PatchDBDataPack data;
    _patchDB.querySelect(nCorners);
    while (_patchDB.getNextRow(data)) {
      // convert to ConvexPatch
      patch.setTopologyStr(data.topology);
      patch.setSingularityStr(data.singularities);
      patch.setNumberOfSingularities(data.nSingularities);
      patch.setNumberOfCorners(data.nCorners);
      patch.setNumberOfFaces(data.nFaces);
      patch.setNumberOfVertices(data.nVertices);
      patch.setNumberOfOccurrences(data.nOccurrences);
      patch.setNumberOfConcaveCorners(data.nConcaveCorners);
      patch.setNumberOfRotations(0);
      patch.setBoundaryFromStr(data.boundaryIDs);
      patch.setExpandedTopologyStr(data.topology);
      patch.reconstructPatch();
      // push into cache
#if __cplusplus >= 201103L || _MSC_VER >= 1800
      _cache[nCorners].push_back(std::move(patch));
#else
      _cache[nCorners].push_back(patch);
#endif
    }
    _patchDB.finalize();

    // patches with n number of corners cached
    _cached[nCorners] = true;
  }
    
    
    ///dxy add
    void myload(const n_corners_type nCorners) const {
        // check DB
        if (!_patchDB.isInitialized()) {
            std::cout << "DB not initialized. Impossible to cache." << std::endl;
            return;
        }
        
        // resize
        if (nCorners >= _cache.size()) {
            _cache.resize(nCorners + 1);
            _cached.resize(nCorners + 1, false);
        } else {
            if (_cached[nCorners]) {
                return;
            }
            _cache[nCorners].clear();
            _cached[nCorners] = false;
        }
        
        // cache the patches
        Patch patch;
        PatchDBDataPack data;
        _patchDB.querySelect(nCorners);
        while (_patchDB.getNextRow(data)) {
            // convert to ConvexPatch
            patch.setTopologyStr(data.topology);
            patch.setSingularityStr(data.singularities);
            patch.setNumberOfSingularities(data.nSingularities);
            patch.setNumberOfCorners(data.nCorners);
            patch.setNumberOfFaces(data.nFaces);
            patch.setNumberOfVertices(data.nVertices);
            patch.setNumberOfOccurrences(data.nOccurrences);
            patch.setNumberOfConcaveCorners(data.nConcaveCorners);
            patch.setNumberOfRotations(0);
            patch.setBoundaryFromStr(data.boundaryIDs);
            patch.setExpandedTopologyStr(data.topology);
            patch.reconstructPatch();
            // push into cache
#if __cplusplus >= 201103L || _MSC_VER >= 1800
            _cache[nCorners].push_back(std::move(patch));
#else
            _cache[nCorners].push_back(patch);
#endif
        }
        _patchDB.finalize();
        
        // patches with n number of corners cached
        _cached[nCorners] = true;
    }
    /// dxy add end

  /**
   * @brief load caches all patches.
   */
  void load() {
    // check DB
    if (!_patchDB.isInitialized()) {
      std::cout << "DB not initialized. Impossible to cache." << std::endl;
      return;
    }

    // clear current cache
    _cache.clear();
    _cached.clear();

    // get the corner numbers
    std::list<n_corners_type> cornerNumbers;
    _patchDB.getCornerNumbers(cornerNumbers);

    // resize and load
    _cache.resize(cornerNumbers.back() + 1);
    _cached.resize(_cache.size(), false);
    std::list<n_corners_type>::iterator cIt;
    for (cIt = cornerNumbers.begin(); cIt != cornerNumbers.end(); ++cIt)
      load(*cIt);
  }

  /**
   * @brief maxNumberOfCorners gives the number of buckets minus one.
   * @return The maximum number of corners.
   */
  n_corners_type maxNumberOfCorners() const {
    return _cache.size();
  }

  /**
   * @brief numberOfCachedPatches gives the number of patches with nCorners corners which are cached.
   * @param nCorners The number of corners.
   * @return how many cached patches.
   */
  count_type numberOfCachedPatches(const n_corners_type nCorners) const {
    if (nCorners >= _cache.size())
      return 0;
    return _cache[nCorners].size();
  }

  /**
   * @brief numberOfCachedPatches gives the total number of cached patches.
   * @return how many cached patches.
   */
  count_type numberOfCachedPatches() const {
    count_type tot = 0;
    for (size_t i = 0; i < _cache.size(); ++i)
      tot += _cache[i].size();
    return tot;
  }

  /**
   * @brief at gives constant access to the Patch at nCorners and patchIndex, or throw the std::out_of_range exeption.
   * @param nCorners The number of corners.
   * @param patchIndex The index of the patch among those with nCorners.
   * @return A constant reference to the requested patch.
   */
  const Patch & at(const n_corners_type nCorners, const count_type patchIndex) const {
    if (nCorners >= _cached.size() || !_cached[nCorners])
      load(nCorners);
    return _cache.at(nCorners).at(patchIndex);
  }

private:
  mutable PatchDBInterface                  _patchDB;   ///< Database.
  mutable std::deque< std::deque<Patch> >   _cache;     ///< Cache to store patches by number of corners.
  mutable std::vector<bool>                 _cached;    ///< true for cached patches.
};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // PATCH_DB_SERVER_H
