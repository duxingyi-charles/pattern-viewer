#ifndef THREAD_SAFE_PATCH_LEARNER_H
#define THREAD_SAFE_PATCH_LEARNER_H

#include <mutex>
#include "../patch_learner.h"

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief The TS_PatchLearner class is a thread-safe version of PatchLearner.
 */
template < typename PolyMeshType >
class TS_PatchLearner : public PatchLearner<PolyMeshType> {
public:
  typedef pl::PatchLearner<PolyMeshType>      PatchLearner;
  typedef typename PatchLearner::Patch        Patch;

  /**
   * @brief TS_PatchLearner Default constructor.
   */
  TS_PatchLearner() : PatchLearner() { }

  /**
   * @brief TS_PatchLearner Constructor setting the DB filename.
   * @param filename
   */
  TS_PatchLearner(const std::string &filename) : PatchLearner(filename) { }

  /**
   * @brief TS_PatchLearner Copy constructor.
   * @param ts_pl The TS_PatchLearner to copy.
   */
  TS_PatchLearner(const TS_PatchLearner &ts_pl) {
    std::lock_guard<std::mutex> ts_pl_lock(ts_pl._mutex);
    PatchLearner::operator =(ts_pl);
  }

  /**
   * @brief TS_PatchLearner Move constructor.
   * @param ts_pl The TS_PatchLearner to move.
   */
  TS_PatchLearner(TS_PatchLearner&& ts_pl) {
    std::lock_guard<std::mutex> ts_pl_lock(ts_pl._mutex);
    PatchLearner::operator =(std::move(ts_pl));
  }

  /**
   * @brief operator = Copy assigning operator.
   * @param ts_pl The TS_PatchLearner to copy.
   * @return A reference to this.
   */
  TS_PatchLearner & operator = (const TS_PatchLearner &ts_pl) {
    if (&ts_pl != this) {
      std::lock_guard<std::mutex> ts_pl_lock(ts_pl._mutex);
      std::lock_guard<std::mutex> lock(_mutex);
      PatchLearner::operator =(ts_pl);
    }
    return *this;
  }

  /**
   * @brief operator = Move assigning operator.
   * @param ts_pl The TS_PatchLearner to move.
   * @return A reference to this.
   */
  TS_PatchLearner & operator = (TS_PatchLearner&& ts_pl) {
    if (&ts_pl != this) {
      std::lock_guard<std::mutex> ts_pl_lock(ts_pl._mutex);
      std::lock_guard<std::mutex> lock(_mutex);
      PatchLearner::operator =(std::move(ts_pl));
    }
    return *this;
  }
  /**
   * @brief getNumberOfPatches gives the number of template patches found.
   * @return The number of template patches.
   */
  count_type getNumberOfPatches() const {
    std::lock_guard<std::mutex> lock(_mutex);
    return PatchLearner::getNumberOfPatches();
  }

  /**
   * @brief getNumberOfNewPatches gives the number of template patches added to the DB.
   * @return The number of new template patches.
   */
  count_type getNumberOfNewPatches() const {
    std::lock_guard<std::mutex> lock(_mutex);
    return PatchLearner::getNumberOfNewPatches();
  }

  /**
   * @brief setPatchDBFilename sets the DB filename.
   * @param dbFilename
   */
  void setPatchDBFilename(const string &dbFilename) {
    std::lock_guard<std::mutex> lock(_mutex);
    PatchLearner::setPatchDBFilename(dbFilename);
  }

  /**
   * @brief setPolygonDBFilename sets the Poltgon DB filename.
   * @param dbFilename
   */
  void setPolygonDBFilename(const std::string &dbFilename) {
    std::lock_guard<std::mutex> lock(_mutex);
    PatchLearner::setPolygonDBFilename(dbFilename);
  }

  /**
   * @brief setMeshFilename sets the current mesh filename.
   * @param meshFilename
   */
  void setMeshFilename(const std::string &meshFilename) {
    std::lock_guard<std::mutex> lock(_mutex);
    PatchLearner::setMeshFilename(meshFilename);
  }

  /**
   * @brief initialize opens the DB and resets temporary data structures.
   */
  void initialize() {
    std::lock_guard<std::mutex> lock(_mutex);
    PatchLearner::initialize();
  }

  /**
   * @brief isInitialized returns true if a DB connection has been already created.
   * @return
   */
  bool isInitialized() const {
    std::lock_guard<std::mutex> lock(_mutex);
    return PatchLearner::isInitialized();
  }

  /**
   * @brief setInsertBufferSize sets the maximum size of the INSERT map.
   * @param bufferSize
   */
  void setInsertBufferSize(const count_type bufferSize = 0) {
    std::lock_guard<std::mutex> lock(_mutex);
    PatchLearner::setInsertBufferSize(bufferSize);
  }

  /**
   * @brief setUpdateBufferSize sets the maximum size of the UPDATE map.
   * @param bufferSize
   */
  void setUpdateBufferSize(const count_type bufferSize = 0) {
    std::lock_guard<std::mutex> lock(_mutex);
    PatchLearner::setUpdateBufferSize(bufferSize);
  }

  /**
   * @brief setPolygonCounterSize sets the maximum size of the counter for each topology key.
   * @param counterSize
   */
  void setPolygonCounterSize(const count_type counterSize = std::numeric_limits<count_type>::max()) {
    std::lock_guard<std::mutex> lock(_mutex);
    PatchLearner::setPolygonCounterSize(counterSize);
  }

  /**
   * @brief acquirePatch stores a patch and updates its probability.
   * @param patch
   */
  void acquirePatch(Patch &patch) {
    std::lock_guard<std::mutex> lock(_mutex);
    PatchLearner::acquirePatch(patch);
  }

  /**
   * @brief bindTopologyToPolygon stores a link between a patch's topology and a polygon where the patch has been found.
   * @param topology The patch's topology string.
   * @param polygonStr The polygon's corners' Pos in string format.
   * @param nFaces The number of faces in the polygon used for sorting.
   */
  void bindTopologyToPolygon(const std::string &topology, const std::string &polygonStr, const num_type nFaces) {
    std::lock_guard<std::mutex> lock(_mutex);
    PatchLearner::bindTopologyToPolygon(topology, polygonStr, nFaces);
  }

  /**
   * @brief finalize writes all the in-memory patches into the DB and closes the DB connection.
   */
  void finalize() {
    std::lock_guard<std::mutex> lock(_mutex);
    PatchLearner::finalize();
  }

private:
  mutable std::mutex          _mutex;         ///< Mutual exclusion.
};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // THREAD_SAFE_PATCH_LEARNER_H
