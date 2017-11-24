#ifndef THREAD_SAFE_POLYGON_REGISTER_H
#define THREAD_SAFE_POLYGON_REGISTER_H

#include <mutex>
#include "../polygon_register.h"

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief The TS_PolygonRegister class is a thread-safe version of PolygonRegister.
 */
template < typename PolyMeshType >
class TS_PolygonRegister : public PolygonRegister<PolyMeshType> {
public:
  typedef pl::PolygonRegister<PolyMeshType>     PolygonRegister;

  /**
   * @brief TS_PolygonRegister Default constructor.
   */
  TS_PolygonRegister() : PolygonRegister() { }

  /**
   * @brief TS_PolygonRegister Constructor which set the database filename specified by filename.
   * @param filename
   */
  TS_PolygonRegister(const std::string &filename) : PolygonRegister(filename) { }

  /**
   * @brief TS_PolygonRegister Copy constructor.
   * @param ts_pr The TS_PolygonRegister to copy.
   */
  TS_PolygonRegister(const TS_PolygonRegister &ts_pr) {
    std::lock_guard<std::mutex> ts_pr_lock(ts_pr._mutex);
    PolygonRegister::operator =(ts_pr);
  }

  /**
   * @brief TS_PolygonRegister Move constructor.
   * @param ts_pr The TS_PolygonRegister to move.
   */
  TS_PolygonRegister(TS_PolygonRegister&& ts_pr) {
    std::lock_guard<std::mutex> ts_pr_lock(ts_pr._mutex);
    PolygonRegister::operator =(std::move(ts_pr));
  }

  /**
   * @brief operator = Copy assignment operator.
   * @param ts_pr The TS_PolygonRegister to copy.
   * @return A reference to this.
   */
  TS_PolygonRegister & operator = (const TS_PolygonRegister &ts_pr) {
    if (&ts_pr != this) {
      std::lock_guard<std::mutex> ts_pr_lock(ts_pr._mutex);
      std::lock_guard<std::mutex> lock(_mutex);
      PolygonRegister::operator =(ts_pr);
    }
    return *this;
  }

  /**
   * @brief operator = Move assignment operator.
   * @param ts_pr The TS_PolygonRegister to move.
   * @return A reference to this.
   */
  TS_PolygonRegister & operator = (TS_PolygonRegister&& ts_pr) {
    if (&ts_pr != this) {
      std::lock_guard<std::mutex> ts_pr_lock(ts_pr._mutex);
      std::lock_guard<std::mutex> lock(_mutex);
      PolygonRegister::operator =(std::move(ts_pr));
    }
    return *this;
  }

  /**
   * @brief getFilename
   * @return
   */
  std::string getFilename() const {
    std::lock_guard<std::mutex> lock(_mutex);
    return PolygonRegister::getFilename();
  }

  /**
   * @brief setFilename
   * @note Call this before Initialize().
   * @param filename
   */
  void setFilename(const std::string &filename) {
    std::lock_guard<std::mutex> lock(_mutex);
    PolygonRegister::setFilename(filename);
  }

  /**
   * @brief initialize just does all the work to make the underlying container ready to use.
   */
  void initialize() {
    std::lock_guard<std::mutex> lock(_mutex);
    PolygonRegister::initialize();
  }

  /**
   * @brief isInitialized returns true if a DB connection has been already created.
   * @return
   */
  bool isInitialized() const {
    std::lock_guard<std::mutex> lock(_mutex);
    return PolygonRegister::isInitialized();
  }

  /**
   * @brief resetDB drops all the tables into the database, if open.
   */
  void resetDB() {
    std::lock_guard<std::mutex> lock(_mutex);
    PolygonRegister::resetDB();
  }

  /**
   * @brief setBufferSize sets the maximum size of the buffer.
   * @param bufferSize
   */
  void setBufferSize(const count_type bufferSize = std::numeric_limits<count_type>::max()) {
    std::lock_guard<std::mutex> lock(_mutex);
    PolygonRegister::setBufferSize(bufferSize);
  }

  /**
   * @brief isRegistered checks if a given polygon (represented by its boundary) has been already registered.
   * If it isn't, then this method register the polygon.
   * @param key is the unique key obtained with BuildKey().
   * @return true if already registered, false otherwise.
   */
  bool isRegistered(const std::string &key) {
    std::lock_guard<std::mutex> lock(_mutex);
    return PolygonRegister::isRegistered(key);
  }

  /**
   * @brief finalize just does all the work to make the underlying container consistent when it is not used anymore.
   */
  void finalize() {
    std::lock_guard<std::mutex> lock(_mutex);
    PolygonRegister::finalize();
  }

private:
  mutable std::mutex            _mutex;     ///< Mutual exclusion.
};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // THREAD_SAFE_POLYGON_REGISTER_H
