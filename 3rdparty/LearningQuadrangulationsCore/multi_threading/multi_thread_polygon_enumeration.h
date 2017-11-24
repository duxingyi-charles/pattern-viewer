#ifndef MULTI_THREAD_POLYGON_ENUMERATION_H
#define MULTI_THREAD_POLYGON_ENUMERATION_H

#include <atomic>
#include <mutex>
#include <condition_variable>
#include <utility>
#include <thread>
#include <memory>
#include <algorithm>
#include "../polygon_enumeration.h"
#include "time_support.h"

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief The MT_PolygonEnumerationNotifiable class provides an interface for any class inheriting this one which want
 * to receive notifications on the ongoing process and/or its completion. See MT_PolygonEnumerator.
 */
class MT_PolygonEnumerationNotifiable {
public:
  /**
   * @brief The NotificationType enum defines the type of event to be sent/received.
   */
  enum NotificationType {
    EMPTY,          ///< Empty notification.
    TREE_BUILT,     ///< Notification on new tree built.
    FLUSH_DB,       ///< Notification on flushing the patches into the DB.
    COMPLETE        ///< Notification on process complete.
  };

  /**
   * @brief notify sends a notification to the receiver.
   * @param notificationType Defines the type of event.
   * @param nTrees The number of trees currently built.
   * @param nPolygons The number of polygons currently enumerated.
   * @param duration The time elapsed since the beginning of the process.
   */
  virtual void notify(NotificationType notificationType,
                      num_type nTrees, count_type nPolygons,
                      std::chrono::nanoseconds duration) = 0;
};



/**
 * @brief The MT_PolygonEnumerator class is a multi-thread version of PolygonEnumerator.
 */
template < typename PolyMeshType >
class MT_PolygonEnumerator : public PolygonEnumerator<PolyMeshType> {
public:
  typedef pl::PolygonEnumerator<PolyMeshType>                 PolygonEnumerator;
  typedef typename PolygonEnumerator::PolygonTree             PolygonTree;
  typedef typename PolygonEnumerator::CornersType             CornersType;
  typedef vcg::tri::pl::PolygonRegister<PolyMeshType>         PolygonRegister;
  typedef vcg::tri::pl::PatchLearner<PolyMeshType>            PatchLearner;
  typedef std::shared_ptr<MT_PolygonEnumerationNotifiable>    MT_Notifiable_Ptr;
  typedef MT_PolygonEnumerationNotifiable::NotificationType   NotificationType;

  /**
   * @brief MT_PolygonEnumerator Default constructor.
   */
  MT_PolygonEnumerator()
      : PolygonEnumerator(), _nTrees(0), _nPolygons(0), _nThreads(std::thread::hardware_concurrency()),
        _stopThreads(false), _nCompleteThreads(_nThreads.load()), _isExecuting(false), _waitTime(std::chrono::seconds(1)) { }

  /**
   * @brief MT_PolygonEnumerator Constructor initializing the size of the forest.
   * @param size The number of trees in the forest.
   */
  MT_PolygonEnumerator(const num_type size)
      : PolygonEnumerator(size), _nTrees(0), _nPolygons(0), _nThreads(std::thread::hardware_concurrency()),
        _stopThreads(false), _nCompleteThreads(_nThreads.load()), _isExecuting(false), _waitTime(std::chrono::seconds(1)) { }

  /**
   * @brief MT_PolygonEnumerator Copy constructor.
   * @param mt_enumerator The MT_PolygonEnumerator to copy.
   */
  MT_PolygonEnumerator(const MT_PolygonEnumerator &mt_enumerator) {
    *this = mt_enumerator;
  }

  /**
   * @brief MT_PolygonEnumerator Move constructor.
   * @param mt_enumerator The MT_PolygonEnumerator to move.
   */
  MT_PolygonEnumerator(MT_PolygonEnumerator&& mt_enumerator) {
    *this = std::move(mt_enumerator);
  }

  /**
   * @brief MT_PolygonEnumerator Destructor.
   */
  ~MT_PolygonEnumerator() {
    terminateAndWaitForThreads();
  }

  /**
   * @brief operator = Copy assignment operator.
   * @param enumerator The MT_PolygonEnumerator to copy.
   * @return A reference to this.
   */
  MT_PolygonEnumerator & operator = (const MT_PolygonEnumerator &mt_enumerator) {
    if (&mt_enumerator != this) {
      mt_enumerator.waitForThreads();
      terminateAndWaitForThreads();
      std::lock_guard<std::mutex> lockR(mt_enumerator._mutex);
      std::lock_guard<std::mutex> lockL(_mutex);
      PolygonEnumerator::operator =(mt_enumerator);
      _faceQueue = mt_enumerator._faceQueue;
      _nTrees = mt_enumerator._nTrees.load();
      _nPolygons = mt_enumerator._nPolygons.load();
      _nThreads = mt_enumerator._nThreads.load();
      _stopThreads = true;
      _nCompleteThreads = _nThreads.load();
      _isExecuting = false;
      _waitTime = mt_enumerator._waitTime;
      _notifiableSet = mt_enumerator._notifiableSet;
    }
    return *this;
  }

  /**
   * @brief operator = Move assignment operator.
   * @param mt_enumerator The MT_PolygonEnumerator to move.
   * @return A reference to this.
   */
  MT_PolygonEnumerator & operator = (MT_PolygonEnumerator&& mt_enumerator) {
    if (&mt_enumerator != this) {
      mt_enumerator.waitForThreads();
      terminateAndWaitForThreads();
      std::lock_guard<std::mutex> lockR(_mutex);
      std::lock_guard<std::mutex> lockL(_mutex);
      PolygonEnumerator::operator =(std::move(mt_enumerator));
      _faceQueue = std::move(mt_enumerator._faceQueue);
      _nTrees = mt_enumerator._nTrees.load();
      _nPolygons = mt_enumerator._nPolygons.load();
      _nThreads = mt_enumerator._nThreads.load();
      _stopThreads = true;
      _nCompleteThreads = _nThreads.load();
      _isExecuting = false;
      _waitTime = mt_enumerator._waitTime;
      _notifiableSet = std::move(mt_enumerator._notifiableSet);
      mt_enumerator._nTrees = 0;
      mt_enumerator._nPolygons = 0;
      mt_enumerator._nThreads = std::thread::hardware_concurrency();
      mt_enumerator._stopThreads = true;
      mt_enumerator._nCompleteThreads = mt_enumerator._nThreads.load();
      mt_enumerator._waitTime = std::chrono::seconds(1);
      mt_enumerator._notifiableSet.clear();
    }
    return *this;
  }

  /**
   * @brief size gives the number of PolygonTrees in the forest.
   * @return The size of the forest.
   */
  num_type size() const {
    std::lock_guard<std::mutex> lock(_mutex);
    return _forest.size();
  }

  /**
   * @brief resize changes the size of the forest.
   * @param size The new size.
   */
  void resize(const num_type size) {
    terminateAndWaitForThreads();
    std::lock_guard<std::mutex> lock(_mutex);
    _nTrees = 0;
    _nPolygons = 0;
    _forest.resize(size);
  }

  /**
   * @brief operator [] gives access to the tree at index.
   * @note This is thread-safe but the returned object is not. Use carefully.
   * @param index The index of the tree.
   * @return A reference to the tree at index.
   */
  PolygonTree & operator [] (const num_type index) {
    static PolygonTree static_tree;
    std::lock_guard<std::mutex> lock(_mutex);
    if (isExecuting())
      return static_tree;
    return _forest[index];
  }

  /**
   * @brief operator[] gives access to the tree at index.
   * @param index The index of the tree.
   * @return A const reference to the tree at index.
   */
  const PolygonTree & operator [] (const num_type index) const {
    static const PolygonTree static_tree_const;
    std::lock_guard<std::mutex> lock(_mutex);
    if (isExecuting())
      return static_tree_const;
    return _forest[index];
  }

  /**
   * @brief numberOfPolygons gives the total number of enumerated polygons.
   * @return The total number of enumerated polygons.
   */
  count_type numberOfPolygons() const {
    return _nPolygons;
  }

  /**
   * @brief numberOfTrees gives the number of built PolygonTrees.
   * @return The number of built PolygonTrees.
   */
  num_type numberOfTrees() const {
    return _nTrees;
  }

  /**
   * @brief maxNumberOfSides gives the maximum number of sides of any enumerated polygon over the whole mesh.
   * @return The maximum number of sides of any enumerated polygon.
   */
  n_corners_type maxNumberOfSides() const {
    std::lock_guard<std::mutex> lock(_mutex);
    if (isExecuting())
      return 0;
    return _forest.maxNumberOfSides();
  }

  /**
   * @brief maxNumberOfSingularities gives the maximum number of singularities of any enumerated polygon over the whole mesh.
   * @return The maximum number of singularities of any enumerated polygon.
   */
  n_corners_type maxNumberOfSingularities() const {
    std::lock_guard<std::mutex> lock(_mutex);
    if (isExecuting())
      return 0;
    return _forest.maxNumberOfSingularities();
  }

  /**
   * @brief isExecuting
   * @return
   */
  bool isExecuting() const {
    return _isExecuting;
  }

  /**
   * @brief numberOfThreads return the number of threads currently sharing the computation.
   * @return The current number of threads.
   */
  u_small_int_type numberOfThreads() const {
    return _nThreads;
  }

  /**
   * @brief setNumberOfThreads sets how many threads will share the computation. Each thread will compute
   * at most size() / nThreads PolygonTrees.
   * @param nThreads The number of threads to set.
   * @return true if the number can be set, false otherwise (i.e. because isExecuting() is true).
   */
  bool setNumberOfThreads(const u_small_int_type nThreads = std::thread::hardware_concurrency()) {
    std::lock_guard<std::mutex> lock(_mutex);
    if (isExecuting())
      return false;
    if (nThreads > 0)
      _nThreads = nThreads;
    else
      _nThreads = 1;
    return true;
  }

  /**
   * @brief waitTime gives the current time to wait before each notification.
   * @return The current time to wait before each notification.
   */
  std::chrono::nanoseconds waitTime() const {
    std::lock_guard<std::mutex> lock(_mutex);
    return _waitTime;
  }

  /**
   * @brief setWaitTime sets the time to wait before each notification.
   * @param duration The new time to wait before each notification.
   */
  template < typename Rep, typename Period >
  void setWaitTime(const std::chrono::duration<Rep,Period> &duration) {
    std::lock_guard<std::mutex> lock(_mutex);
    if (!isExecuting())
      _waitTime = std::chrono::duration_cast<std::chrono::nanoseconds>(duration);
  }

  /**
   * @brief registerNotifiable registers a MT_PolygonEnumerationNotifiable object to receive
   * notifications.
   * @note It can not be registered during execution.
   * @param notifiable_ptr The MT_Notifiable_Ptr object to register.
   * @return true if successfully registered, false otherwise.
   */
  bool registerNotifiable(MT_Notifiable_Ptr &notifiable_ptr) {
    std::lock_guard<std::mutex> lock(_mutex);
    if (isExecuting())
      return false;
    MT_Notifiable_Set_Iterator it = _notifiableSet.find(notifiable_ptr);
    if (it == _notifiableSet.end())
      _notifiableSet.insert(notifiable_ptr);
    return true;
  }

  /**
   * @brief deregisterNotifiable deregisters a MT_PolygonEnumerationNotifiable object to stop
   * receiving notifications.
   * @note It can not be deregistered during execution.
   * @param notifiable_ptr The MT_Notifiable_Ptr object to deregister.
   * @return true if successfully deregistered, false otherwise.
   */
  bool deregisterNotifiable(MT_Notifiable_Ptr &notifiable_ptr) {
    std::lock_guard<std::mutex> lock(_mutex);
    if (isExecuting())
      return false;
    MT_Notifiable_Set_Iterator it = _notifiableSet.find(notifiable_ptr);
    if (it != _notifiableSet.end())
      _notifiableSet.erase(it);
    return true;
  }

  /**
   * @brief enumerate finds all the polygons in mesh with at most maxNFaces.
   * @note This is a multi-threaded version: it launches _nThreads threads and soon returns. To know when it has finished
   * to enumerate call isExecuting. A better solution is to call registerNotifiable() before this.
   * @param mesh The input mesh.
   * @param polygonRegister The register which filters out the duplicated polygons.
   * @param patchLearner The acquisition engine.
   * @param minNFaces The minimum number of per-polygon faces for counting and template storing.
   * @param maxNFaces The maximum number of faces each convex polygon is allowed to have.
   * @param cornersType The type of corners allowed.
   * @param bindPolygon true if the template patch must be linked to the polygon whre it's been found, false otherwise.
   * @param collapsePolygon true if the polygons must be simplified before acquisition, false otherwise.
   * @param keepInMemory true if all the enumerated polygons must be kept in memory, false otherwise.
   */
  void enumerate(PolyMeshType &mesh,
                 PolygonRegister *polygonRegister,
                 PatchLearner &patchLearner,
                 const num_type minNFaces = 0,
                 const num_type maxNFaces = std::numeric_limits<num_type>::max(),
                 const CornersType cornersType = CornersType::CONVEX,
                 const bool bindPolygon = false,
                 const bool collapsePolygon = true,
                 const bool keepInMemory = false) {
    terminateAndWaitForThreads();
    std::lock_guard<std::mutex> lock(_mutex);

    // check if the forest size is the same
    if (_forest.size() != (num_type)mesh.face.size())
      _forest.resize(mesh.face.size());

    // set executing flag
    _isExecuting = true;

    // enumerate
    _threadsManager = std::thread(&MT_PolygonEnumerator::manageThreads, this,
                                  std::ref(mesh), polygonRegister, std::ref(patchLearner), minNFaces,
                                  maxNFaces, cornersType, bindPolygon, collapsePolygon, keepInMemory);
  }

  /**
   * @brief terminateAndWaitForThreads forces to stop any threads (by setting a flag) and then waits for their termination.
   */
  void terminateAndWaitForThreads() {
    _stopThreads = true;
    waitForThreads();
  }

  /**
   * @brief waitForThreads guarantees that all threads finish before doing anything else.
   */
  void waitForThreads() {
    if (_threadsManager.joinable())
      _threadsManager.join();
    if (isExecuting()) {
      std::cout << "Warning. Threads are complete, but flag is still true." << std::endl;
      _isExecuting = false;
    }
  }

private:
  /**
   * @brief manageThreads spans the enumeration process across _nThreads threads and manages their progresses
   * and completion notifications.
   * @param mesh The input mesh.
   * @param polygonRegister The register which filters out the duplicated polygons.
   * @param patchLearner The acquisition engine.
   * @param minNFaces The minimum number of per-polygon faces for counting and template storing.
   * @param maxNFaces The maximum number of faces each polygon is allowed to have.
   * @param cornersType The type of corners allowed.
   * @param bindPolygon true if the template patch must be linked to the polygon whre it's been found, false otherwise.
   * @param collapsePolygon true if the polygons must be simplified before acquisition, false otherwise.
   * @param keepInMemory true if all the enumerated polygons must be kept in memory, false otherwise.
   */
  void manageThreads(PolyMeshType &mesh,
                     PolygonRegister *polygonRegister,
                     PatchLearner &patchLearner,
                     const num_type minNFaces,
                     num_type maxNFaces,
                     const CornersType cornersType,
                     const bool bindPolygon,
                     const bool collapsePolygon,
                     const bool keepInMemory) {
    std::mutex managerMutex;
    std::unique_lock<std::mutex> lock(managerMutex);
    std::chrono::steady_clock::time_point currentTime;
    std::chrono::nanoseconds elapsedTime;
    // initial time point
    std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();

    // initialize
    _nTrees = 0;
    _nPolygons = 0;
    _nCompleteThreads = 0;
    _stopThreads = false;

    // initialize the face queue
    num_type _size = _forest.size();
    _faceQueue.init(_size);

    // throw nThreads threads
    std::list<std::thread> threads;
    num_type range = _size;
    if (_nThreads < range) {
      range += _size % _nThreads;
      range /= _nThreads;
    }
    for (num_type t = 0; t < _size; t += range)
      threads.push_back(std::thread(&MT_PolygonEnumerator::enumerationWorker, this, std::ref(mesh), polygonRegister, std::ref(patchLearner),
                                    minNFaces, maxNFaces, cornersType, bindPolygon, collapsePolygon, keepInMemory));

    // wait and notify progresses
    if (!_notifiableSet.empty()) {
      while (_nCompleteThreads < threads.size()) {
        _threadsCV.wait_for(lock, _waitTime);

        // current time
        currentTime = std::chrono::steady_clock::now();

        // elapsed time
        elapsedTime = std::chrono::duration_cast<std::chrono::nanoseconds>(currentTime - startTime);

        // notify progress
        for (MT_Notifiable_Ptr notifiable_ptr : _notifiableSet)
          notifiable_ptr->notify(NotificationType::TREE_BUILT, _nTrees, _nPolygons, elapsedTime);
      }
    }

    // wait for threads to finish and delete them
    for (std::thread & th : threads)
      th.join();
    threads.clear();
    _nCompleteThreads = _nThreads.load();

    // flush
    _nCompleteThreads = _nThreads - 1;
    threads.push_back(std::thread(&MT_PolygonEnumerator::finalizerWorker, this,
                                  polygonRegister, std::ref(patchLearner)));

    // wait and notify progresses
    if (!_notifiableSet.empty()) {
      while (_nCompleteThreads < _nThreads) {
        _threadsCV.wait_for(lock, _waitTime);

        // current time
        currentTime = std::chrono::steady_clock::now();

        // elapsed time
        elapsedTime = std::chrono::duration_cast<std::chrono::nanoseconds>(currentTime - startTime);

        // notify progress
        for (MT_Notifiable_Ptr notifiable_ptr : _notifiableSet)
          notifiable_ptr->notify(NotificationType::FLUSH_DB, _nTrees, _nPolygons, elapsedTime);
      }
    }

    // wait for threads to finish and delete them
    for (std::thread & th : threads)
      th.join();
    threads.clear();

    // current time
    currentTime = std::chrono::steady_clock::now();

    // elapsed time
    elapsedTime = std::chrono::duration_cast<std::chrono::nanoseconds>(currentTime - startTime);

    // end executing
    _isExecuting = false;

    // notify completion
    for (MT_Notifiable_Ptr notifiable_ptr : _notifiableSet)
      notifiable_ptr->notify(NotificationType::COMPLETE, _nTrees, _nPolygons, elapsedTime);

    // notify
    _threadsCV.notify_all();
  }

  /**
   * @brief enumerationWorker builds PolygonTree from fromFace to toFace indices. Call this in a new thread.
   * @param mesh The input mesh. It is copied before the enumeration, and after that the tree mesh reference is changed back to mesh.
   * @param polygonRegister The register which avoids duplicate polygons.
   * @param patchLearner The patch acquisitor.
   * @param minNFaces The minimum number of per-polygon faces for counting and template storing.
   * @param maxNFaces The maximum number of per-polygon faces.
   * @param cornersType The type of corners allowed.
   * @param bindPolygon true if the template patch must be linked to the polygon whre it's been found, false otherwise.
   * @param collapsePolygon true if the polygon must be simplified before acquisition, false otherwise.
   * @param keepInMemory true if the nodes must be kept in memory, false otherwise.
   */
  void enumerationWorker(PolyMeshType &mesh,
                         PolygonRegister *polygonRegister,
                         PatchLearner &patchLearner,
                         const num_type minNFaces,
                         num_type maxNFaces,
                         const CornersType cornersType,
                         const bool bindPolygon,
                         const bool collapsePolygon,
                         const bool keepInMemory) {
    // make a copy of the mesh
    PolyMeshType meshCopy;
    vcg::tri::Append<PolyMeshType,PolyMeshType>::MeshCopy(meshCopy, mesh);
    vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(meshCopy);

    // enumerate
    for (num_type f = _faceQueue.front(); f >= 0 && !_stopThreads; f = _faceQueue.front()) {
      // build the tree
      _forest.buildPolygonTreeAt(meshCopy, f, polygonRegister, patchLearner, minNFaces, maxNFaces,
                                 cornersType, bindPolygon, collapsePolygon, keepInMemory);
      // change the reference mesh from meshCopy to mesh
      _forest[f].changeReferencedMesh(meshCopy, mesh);
      // update the total number of enumerated convex polygons
      _nPolygons += _forest[f].numberOfPolygons();
      // increment the number of built trees
      ++_nTrees;
      // notify for progress
      _threadsCV.notify_all();
    }
    // notify the manager that this worker thread has finished
    ++_nCompleteThreads;
    _threadsCV.notify_all();
  }

  /**
   * @brief finalizerWorker flushes all the learnt patches into the DB and finalizes it. Call this in another thread.
   * @param polygonRegister The register which avoids duplicate polygons.
   * @param patchLearner The patch acquisitor.
   */
  void finalizerWorker(PolygonRegister *polygonRegister,
                       PatchLearner &patchLearner) {
    patchLearner.finalize();
    if (polygonRegister)
      polygonRegister->finalize();
    // notify the manager that this worker thread has finished
    ++_nCompleteThreads;
    _threadsCV.notify_all();
  }

  /**
   * @brief The FaceQueue class is a thread-safe queue of indices.
   */
  class FaceQueue {
  public:
    /**
     * @brief FaceQueue Default constructor.
     */
    FaceQueue() { }


    /**
     * @brief FaceQueue Constructor initializer.
     * @param size The initial size.
     */
    FaceQueue(const num_type size) {
      init(size);
    }

    /**
     * @brief FaceQueue Copy constructor.
     * @param fq The FaceQueue to copy.
     */
    FaceQueue(const FaceQueue &fq) {
      *this = fq;
    }

    /**
     * @brief FaceQueue Move constructor.
     * @param fq The FaceQueue to move.
     */
    FaceQueue(FaceQueue&& fq) {
      *this = fq;
    }

    /**
     * @brief operator = Copy assignment.
     * @param fq The FaceQueue to copy.
     * @return A reference to this.
     */
    FaceQueue & operator = (const FaceQueue &fq) {
      if (this != &fq) {
        std::lock_guard<std::mutex> lock_s(_mutex);
        std::lock_guard<std::mutex> lock_r(fq._mutex);
        _indices = fq._indices;
      }
      return *this;
    }

    /**
     * @brief operator = Move assignment.
     * @param fq The FaceQueue to move.
     * @return A reference to this.
     */
    FaceQueue & operator = (FaceQueue&& fq) {
      if (this != &fq) {
        std::lock_guard<std::mutex> lock_s(_mutex);
        std::lock_guard<std::mutex> lock_r(fq._mutex);
        _indices = fq._indices;
      }
      return *this;
    }

    /**
     * @brief clear empties the queue.
     */
    void clear() {
      std::lock_guard<std::mutex> lock(_mutex);
      while(!_indices.empty())
        _indices.pop();
    }

    /**
     * @brief init sets a new size for the queue, initializing it.
     * @param size The new size.
     */
    void init(const num_type size) {
      std::lock_guard<std::mutex> lock(_mutex);
      while((num_type)_indices.size() > size)
        _indices.pop();
      for (num_type i = (num_type)_indices.size(); i < size; ++i)
        _indices.push(i);
    }

    /**
     * @brief front gives and pops the index at the top of the queue.
     * @return The top index in the queue.
     */
    num_type front() {
      std::lock_guard<std::mutex> lock(_mutex);
      if (_indices.empty())
        return -1;
      num_type i = _indices.front();
      _indices.pop();
      return i;
    }

  private:
    std::queue<num_type>  _indices;     ///< The queue of indices of faces.
    mutable std::mutex    _mutex;       ///< Mutual exclusion.
  };

  using PolygonEnumerator::_forest;
  FaceQueue                         _faceQueue;                 ///< Thread-safe queue of face indices for the enumeration workers.
  std::atomic<num_type>             _nTrees;                    ///< Current number of built trees.
  std::atomic<count_type>           _nPolygons;                 ///< Current number of polygons.
  std::atomic<u_small_int_type>     _nThreads;                  ///< Number of threads sharing the enumeration process.
  std::atomic<bool>                 _stopThreads;               ///< Flag used to stop threads before they finish.
  std::atomic<u_small_int_type>     _nCompleteThreads;          ///< Number of threads which have finished computing.
  std::atomic<bool>                 _isExecuting;               ///< Computation is ongoing or complete.
  std::chrono::nanoseconds          _waitTime;                  ///< Time to wait before notifying.
  mutable std::mutex                _mutex;                     ///< Mutual exclusion.
  std::thread                       _threadsManager;            ///< Main thread which manages computing threads and emitting notifications.
  std::condition_variable           _threadsCV;                 ///< Condition variable used to wait the threads sharing the enumeration process.

  using MT_Notifiable_Set           = std::unordered_set<MT_Notifiable_Ptr>;    ///< Type alias.
  using MT_Notifiable_Set_Iterator  = MT_Notifiable_Set::iterator;              ///< Type alias.
  MT_Notifiable_Set                 _notifiableSet;             ///< Set of MT_Notifiable objects to notify.
};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // MULTI_THREAD_POLYGON_ENUMERATION_H
