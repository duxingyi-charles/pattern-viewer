#ifndef MULTI_THREAD_PATCH_PRIORITY_LOOKUP_H
#define MULTI_THREAD_PATCH_PRIORITY_LOOKUP_H

#include <atomic>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <algorithm>
#include "../patch_priority_lookup.h"
#include "thread_safe_predicate.h"
#include "time_support.h"

///dxy test
//#include "db_test.h"

///

namespace vcg {
    namespace tri {
        
        // namespace pl: patch learning
        namespace pl {
            
            /**
             * @brief The MT_PatchLookupNotifiable class provides an interface for any class inheriting this one which want
             * to receive notifications on the ongoing process and/or its completion. See MT_PatchPriorityLookup.
             */
            class MT_PatchPriorityLookupNotifiable {
            public:
                /**
                 * @brief The NotificationType enum defines the type of event to be sent/received.
                 */
                enum NotificationType {
                    EMPTY,                  ///< Empty notification.
                    STRUCTURE_FILTER,       ///< Notification on structure filtering.
                    FLOW_FILTER,            ///< Notification on flow-based filtering.
                    FLUSHING,               ///< Notification on flushing the patch array.
                    ILP_FILTER,             ///< Notification on ILP and uniqueness filtering.
                    SORT,                   ///< Notification on sorting.
                    COMPLETE                ///< Notification on process complete.
                };
                
                /**
                 * @brief notify sends a notification to the receiver.
                 * @param notificationType Defines the type of event.
                 * @param nTemplatePatches The number of template patches that passed the filters.
                 * @param nPatches The number of expanded patches.
                 * @param duration The time elapsed since the beginning of the process.
                 */
                virtual void notify(NotificationType notificationType,
                                    count_type nTemplatePatches, count_type nPatches,
                                    std::chrono::nanoseconds duration) = 0;
            };
            
            
            
            /**
             * @brief The MT_PatchPriorityLookup class is a multi-thread version of PatchPriorityLookup.
             */
            template < typename PolyMeshType >
            class MT_PatchPriorityLookup : public PatchPriorityLookup<PolyMeshType> {
            private:
                /**
                 * @brief The JobType enum defines the job type to do in other threads.
                 */
                enum JobType {
                    NONE,               ///< Nothing to do.
                    FILTER_STRUCTURE,   ///< Filter by perimeter and singularities and concavities.
                    FILTER_FLOWS,       ///< Filter (further) by flows.
                    FILTER_ILP,         ///< Filter (further) by solving ILPs.
                    SORT                ///< Work type for sorting patches.
                };
                
            public:
                typedef pl::PatchPriorityLookup<PolyMeshType>                       PatchPriorityLookup;
                typedef typename PatchPriorityLookup::PatchDBServer                 PatchDBServer;
                typedef typename PatchPriorityLookup::PatchDBServerPointer          PatchDBServerPointer;
                typedef typename PatchPriorityLookup::Patch                         Patch;
                typedef typename PatchPriorityLookup::PosType                       PosType;
                typedef typename PatchPriorityLookup::FacePointer                   FacePointer;
                typedef typename PatchPriorityLookup::VertexPointer                 VertexPointer;
                typedef typename PatchPriorityLookup::PatchContainer                PatchContainer;
                typedef typename PatchPriorityLookup::CompareCriteria               CompareCriteria;
                typedef typename PatchPriorityLookup::PolychordSupport              PolychordSupport;
                typedef typename PatchPriorityLookup::FlowMap                       FlowMap;
                typedef typename PatchPriorityLookup::FlowMapConstIterator          FlowMapConstIterator;
                typedef typename PatchPriorityLookup::PatchGeometrizator            PatchGeometrizator;
                typedef typename PatchPriorityLookup::PatchGeometrizatorPointer     PatchGeometrizatorPointer;
                typedef std::shared_ptr<MT_PatchPriorityLookupNotifiable>           MT_Notifiable_Ptr;
                typedef MT_PatchPriorityLookupNotifiable::NotificationType          MT_NotificationType;
                
                /**
                 * @brief MT_PatchPriorityLookup Default constructor.
                 */
                MT_PatchPriorityLookup()
                : PatchPriorityLookup(), _nThreads(std::thread::hardware_concurrency()), _nCompleteThreads(_nThreads.load()), _isExecuting(false),
                _stopped(false), _waitTime(std::chrono::seconds(1)), _printMsg(false), _nCorners(0), _nPatches(0), _nTemplatePatches(0) { }
                
                /**
                 * @brief MT_PatchPriorityLookup Constructor with a DB server parameter.
                 * @param patchDBServer A pointer to the DB server.
                 */
                MT_PatchPriorityLookup(const PatchDBServerPointer &patchDBServer)
                : PatchPriorityLookup(patchDBServer), _nThreads(std::thread::hardware_concurrency()), _nCompleteThreads(_nThreads.load()), _isExecuting(false),
                _stopped(false), _waitTime(std::chrono::seconds(1)), _printMsg(false), _nCorners(0), _nTemplatePatches(0), _nPatches(0) { }
                
                /**
                 * @brief MT_PatchPriorityLookup Copy constructor.
                 * @param mt_lookup The MT_PatchPriorityLookup to copy.
                 */
                MT_PatchPriorityLookup(const MT_PatchPriorityLookup &mt_lookup) {
                    *this = mt_lookup;
                }
                
                /**
                 * @brief MT_PatchPriorityLookup Move constructor.
                 * @param mt_lookup The MT_PatchPriorityLookup to move.
                 */
                MT_PatchPriorityLookup(MT_PatchPriorityLookup&& mt_lookup) {
                    *this = mt_lookup;
                }
                
                /**
                 * @brief MT_PatchPriorityLookup Destructor.
                 */
                ~MT_PatchPriorityLookup() {
                    clear();
                }
                
                /**
                 * @brief operator = Copy assignment operator.
                 * @param mt_lookup The MT_PatchPriorityLookup to copy.
                 * @return A reference to this.
                 */
                MT_PatchPriorityLookup & operator = (const MT_PatchPriorityLookup &mt_lookup) {
                    if (&mt_lookup != this) {
                        mt_lookup.waitForThreads();
                        clear();
                        std::lock_guard<std::mutex> lockR(mt_lookup._mutex);
                        std::lock_guard<std::mutex> lockL(_mutex);
                        PatchPriorityLookup::operator =(mt_lookup);
                        _predicate = mt_lookup._predicate;
                        _nThreads = mt_lookup._nThreads.load();
                        _nCompleteThreads = _nThreads.load();
                        _isExecuting = false;
                        _stopped = mt_lookup._stopped.load();
                        _waitTime = mt_lookup._waitTime;
                        _printMsg = mt_lookup._printMsg.load();
                        _nCorners = mt_lookup._nCorners;
                        _boundaryLength = mt_lookup._boundaryLength;
                        _flows = mt_lookup._flows;
                        _geometrizator = mt_lookup._geometrizator;
                        _perimeterPatches = mt_lookup._perimeterPatches;
                        _flowPatches = mt_lookup._flowPatches;
                        _nPatches = mt_lookup._nPatches.load();
                        _nTemplatePatches = mt_lookup._nTemplatePatches.load();
                        _notifiableSet = mt_lookup._notifiableSet;
                    }
                    return *this;
                }
                
                /**
                 * @brief operator = Move assignment operator.
                 * @param mt_lookup The MT_PatchPriorityLookup to move.
                 * @return A reference to this.
                 */
                MT_PatchPriorityLookup & operator = (MT_PatchPriorityLookup&& mt_lookup) {
                    if (&mt_lookup != this) {
                        mt_lookup.waitForThreads();
                        clear();
                        std::lock_guard<std::mutex> lockR(mt_lookup._mutex);
                        std::lock_guard<std::mutex> lockL(_mutex);
                        PatchPriorityLookup::operator =(std::move(mt_lookup));
                        _predicate = std::move(mt_lookup._predicate);
                        _nThreads = mt_lookup._nThreads.load();
                        _nCompleteThreads = _nThreads.load();
                        _isExecuting = false;
                        _stopped = mt_lookup._stopped.load();
                        _waitTime = mt_lookup._waitTime;
                        _printMsg = mt_lookup._printMsg.load();
                        _nCorners = mt_lookup._nCorners;
                        _boundaryLength = std::move(mt_lookup._boundaryLength);
                        _flows = std::move(mt_lookup._flows);
                        _geometrizator = mt_lookup._geometrizator;
                        _perimeterPatches = std::move(mt_lookup._perimeterPatches);
                        _flowPatches = std::move(mt_lookup._flowPatches);
                        _nPatches = mt_lookup._nPatches.load();
                        _nTemplatePatches = mt_lookup._nTemplatePatches.load();
                        _notifiableSet = std::move(mt_lookup._notifiableSet);
                        mt_lookup._nThreads = std::thread::hardware_concurrency();
                        mt_lookup._nCompleteThreads = mt_lookup._nThreads.load();
                        mt_lookup._waitTime = std::chrono::seconds(1);
                        mt_lookup._nCorners = 0;
                        mt_lookup._boundaryLength.clear();
                        mt_lookup._flows.clear();
                        mt_lookup._perimeterPatches.clear();
                        mt_lookup._flowPatches.clear();
                        mt_lookup._nPatches = 0;
                        mt_lookup._nTemplatePatches = 0;
                        mt_lookup._notifiableSet.clear();
                    }
                    return *this;
                }
                
                /**
                 * @brief setPatchDBServer changes the patch DB server.
                 * @param patchDBServer
                 */
                void setPatchDBServer(const PatchDBServerPointer &patchDBServer) {
                    std::lock_guard<std::mutex> lock(_mutex);
                    if (isExecuting())
                        return;
                    PatchPriorityLookup::setPatchDBServer(patchDBServer);
                    reset();
                }
                
                /**
                 * @brief getSortCriteria gets the current criteria for comparing and sorting patches.
                 * @return The current sequence of criteria.
                 */
                const std::list<CompareCriteria> & getSortCriteria() const {
                    std::lock_guard<std::mutex> lock(_mutex);
                    return PatchPriorityLookup::getSortCriteria();
                }
                
                /**
                 * @brief setSortCriteria sets the criteria for comparing and sorting patches.
                 * @param criteria The sequence of criteria to copy.
                 */
                void setSortCriteria(const std::list<CompareCriteria> &criteria) {
                    std::lock_guard<std::mutex> lock(_mutex);
                    if (isExecuting())
                        return;
                    PatchPriorityLookup::setSortCriteria(criteria);
                }
                
#if __cplusplus >= 201103L || _MSC_VER >= 1800
                /**
                 * @brief setSortCriteria sets the criteria for comparing and sorting patches.
                 * @param criteria The sequence of criteria to move.
                 */
                void setSortCriteria(std::list<CompareCriteria>&& criteria) {
                    std::lock_guard<std::mutex> lock(_mutex);
                    if (isExecuting())
                        return;
                    PatchPriorityLookup::setSortCriteria(std::move(criteria));
                }
#endif
                
                /**
                 * @brief setMinValence See _minValence.
                 * @param minValence
                 */
                void setMinValence(const valence_type minValence = 2) {
                    std::lock_guard<std::mutex> lock(_mutex);
                    if (isExecuting())
                        return;
                    PatchPriorityLookup::setMinValence(minValence);
                    reset();
                }
                
                /**
                 * @brief getMinValence See _minValence.
                 * @return
                 */
                valence_type getMinValence() const {
                    std::lock_guard<std::mutex> lock(_mutex);
                    return PatchPriorityLookup::getMinValence();
                }
                
                /**
                 * @brief setMaxValence See _maxValence.
                 * @param maxValence
                 */
                void setMaxValence(const valence_type maxValence = std::numeric_limits<valence_type>::max()) {
                    std::lock_guard<std::mutex> lock(_mutex);
                    if (isExecuting())
                        return;
                    PatchPriorityLookup::setMaxValence(maxValence);
                    reset();
                }
                
                /**
                 * @brief getMaxValence See _maxValence.
                 * @return
                 */
                virtual valence_type getMaxValence() const {
                    std::lock_guard<std::mutex> lock(_mutex);
                    return PatchPriorityLookup::getMaxValence();
                }
                
                /**
                 * @brief setMaxSigularities See _maxSingularities.
                 * @param maxSingularities
                 */
                void setMaxSingularities(const num_type maxSingularities = std::numeric_limits<num_type>::max()) {
                    std::lock_guard<std::mutex> lock(_mutex);
                    if (isExecuting())
                        return;
                    PatchPriorityLookup::setMaxSingularities(maxSingularities);
                    reset();
                }
                
                /**
                 * @brief getMaxSingularities See _maxSingularities.
                 * @return
                 */
                num_type getMaxSingularities() const {
                    std::lock_guard<std::mutex> lock(_mutex);
                    return PatchPriorityLookup::getMaxSingularities();
                }
                
                /**
                 * @brief setMaxConcavities Patches with more than this number are discarded.
                 * @param maxConcavities
                 */
                void setMaxConcavities(const n_corners_type maxConcavities = std::numeric_limits<n_corners_type>::max()) {
                    std::lock_guard<std::mutex> lock(_mutex);
                    if (isExecuting())
                        return;
                    PatchPriorityLookup::setMaxConcavities(maxConcavities);
                    reset();
                }
                
                /**
                 * @brief getMaxConcavities Patches with more than this number are discarded.
                 * @return
                 */
                n_corners_type getMaxConcavities() const {
                    std::lock_guard<std::mutex> lock(_mutex);
                    return PatchPriorityLookup::getMaxConcavities();
                }
                
                /**
                 * @brief setMinSizeThresholdILP See _minSizeThresholdILP.
                 * @param th
                 */
                void setMinSizeThresholdILP(const num_type th) {
                    std::lock_guard<std::mutex> lock(_mutex);
                    if (isExecuting())
                        return;
                    PatchPriorityLookup::setMinSizeThresholdILP(th);
                    reset();
                }
                
                /**
                 * @brief getMinSizeThresholdILP See _minSizeThresholdILP.
                 * @return
                 */
                num_type getMinSizeThresholdILP() const {
                    std::lock_guard<std::mutex> lock(_mutex);
                    return PatchPriorityLookup::getMinSizeThresholdILP();
                }
                
                /**
                 * @brief setFlowConstraintsILPtolerance See _flowConstraintsILPtolerance.
                 * @param tolerance
                 */
                void setFlowConstraintsILPtolerance(const u_small_int_type tolerance) {
                    std::lock_guard<std::mutex> lock(_mutex);
                    if (isExecuting())
                        return;
                    PatchPriorityLookup::setFlowConstraintsILPtolerance(tolerance);
                    reset();
                }
                
                /**
                 * @brief getFlowConstraintsILPtolerance See _flowConstraintsILPtolerance.
                 * @return
                 */
                u_small_int_type getFlowConstraintsILPtolerance() const {
                    std::lock_guard<std::mutex> lock(_mutex);
                    return PatchPriorityLookup::getFlowConstraintsILPtolerance();
                }
                
                /**
                 * @brief numberOfPatches gives the number of patches satisfying the given boundary constraints.
                 * @return The number of patches.
                 */
                count_type numberOfPatches() const {
                    return _nPatches;
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
                 * @brief setNumberOfThreads sets how many threads will share the computation.
                 * @param nThreads The number of threads to set.
                 * @return true if sucessfully set, false otherwise (when computing it is not possible to change the number).
                 */
                bool setNumberOfThreads(const u_small_int_type nThreads = std::thread::hardware_concurrency()) {
                    std::lock_guard<std::mutex> lock(_mutex);
                    if (isExecuting())
                        return false;
                    if (nThreads > 0)
                        _nThreads = nThreads;
                    else
                        _nThreads = 1;
                    _nCompleteThreads = _nThreads.load();
                    return true;
                }
                
                /**
                 * @brief waitTime gives the current time to wait for notifications.
                 * @return The current time allowed for processing.
                 */
                std::chrono::nanoseconds waitTime() const {
                    std::lock_guard<std::mutex> lock(_mutex);
                    return _waitTime;
                }
                
                /**
                 * @brief setWaitTime sets the time to wait for notifications.
                 * @param duration The new time allowed for processing.
                 * @return true if sucessfully set, false otherwise (when computing it is not possible to change the duration).
                 */
                template < typename Rep, typename Period >
                bool setWaitTime(const std::chrono::duration<Rep,Period> &duration) {
                    std::lock_guard<std::mutex> lock(_mutex);
                    if (isExecuting())
                        return false;
                    _waitTime = std::chrono::duration_cast<std::chrono::nanoseconds>(duration);
                    return true;
                }
                
                /**
                 * @brief showStatusMessages sets the flag for printing messages abour the ongoing lookup process.
                 * @param show true to print messages, false otherwise.
                 */
                void showStatusMessages(const bool show = false) {
                    _printMsg = show;
                }
                
                /**
                 * @brief registerNotifiable registers a MT_PatchPriorityLookupNotifiable object to receive
                 * notifications.
                 * @note It can not be registered during execution.
                 * @param notifiable_ptr The MT_Notifiable_Ptr object to register.
                 * @return true if successfully registered, false otherwise.
                 */
                bool registerNotifiable(const MT_Notifiable_Ptr &notifiable_ptr) {
                    std::lock_guard<std::mutex> lock(_mutex);
                    if (isExecuting())
                        return false;
                    MT_Notifiable_Set_Iterator it = _notifiableSet.find(notifiable_ptr);
                    if (it == _notifiableSet.end())
                        _notifiableSet.insert(notifiable_ptr);
                    return true;
                }
                
                /**
                 * @brief deregisterNotifiable deregisters a MT_PatchPriorityLookupNotifiable object to stop
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
                 * @brief findPatchesFromBoundaryLength makes a new query to find Patches with nCorners and boundaryLength and flows.
                 * @param nCorners The number of corners (it may be 0).
                 * @param boundaryLength The boundary's sides' length.
                 * @param flows A list of flows to check (patches without such flows, if any, are discarded).
                 * @param geometrizator The responsible of assigning geometry to patches.
                 * @param predicate A predicate whose evaluation allows for the continuation of the extraction of patch.
                 */
                void findPatchesFromBoundaryLength(const n_corners_type nCorners, const std::vector<num_type> &boundaryLength,
                                                   const FlowMap &flows = FlowMap(),
                                                   const PatchGeometrizatorPointer &geometrizator = PatchGeometrizatorPointer(0),
                                                   const PredicatePointer &predicate = PredicatePointer(0)) {
                    // terminate and wait any previous computation
                    terminateAndWaitForThreads();
                    std::lock_guard<std::mutex> lock(_mutex);
                    if (_stopped) {
                        reset();
                        _stopped = false;
                    }
                    // set the job
                    JobType jobType = NONE;
                    if (nCorners != _nCorners || boundaryLength != _boundaryLength) {
                        _nCorners = nCorners;
                        _boundaryLength = boundaryLength;
                        _flows = flows;
                        jobType = FILTER_STRUCTURE;
                    } else if (flows != _flows) {
                        _flows = flows;
                        jobType = FILTER_FLOWS;
                    } else {
                        jobType = SORT;
                    }
                    // reset
                    _isExecuting = true;
                    _geometrizator = geometrizator;
                    _predicate.setPredicate(predicate);
                    // check valence range
                    if (_maxValence < _minValence)
                        std::swap(_minValence, _maxValence);
                    // lookup
                    _threadsManager = std::thread(&MT_PatchPriorityLookup::manageThreads, this, jobType);
                }
                
                /**
                 * @brief sortPatches sorts the previously found patches.
                 * @param flows A list of flows to check (patches without such flows, if any, are discarded).
                 * @param geometrizator An object responsible for assigning geometry.
                 * @param predicate A predicate whose evaluation allows for the continuation of the extraction of patch.
                 */
                void sortPatches(const FlowMap &flows,
                                 const PatchGeometrizatorPointer &geometrizator = PatchGeometrizatorPointer(0),
                                 const PredicatePointer &predicate = PredicatePointer(0)) {
                    // wait any previous computation
                    terminateAndWaitForThreads();
                    std::lock_guard<std::mutex> lock(_mutex);
                    // reset
                    _isExecuting = true;
                    _flows = flows;
                    _geometrizator = geometrizator;
                    _predicate.setPredicate(predicate);
                    // lookup
                    _threadsManager = std::thread(&MT_PatchPriorityLookup::manageThreads, this, SORT);
                }
                
                /**
                 * @brief at gives access to the patch at index.
                 * @note Call this method when not executing, otherwise no valid patch is returned.
                 * @param index The index of the patch.
                 * @return A reference to the the requested ConvexPatch. A reference to a NULL Patch is returned
                 * in case of a out-of-bounds index.
                 */
                Patch & at(const count_type index) {
                    // a static object to return in case of error
                    static Patch staticPatch;
                    staticPatch.reset();
                    
                    if (isExecuting())
                        return staticPatch;
                    
                    std::lock_guard<std::mutex> lock(_mutex);
                    // if the template patch at tIndex has not been found, return nothing
                    if (index >= _patches.size())
                        return staticPatch;
                    
                    // else return the requested patches
                    return _patches.at(index);
                }
                
                /**
                 * @brief at gives access to the patch at index.
                 * @note This method can be used during execution. It is useful for dynamic updating on-screen results.
                 * @param index The index of the patch.
                 * @param patch The output Patch. In case of out-of-bounds, it is set to NULL.
                 */
                void at(const count_type index, Patch &patch) const {
                    patch.reset();
                    std::lock_guard<std::mutex> lock(_mutex);
                    if (index < _patches.size())
                        patch = _patches.at(index);
                }
                
                /**
                 * @brief at gives access to the range of patches at [index, index+offset-1].
                 * @param index The starting index of the range.
                 * @param offset The number of patches in the range.
                 * @param patches The output container of patches
                 */
                void at(const count_type index, const count_type offset, PatchContainer &patches) const {
                    patches.clear();
                    std::lock_guard<std::mutex> lock(_mutex);
                    for (count_type i = 0; i < offset && index + i < _patches.size(); ++i)
                        patches.push_back(_patches.at(index + i));
                }
                
                /**
                 * @brief clear empties the buffer of patches.
                 */
                void clear() {
                    terminateAndWaitForThreads();
                    std::lock_guard<std::mutex> lock(_mutex);
                    reset();
                }
                
                /**
                 * @brief sendInterrupt tells threads to interrupt their execution.
                 */
                void sendInterrupt() {
                    _predicate.disable();
                }
                
                /**
                 * @brief stop interrupts the execution.
                 */
                void stop() {
                    terminateAndWaitForThreads();
                    _stopped = true;
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
                typedef typename PatchPriorityLookup::TopologySetType           TopologySetType;
                typedef typename PatchPriorityLookup::PatchIndexRotationType    PatchIndexRotationType;
                typedef typename PatchPriorityLookup::ILPArrayType              ILPArrayType;
                typedef typename PatchPriorityLookup::ILPMatrixType             ILPMatrixType;
                using PatchPriorityLookup::solveILP;
                using PatchPriorityLookup::compareSidesLength;
                using PatchPriorityLookup::isUnique;
                using PatchPriorityLookup::addNewPatch;
                
                /**
                 * @brief terminateAndWaitForThreads forces to stop any threads (by setting a flag) and then waits for their termination.
                 */
                void terminateAndWaitForThreads() {
                    sendInterrupt();
                    waitForThreads();
                }
                
                /**
                 * @brief reset deletes any current state.
                 */
                void reset() {
                    _nCorners = 0;
                    _boundaryLength.clear();
                    _flows.clear();
                    _perimeterPatches.clear();
                    _flowPatches.clear();
                    _nPatches = 0;
                    _nTemplatePatches = 0;
                }
                
                /**
                 * @brief manageThreads spans the finding process across _nThreads threads and manages their progresses
                 * and completion notifications.
                 * @param jobType What job to do.
                 */
                void manageThreads(const JobType jobType) {
                    // some variables
                    std::list<std::list<PatchIndexRotationType>> tmpTemplatePatchLists;
                    typename std::list<PatchIndexRotationType>::const_iterator it, nextIt;
                    std::list<std::thread> threads;
                    std::chrono::steady_clock::time_point currentTime;
                    std::chrono::nanoseconds elapsedTime;
                    count_type tot = 0, range = 0, t = 0;
                    num_type perimeter = std::accumulate(_boundaryLength.begin(), _boundaryLength.end(), 0);
                    std::mutex managerMutex;
                    std::unique_lock<std::mutex> lock(managerMutex);
                    std::string durationsStr;
                    // initial time point
                    std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
                    
                    _nPatches = 0;
                    _predicate.enable();
                    
                    // which job?
                    switch (jobType) {
                            // the job is to look up patches
                        case FILTER_STRUCTURE:
                            tot = _patchDBServer->numberOfCachedPatches(_nCorners);
                            _nTemplatePatches = tot * _boundaryLength.size();
                            if (_printMsg) {
                                std::cout << tot << " initial templates." << std::endl;
                                std::cout << _nTemplatePatches << " initial rotate templates." << std::endl;
                            }
                            
                            // if no rows
                            if (tot == 0) {
                                // than conclude
                                _nCompleteThreads = _nThreads.load();
                                _isExecuting = false;
                                // notify
                                _threadsCV.notify_all();
                                // current time
                                currentTime = std::chrono::steady_clock::now();
                                // elapsed time
                                elapsedTime = std::chrono::duration_cast<std::chrono::nanoseconds>(currentTime - startTime);
                                // notify completion
                                for (MT_Notifiable_Ptr notifiable_ptr : _notifiableSet)
                                    notifiable_ptr->notify(MT_NotificationType::COMPLETE, _nTemplatePatches, _nPatches, elapsedTime);
                                return;
                            }
                            
                            // first filter: perimeter
                            _nCompleteThreads = 0;
                            range = tot;
                            if (_nThreads < range) {
                                range += tot % _nThreads;
                                range /= _nThreads;
                            }
                            for (t = 0; t < tot; t += range) {
                                tmpTemplatePatchLists.push_back(std::list<PatchIndexRotationType>());
                                threads.push_back(std::thread(&MT_PatchPriorityLookup::structureFilterWorker, this,
                                                              std::ref(tmpTemplatePatchLists.back()), perimeter, t, t + range));
                            }
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
                                        notifiable_ptr->notify(MT_NotificationType::STRUCTURE_FILTER, _nTemplatePatches, _nPatches, elapsedTime);
                                }
                            }
                            // wait for threads to finish and delete them
                            for (std::thread & th : threads)
                                th.join();
                            threads.clear();
                            _perimeterPatches.clear();
                            for (std::list<PatchIndexRotationType> & list : tmpTemplatePatchLists)
                                _perimeterPatches.splice(_perimeterPatches.end(), list);
                            tmpTemplatePatchLists.clear();
                            _nCompleteThreads = _nThreads.load();
                            
                            if (_printMsg) {
                                // elapsed time
                                elapsedTime = std::chrono::duration_cast<std::chrono::nanoseconds>(currentTime - startTime);
                                TimeSupport::TimeToString(elapsedTime, durationsStr);
                                std::cout << _nTemplatePatches << " Structure-filtered rotated templates after " << durationsStr << std::endl;
                            }
                            
                            
                        case FILTER_FLOWS:
                            // second filter: flows
                            tot = _perimeterPatches.size();
                            _nTemplatePatches = tot;
                            if (!_flows.empty() && tot > 0) {
                                _nCompleteThreads = 0;
                                range = tot;
                                if (_nThreads <= range) {
                                    range += tot % _nThreads;
                                    range /= _nThreads;
                                }
                                it = _perimeterPatches.begin();
                                nextIt = it;
                                std::advance(nextIt, range);
                                for (t = 0; t < tot; t += range, it = nextIt, std::advance(nextIt, std::min(range, tot - t))) {
                                    tmpTemplatePatchLists.push_back(std::list<PatchIndexRotationType>());
                                    threads.push_back(std::thread(&MT_PatchPriorityLookup::flowFilterWorker, this,
                                                                  std::ref(tmpTemplatePatchLists.back()), it, nextIt));
                                }
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
                                            notifiable_ptr->notify(MT_NotificationType::FLOW_FILTER, _nTemplatePatches, _nPatches, elapsedTime);
                                    }
                                }
                                // wait for threads to finish and delete them
                                for (std::thread & th : threads)
                                    th.join();
                                threads.clear();
                                _flowPatches.clear();
                                for (std::list<PatchIndexRotationType> &list : tmpTemplatePatchLists)
                                    _flowPatches.splice(_flowPatches.end(), list);
                                tmpTemplatePatchLists.clear();
                                _nCompleteThreads = _nThreads.load();
                                
                                if (_printMsg) {
                                    // elapsed time
                                    elapsedTime = std::chrono::duration_cast<std::chrono::nanoseconds>(currentTime - startTime);
                                    TimeSupport::TimeToString(elapsedTime, durationsStr);
                                    std::cout << _nTemplatePatches << " Flow-filtered rotated templates after " << durationsStr << std::endl;
                                }
                            } else {
                                _flowPatches = _perimeterPatches;
                            }
                            
                            
                        case FILTER_ILP:
                            // clear the patch array
                            tot = _patches.size();
                            if (tot > 0) {
                                _nPatches = tot;
                                range = tot;
                                if (_nThreads <= range) {
                                    range += tot % _nThreads;
                                    range /= _nThreads;
                                }
                                _nCompleteThreads = 0;
                                for (t = 0; t < tot; t += range)
                                    threads.push_back(std::thread(&MT_PatchPriorityLookup::resetPatchesWorker, this,
                                                                  t, t + range));
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
                                            notifiable_ptr->notify(MT_NotificationType::FLUSHING, _nTemplatePatches, _nPatches, elapsedTime);
                                    }
                                }
                                // wait for threads to finish and delete them
                                for (std::thread & th : threads)
                                    th.join();
                                threads.clear();
                                _patches.clear();
                                _nPatches = 0;
                            }
                            // last filter: ILP and uniqueness
                            tot = _flowPatches.size();
                            if (tot > 0) {
                                range = tot;
                                if (_nThreads <= range) {
                                    range += tot % _nThreads;
                                    range /= _nThreads;
                                }
                                // ilp filter starts here
                                _nCompleteThreads = 0;
                                it = _flowPatches.begin();
                                nextIt = it;
                                std::advance(nextIt, range);
                                for (t = 0; t < tot; t += range, it = nextIt, std::advance(nextIt, std::min(range, tot - t)))
                                    threads.push_back(std::thread(&MT_PatchPriorityLookup::ilpFilterWorker, this,
                                                                  perimeter, it, nextIt));
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
                                            notifiable_ptr->notify(MT_NotificationType::ILP_FILTER, _nTemplatePatches, _nPatches, elapsedTime);
                                    }
                                }
                                // wait for threads to finish and delete them
                                for (std::thread & th : threads)
                                    th.join();
                                threads.clear();
                                _uniquePatches.clear();
                                _nCompleteThreads = _nThreads.load();
                            }
                            
                            if (_printMsg) {
                                // elapsed time
                                elapsedTime = std::chrono::duration_cast<std::chrono::nanoseconds>(currentTime - startTime);
                                TimeSupport::TimeToString(elapsedTime, durationsStr);
                                std::cout << _nTemplatePatches << " ILP-filtered rotated templates after " << durationsStr << std::endl << std::endl;
                            }
                            break;
                            
                            
                        case SORT:
                            // clear the patch array
                            tot = _patches.size();
                            if (tot > 0) {
                                PatchContainer tmpPatches(std::move(_patches));
                                _patches.clear();
                                range = tot;
                                if (_nThreads <= range) {
                                    range += tot % _nThreads;
                                    range /= _nThreads;
                                }
                                _nCompleteThreads = 0;
                                _nPatches = 0;
                                // thow the threads
                                for (t = 0; t < tot; t += range)
                                    threads.push_back(std::thread(&MT_PatchPriorityLookup::sortWorker, this, std::ref(tmpPatches), t, t + range));
                                // wait and notify progresses
                                if (!_notifiableSet.empty()) {
                                    while (_nCompleteThreads < threads.size()) {
                                        _threadsCV.wait_for(lock, _waitTime);
                                        // current time
                                        currentTime = std::chrono::steady_clock::now();
                                        
                                        // elapsed time
                                        elapsedTime = std::chrono::duration_cast<std::chrono::nanoseconds>(currentTime - startTime);
                                        
                                        // notify sorting
                                        for (MT_Notifiable_Ptr notifiable_ptr : _notifiableSet)
                                            notifiable_ptr->notify(MT_NotificationType::SORT, _nTemplatePatches, _nPatches, elapsedTime);
                                    }
                                }
                                // wait for threads to finish and delete them
                                for (std::thread & th : threads)
                                    th.join();
                                threads.clear();
                            }
                            
                            if (_printMsg) {
                                // elapsed time
                                elapsedTime = std::chrono::duration_cast<std::chrono::nanoseconds>(currentTime - startTime);
                                TimeSupport::TimeToString(elapsedTime, durationsStr);
                                std::cout << _nPatches << " patches sorted after " << durationsStr << std::endl << std::endl;
                            }
                            break;
                            
                        default:
                            break;
                    }
                    
                    // current time
                    currentTime = std::chrono::steady_clock::now();
                    
                    // elapsed time
                    elapsedTime = std::chrono::duration_cast<std::chrono::nanoseconds>(currentTime - startTime);
                    
                    // notify completion
                    for (MT_Notifiable_Ptr notifiable_ptr : _notifiableSet)
                        notifiable_ptr->notify(MT_NotificationType::COMPLETE, _nTemplatePatches, _nPatches, elapsedTime);
                    
                    // end executing
                    _isExecuting = false;
                    
                    // notify
                    _threadsCV.notify_all();
                }
                
                /**
                 * @brief structureFilterWorker lists all the patches (and their rotations) that pass the perimeter and singularities and concavities test.
                 * @param outPatches The list of passed patches.
                 * @param startIndex The patches starting index.
                 * @param endIndex The patches ending index (not included).
                 */
                void structureFilterWorker(std::list<PatchIndexRotationType> &outPatches,
                                           const num_type perimeter,
                                           const count_type startIndex = 0,
                                           const count_type endIndex = std::numeric_limits<count_type>::max()) {
                    outPatches.clear();
                    bool firstCheck = true;
                    if (perimeter < 4)
                        firstCheck = false;
                    else if (perimeter % 2 != 0)
                        firstCheck = false;
                    if (!_patchDBServer)
                        firstCheck = false;
                    if (!firstCheck) {
                        _nTemplatePatches = 0;
                        // complete! notify
                        ++_nCompleteThreads;
                        _threadsCV.notify_all();
                        return;
                    }
                    std::vector< std::pair<valence_type,num_type> > singularityVec;
                    valence_type minValence, maxValence;
                    // for each patch
                    for (count_type index = startIndex; index < endIndex && index < _patchDBServer->numberOfCachedPatches(_nCorners); ++index) {
                        _patchDBServer->at(_nCorners, index).getSingularityVec(singularityVec);
                        if (singularityVec.empty()) {
                            minValence = _minValence;
                            maxValence = _maxValence;
                        } else {
                            minValence = singularityVec.front().first;
                            maxValence = singularityVec.back().first;
                        }
                        if (_patchDBServer->at(_nCorners, index).getNumberOfSingularities() <= _maxSingularities &&
                            minValence >= _minValence && maxValence <= _maxValence &&
                            _patchDBServer->at(_nCorners, index).getNumberOfConcaveCorners() <= _maxConcavities &&
                            _patchDBServer->at(_nCorners, index).getPerimeter() <= perimeter) {
                            // for each rotation, only those patches that pass the filter are added
                            for (n_corners_type nRotations = 0; nRotations < _boundaryLength.size(); ++nRotations) {
                                if (!_predicate()) {
                                    outPatches.clear();
                                    // complete! notify
                                    ++_nCompleteThreads;
                                    _threadsCV.notify_all();
                                    return;
                                }
                                if (compareSidesLength(_boundaryLength, _patchDBServer->at(_nCorners, index).getBoundaryLenVec(), nRotations))
                                    outPatches.push_back(PatchIndexRotationType(index, nRotations));
                                else
                                    --_nTemplatePatches;
                                // notify
                                _threadsCV.notify_all();
                            }
                        } else {
                            _nTemplatePatches -= _boundaryLength.size();
                        }
                    }
                    
                    // complete! notify
                    ++_nCompleteThreads;
                    _threadsCV.notify_one();
                }
                
                /**
                 * @brief flowFilterWorker lists all patches (and their rotations) that pass the flow test.
                 * The flow test is passed if there are the polychords from/to the sides specified in
                 * flows, or if there is any loop when flows specifies it.
                 * @param outPatches The list of passed patches.
                 * @param firstIt The first patch to check.
                 * @param lastIt The last patch to check (not included).
                 * @param flows The list of flows to check against.
                 * @param predicate A predicate whose evaluation allows for the continuation of the process.
                 */
                void flowFilterWorker(std::list<PatchIndexRotationType> &outPatches,
                                      typename std::list<PatchIndexRotationType>::const_iterator firstIt,
                                      typename std::list<PatchIndexRotationType>::const_iterator lastIt) {
                    outPatches.clear();
                    
                    if (!_patchDBServer) {
                        outPatches.clear();
                        _nTemplatePatches = 0;
                        // complete! notify
                        ++_nCompleteThreads;
                        _threadsCV.notify_all();
                        return;
                    }
                    
                    FlowMapConstIterator flowsIt;
                    n_corners_type startSide, endSide, nSides = _boundaryLength.size();
                    var_label_type boundaryID = 0;
                    PolyMeshType mesh;
                    PosType startPos;
                    std::deque<PosType> corners;
                    n_corners_type nc = 0;
                    std::deque<PosType> polychords;
                    std::vector<bool> isFlowPassed(_flows.size(), false);
                    size_t fInd = 0;
                    bool pass = true;
                    
                    // for each patch
                    for(; firstIt != lastIt; ++firstIt) {
                        if (!_predicate()) {
                            outPatches.clear();
                            // complete! notify
                            ++_nCompleteThreads;
                            _threadsCV.notify_all();
                            return;
                        }
                        // get polychords
                        _patchDBServer->at(_nCorners, firstIt->first).reconstructPatch(mesh, startPos);
                        startPos.FlipE();
                        while (!startPos.IsBorder()) {
                            startPos.FlipF();
                            startPos.FlipE();
                        }
                        startPos.FlipV();
                        nc = PolychordSupport::FindPolygonCorners(mesh, corners, false, false, startPos);
                        assert(nc == _nCorners);
                        // for each flow
                        fInd = 0;
                        for (flowsIt = _flows.begin(); flowsIt != _flows.end(); ++flowsIt) {
                            if (!_predicate()) {
                                outPatches.clear();
                                // complete! notify
                                ++_nCompleteThreads;
                                _threadsCV.notify_all();
                                return;
                            }
                            isFlowPassed[fInd] = false;
                            // border-to-border case
                            if (flowsIt->first.startSide < nSides && flowsIt->first.endSide < nSides) {
                                startSide = (flowsIt->first.startSide + firstIt->second) % nSides;
                                endSide = (flowsIt->first.endSide + firstIt->second) % nSides;
                                for (size_t edge = 0; edge < _patchDBServer->at(_nCorners, firstIt->first).getBoundaryIDsVec()[startSide].size(); ++edge) {
                                    boundaryID = _patchDBServer->at(_nCorners, firstIt->first).getBoundaryIDsVec()[startSide][edge];
                                    if ((_patchDBServer->at(_nCorners, firstIt->first).getBoundaryIDsEnds()[boundaryID].first == startSide
                                         && _patchDBServer->at(_nCorners, firstIt->first).getBoundaryIDsEnds()[boundaryID].second == endSide) ||
                                        (_patchDBServer->at(_nCorners, firstIt->first).getBoundaryIDsEnds()[boundaryID].second == startSide
                                         && _patchDBServer->at(_nCorners, firstIt->first).getBoundaryIDsEnds()[boundaryID].first == endSide)) {
                                            isFlowPassed[fInd] = true;
                                            break;
                                        }
                                }
                            } else {    // loop case
                                PolychordSupport::FindPolychordsWithTopology(mesh, corners, flowsIt->first, polychords);
                                if (polychords.size() >= flowsIt->second.size())
                                    isFlowPassed[fInd] = true;
                            }
                            ++fInd;
                        }
                        // if pass, push out
                        pass = true;
                        for (fInd = 0; fInd < isFlowPassed.size(); ++fInd)
                            pass = pass && isFlowPassed[fInd];
                        if (pass)
                            outPatches.push_back(*firstIt);
                        else
                            --_nTemplatePatches;
                        // notify
                        _threadsCV.notify_all();
                    }
                    
                    // complete! notify
                    ++_nCompleteThreads;
                    _threadsCV.notify_all();
                }
                
                /**
                 * @brief resetPatchesWorker resets the patches from startIndex to endIndex (not included).
                 * @param startIndex The first patch to clear.
                 * @param endIndex The last patch to clear (not included).
                 */
                void resetPatchesWorker(const count_type startIndex = 0,
                                        const count_type endIndex = std::numeric_limits<count_type>::max()) {
                    for (count_type index = startIndex; index < endIndex && index < _patches.size(); ++index) {
                        if (!_predicate()) {
                            // complete! notify
                            ++_nCompleteThreads;
                            _threadsCV.notify_all();
                            return;
                        }
                        // clear patch at index
                        _patches.at(index).reset();
                        --_nPatches;
                        _threadsCV.notify_all();
                    }
                    // complete! notify
                    ++_nCompleteThreads;
                    _threadsCV.notify_all();
                }
                
                /**
                 * @brief ilpFilterWorker adds all the patches that pass the ILP filter and that are unique.
                 * @param perimeter It is equal to the sum of the boundary's sides' length.
                 * @param firstIt The first patch to check.
                 * @param lastIt The last patch to check (not included).
                 */
                void ilpFilterWorker(const num_type perimeter,
                                     typename std::list<PatchIndexRotationType>::const_iterator firstIt,
                                     typename std::list<PatchIndexRotationType>::const_iterator lastIt) {
                    if (!_patchDBServer) {
                        _nTemplatePatches = 0;
                        // complete! notify
                        ++_nCompleteThreads;
                        _threadsCV.notify_all();
                        return;
                    }
                    
                    ILPMatrixType solutions;
                    Patch patch;
                    std::vector<num_type> boundaryVarVec;
                    size_t nDiscarded = 0;
                    std::deque<PosType> corners;
                    PosType startPos;
                    n_corners_type nc = 0;
                    std::deque<PosType> polychords;
                    FlowMapConstIterator flowsIt;
                    bool pass = true;
                    for (; firstIt != lastIt; ++firstIt) {
                        if (!_predicate()) {
                            // complete! notify
                            ++_nCompleteThreads;
                            _threadsCV.notify_all();
                            return;
                        }
                        // solve
                        solveILP(_nCorners, _boundaryLength, *firstIt, _flows, solutions);
                        if (!solutions.empty()) {
                            ///dxy test
                            std::cout << "there are " << solutions.size() << " solutions." << std::endl;
                            ///
                            // init patch
                            patch = _patchDBServer->at(_nCorners, firstIt->first);
                            patch.setNumberOfRotations(firstIt->second);
                            ///dxy add: save old LenVec
                            auto oldLenVec = patch.getBoundaryLenVec();
                            ///
                            patch.setBoundaryLenVec(_boundaryLength);
                            patch.setPerimeter(perimeter);
                            boundaryVarVec = patch.getBoundaryVarVec();
                            nDiscarded = solutions.size();
                            // for each solution
                            for (size_t i = 0; i < solutions.size(); ++i) {
                                if (!_predicate()) {
                                    // complete! notify
                                    ++_nCompleteThreads;
                                    _threadsCV.notify_all();
                                    return;
                                }
                                // set current solution
                                for (size_t v = 0; v < boundaryVarVec.size(); ++v)
                                    boundaryVarVec[v] = (num_type)solutions[i][v];
                                patch.setBoundaryVarVec(boundaryVarVec);
                                ///dxy add: recover oldLenVec
                                patch.setBoundaryLenVec(oldLenVec);
                                ///
                                patch.buildExpandedTopologyStr();
                                ///dxy test: save expanded patch
//                                DB_Test<PolyMesh> dbtest;
//                                auto patch_copy = patch;
//                                dbtest.geometrize(patch_copy, 90);
//                                std::string mesh_name = std::to_string(_nCorners)
//                                                        + "_" + std::to_string(firstIt->first)
//                                                        + "_r" + std::to_string(firstIt->second)
//                                                        + "_s" + std::to_string(i) + ".obj";
//                                dbtest.save_mesh(patch_copy.getMesh(), mesh_name);
                                ///
                                ///dxy add: recover newLenVec
                                patch.setBoundaryLenVec(_boundaryLength);
                                ///
                                // check flows
                                pass = true;
                                startPos = patch.getStartCorner();
                                startPos.FlipE();
                                while (!startPos.IsBorder()) {
                                    startPos.FlipF();
                                    startPos.FlipE();
                                }
                                startPos.FlipV();
                                nc = PolychordSupport::FindPolygonCorners(patch.getMesh(), corners, false, false, startPos);  //!!!
                                assert(nc == _nCorners);
                                for (flowsIt = _flows.begin(); flowsIt != _flows.end(); ++flowsIt) {
                                    PolychordSupport::FindPolychordsWithTopology(patch.getMesh(), corners, flowsIt->first, polychords);
                                    if (polychords.size() < flowsIt->second.size()) {
                                        pass = false;
                                        break;
                                    }
                                }
                                // check uniqueness
                                if (pass) {
                                    _mutex.lock();
                                    pass = isUnique(patch.getExpandedTopologyStr());
                                    _mutex.unlock();
                                }
                                // geometrize and check
                                if (pass && _geometrizator)
                                    pass = _geometrizator->geometrizePatch(patch, _flows);
                                // if new, add it
                                if (pass) {
                                    _mutex.lock();
                                    addNewPatch(patch);
                                    _mutex.unlock();
                                    ++_nPatches;
                                } else {
                                    pass = false;
                                    --nDiscarded;
                                }
                                // notify
                                if (pass)
                                    _threadsCV.notify_all();
#ifdef SPLIT_POLYCHORD_LOOPS
                                if (pass) {
                                    flowsIt = _flows.find(FlowTopology());
                                    if (flowsIt != _flows.end() && flowsIt->second.size() == 1) {
                                        Patch splitLoopPatch = patch;
                                        startPos = splitLoopPatch.getStartCorner();
                                        startPos.FlipE();
                                        while (!startPos.IsBorder()) {
                                            startPos.FlipF();
                                            startPos.FlipE();
                                        }
                                        startPos.FlipV();
                                        nc = PolychordSupport::FindPolygonCorners(splitLoopPatch.getMesh(), corners, false, false, startPos);
                                        assert(nc == _nCorners);
                                        PolychordSupport::FindPolychordsWithTopology(splitLoopPatch.getMesh(), corners, FlowTopology(), polychords);
                                        if (polychords.size() == 1) {
                                            std::vector<FacePointer *> facesToUpdate;
                                            std::vector<VertexPointer *> verticesToUpdate;
                                            facesToUpdate.push_back(&splitLoopPatch.getStartCorner().F());
                                            verticesToUpdate.push_back(&splitLoopPatch.getStartCorner().V());
                                            for (size_t i = 0; i < polychords.size(); ++i) {
                                                facesToUpdate.push_back(&polychords[i].F());
                                                verticesToUpdate.push_back(&polychords[i].V());
                                            }
                                            for (size_t i = 0; i < polychords.size(); ++i)
                                                vcg::tri::PolychordCollapse<PolyMeshType>::SplitPolychord(splitLoopPatch.getMesh(), polychords[i],
                                                                                                          (splitLoopPatch.getMesh().FN() > 40 ? 3 : 2),
                                                                                                          facesToUpdate, verticesToUpdate);
                                            std::string expandedTopology;
                                            Patch::BuildTopologyStr(splitLoopPatch.getMesh(), splitLoopPatch.getStartCorner(), expandedTopology);
                                            splitLoopPatch.setExpandedTopologyStr(expandedTopology);
                                            // check uniqueness
                                            if (pass) {
                                                _mutex.lock();
                                                pass = isUnique(patch.getExpandedTopologyStr());
                                                _mutex.unlock();
                                            }
                                            // geometrize and check
                                            if (pass && _geometrizator)
                                                pass = _geometrizator->geometrizePatch(patch, _flows);
                                            // if new, add it
                                            if (pass) {
                                                _mutex.lock();
                                                addNewPatch(patch);
                                                _mutex.unlock();
                                                ++_nPatches;
                                            }
                                            // notify
                                            if (pass)
                                                _threadsCV.notify_all();
                                        }
                                    }
                                }
#endif
                            }
                            if (nDiscarded == 0)
                                --_nTemplatePatches;
                            else // notify
                                _threadsCV.notify_all();
                        } else {
                            --_nTemplatePatches;
                            // notify
                            _threadsCV.notify_all();
                        }
                    }
                    
                    // complete! notify
                    ++_nCompleteThreads;
                    _threadsCV.notify_all();
                }
                
                /**
                 * @brief sortWorker sorts the patches (previously found) by storing them into _orderedPatches.
                 */
                void sortWorker(PatchContainer &tmpPatches,
                                const count_type startIndex = 0,
                                const count_type endIndex = std::numeric_limits<count_type>::max()) {
                    for (count_type index = startIndex; index < endIndex && index < tmpPatches.size(); ++index) {
                        if (!_predicate()) {
                            // complete! notify
                            ++_nCompleteThreads;
                            _threadsCV.notify_all();
                            return;
                        }
                        // compute new geometry
                        if (_geometrizator)
                            _geometrizator->geometrizePatch(tmpPatches[index], _flows);
                        // add in new position
                        _mutex.lock();
                        addNewPatch(std::move(tmpPatches[index]));
                        _mutex.unlock();
                        ++_nPatches;
                        _threadsCV.notify_all();
                    }
                    ++_nCompleteThreads;
                    // notify
                    _threadsCV.notify_all();
                }
                
                using PatchPriorityLookup::_patchDBServer;
                using PatchPriorityLookup::_patches;
                using PatchPriorityLookup::_uniquePatches;
                using PatchPriorityLookup::_minValence;
                using PatchPriorityLookup::_maxValence;
                using PatchPriorityLookup::_maxSingularities;
                using PatchPriorityLookup::_maxConcavities;
                
                TS_PredicateWrapper               _predicate;                 ///< Thread-safe predicate able to interrupt the process.
                std::atomic<u_small_int_type>     _nThreads;                  ///< Number of threads sharing the lookup process.
                std::atomic<u_small_int_type>     _nCompleteThreads;          ///< Number of threads which have finished computing.
                std::atomic<bool>                 _isExecuting;               ///< Computation is ongoing or complete.
                std::atomic<bool>                 _stopped;                   ///< Computation interrupted.
                std::chrono::nanoseconds          _waitTime;                  ///< Time to wait before notify.
                std::atomic<bool>                 _printMsg;                  ///< true for printing status messages.
                mutable std::mutex                _mutex;                     ///< Mutual exclusion.
                std::thread                       _threadsManager;            ///< Main thread which manages computing threads and emitting notifications.
                std::condition_variable           _threadsCV;                 ///< Condition variable used to wait the threads sharing the lookup process.
                n_corners_type                    _nCorners;                  ///< Current query's number of corners.
                std::vector<num_type>             _boundaryLength;            ///< Current query's boundary lengths.
                FlowMap                           _flows;                     ///< List of flows to filter patches.
                PatchGeometrizatorPointer         _geometrizator;             ///< The responsible of giving geometry.
                std::list<PatchIndexRotationType> _perimeterPatches;          ///< List of templates currently passed the perimeter filter.
                std::list<PatchIndexRotationType> _flowPatches;               ///< List of templates currently passed the flow filter.
                std::atomic<count_type>           _nPatches;                  ///< Current number of extracted patches.
                std::atomic<count_type>           _nTemplatePatches;          ///< Current number of filtered template patches.
                
                using MT_Notifiable_Set           = std::unordered_set<MT_Notifiable_Ptr>;  ///< Type alias.
                using MT_Notifiable_Set_Iterator  = MT_Notifiable_Set::iterator;            ///< Type alias.
                MT_Notifiable_Set                 _notifiableSet;             ///< Set of MT_Notifiable objects to notify.
            };
            
        } // end namespace pl
        
    } // end namespace tri
} // end namespace vcg

#endif // MULTI_THREAD_PATCH_PRIORITY_LOOKUP_H
