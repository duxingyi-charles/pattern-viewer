#ifndef THREAD_SAFE_PREDICATE_H
#define THREAD_SAFE_PREDICATE_H

#include <mutex>
#include "../predicate.h"

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

class TS_PredicateWrapper : public Predicate {
public:
  /**
   * @brief TS_PredicateWrapper Default constructor.
   */
  TS_PredicateWrapper() : _predicate(PredicatePointer(0)), _enabled(true) { }

  /**
   * @brief TS_PredicateWrapper Copy initializer constructor.
   * @param predicate The Predicate to wrap.
   */
  TS_PredicateWrapper(const PredicatePointer &predicate) : _predicate(predicate), _enabled(true) { }

  /**
   * @brief TS_PredicateWrapper Copy constructor.
   * @param ts_pred The TS_PredicateWrapper to copy.
   */
  TS_PredicateWrapper(const TS_PredicateWrapper &ts_pred) {
    std::lock_guard<std::mutex> rhs_lock(ts_pred._mutex);
    _predicate = ts_pred._predicate;
    _enabled = ts_pred._enabled.load();
  }

  /**
   * @brief TS_PredicateWrapper Move constructor.
   * @param ts_pred The TS_PredicateWrapper to move.
   */
  TS_PredicateWrapper(TS_PredicateWrapper&& ts_pred) {
    std::lock_guard<std::mutex> rhs_lock(ts_pred._mutex);
    _predicate = ts_pred._predicate;
    _enabled = ts_pred._enabled.load();
    ts_pred._predicate = PredicatePointer(0);
  }

  /**
   * @brief operator = Copy assignement operator.
   * @param ts_pred The TS_PredicateWrapper to copy.
   * @return A reference to this.
   */
  TS_PredicateWrapper & operator = (const TS_PredicateWrapper &ts_pred) {
    if (this != &ts_pred) {
      std::lock_guard<std::mutex> rhs_lock(ts_pred._mutex);
      std::lock_guard<std::mutex> lock(_mutex);
      _predicate = ts_pred._predicate;
      _enabled = ts_pred._enabled.load();
    }
    return *this;
  }

  /**
   * @brief operator = Move assignment operator.
   * @param ts_pred The TS_PredicateWrapper to move.
   * @return A reference to this.
   */
  TS_PredicateWrapper & operator = (TS_PredicateWrapper&& ts_pred) {
    if (this != &ts_pred) {
      std::lock_guard<std::mutex> rhs_lock(ts_pred._mutex);
      std::lock_guard<std::mutex> lock(_mutex);
      _predicate = ts_pred._predicate;
      _enabled = ts_pred._enabled.load();
      ts_pred._predicate = PredicatePointer(0);
    }
    return *this;
  }

  /**
   * @brief setPredicate wraps a new predicate object.
   * @param predicate The Predicate to wrap.
   */
  void setPredicate(const PredicatePointer &predicate) {
    std::lock_guard<std::mutex> lock(_mutex);
    _predicate = predicate;
  }

  /**
   * @brief enable enables the wrapped predicate.
   */
  void enable() {
    _enabled = true;
  }

  /**
   * @brief disable disables the wrapper predicate.
   */
  void disable() {
    _enabled = false;
  }

  /**
   * @brief operator () evaluates the predicate and returns its value. Thread-safe version.
   * @return The evaluation of this predicate.
   */
  bool operator () () const {
    std::lock_guard<std::mutex> lock(_mutex);
    return _enabled && (!_predicate || _predicate->operator()());
  }

private:
  PredicatePointer      _predicate;   ///< The Predicate to wrap.
  mutable std::mutex    _mutex;       ///< Mutual exclusion.
  std::atomic<bool>     _enabled;     ///< Is the wrapped predicate enabled or disabled?
};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // THREAD_SAFE_PREDICATE_H
