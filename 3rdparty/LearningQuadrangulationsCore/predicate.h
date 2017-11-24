#ifndef PREDICATE_H
#define PREDICATE_H

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  #include <memory>
#endif

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief The Predicate class represents a predicate which always evaluates to true.
 * Inherit it and override operator() to change its behaviour.
 */
class Predicate {
public:
  /**
   * @brief Predicate Default constructor.
   */
  Predicate() { }

  /**
   * @brief ~Predicate Default destructor.
   */
  virtual ~Predicate() { }

  /**
   * @brief operator () evaluates the predicate and returns its value.
   * @return The evaluation of this predicate; always true.
   */
  virtual bool operator () () const {
    return true;
  }
};

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  typedef std::shared_ptr<Predicate>  PredicatePointer;
#else
  typedef Predicate *                 PredicatePointer;
#endif

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // PREDICATE_H
