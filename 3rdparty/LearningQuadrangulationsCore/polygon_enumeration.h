#ifndef POLYGON_ENUMERATION_H
#define POLYGON_ENUMERATION_H

#include "polygon_forest.h"

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief The PolygonEnumerator class provides a simplified interface to PolygonForest
 * for the enumeration of polygons.
 */
template < typename PolyMeshType >
class PolygonEnumerator {
public:
  typedef pl::PolygonForest<PolyMeshType>         PolygonForest;
  typedef typename PolygonForest::PolygonTree     PolygonTree;
  typedef typename PolygonForest::CornersType     CornersType;
  typedef pl::PolygonRegister<PolyMeshType>       PolygonRegister;
  typedef pl::PatchLearner<PolyMeshType>          PatchLearner;

  /**
   * @brief PolygonEnumerator Default constructor.
   */
  PolygonEnumerator() { }

  /**
   * @brief PolygonEnumerator Constructor initializing the size of the forest.
   * @param size The number of trees in the forest.
   */
  PolygonEnumerator(const num_type size) : _forest(size) { }

  /**
   * @brief PolygonEnumerator Copy constructor.
   * @param enumerator The PolygonEnumerator to copy.
   */
  PolygonEnumerator(const PolygonEnumerator &enumerator) : _forest(enumerator._forest) { }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief PolygonEnumerator Move contructor.
   * @param enumerator The PolygonEnumerator to move.
   */
  PolygonEnumerator(PolygonEnumerator&& enumerator) : _forest(std::move(enumerator._forest)) {
    enumerator.resize(0);
  }
#endif

  /**
   * @brief ~PolygonEnumerator Default destructor.
   */
  virtual ~PolygonEnumerator() { }

  /**
   * @brief operator = Copy assignment operator.
   * @param enumerator The PolygonEnumerator to copy.
   * @return A reference to this.
   */
  PolygonEnumerator & operator = (const PolygonEnumerator &enumerator) {
    if (&enumerator != this)
      _forest = enumerator._forest;
    return *this;
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief operator = Move assignment operator.
   * @param enumerator The PolygonEnumerator to move.
   * @return A reference to this.
   */
  PolygonEnumerator & operator = (PolygonEnumerator&& enumerator) {
    if (&enumerator != this) {
      _forest = std::move(enumerator._forest);
      enumerator.resize(0);
    }
    return *this;
  }
#endif

  /**
   * @brief size gives the number of PolygonTrees in the forest.
   * @return The size of the forest.
   */
  virtual num_type size() const {
    return _forest.size();
  }

  /**
   * @brief resize changes the size of the forest.
   * @param size The new size.
   */
  virtual void resize(const num_type size) {
    _forest.resize(size);
  }

  /**
   * @brief operator [] gives access to the tree at index.
   * @param index The index of the tree.
   * @return A reference to the tree at index.
   */
  virtual PolygonTree & operator [] (const num_type index) {
    return _forest[index];
  }

  /**
   * @brief operator [] gives access to the tree at index.
   * @param index The index of the tree.
   * @return A const reference to the tree at index.
   */
  virtual const PolygonTree & operator [] (const num_type index) const {
    return _forest[index];
  }

  /**
   * @brief numberOfPolygons gives the total number of enumerated polygons.
   * @return The total number of enumerated polygons.
   */
  virtual count_type numberOfPolygons() const {
    return _forest.numberOfPolygons();
  }

  /**
   * @brief maxNumberOfSides gives the maximum number of sides of any enumerated polygon over the whole mesh.
   * @return The maximum number of sides of any enumerated polygon.
   */
  virtual n_corners_type maxNumberOfSides() const {
    return _forest.maxNumberOfSides();
  }

  /**
   * @brief maxNumberOfSingularities gives the maximum number of singularities of any enumerated polygon over the whole mesh.
   * @return The maximum number of singularities of any enumerated polygon.
   */
  virtual n_corners_type maxNumberOfSingularities() const {
    return _forest.maxNumberOfSingularities();
  }

  /**
   * @brief enumerate finds all the polygons in mesh with at most maxNFaces.
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
  virtual void enumerate(PolyMeshType &mesh,
                         PolygonRegister *polygonRegister,
                         PatchLearner &patchLearner,
                         const num_type minNFaces = 0,
                         const num_type maxNFaces = std::numeric_limits<num_type>::max(),
                         const CornersType cornersType = CornersType::CONVEX,
                         const bool bindPolygon = false,
                         const bool collapsePolygon = true,
                         const bool keepInMemory = false) {
    // check if the forest size is the same
    if (_forest.size() != (num_type)mesh.face.size())
      _forest.resize(mesh.face.size());

    // enumerate
    for (num_type i = 0; i < _forest.size(); ++i)
      _forest.buildPolygonTreeAt(mesh, i, polygonRegister, patchLearner, minNFaces, maxNFaces, cornersType, bindPolygon, collapsePolygon, keepInMemory);
  }

protected:
  PolygonForest         _forest;      ///< The set of trees of convex polygons for each face of the mesh.
};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // POLYGON_ENUMERATION_H
