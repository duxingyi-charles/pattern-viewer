#ifndef POLYGON_FOREST_H
#define POLYGON_FOREST_H

#include "polygon_tree.h"

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief The PolygonForest class is able to build and store a forest of PolygonTree objects, one for each face of a mesh.
 */
template < typename PolyMeshType >
class PolygonForest {
public:
  typedef pl::PolygonTree<PolyMeshType>               PolygonTree;
  typedef typename PolygonTree::CornersType           CornersType;
  typedef pl::PolygonRegister<PolyMeshType>           PolygonRegister;
  typedef pl::PatchLearner<PolyMeshType>              PatchLearner;

  /**
   * @brief PolygonForest Default constructor.
   */
  PolygonForest() { }

  /**
   * @brief PolygonForest Constructor initializing the size of the forest.
   * @param size The size used to resize the forest.
   */
  PolygonForest(const num_type size) : _trees(size) { }

  /**
   * @brief PolygonForest Copy constructor.
   * @param forest The forest to be copied.
   */
  PolygonForest(const PolygonForest &forest) : _trees(forest._trees) { }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief PolygonForest Move constructor.
   * @param forest The PolygonForest to move.
   */
  PolygonForest(PolygonForest&& forest) : _trees(std::move(forest._trees)) {
    forest.resize(0);
  }
#endif

  /**
   * @brief operator = Assignment operator. Copies the trees from forest to this.
   * @param forest The forest to be copied.
   * @return A reference to this.
   */
  PolygonForest & operator = (const PolygonForest &forest) {
    if (&forest != this)
      _trees = forest._trees;
    return *this;
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief operator = Move assignment operator.
   * @param forest the PolygonForest to move.
   * @return A reference to this.
   */
  PolygonForest & operator = (PolygonForest&& forest) {
    if (&forest != this) {
      _trees = std::move(forest._trees);
      forest.resize(0);
    }
    return *this;
  }
#endif

  /**
   * @brief size gives the number of PolygonTrees in this forest.
   * @return The number of trees.
   */
  num_type size() const {
    return (num_type)_trees.size();
  }

  /**
   * @brief resize changes the size of this forest.
   * @param size The new size.
   */
  void resize(num_type size) {
    _trees.clear();
    _trees.resize(size);
  }

  /**
   * @brief operator [] gives access to the tree at index.
   * @param index The index of the tree to be accessed.
   * @return A reference to the tree at index.
   */
  PolygonTree & operator [] (const num_type index) {
    return _trees.at(index);
  }

  /**
   * @brief operator [] gives access to the tree at index.
   * @param index The index of the tree to be accessed.
   * @return A const reference to the tree at index.
   */
  const PolygonTree & operator [] (const num_type index) const {
    return _trees.at(index);
  }

  /**
   * @brief buildPolygonTreeAt builds a tree from face at faceIndex.
   * It results into a tree where each node is a polygon, and each child of a node is a polygon built
   * by expanding one of the sides of the node's polygon.
   * @param mesh The input mesh.
   * @param faceIndex The starting face (a.k.a. the root polygon).
   * @param polygonRegister The register which avoids duplicate polygons.
   * @param patchLearner The patch acquisitor.
   * @param minNFaces The minimum number of per-polygon faces for counting and template storing.
   * @param maxNFaces The maximum number of per-polygon faces.
   * @param cornersType The type of corners allowed.
   * @param bindPolygon true if the template patch must be linked to the polygon whre it's been found, false otherwise.
   * @param collapsePolygon true if the polygon must be simplified before acquisition, false otherwise.
   * @param keepInMemory true if the nodes must be kept in memory, false otherwise.
   */
  void buildPolygonTreeAt(PolyMeshType &mesh,
                          const num_type faceIndex,
                          PolygonRegister *polygonRegister,
                          PatchLearner &patchLearner,
                          const num_type minNFaces = 0,
                          const num_type maxNFaces = std::numeric_limits<num_type>::max(),
                          const CornersType cornersType = CornersType::CONVEX,
                          const bool bindPolygon = false,
                          const bool collapsePolygon = true,
                          const bool keepInMemory = false) {
    _trees.at(faceIndex).buildPolygonsTree(mesh, faceIndex, polygonRegister, patchLearner, minNFaces, maxNFaces,
                                           cornersType, bindPolygon, collapsePolygon, keepInMemory);
  }

  /**
   * @brief numberOfPolygons gives the total number of enumerated polygons.
   * @return The total number of enumerated polygons.
   */
  count_type numberOfPolygons() const {
    count_type sum = 0;
    for (size_t i = 0; i < _trees.size(); ++i)
      sum += _trees.at(i).numberOfPolygons();
    return sum;
  }

  /**
   * @brief maxNumberOfSides gives the maximum number of sides of any enumerated polygon over the whole mesh.
   * @return The maximum number of sides of any enumerated polygon.
   */
  n_corners_type maxNumberOfSides() const {
    n_corners_type maxSides = 0;
    for (size_t i = 0; i < _trees.size(); ++i)
      if (_trees.at(i).maxNumberOfSides() > maxSides)
        maxSides = _trees.at(i).maxNumberOfSides();
    return maxSides;
  }

  /**
   * @brief maxNumberOfSingularities gives the maximum number of singularities of any enumerated polygon over the whole mesh.
   * @return The maximum number of singularities of any enumerated polygon.
   */
  n_corners_type maxNumberOfSingularities() const {
    n_corners_type maxSingularities = 0;
    for (size_t i = 0; i < _trees.size(); ++i)
      if (_trees.at(i).maxNumberOfSingularities() > maxSingularities)
        maxSingularities = _trees.at(i).maxNumberOfSingularities();
    return maxSingularities;
  }

private:
  std::vector<PolygonTree>      _trees;     ///< A forest consists of a set of trees.
};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // POLYGON_FOREST_H
