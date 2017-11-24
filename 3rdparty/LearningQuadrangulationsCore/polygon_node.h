#ifndef POLYGON_NODE_H
#define POLYGON_NODE_H

#include <vector>
#include <list>
#if __cplusplus >= 201103L || _MSC_VER >= 1800
  #include <unordered_set>
#else
  #include <set>
#endif
#include <vcg/complex/complex.h>
#include <vcg/simplex/face/pos.h>
#include "default_types.h"
#include "mesh_support.h"
#include "polychord_support.h"

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief The PolygonNode class represents a polygon, storing information like corners and pointers to other nodes
 * representing polygons obtained by expansion.
 */
template < typename PolyMeshType >
class PolygonNode {
public:
  typedef typename PolyMeshType::VertexType             VertexType;
  typedef typename PolyMeshType::VertexPointer          VertexPointer;
  typedef typename PolyMeshType::FaceType               FaceType;
  typedef typename PolyMeshType::FacePointer            FacePointer;
  typedef PolygonNode<PolyMeshType> *                   PolygonNodePointer;
  typedef pl::PolychordSupport<PolyMeshType>            PolychordSupport;
  typedef typename PolychordSupport::PosType            PosType;
  typedef pl::MeshSupport<PolyMeshType>                 MeshSupport;
  typedef typename MeshSupport::FaceSet                 FaceSet;
  typedef typename MeshSupport::FaceSetIterator         FaceSetIterator;

  /**
   * @brief PolygonNode default constructor.
   */
  PolygonNode() : _nCorners(0),
                  _nFaces(0),
                  _nVertices(0),
                  _nEdges(0),
                  _nSingularities(0),
                  _parent(0),
                  _parentSide(std::numeric_limits<n_corners_type>::max()) { }

  /**
   * @brief PolygonNode constructs a new polygon node from a mesh face.
   * @param mesh The input mesh.
   * @param f The index of the face.
   */
  PolygonNode(PolyMeshType &mesh, const num_type f) {
    _nCorners = (n_corners_type)mesh.face.at(f).VN();
    assert(_nCorners > 0);
    _corners.resize(_nCorners);
    PosType PosType;
    for (size_t c = 0; c < _corners.size(); ++c) {
      PosType.Set(&mesh.face.at(f), c, mesh.face.at(f).V(c));
      _corners[c] = PosType;
    }
    _nFaces = 1;
    _nVertices = _nCorners;
    _nEdges = _nCorners;
    _nSingularities = 0;
    _faces.insert((num_type)f);
    _parent = 0;
    _parentSide = std::numeric_limits<n_corners_type>::max();
    _children.resize(_nCorners, 0);
  }

  /**
   * @brief PolygonNode Safe copy constructor.
   * @note The references to corners, parent and children are set null. These must be reset after the copy.
   * @param node The PolygonNode to be copied.
   */
  PolygonNode(const PolygonNode &node) : _nCorners(node._nCorners),
                                         _corners(node._corners),
                                         _nFaces(node._nFaces),
                                         _nVertices(node._nVertices),
                                         _nEdges(node._nEdges),
                                         _nSingularities(node._nSingularities),
                                         _faces(node._faces),
                                         _parent(0),
                                         _parentSide(std::numeric_limits<n_corners_type>::max()),
                                         _children(node._children.size(), 0) { }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief PolygonNode Move constructor.
   * @param node The PolygonNode to move.
   */
  PolygonNode(PolygonNode&& node) : _nCorners(node._nCorners),
                                    _corners(std::move(node._corners)),
                                    _nFaces(node._nFaces),
                                    _nVertices(node._nVertices),
                                    _nEdges(node._nEdges),
                                    _nSingularities(node._nSingularities),
                                    _faces(std::move(node._faces)),
                                    _parent(node._parent),
                                    _parentSide(node._parentSide),
                                    _children(std::move(node._children)) {
    node._nCorners = 0;
    node._corners.clear();
    node._nFaces = 0;
    node._nVertices = 0;
    node._nEdges = 0;
    node._nSingularities = 0;
    node._faces.clear();
    node._parent = 0;
    node._parentSide = std::numeric_limits<n_corners_type>::max();
    node._children.clear();
  }
#endif

  /**
   * @brief ~PolygonNode Destructor.
   * @note Parent and children are not destructed.
   */
  ~PolygonNode() {
    _nCorners = 0;
    _corners.clear();
    _nFaces = 0;
    _nVertices = 0;
    _nEdges = 0;
    _nSingularities = 0;
    _faces.clear();
    _parent = 0;
    _parentSide = std::numeric_limits<n_corners_type>::max();
    _children.clear();
  }

  /**
   * @brief operator = Assignment operator. See the copy constructor.
   * @param node The PolygonNode to be copied.
   * @return A reference to this.
   */
  PolygonNode & operator = (const PolygonNode &node) {
    if (&node != this) {
      _nCorners = node._nCorners;
      _corners = node._corners;
      _nFaces = node._nFaces;
      _nVertices = node._nVertices;
      _nEdges = node._nEdges;
      _nSingularities = node._nSingularities;
      _faces = node._faces;
      _parent = 0;
      _parentSide = std::numeric_limits<n_corners_type>::max();
      _children.clear();
      _children.resize(node._children.size(), 0);
    }
    return *this;
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief operator = Move assignment operator.
   * @param node The PolygonNode to move.
   * @return A reference to this.
   */
  PolygonNode & operator = (PolygonNode&& node) {
    if (&node != this) {
      _nCorners = node._nCorners;
      _corners = std::move(node._corners);
      _nFaces = node._nFaces;
      _nVertices = node._nVertices;
      _nEdges = node._nEdges;
      _nSingularities = node._nSingularities;
      _faces = std::move(node._faces);
      _parent = node._parent;
      _parentSide = node._parentSide;
      _children = std::move(node._children);
      node._nCorners = 0;
      node._corners.clear();
      node._nFaces = 0;
      node._nVertices = 0;
      node._nEdges = 0;
      node._nSingularities = 0;
      node._faces.clear();
      node._parent = 0;
      node._parentSide = std::numeric_limits<n_corners_type>::max();
      node._children.clear();
    }
    return *this;
  }
#endif

  /**
   * @brief numberOfCorners see _nCorners.
   * @return
   */
  n_corners_type numberOfCorners() const {
    return _nCorners;
  }

  /**
   * @brief numberOfCorners see _nCorners.
   * @return
   */
  n_corners_type & numberOfCorners() {
    return _nCorners;
  }

  /**
   * @brief numberOfSides returns the number of sides which is alwayas at least 1. See _nCorners and _corners.
   * @return
   */
  n_corners_type numberOfSides() const {
    return (n_corners_type)_corners.size();
  }

  /**
   * @brief corners returns the corners' position. See _corners.
   * @param corners
   */
  void corners(std::deque<PosType> &corners) const {
    corners = _corners;
  }

  /**
   * @brief corners returns the corners' position. See _corners.
   * @return
   */
  std::deque< vcg::face::Pos<FaceType> > & corners() {
    return _corners;
  }

  /**
   * @brief corners returns the corners' position. See _corners.
   * @return
   */
  const std::deque<PosType> & corners() const {
    return _corners;
  }

  /**
   * @brief corner returns the starting corner of the requested side (or the first vertex, if _nCorners=0). See _corners.
   * @param side
   * @param corner
   */
  void corner(const n_corners_type side, PosType &corner) const {
    corner = _corners.at(side);
  }

  /**
   * @brief corner returns the starting corner of the requested side (or the first vertex, if _nCorners=0). See _corners.
   * @param side
   * @return
   */
  PosType & corner(const n_corners_type side) {
    return _corners.at(side);
  }

  /**
   * @brief corner returns the starting corner of the requested side (or the first vertex, if _nCorners=0). See _corners.
   * @param side
   * @return
   */
  const PosType & corner(const n_corners_type side) const {
    return _corners.at(side);
  }

  /**
   * @brief numberOfFaces see _nFaces.
   * @return
   */
  num_type numberOfFaces() const {
    return _nFaces;
  }

  /**
   * @brief numberOfFaces see _nFaces.
   * @return
   */
  num_type & numberOfFaces() {
    return _nFaces;
  }

  /**
   * @brief numberOfVertices see _nVertices.
   * @return
   */
  num_type numberOfVertices() const {
    return _nVertices;
  }

  /**
   * @brief numberOfVertices see _nVertices.
   * @return
   */
  num_type & numberOfVertices() {
    return _nVertices;
  }

  /**
   * @brief numberOfEdges see _nEdges.
   * @return
   */
  num_type numberOfEdges() const {
    return _nEdges;
  }

  /**
   * @brief numberOfEdges see _nEdges.
   * @return
   */
  num_type & numberOfEdges() {
    return _nEdges;
  }

  /**
   * @brief numberOfSingularities see _nSingularities.
   * @return
   */
  num_type numberOfSingularities() const {
    return _nSingularities;
  }

  /**
   * @brief numberOfSingularities see _nSingularities.
   * @return
   */
  num_type & numberOfSingularities() {
    return _nSingularities;
  }

  /**
   * @brief faces see _faces.
   * @return
   */
  const FaceSet & faces() const {
    return _faces;
  }

  /**
   * @brief faces see _faces.
   * @return
   */
  FaceSet & faces() {
    return _faces;
  }

  /**
   * @brief allFaces returns a set of all the faces composing this polygon (from this node to the root ancestor). See _faces.
   * @param allFaces
   */
  void allFaces(FaceSet &allFaces) const {
    allFaces = _faces;
    PolygonNodePointer node = _parent;
    while (node != 0) {
      allFaces.insert(node->_faces.begin(), node->_faces.end());
      node = node->_parent;
    }
  }

  /**
   * @brief parent see _parent.
   * @return
   */
  PolygonNodePointer parent() const {
    return _parent;
  }

  /**
   * @brief parent see _parent.
   * @return
   */
  PolygonNodePointer & parent() {
    return _parent;
  }

  /**
   * @brief parentSide see _parentSide.
   * @return
   */
  n_corners_type parentSide() const {
    return _parentSide;
  }

  /**
   * @brief parentSide see _parentSide.
   * @return
   */
  n_corners_type & parentSide() {
    return _parentSide;
  }

  /**
   * @brief children see _children.
   * @param children
   */
  void children(std::deque<PolygonNodePointer> & children) const {
    children = _children;
  }

  /**
   * @brief children see _children.
   * @return
   */
  std::deque<PolygonNodePointer> & children() {
    return _children;
  }

  /**
   * @brief children see _children.
   * @return
   */
  const std::deque<PolygonNodePointer> & children() const {
    return _children;
  }

  /**
   * @brief child returns the child node expanded from the requested side. See _children.
   * @param side
   */
  PolygonNodePointer child(const n_corners_type side) const {
    return _children.at(side);
  }

  /**
   * @brief child returns the child node expanded from the requested side. See _children.
   * @param side
   */
  PolygonNodePointer & child(const n_corners_type side) {
    return _children.at(side);
  }

  /**
   * @brief detachFromChildren sets its children as roots (_parent = 0) and copies _faces into those of the children.
   */
  void detachFromChildren() {
    // for each child
    for (size_t c = 0; c < _children.size(); ++c)
      if (_children[c] != 0) {
        // copy this convex polygon faces into this child faces
        _children[c]->_faces.insert(_faces.begin(), _faces.end());
        // detach the child
        _children[c]->_parent = 0;
        // unlink
        _children[c] = 0;
      }
  }

  /**
   * @brief expandPolygon constructs a new polygon by expanding this from a selected side.
   * @param mesh The input mesh.
   * @param side The side to expand from.
   * @param convexOnly true if convex only corners must be found, false otherwise.
   * @param allEdges true if all edges must be considered as corners, false otherwise.
   * @param singularityVec A vector of pairs <valence,#vertices> in order of valence.
   * @param maxNFaces If the new polygon has more than this threshold, it is discarded.
   * @return true if it has been successfully expanded, false otherwise.
   */
  bool expandPolygon(PolyMeshType &mesh, const n_corners_type side, const bool convexOnly, const bool allEdges,
                     std::vector< std::pair<valence_type, num_type> > &singularityVec,
                     const num_type maxNFaces = std::numeric_limits<num_type>::max()) {
    // requirements
    vcg::tri::RequireFFAdjacency(mesh);
    vcg::tri::RequirePerFaceMark(mesh);
    vcg::tri::RequirePerVertexMark(mesh);
    // check input
    if (side >= numberOfSides()) {
      std::cout << "Wrong side (" << side << " of " << numberOfSides() << "). Skipped." << std::endl;
      return false;
    }
    if (_children[side] != 0) {
      std::cout << "Side " << side << " of " << numberOfSides() << " already expanded. Skipped." << std::endl;
      return true;
    }

    // some variables
    PolygonNodePointer child = 0;
    PosType pos, runPos, firstPos;
    num_type X = 0;
    FaceSet allFaces;
    FaceSetIterator fIt;
    std::deque<PosType> corners;
    num_type nCorners = -1;
    size_t firstCorner = 0;

    // build a new node
    try {
      child = new PolygonNode<PolyMeshType>();
    } catch(std::bad_alloc ba) {
      std::cout << "Bad allocation: " << ba.what() << std::endl;
      child = 0;
    }
    if (child == 0) {
      std::cout << "Error allocating memory for a new node. Skipped." << std::endl;
      return false;
    }
    child->_parent = this;
    child->_parentSide = side;
    child->_nFaces = _nFaces;
    child->_nVertices = _nVertices;
    child->_nEdges = _nEdges;
    _children[side] = child;

    // increment the mark!!
    MeshSupport::IncrementMark(mesh);

    // mark the current polygon
    markAllFacesAndVertices(mesh);

    // scan side edge by edge
    pos = _corners[side];
    do {
      // add the next external face to the set
      if (!pos.IsBorder() && !vcg::tri::IsMarked(mesh, pos.FFlip()))
        child->_faces.insert(vcg::tri::Index(mesh, pos.FFlip()));
      else {
        delete child;
        _children[side] = child = 0;
        return false;
      }

      // go on next edge
      pos.FlipV();
      pos.FlipE();
      while (!pos.IsBorder() && vcg::tri::IsMarked(mesh, pos.FFlip())) {
        pos.FlipF();
        pos.FlipE();
      }
    } while (pos != _corners[(side+1)%_corners.size()]);
    // if the whole side is on border, it's not possible to expand
    if (child->_faces.empty()) {
      delete child;
      _children[side] = child = 0;
      return false;
    }

    // mark the new added faces
    for (fIt = child->_faces.begin(); fIt != child->_faces.end(); ++fIt) {
      // for each vertex/edge of the new face
      for (int j = 0; j < mesh.face[*fIt].VN(); ++j) {
        // if this vertex is not marked, increment this polygon's number of vertices
        if (!vcg::tri::IsMarked(mesh, mesh.face[*fIt].V(j))) {
          ++child->_nVertices;
          // and mark it
          vcg::tri::Mark(mesh, mesh.face[*fIt].V(j));
        }
        // if it's on border, or its adjacent face isn't marked, increment this polygon's number of edges
        if (vcg::face::IsBorder(mesh.face[*fIt], j) || !vcg::tri::IsMarked(mesh, mesh.face[*fIt].cFFp(j)))
          ++child->_nEdges;
      }

      // increment this polygon's number of faces
      ++child->_nFaces;

      // check if the number of faces is below the maximum allowed
      if (child->_nFaces > maxNFaces) {
        delete child;
        _children[side] = child = 0;
        return false;
      }

      // finally, mark this new face
      vcg::tri::Mark(mesh, &mesh.face[*fIt]);
    }

    do {
      // compute the Euler characteristic and check if it is = 1, otherwise return (not disk-like)
      X = child->_nVertices - child->_nEdges + child->_nFaces;
      if (X != 1) {
        delete child;
        _children[side] = child = 0;
        return false;
      }

      // find corners and check if convex
      pos = _corners[side];
      pos.FlipF();
      pos.FlipE();
      nCorners = PolychordSupport::FindPolygonCorners(mesh, corners, convexOnly, allEdges, pos);

      // if no boundary is found, this polygon is closed
      if (corners.empty()) {
        delete child;
        _children[side] = child = 0;
        return false;
      }

      if (convexOnly && nCorners < 0) {
        // get the concave corner
        runPos = corners.back();

        // if the concave corner is on border, it's not possible to expand
        if (runPos.IsBorder()) {
          delete child;
          _children[side] = child = 0;
          return false;
        }

        // add the next external face
        child->_faces.insert(vcg::tri::Index(mesh, runPos.FFlip()));

        // for each vertex/edge of the new face
        for (num_type j = 0; j < (num_type)runPos.FFlip()->VN(); ++j) {
          // if this vertex is not marked, increment this polygon's number of vertexes
          if (!vcg::tri::IsMarked(mesh, runPos.FFlip()->V(j))) {
            ++child->_nVertices;
            // and mark it
            vcg::tri::Mark(mesh, runPos.FFlip()->V(j));
          }
          // if it's on border, or its adjacent face isn't marked, increment this polygon's number of edges
          if (vcg::face::IsBorder(*runPos.FFlip(), j) || !vcg::tri::IsMarked(mesh, runPos.FFlip()->cFFp(j)))
            ++child->_nEdges;
        }

        // increment this polygon's number of faces
        ++child->_nFaces;

        // check if the number of faces is below the maximum allowed
        if (child->_nFaces > maxNFaces) {
          delete child;
          _children[side] = child = 0;
          return false;
        }

        // finally, mark this new face
        vcg::tri::Mark(mesh, runPos.FFlip());
      }
    } while (convexOnly && nCorners < 0);

    // set the number of corners
    assert(nCorners >= 0);
    child->_nCorners = nCorners;

    // set the corners
    child->_corners.resize(corners.size());
    child->_children.resize(corners.size(), 0);
    // find the first corner (in order of vertex index)
    if (child->_nCorners > 0) {
      firstCorner = 0;
      for (size_t i = 1; i < corners.size(); ++i)
        if (corners[i].V() < corners[firstCorner].V())
          firstCorner = i;
      for (size_t i = 0; i < corners.size(); ++i)
        child->_corners[i] = corners[(firstCorner + i) % corners.size()];
    } else {
      firstPos = corners.front();
      pos = runPos = firstPos;
      do {
        if (runPos.V() < firstPos.V())
          firstPos = runPos;
        runPos.FlipV();
        runPos.FlipE();
        runPos.FlipF();
        runPos.FlipE();
      } while (runPos != pos);
      child->_corners.front() = firstPos;
    }

    // count singularities
    child->allFaces(allFaces);
    child->_nSingularities = MeshSupport::CountInternalSingularities(mesh, allFaces, singularityVec);

    // return
    return true;
  }



private:
  /**
   * @brief markAllFacesAndVertices marks all the faces and vertices of this polygon (from this node to the root ancestor). See _faces.
   * @param mesh The mesh the faces belong to.
   */
  void markAllFacesAndVertices(PolyMeshType &mesh) {
    FaceSetIterator fIt;
    PolygonNodePointer node = this;
    while (node != 0) {
      // mark all the added faces by the _current node
      for (fIt = node->_faces.begin(); fIt != node->_faces.end(); ++fIt) {
        vcg::tri::Mark(mesh, &mesh.face[*fIt]);
        // mark all its vertexes
        for (int j = 0; j < mesh.face[*fIt].VN(); ++j)
          vcg::tri::Mark(mesh, mesh.face[*fIt].V(j));
      }
      // go to parent
      node = node->_parent;
    }
  }

  // this polygon node data
  n_corners_type                    _nCorners;          ///< Number of corners. When it is 0, there is 1 side and thus a starting vertex in _corners[0].
  std::deque<PosType>               _corners;           ///< Positions of corners. If _nCorners=0, there will be a single position of a border vertex.
  num_type                          _nFaces;            ///< Number of quad faces.
  num_type                          _nVertices;         ///< Number of vertices.
  num_type                          _nEdges;            ///< Number of edges.
  num_type                          _nSingularities;    ///< Number of singularities.
  FaceSet                           _faces;             ///< Indices of the faces of this node (eventually only those added w.r.t. its parent).
  // other polygons' references
  PolygonNodePointer                _parent;            ///< The parent node.
  n_corners_type                    _parentSide;        ///< The parent's side this polygon has been expanded from.
  std::deque<PolygonNodePointer>    _children;          ///< Children nodes (one for each side).
};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // POLYGON_NODE_H
