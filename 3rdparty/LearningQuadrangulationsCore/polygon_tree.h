#ifndef POLYGON_TREE_H
#define POLYGON_TREE_H

#include <stack>
#include "polygon_node.h"
#include "polygon_register.h"
#include "patch.h"
#include "patch_learner.h"

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief The PolygonTree class is responsible of enumerating all the polygons that can be built
 * starting from a given face.
 */
template < typename PolyMeshType >
class PolygonTree {
public:
  typedef typename PolyMeshType::VertexType                                   VertexType;
  typedef typename PolyMeshType::VertexPointer                                VertexPointer;
  typedef typename PolyMeshType::FaceType                                     FaceType;
  typedef typename PolyMeshType::FacePointer                                  FacePointer;
  typedef pl::PolygonNode<PolyMeshType>                                       PolygonNode;
  typedef typename PolygonNode::FaceSet                                       FaceSet;
  typedef typename PolygonNode::FaceSetIterator                               FaceSetIterator;
  typedef typename PolygonNode::PolygonNodePointer                            PolygonNodePointer;
  typedef typename PolygonNode::PosType                                       PosType;
  typedef pl::PolygonRegister<PolyMeshType>                                   PolygonRegister;
  typedef pl::Patch<PolyMeshType>                                             Patch;
  typedef pl::PatchLearner<PolyMeshType>                                      PatchLearner;

  /**
   * @brief The CornersType enum defines the type of corners allowed by the enumeration.
   */
  enum CornersType {
    CONVEX      = 0x1,                  ///< convex only
    CONCAVE     = 0x2,                  ///< concave only
    REGULAR     = 0x4,                  ///< regular only (not working and not useful at all)
    IRREGULAR   = CONVEX | CONCAVE,     ///< both convex and concave
    NOTCONVEX   = CONCAVE | REGULAR,    ///< special concave made by edge expansion (instead of by whole side)
    ANY         = IRREGULAR | REGULAR   ///< everything
  };

  /**
   * @brief PolygonTree default constructor. No nodes, no root.
   */
  PolygonTree() : _root(0), _current(0), _numberOfPolygons(0), _maxNSides(0), _maxNSingularities(0) { }

  /**
   * @brief PolygonTree constructs a tree of polygons starting from the given face.
   * It results into a tree where each node is a polygon, and each child of a node is a polygon built
   * by expanding one of the sides of the node's polygon.
   * @param mesh The input mesh.
   * @param faceIndex The starting face (a.k.a. the root polygon).
   * @param polygonRegister The register which avoids duplicate polygons.
   * @param learner The patch acquisitor.
   * @param minNFaces The minimum number of per-polygon faces for counting and template storing.
   * @param maxNFaces The maximum number of per-polygon faces.
   * @param cornersType The type of corners allowed.
   * @param bindPolygon true if the template patch must be linked to the polygon whre it's been found, false otherwise.
   * @param collapsePolygon true if the polygon must be simplified before acquisition, false otherwise.
   * @param keepInMemory true if the nodes must be kept in memory, false otherwise.
   */
  PolygonTree(PolyMeshType &mesh,
              const num_type faceIndex,
              PolygonRegister *polygonRegister,
              PatchLearner &learner,
              const num_type minNFaces = 0,
              const num_type maxNFaces = std::numeric_limits<num_type>::max(),
              const CornersType cornersType = CONVEX,
              const bool bindPolygon = false,
              const bool collapsePolygon = true,
              const bool keepInMemory = false) : _root(0), _current(0), _numberOfPolygons(0), _maxNSides(0), _maxNSingularities(0) {
    buildPolygonsTree(mesh, faceIndex, polygonRegister, learner, minNFaces, maxNFaces, cornersType, bindPolygon, collapsePolygon, keepInMemory);
  }

  /**
   * @brief PolygonTree Copy constructor.
   * @param tree The PolygonTree to be copied.
   */
  PolygonTree(const PolygonTree &tree) {
    copyPolygonTree(tree);
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief PolygonTree Move constructor.
   * @param tree The PolygonTree to move.
   */
  PolygonTree(PolygonTree&& tree) {
    _root = tree._root;
    _current = tree._current;
    _numberOfPolygons = tree._numberOfPolygons;
    _maxNSides = tree._maxNSides;
    _maxNSingularities = tree._maxNSingularities;
    tree._root = 0;
    tree._current = 0;
    tree._numberOfPolygons = 0;
    tree._maxNSides = 0;
    tree._maxNSingularities = 0;
  }
#endif

  /**
   * @brief operator = Assignment operator.
   * @param tree The PolygonTree to be copied.
   * @return A reference to this.
   */
  PolygonTree & operator = (const PolygonTree &rhs_tree) {
    if (&rhs_tree != this)
      copyPolygonTree(rhs_tree);
    return *this;
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief operator = Move assignment operator.
   * @param tree The PolygonTree to move.
   * @return A reference to this.
   */
  PolygonTree & operator = (PolygonTree&& tree) {
    if (&tree != this) {
      deAllocate();
      _root = tree._root;
      _current = tree._current;
      _numberOfPolygons = tree._numberOfPolygons;
      _maxNSides = tree._maxNSides;
      _maxNSingularities = tree._maxNSingularities;
      tree._root = 0;
      tree._current = 0;
      tree._numberOfPolygons = 0;
      tree._maxNSides = 0;
      tree._maxNSingularities = 0;
    }
    return *this;
  }
#endif

  /**
   * @brief PolygonTree destructor.
   */
  ~PolygonTree() {
    deAllocate();
  }

  /**
   * @brief deAllocate
   */
  void deAllocate() {
    if (_root != 0) {
      std::stack<PolygonNodePointer> stack;
      stack.push(_root);
      while (!stack.empty()) {
        // pick the top polygon from the stack
        _current = stack.top();
        // pop from the stack
        stack.pop();
        // for each child
        for (n_corners_type i = 0; i < _current->numberOfSides(); ++i) {
          // if a new convex polygon has been built, push it on queue
          if (_current->child(i) != 0)
            stack.push(_current->child(i));
        }
        // delete this node
        delete _current;
      }
    }
    _root = _current = 0;
    _numberOfPolygons = 0;
    _maxNSides = 0;
    _maxNSingularities = 0;
  }

  /**
   * @brief changeReferencedMesh changes pointers of corners' Pos from oldMesh to newMesh. It may be useful when copying trees.
   * @param oldMesh
   * @param newMesh
   */
  void changeReferencedMesh(const PolyMeshType &oldMesh, PolyMeshType &newMesh) {
    // requirements
    vcg::tri::RequireFFAdjacency(oldMesh);
    vcg::tri::RequireFFAdjacency(newMesh);

    // general checks
    if (oldMesh.FN() != newMesh.FN() || oldMesh.VN() != newMesh.VN()
        || oldMesh.face.size() != newMesh.face.size() || oldMesh.vert.size() != newMesh.vert.size()) {
      std::cout << "Impossible to change refereced mesh of ConvexPolygonTree. Meshes are not equal. Nothing changed." << std::endl;
      return;
    }
    for (size_t f = 0; f < oldMesh.face.size(); ++f) {
      if ((vcg::tri::HasPerFaceFlags(oldMesh) && !vcg::tri::HasPerFaceFlags(newMesh))
          || (!vcg::tri::HasPerFaceFlags(oldMesh) && vcg::tri::HasPerFaceFlags(newMesh))
          || (vcg::tri::HasPerFaceFlags(oldMesh) && oldMesh.face[f].IsD() && !newMesh.face[f].IsD())
          || (vcg::tri::HasPerFaceFlags(oldMesh) && !oldMesh.face[f].IsD() && newMesh.face[f].IsD())
          || oldMesh.face[f].VN() != newMesh.face[f].VN()) {
        std::cout << "Impossible to change refereced mesh of ConvexPolygonTree. Meshes are not equal. Nothing changed." << std::endl;
        return;
      }
      for (int v = 0; v < oldMesh.face[f].VN(); ++v)
        if (vcg::tri::Index(oldMesh, oldMesh.face[f].cFFp(v)) != vcg::tri::Index(newMesh, newMesh.face[f].cFFp(v))
            || oldMesh.face[f].cFFi(v) != newMesh.face[f].cFFi(v)
            || vcg::tri::Index(oldMesh, oldMesh.face[f].cV(v)) != vcg::tri::Index(newMesh, newMesh.face[f].cV(v))) {
          std::cout << "Impossible to change refereced mesh of ConvexPolygonTree. Meshes are not equal. Nothing changed." << std::endl;
          return;
        }
    }

    // if empty, return
    if (_root == 0)
      return;

    // update corners references
    std::stack<PolygonNodePointer> stack;
    PolygonNodePointer current = 0;
    stack.push(_root);
    while (!stack.empty()) {
      current = stack.top();
      stack.pop();
      for (n_corners_type c = 0; c < current->numberOfSides(); ++c)
        current->corner(c).Set(&newMesh.face[vcg::tri::Index(oldMesh, current->corner(c).F())],
                                current->corner(c).E(),
                               &newMesh.vert[vcg::tri::Index(oldMesh, current->corner(c).V())]);
      for (n_corners_type i = 0; i < current->numberOfSides(); ++i)
        if (current->child(i) != 0)
          stack.push(current->child(i));
    }
  }

  /**
   * @brief buildPolygonsTree enumerates all the polygons starting from a given face.
   * It results into a tree where each node is a polygon, and each child of a node is a polygon built
   * by expanding one of the sides of the node's polygon.
   * @param mesh The input mesh.
   * @param faceIndex The starting face (a.k.a. the root polygon).
   * @param polygonRegister The register which avoids duplicate polygons.
   * @param learner The patch acquisitor.
   * @param minNFaces The minimum number of per-polygon faces for counting and template storing.
   * @param maxNFaces The maximum number of per-polygon faces.
   * @param cornersType The type of corners allowed.
   * @param bindPolygon true if the template patch must be linked to the polygon whre it's been found, false otherwise.
   * @param collapsePolygon true if the polygon must be simplified before acquisition, false otherwise.
   * @param keepInMemory true if the nodes must be kept in memory, false otherwise.
   */
  void buildPolygonsTree(PolyMeshType &mesh, const num_type faceIndex,
                         PolygonRegister *polygonRegister,
                         PatchLearner &learner,
                         const num_type minNfaces = 0,
                         const num_type maxNFaces = std::numeric_limits<num_type>::max(),
                         const CornersType cornersType = CONVEX,
                         const bool bindPolygon = false,
                         const bool collapsePolygon = true,
                         const bool keepInMemory = false) {
    // requirements
    vcg::tri::RequireFFAdjacency(mesh);
    vcg::tri::RequirePerFaceMark(mesh);
    vcg::tri::RequirePerVertexMark(mesh);

    // first clear this tree
    deAllocate();

    // general checks
    if (faceIndex >= (num_type)mesh.face.size()) {
      std::cout << "Wrong face index: " << faceIndex << ". Mesh has " << mesh.face.size() << " faces." << std::endl;
      return;
    }

    if (maxNFaces <= 0)
      return;

    if (vcg::tri::HasPerFaceFlags(mesh) && mesh.face[faceIndex].IsD())
      return;

    // some variables
    bool found = false;
    std::vector< std::pair<valence_type, num_type> > singularityVec;
    std::stack<PolygonNodePointer> stack;
    std::string key;

    // create a new polygon (_root node) with mesh.face[faceIndex].VN() sides and without parent
    try {
      _root = new PolygonNode(mesh, faceIndex);
    } catch(std::bad_alloc ba) {
      std::cout << "Bad allocation: " << ba.what() << std::endl;
      _root = 0;
    }
    if (_root == 0) {
      std::cout << "Error allocating memory for the root of a new tree. Tree skipped." << std::endl;
      return;
    }
    // update the maximum number of sides found
    _maxNSides = _root->numberOfSides();
    // filter
    // check if already registered, else registered it (it SHOULD not be already registered)
    if (polygonRegister) {
      PolygonRegister::BuildKey(mesh, _root->corners(), key);
      if (polygonRegister->isRegistered(key)) {
        deAllocate();
        return;
      }
    }
    // acquire the patch
    singularityVec.clear();
    AcquirePolygon(mesh, _root, learner, singularityVec, cornersType, bindPolygon, collapsePolygon);

    // count this first polygon
    if (minNfaces < 1)
      _numberOfPolygons = 1;

    // now build the rest of the tree by ENUMERATING the polygons
    stack.push(_root);
    // manage memory
    if (!keepInMemory)
      _root = 0;
    while (!stack.empty()) {
      // pick the top polygon from the stack
      _current = stack.top();
      // pop from the stack
      stack.pop();
      // for each side, try to expand
      for (n_corners_type i = 0; i < _current->children().size(); ++i) {
        // expand side i
        found = _current->expandPolygon(mesh, i, (cornersType == CONVEX), (cornersType & REGULAR), singularityVec, maxNFaces);
        // if found, register if not yet
        if (found && polygonRegister) {
          // check if already registered, else register it
          PolygonRegister::BuildKey(mesh, _current->child(i)->corners(), key);
          found = !polygonRegister->isRegistered(key);
          // if already registered, deallocate child at i
          if (!found) {
            delete _current->child(i);
            _current->child(i) = 0;
          }
        }
        // if a new convex polygon has been built and registered, push it on stack
        if (found) {
          if (_current->child(i)->numberOfFaces() >= minNfaces) {
            // update the maximum number of sides
            if (_current->child(i)->numberOfSides() > _maxNSides)
              _maxNSides = _current->child(i)->numberOfSides();
            // update the maximum number of singularities
            if (_current->child(i)->numberOfSingularities() > _maxNSingularities)
              _maxNSingularities = _current->child(i)->numberOfSingularities();
            // increment the number of polygons
            ++_numberOfPolygons;
            // acquire the patch
            AcquirePolygon(mesh, _current->child(i), learner, singularityVec, cornersType, bindPolygon, collapsePolygon);
          }
          // push on stack
          stack.push(_current->child(i));
        }
      }
      // manage memory
      if (!keepInMemory) {
        _current->detachFromChildren();
        delete _current;
        _current = 0;
      }
    }
    // always start from root when navigating
    _current = _root;
  }

  /**
   * @brief AcquirePolygon builds a Patch and lets a PatchLearner acquire it.
   * @param mesh The input mesh.
   * @param polygonNode The polygon to set as a Patch and acquire.
   * @param learner The convex patch acquisitor.
   * @param singularityVec A vector of pairs <valence_V,number_of_vertices_with_valence_V>.
   * @param cornersType The type of corners allowed.
   * @param bindPolygon true if the template patch must be linked to the polygon whre it's been found, false otherwise.
   * @param collapsePolygon true if the polygon must be simplified before acquisition, false otherwise.
   */
  static void AcquirePolygon(PolyMeshType &mesh, PolygonNodePointer polygonNode,
                             PatchLearner &learner,
                             const std::vector< std::pair<valence_type,num_type> > &singularityVec,
                             const CornersType cornersType = CONVEX,
                             const bool bindPolygon = false,
                             const bool collapsePolygon = true) {
    // check input
    if (polygonNode == 0)
      return;
    if (polygonNode->numberOfFaces() == 0)
      return;
    if (mesh.IsEmpty())
      return;

    // some variables
    Patch patch;
    std::string polygonStr;

    if (collapsePolygon) {
      PolyMeshType patchMesh;
      FaceSet allFaces;
      FaceSetIterator fIt;
      // get all faces
      polygonNode->allFaces(allFaces);
      // select faces in the mesh
      vcg::tri::UpdateFlags<PolyMeshType>::Clear(mesh);
      for (fIt = allFaces.begin(); fIt != allFaces.end(); ++fIt) {
        mesh.face[*fIt].SetS();
        for (int j = 0; j < mesh.face[*fIt].VN(); ++j)
          mesh.face[*fIt].V(j)->SetS();
      }
      // copy the patch into a new mesh
      vcg::tri::Append<PolyMeshType,PolyMeshType>::MeshCopy(patchMesh, mesh, true, true);
      vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(patchMesh);
      // clear selection
      for (fIt = allFaces.begin(); fIt != allFaces.end(); ++fIt) {
        mesh.face[*fIt].ClearS();
        for (int j = 0; j < mesh.face[*fIt].VN(); ++j)
          mesh.face[*fIt].V(j)->ClearS();
      }
      // collapse it
      vcg::tri::PolychordCollapse<PolyMeshType>::CollapseAllPolychords(patchMesh);
      vcg::tri::Allocator<PolyMeshType>::CompactFaceVector(patchMesh);
      vcg::tri::Allocator<PolyMeshType>::CompactVertexVector(patchMesh);
      vcg::tri::InitFaceIMark(patchMesh);
      vcg::tri::InitVertexIMark(patchMesh);
      vcg::tri::IMark(patchMesh) = 0;
      // find corners
      num_type nCorners = 0;
      PosType startCornerPos;
      startCornerPos.Set(&patchMesh.face[0], 0, patchMesh.face[0].V(0));
      std::deque<PosType> corners;
      nCorners = PolychordSupport<PolyMeshType>::FindPolygonCorners(patchMesh, corners, (cornersType == CONVEX), false, startCornerPos);
      // set patch
      if (!corners.empty()) {
        patch.setPatch(patchMesh, corners, nCorners, patchMesh.FN(), patchMesh.VN(), singularityVec);
      } else {
        std::cout << "Error. Corners not found for collapsed patch with " << polygonNode->numberOfCorners()
                  << " corners (returned " << nCorners << "). Not acquired." << std::endl;
        return;
      }
    } else {
      patch.setPatch(mesh, polygonNode->corners(), polygonNode->numberOfCorners(),
                     polygonNode->numberOfFaces(), polygonNode->numberOfVertices(),
                     singularityVec);
    }
    // acquire!!
    if (!patch.isNull() && (cornersType & CONVEX || patch.getNumberOfConcaveCorners() > 0)) {
      learner.acquirePatch(patch);
      if (bindPolygon) {
        PolygonDBInterface::PolygonCornersToString(mesh, polygonNode->corners(), polygonStr);
        learner.bindTopologyToPolygon(patch.getTopologyStr(), polygonStr, polygonNode->numberOfFaces());
      }
    }
  }

  /**
   * @brief numberOfPolygons gets how many polygons have been enumerated.
   * @return The number of enumerated polygons.
   */
  count_type numberOfPolygons() const {
    return _numberOfPolygons;
  }

  /**
   * @brief isEmpty returns if this tree is empty.
   * @return true if empty, false otherwise.
   */
  bool isEmpty() const {
    return _root == 0;
  }

  /**
   * @brief selectRoot sets the current node to the root.
   */
  void selectRoot() {
    _current = _root;
  }

  /**
   * @brief maxNumberOfSides gives the maximum number of sides that any enumerated polygon has.
   * @return The maximum number of sides of any enumerated polygon.
   */
  n_corners_type maxNumberOfSides() const {
    return _maxNSides;
  }

  /**
   * @brief numberOfCorners gives the number of corners of the current polygon, which may be 0.
   * @return The number of corners of the current polygon.
   */
  n_corners_type numberOfCorners() const {
    if (_current != 0)
      return _current->numberOfCorners();
    return 0;
  }

  /**
   * @brief numberOfSides gives the number of sides of the current polygon, which is always greater than 0.
   * @return The number of sides of the current polygon.
   */
  n_corners_type numberOfSides() const {
    if (_current != 0)
      return _current->numberOfSides();
    return 0;
  }

  /**
   * @brief maxNumberOfSingularities gives the maximum number of singularities inside any enumerated polygon.
   * @return The maximum number of singularities inside any enumerated polygon.
   */
  num_type maxNumberOfSingularities() const {
    return _maxNSingularities;
  }

  /**
   * @brief corners gives the corners of the current polygon as a vector of vcg::face::Pos.
   * @param corners An output vector of vcg::face::Pos.
   * @return true if the current polygon exists, false otherwise.
   */
  bool corners(std::deque<PosType> &corners) const {
    if (_current != 0) {
      _current->corners(corners);
      return true;
    }
    corners.clear();
    return false;
  }

  /**
   * @brief faces gives the set of (the pointers of) the faces composing the current polygon.
   * @param faces A set of faces. See PolygonNode::FaceSet.
   * @return true if the current polygon exists, false otherwise.
   */
  bool faces(FaceSet &faces) const {
    if (_current != 0) {
      _current->allFaces(faces);
      return true;
    }
    return false;
  }

  /**
   * @brief isRoot checks if the current node is the root.
   * @return true if it is the root, false otherwise.
   */
  bool isRoot() const {
    return _current == _root;
  }

  /**
   * @brief isParentNull checks if the parent of the current node is NULL.
   * @return true if it is NULL, false otherwise.
   */
  bool isParentNull() const {
    return _current == 0 || _current->parent() == 0;
  }
  /**
   * @brief selectParent sets the current node to its parent, if not null.
   * @return true if successful, false otherwise.
   */
  bool selectParent() {
    if (_current != 0 && _current->parent() != 0) {
      _current = _current->parent();
      return true;
    }
    return false;
  }

  /**
   * @brief parentSide gives the index of the parent (polygon) side this polygon has been expanded from.
   * @return The parent side index this polygon has been expanded from.
   */
  n_corners_type parentSide() const {
    if (_current != 0)
      return _current->parentSide();
    return std::numeric_limits<n_corners_type>::max();
  }

  /**
   * @brief isChildNull checks if the child at side is NULL.
   * @param side The index of the side the child could have been expanded from.
   * @return true if it is NULL, false otherwise.
   */
  bool isChildNull(const n_corners_type side) const {
    return _current == 0 || side >= _current->numberOfSides() || _current->child(side) == 0;
  }

  /**
   * @brief selectChild sets the current node to the child at side "side".
   * @param side The index of the side the child could have been expanded from.
   * @return true if that child exists, false otherwise.
   */
  bool selectChild(const n_corners_type side) {
    if (_current != 0 && side < _current->numberOfSides() && _current->child(side) != 0) {
      _current = _current->child(side);
      return true;
    }
    return false;
  }

private:
  /**
   * @brief copyPolygonTree builds a new PolygonTree by copying data from tree.
   * @param tree The PolygonTree to be copied.
   */
  void copyPolygonTree(const PolygonTree &tree) {
    // first de-allocate this tree
    deAllocate();
    // copy data
    _root = 0;
    _current = 0;
    _numberOfPolygons = tree._numberOfPolygons;
    _maxNSides = tree._maxNSides;
    _maxNSingularities = tree._maxNSingularities;
    // copy nodes
    if (tree._root == 0)
      return;
    try {
      _root = new PolygonNode(*tree._root);
    } catch(std::bad_alloc ba) {
      std::cout << "Bad allocation: " << ba.what() << " - Copy skipped." << std::endl;
      deAllocate();
      return;
    }
    if (_root == 0) {
      std::cout << "Error allocating memory for a new node. Copy skipped." << std::endl;
      deAllocate();
      return;
    }
    std::stack< std::pair<PolygonNodePointer, PolygonNodePointer> > stack;
    PolygonNodePointer rhs_current = 0;
    stack.push(std::make_pair(_root, tree._root));
    while (!stack.empty()) {
      _current = stack.top().first;
      rhs_current = stack.top().second;
      stack.pop();
      for (num_type i = 0; i < rhs_current->numberOfSides(); ++i)
        if (rhs_current->child(i) != 0) {
          try {
            _current->child(i) = new PolygonNode(*rhs_current->child(i));
          } catch(std::bad_alloc ba) {
            std::cout << "Bad allocation: " << ba.what() << " - Copy skipped." << std::endl;
            deAllocate();
            return;
          }
          if (_current->child(i) == 0) {
            std::cout << "Error allocating memory for a new node. Copy skipped." << std::endl;
            deAllocate();
            return;
          }
          stack.push(std::make_pair(_current->child(i), rhs_current->child(i)));
        }
    }
    _current = _root;
  }

  PolygonNodePointer      _root;                ///< Since this is a tree, we need a root node
  PolygonNodePointer      _current;             ///< Node pointer used to traverse the tree.
  count_type              _numberOfPolygons;    ///< Number of polygons.
  n_corners_type          _maxNSides;           ///< The maximum number of sides of any enumerated polygon.
  num_type                _maxNSingularities;   ///< The maximum number of singular vertexes found in all the polygons.
};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // POLYGON_TREE_H
