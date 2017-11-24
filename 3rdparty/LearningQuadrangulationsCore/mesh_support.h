#ifndef MESH_SUPPORT_H
#define MESH_SUPPORT_H

#include <deque>
#include <vector>
#include <list>
#include <utility>
#if __cplusplus >= 201103L || _MSC_VER >= 1800
  #include <unordered_set>
  #include <unordered_map>
#else
  #include <set>
  #include <map>
#endif
#include <map>
#include <limits>
#include <vcg/complex/complex.h>
#include <vcg/simplex/face/pos.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include "default_types.h"

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief The MeshSupport class provides some utility tools for meshes.
 */
template < typename PolyMeshType >
class MeshSupport {
public:
  typedef typename PolyMeshType::VertexPointer                        VertexPointer;
  typedef typename PolyMeshType::FaceType                             FaceType;
  typedef typename PolyMeshType::FacePointer                          FacePointer;
  typedef vcg::face::Pos<FaceType>                                    PosType;
#if __cplusplus >= 201103L || _MSC_VER >= 1800
  typedef std::unordered_set<num_type>                                FaceSet;
#else
  typedef std::set<num_type>                                          FaceSet;
#endif
  typedef typename FaceSet::iterator                                  FaceSetIterator;

  typedef std::pair<size_t, int>                                      HalfEdge;
#if __cplusplus >= 201103L || _MSC_VER >= 1800
  struct HalfEdgeHash {
    size_t operator() (const HalfEdge &he) const {
      return he.first << 2 | (std::size_t)(he.second & 3);
    }
  };
  typedef std::unordered_set<HalfEdge, HalfEdgeHash>                  HalfEdgeSet;
#else
  typedef std::set<HalfEdge>                                          HalfEdgeSet;
#endif
  typedef typename HalfEdgeSet::iterator                              HalfEdgeSetIterator;

  /**
   * @brief IncrementMark increments and returns the next mark to use (modulo std::numeric_limits<int>::max()).
   * @param mesh The input mesh.
   * @return The next mark to use.
   */
  static int IncrementMark(PolyMeshType &mesh) {
    if (vcg::tri::IMark(mesh) == std::numeric_limits<int>::max()) {
      vcg::tri::InitFaceIMark(mesh);
      vcg::tri::InitVertexIMark(mesh);
      vcg::tri::IMark(mesh) = 0;
    }
    vcg::tri::UnMarkAll(mesh);
    return vcg::tri::IMark(mesh);
  }

  /**
   * @brief MarkSubmeshFromBoundary marks and collects all the faces inside a given boundary bounded by corners.
   * @note The vertices are also marked.
   * @param mesh The inpunt mesh.
   * @param corners The corners of the boundary (counterclockwisely).
   * @param markedFaces The output collected faces.
   */
  static void MarkSubmeshFromBoundary(PolyMeshType &mesh,
                                      const std::deque<PosType> &corners,
                                      FaceSet &markedFaces) {
    markedFaces.clear();

    // check input
    vcg::tri::RequirePerFaceMark(mesh);
    vcg::tri::RequirePerVertexMark(mesh);
    if (mesh.IsEmpty() || corners.empty())
      return;

    // first: increment the mark
    IncrementMark(mesh);

    // second: mark the boundary
    std::queue<FacePointer> fQueue;
    HalfEdgeSet boundarySet;
    PosType runPos = corners.front();
    for (size_t c = 0; c < corners.size(); ++c) {
      do {
        // set this half edge as border
        boundarySet.insert(HalfEdge(vcg::tri::Index(mesh, runPos.F()), runPos.E()));

        // collect this face
        markedFaces.insert((num_type)vcg::tri::Index(mesh, runPos.F()));
        // enqueue this face
        fQueue.push(runPos.F());

        // mark face and vertices
        vcg::tri::Mark(mesh, runPos.F());
        for (int v = 0; v < runPos.F()->VN(); ++v)
          vcg::tri::Mark(mesh, runPos.F()->V(v));

        // go to next edge
        runPos.FlipV();
        runPos.FlipE();
        if (runPos.V() == corners[(c + 1) % corners.size()].V()) {
          while (runPos != corners[(c + 1) % corners.size()]) {
            runPos.FlipF();
            runPos.FlipE();
          }
        } else {
          assert(!runPos.IsBorder());
          runPos.FlipF();
          runPos.FlipE();
        }
      } while (runPos != corners[(c + 1) % corners.size()]);
    }

    // third: explore inside boundary and mark
    FacePointer f = 0;
    while (!fQueue.empty()) {
      f = fQueue.front();
      fQueue.pop();
      for (int e = 0; e < f->VN(); ++e)
        if (!vcg::face::IsBorder(*f, e) &&
            boundarySet.find(HalfEdge(vcg::tri::Index(mesh, f), e)) == boundarySet.end() &&
            !vcg::tri::IsMarked(mesh, f->FFp(e))) {
          // mark face and vertices
          vcg::tri::Mark(mesh, f->FFp(e));
          for (int v = 0; v < f->FFp(e)->VN(); ++v)
            if (!vcg::tri::IsMarked(mesh, f->FFp(e)->V(v)))
              vcg::tri::Mark(mesh, f->FFp(e)->V(v));
          // collect this face
          markedFaces.insert((num_type)vcg::tri::Index(mesh, f->FFp(e)));
          // enqueue this face
          fQueue.push(f->FFp(e));
        }
    }
  }

  /**
   * @brief CountInternalSingularities counts the number of singular vertices (with valence not equal to 4) and,
   * for each singular valence, counts how many singular vertices with such a valence there are inside the mesh,
   * i.e. without counting border vertices.
   * @param mesh The input mesh.
   * @param faces A set of faces the interested piece of mesh is composed by.
   * @param singularityVec The output vector of pairs <valence_type,num_type>.
   * @return The total number of singularities.
   */
  static num_type CountInternalSingularities(const PolyMeshType &mesh,
                                             FaceSet &faces,
                                             std::vector< std::pair<valence_type, num_type> > &singularityVec) {
    FaceSetIterator fIt;
#if __cplusplus >= 201103L || _MSC_VER >= 1800
    std::unordered_map<VertexPointer, valence_type> valences;
    typename std::unordered_map<VertexPointer, valence_type>::iterator valIt;
#else
    std::map<VertexPointer, valence_type> valences;
    typename std::map<VertexPointer, valence_type>::iterator valIt;
#endif
    std::map<valence_type, num_type> singularities;
    std::map<valence_type, num_type>::iterator sIt;
    num_type nSingularities = 0;
    // update the number of singularities
    for (fIt = faces.begin(); fIt != faces.end(); ++fIt)
      // for each vertex/edge of the added face
      for (int j = 0; j < mesh.face[*fIt].VN(); ++j) {
        valIt = valences.find(mesh.face[*fIt].V(j));
        if (valIt != valences.end())
          ++valIt->second;
        else
          valences.insert(std::pair<VertexPointer, valence_type>(mesh.face[*fIt].V(j), 1));
      }
    // un-map the vertexes lying on the boundary of (the new part of) the polygon
    for (fIt = faces.begin(); fIt != faces.end(); ++fIt)
      // for each vertex/edge of the added face
      for (int j = 0; j < mesh.face[*fIt].VN(); ++j)
        if (vcg::face::IsBorder(mesh.face[*fIt], j) || !vcg::tri::IsMarked(mesh, mesh.face[*fIt].cFFp(j))) {
          valIt = valences.find(mesh.face[*fIt].V(j));
          if (valIt != valences.end())
            valences.erase(valIt);
          valIt = valences.find(mesh.face[*fIt].V((j+1)%mesh.face[*fIt].VN()));
          if (valIt != valences.end())
            valences.erase(valIt);
        }
    // update the number of singularities
    for (valIt = valences.begin(); valIt != valences.end(); ++valIt)
      if (valIt->second != 4) {
        ++nSingularities;
        sIt = singularities.find(valIt->second);
        if (sIt != singularities.end())
          ++sIt->second;
        else
          singularities.insert(std::pair<valence_type, num_type>(valIt->second, 1));
      }
    // build the output vector
    singularityVec.resize(singularities.size());
    size_t ind = 0;
    for (sIt = singularities.begin(); sIt != singularities.end(); ++sIt, ++ind)
      singularityVec[ind] = *sIt;

    return nSingularities;
  }

  /**
   * @brief Open Loads a polygonal mesh from file.
   * @param mesh The mesh object into which to extract the mesh.
   * @param filename The filename where the mesh is stored.
   * @return An result code.
   */
  static int OpenMesh(PolyMeshType &mesh, const char *filename) {
    // try to load the mesh from the file
    int err = vcg::tri::io::Importer<PolyMeshType>::Open(mesh, filename);

    // check if successfully loaded
    if (err == 0) {
      // update bounding box
      vcg::tri::UpdateBounding<PolyMeshType>::Box(mesh);
      // update topology
      vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(mesh);
      // update normals
      vcg::tri::UpdateNormal<PolyMeshType>::PerPolygonalFaceNormalized(mesh);
      // update flags
      vcg::tri::UpdateFlags<PolyMeshType>::Clear(mesh);
      // initialize mark
      vcg::tri::InitFaceIMark(mesh);
      vcg::tri::InitVertexIMark(mesh);
      vcg::tri::IMark(mesh) = 0;
    }

    return err;
  }

  /**
   * @brief SaveMesh Writes a polygonal mesh into a file.
   * @param mesh The mesh to write.
   * @param filename The filename where to store the mesh.
   * @return An result code.
   */
  static int SaveMesh(PolyMeshType &mesh, const char *filename) {
    // try to write the mesh into the file
    return vcg::tri::io::Exporter<PolyMeshType>::Save(mesh, filename);
  }
};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // MESH_SUPPORT_H
