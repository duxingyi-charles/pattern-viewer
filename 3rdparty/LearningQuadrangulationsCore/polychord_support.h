#ifndef POLYCHORD_SUPPORT_H
#define POLYCHORD_SUPPORT_H

#include "mesh_support.h"
#include <vcg/complex/algorithms/polygon_polychord_collapse.h>

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief The FlowTopology struct represents a type of flow which starts and ends at sides of a Patch,
 * or is a cycle. It may also have a number of self intersections.
 */
struct FlowTopology {
  n_corners_type    startSide;              ///< The index of the starting side.
  n_corners_type    endSide;                ///< The index of the ending side.

  /**
   * @brief FlowTopology Default constructor.
   */
  FlowTopology() : startSide(std::numeric_limits<n_corners_type>::max()),
                   endSide(std::numeric_limits<n_corners_type>::max()) { }

  /**
   * @brief operator < gives the strick weak order of two FlowTopology objects.
   * @param flowTopology
   * @return
   */
  bool operator < (const FlowTopology &flowTopology) const {
    if (startSide != flowTopology.startSide)
      return startSide < flowTopology.startSide;
    return endSide < flowTopology.endSide;
  }

  /**
   * @brief operator > determines the relational order between this and flowTopology.
   * @param flowTopology The flowTopology to compare with.
   * @return true if this > flowTopology, false otherwise.
   */
  bool operator > (const FlowTopology &flowTopology) const {
    return flowTopology < *this;
  }

  /**
   * @brief operator <= determines the relational order between this and flowTopology.
   * @param flowTopology The flowTopology to compare with.
   * @return true if this <= flowTopology, false otherwise.
   */
  bool operator <= (const FlowTopology &flowTopology) const {
    return !(flowTopology < *this);
  }

  /**
   * @brief operator >= determines the relational order between this and flowTopology.
   * @param flowTopology The flowTopology to compare with.
   * @return true if this >= flowTopology, false otherwise.
   */
  bool operator >= (const FlowTopology &flowTopology) const {
    return !(*this < flowTopology);
  }

  /**
   * @brief operator == determines the relational order between this and flowTopology.
   * @param flowTopology The flowTopology to compare with.
   * @return true if this == flowTopology, false otherwise.
   */
  bool operator == (const FlowTopology &flowTopology) const {
    return !(*this < flowTopology || flowTopology < *this);
  }

  /**
   * @brief operator != determines the relational order between this and flowTopology.
   * @param flowTopology The flowTopology to compare with.
   * @return true if this == flowTopology, false otherwise.
   */
  bool operator != (const FlowTopology &flowTopology) const {
    return *this < flowTopology || flowTopology < *this;
  }

  /**
   * @brief swap inverts startSide with endSide.
   */
  void swap() {
    std::swap(startSide, endSide);
  }
};

/**
 * @brief The FlowPoints struct represents a flow as a sequence of points and its
 * parametric version given by arc lengths of each point with respect to the first one.
 */
template < typename MeshType >
struct FlowPoints {
  typedef typename MeshType::ScalarType             ScalarType;
  typedef std::deque<ScalarType>                    ScalarContainer;
  typedef typename ScalarContainer::iterator        ScalarIterator;
  typedef typename ScalarContainer::const_iterator  ScalarConstIterator;
  typedef typename MeshType::CoordType              PointType;
  typedef std::deque<PointType>                     PointContainer;
  typedef typename PointContainer::iterator         PointIterator;
  typedef typename PointContainer::const_iterator   PointConstIterator;

  num_type            startEdge;      ///< Index of the edge on the start side.
  num_type            endEdge;        ///< Index of the edge on the end side.
  PointContainer      points;         ///< The sequence of points.
  ScalarContainer     arcLengths;     ///< Their arc lengths.
  PointContainer      samplePoints;   ///< The sample points, ideally uniformly distanced in parametric space.

  /**
   * @brief FlowPoints Default constructor.
   */
  FlowPoints() : startEdge(-1), endEdge(-1) { }

  /**
   * @brief FlowPoints Copy constructor.
   * @param flowPoints The FlowPoints to copy.
   */
  FlowPoints(const FlowPoints &flowPoints) : startEdge(flowPoints.startEdge),
                                             endEdge(flowPoints.endEdge),
                                             points(flowPoints.points),
                                             arcLengths(flowPoints.arcLengths),
                                             samplePoints(flowPoints.samplePoints) { }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief FlowPoints Move constructor.
   * @param flowPoints The FlowPoints to move.
   */
  FlowPoints(FlowPoints&& flowPoints) : startEdge(flowPoints.startEdge),
                                        endEdge(flowPoints.endEdge),
                                        points(std::move(flowPoints.points)),
                                        arcLengths(std::move(flowPoints.arcLengths)),
                                        samplePoints(std::move(flowPoints.samplePoints)) { }
#endif

  /**
   * @brief operator = Copy assignment.
   * @param flowPoints The FlowPoints to copy.
   * @return A reference to this.
   */
  FlowPoints & operator = (const FlowPoints &flowPoints) {
    if (this != &flowPoints) {
      startEdge = flowPoints.startEdge;
      endEdge = flowPoints.endEdge;
      points = flowPoints.points;
      arcLengths = flowPoints.arcLengths;
      samplePoints = flowPoints.samplePoints;
    }
    return *this;
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief operator = Move assignment.
   * @param flowPoints The FlowPoints to move.
   * @return A reference to this.
   */
  FlowPoints & operator = (FlowPoints&& flowPoints) {
    if (this != &flowPoints) {
      startEdge = flowPoints.startEdge;
      endEdge = flowPoints.endEdge;
      points = std::move(flowPoints.points);
      arcLengths = std::move(flowPoints.arcLengths);
      samplePoints = std::move(flowPoints.samplePoints);
    }
    return *this;
  }
#endif

  /**
   * @brief operator < gives the strick weak order of two FlowPoints objects.
   * @param flow
   * @return
   */
  bool operator < (const FlowPoints &flow) const {
    if (startEdge != flow.startEdge)
      return startEdge < flow.startEdge;
    if (endEdge != flow.endEdge)
      return endEdge < flow.endEdge;
    if (points != flow.points)
      return points < flow.points;
    if (arcLengths != flow.arcLengths)
      return arcLengths < flow.arcLengths;
    return samplePoints < flow.samplePoints;
  }

  /**
   * @brief operator > determines the relational order between this and flow.
   * @param flow The flow to compare with.
   * @return true if this > flow, false otherwise.
   */
  bool operator > (const FlowPoints &flow) const {
    return flow < *this;
  }

  /**
   * @brief operator <= determines the relational order between this and flow.
   * @param flow The flow to compare with.
   * @return true if this <= flow, false otherwise.
   */
  bool operator <= (const FlowPoints &flow) const {
    return !(flow < *this);
  }

  /**
   * @brief operator >= determines the relational order between this and flow.
   * @param flow The flow to compare with.
   * @return true if this >= flow, false otherwise.
   */
  bool operator >= (const FlowPoints &flow) const {
    return !(*this < flow);
  }

  /**
   * @brief operator == determines the relational order between this and flow.
   * @param flow The flow to compare with.
   * @return true if this == flow, false otherwise.
   */
  bool operator == (const FlowPoints &flow) const {
    return !(*this < flow || flow < *this);
  }

  /**
   * @brief operator != determines the relational order between this and flow.
   * @param flow The flow to compare with.
   * @return true if this == flow, false otherwise.
   */
  bool operator != (const FlowPoints &flow) const {
    return *this < flow || flow < *this;
  }

  /**
   * @brief updateArcLengths computes the parametric space of the points.
   */
  void updateArcLengths() {
    ComputeArcLengths(points, arcLengths);
  }

  /**
   * @brief computeSamples computes samples equally distanced in parametric space.
   * @param nSamples The total number of samples. It must be greater than 1.
   */
  void computeSamples(const count_type nSamples) {
    samplePoints.clear();

    if (points.empty() || points.size() != arcLengths.size()) {
      std::cout << "Points and arc lengths not coherent. Impossible to sample." << std::endl;
      return;
    }

    if (nSamples < 2) {
      std::cout << "The number of samples must be greater than 1. Not sampled." << std::endl;
      return;
    }

    SampleFlowPoints(points, arcLengths.back() / ScalarType(nSamples - 1), samplePoints);
    while (samplePoints.size() > nSamples)
      samplePoints.pop_back();
    assert(samplePoints.size() == nSamples);
  }

  /**
   * @brief invertPoints inverts point and samplePoints so that they are ordered
   * from the end to the beginning.
   * @note After this, always call updateArcLengths().
   */
  void invertPoints() {
    std::swap(startEdge, endEdge);
    PointContainer tmpPoints(points);
    points.assign(tmpPoints.rbegin(), tmpPoints.rend());
    tmpPoints = samplePoints;
    samplePoints.assign(tmpPoints.rbegin(), tmpPoints.rend());
  }

  /**
   * @brief ComputeArcLengths computes the arc lengths of a sequence of points.
   * @param pp Theinput points.
   * @param al The output arc lengths.
   */
  static void ComputeArcLengths(const PointContainer &pp, ScalarContainer &al) {
    al.clear();
    PointType lastPoint;
    PointConstIterator pIt = pp.begin();
    if (pIt != pp.end()) {
      lastPoint = *pIt;
      al.push_back(0);
      ++pIt;
    }
    for (; pIt != pp.end(); lastPoint = *pIt, ++pIt)
      al.push_back(al.back() + (*pIt - lastPoint).Norm());
  }

  /**
   * @brief SampleFlowPoints computes samples at stepLength distance from each other.
   * @param points The input points.
   * @param stepLength The distance between each point.
   * @param samples The output sample points.
   */
  static void SampleFlowPoints(const PointContainer &points, const ScalarType stepLength, PointContainer &samples) {
    samples.clear();
    if (stepLength == ScalarType(0))
      return;
    if (points.empty())
      return;
    ScalarType currentLineLength = ScalarType(0), currentArcLength = ScalarType(0), l = ScalarType(0);
    PointType currentDir(0, 0, 0);
    for (size_t p = 0; p < points.size() - 1; ++p) {
      // set current direction
      currentDir = points[p+1] - points[p];
      currentLineLength = currentDir.Norm();
      while (currentArcLength < currentLineLength) {
        l = currentArcLength / currentLineLength;
        // pick a new sample
        samples.push_back(points[p] * (1.0 - l) + points[p+1] * l);
        currentArcLength += stepLength;
      }
      currentArcLength -= currentLineLength;
    }
    samples.push_back(points.back());
  }
};

/**
 * @brief The PolychordPoints struct represents a polychord starting at startPos as a sequence of
 * points and its parametric version given by arc lengths of each point with respect to the first one.
 */
template < typename PolyMeshType >
struct PolychordPoints : public FlowPoints<PolyMeshType> {
  typedef typename PolyMeshType::FaceType                     FaceType;
  typedef pl::MeshSupport<PolyMeshType>                       MeshSupport;
  typedef typename MeshSupport::PosType                       PosType;
  typedef typename FlowPoints<PolyMeshType>::ScalarContainer  ScalarContainer;
  using FlowPoints<PolyMeshType>::startEdge;
  using FlowPoints<PolyMeshType>::endEdge;
  using FlowPoints<PolyMeshType>::points;
  using FlowPoints<PolyMeshType>::arcLengths;
  using FlowPoints<PolyMeshType>::samplePoints;
  using FlowPoints<PolyMeshType>::computeSamples;
  using FlowPoints<PolyMeshType>::updateArcLengths;

  PosType                     startPos;     ///< The starting Pos of the polychord.
  ScalarContainer             distances;    ///< Distance for each sample point.

  /**
   * @brief PolychordPoints Default constructor.
   */
  PolychordPoints() : FlowPoints<PolyMeshType>() { }

  /**
   * @brief PolychordPoints Copy constructor.
   * @param polychordPoints The PolychordPoints to copy.
   */
  PolychordPoints(const PolychordPoints &polychordPoints) : FlowPoints<PolyMeshType>(polychordPoints),
                                                            startPos(polychordPoints.startPos),
                                                            distances(polychordPoints.distances) { }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief PolychordPoints Move constructor.
   * @param polychordPoints The PolychordPoints to move.
   */
  PolychordPoints(PolychordPoints&& polychordPoints) : FlowPoints<PolyMeshType>(std::move(polychordPoints)),
                                                       startPos(polychordPoints.startPos),
                                                       distances(std::move(polychordPoints.distances)) { }
#endif

  /**
   * @brief operator = Copy assignment.
   * @param polychordPoints The PolychordPoints to copy.
   * @return A reference to this.
   */
  PolychordPoints & operator = (const PolychordPoints &polychordPoints) {
    if (this != &polychordPoints) {
      FlowPoints<PolyMeshType>::operator =(polychordPoints);
      startPos = polychordPoints.startPos;
      distances = polychordPoints.distances;
    }
    return *this;
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief operator = Move assignment.
   * @param polychordPoints The PolychordPoints to move.
   * @return A reference to this.
   */
  PolychordPoints & operator = (PolychordPoints&& polychordPoints) {
    if (this != &polychordPoints) {
      FlowPoints<PolyMeshType>::operator =(std::move(polychordPoints));
      startPos = polychordPoints.startPos;
      distances = std::move(polychordPoints.distances);
    }
    return *this;
  }
#endif

  /**
   * @brief operator < gives the strick weak order of two PolychordPoints objects.
   * @param polychord
   * @return
   */
  bool operator < (const PolychordPoints &polychord) const {
    if (points != polychord.points)
      return points < polychord.points;
    if (arcLengths != polychord.arcLengths)
      return arcLengths < polychord.arcLengths;
    if (samplePoints != polychord.samplePoints)
      return samplePoints < polychord.samplePoints;
    return distances < polychord.distances;
  }

  /**
   * @brief operator > determines the relational order between this and polychord.
   * @param polychord The polychord to compare with.
   * @return true if this > polychord, false otherwise.
   */
  bool operator > (const PolychordPoints &polychord) const {
    return polychord < *this;
  }

  /**
   * @brief operator <= determines the relational order between this and polychord.
   * @param polychord The polychord to compare with.
   * @return true if this <= polychord, false otherwise.
   */
  bool operator <= (const PolychordPoints &polychord) const {
    return !(polychord < *this);
  }

  /**
   * @brief operator >= determines the relational order between this and polychord.
   * @param polychord The polychord to compare with.
   * @return true if this >= polychord, false otherwise.
   */
  bool operator >= (const PolychordPoints &polychord) const {
    return !(*this < polychord);
  }

  /**
   * @brief operator == determines the relational order between this and polychord.
   * @param polychord The polychord to compare with.
   * @return true if this == polychord, false otherwise.
   */
  bool operator == (const PolychordPoints &polychord) const {
    return !(*this < polychord || polychord < *this);
  }

  /**
   * @brief operator != determines the relational order between this and polychord.
   * @param polychord The polychord to compare with.
   * @return true if this == polychord, false otherwise.
   */
  bool operator != (const PolychordPoints &polychord) const {
    return *this < polychord || polychord < *this;
  }

  /**
   * @brief polychordToPoints converts the polychord starting at startPos into a sequence of
   * points, picking them at the center of each crossed edge and at the center of each face.
   */
  void polychordToPoints() {
    distances.clear();
    samplePoints.clear();
    arcLengths.clear();
    points.clear();
    if (startPos.IsNull())
      return;
    PosType runPos = startPos;
    do {
      points.push_back((runPos.V()->P() + runPos.VFlip()->P()) / 2.0);
      runPos.FlipE();
      runPos.FlipV();
      runPos.FlipE();
      points.push_back((points.back() + (runPos.V()->P() + runPos.VFlip()->P()) / 2.0) / 2.0);
      runPos.FlipF();
    } while (runPos != startPos && !runPos.IsBorder());
    points.push_back((runPos.V()->P() + runPos.VFlip()->P()) / 2.0);
    updateArcLengths();
  }

  /**
   * @brief computePolychordStartEndEdges finds the index of the edge on the start side and the one on the end side.
   */
  void computePolychordStartEndEdges() {
    // if null or loop, set -1
    if (startPos.IsNull() || !startPos.IsBorder()) {
      startEdge = -1;
      endEdge = -1;
      return;
    }

    // find start edge
    vcg::face::Pos<FaceType> runPos = startPos;
    startEdge = 0;
    do {
      runPos.FlipE();
      if (!runPos.IsBorder()) {
        runPos.FlipF();
        runPos.FlipE();
        if (!runPos.IsBorder())
          break;
        else {
          runPos.FlipV();
          ++startEdge;
        }
      }
    } while (!runPos.IsBorder());

    // for the end edge, first go to the other side of the polychord
    runPos = startPos;
    runPos.FlipV();
    do {
      runPos.FlipE();
      runPos.FlipV();
      runPos.FlipE();
      runPos.FlipF();
    } while (!runPos.IsBorder());
    // then find end edge
    endEdge = 0;
    do {
      runPos.FlipE();
      if (!runPos.IsBorder()) {
        runPos.FlipF();
        runPos.FlipE();
        if (!runPos.IsBorder())
          break;
        else {
          runPos.FlipV();
          ++endEdge;
        }
      }
    } while (!runPos.IsBorder());
  }
};

/**
 * @brief The PolychordSupport class provides methods concerning polychords.
 */
template < typename PolyMeshType >
class PolychordSupport {
public:
  typedef typename PolyMeshType::FaceType             FaceType;
  typedef typename PolyMeshType::FacePointer          FacePointer;
  typedef pl::FlowPoints<PolyMeshType>                FlowPoints;
  typedef typename FlowPoints::ScalarType             ScalarType;
  typedef typename FlowPoints::ScalarContainer        ScalarContainer;
  typedef typename FlowPoints::ScalarIterator         ScalarIterator;
  typedef typename FlowPoints::ScalarConstIterator    ScalarConstIterator;
  typedef typename FlowPoints::PointType              PointType;
  typedef typename FlowPoints::PointContainer         PointContainer;
  typedef typename FlowPoints::PointIterator          PointIterator;
  typedef typename FlowPoints::PointConstIterator     PointConstIterator;
  typedef std::deque<FlowPoints>                      FlowContainer;
  typedef typename FlowContainer::iterator            FlowContainerIterator;
  typedef typename FlowContainer::const_iterator      FlowContainerConstIterator;
  typedef std::map<FlowTopology, FlowContainer>       FlowMap;
  typedef typename FlowMap::iterator                  FlowMapIterator;
  typedef typename FlowMap::const_iterator            FlowMapConstIterator;
  typedef pl::PolychordPoints<PolyMeshType>           PolychordPoints;
  typedef typename PolychordPoints::MeshSupport       MeshSupport;
  typedef typename PolychordPoints::PosType           PosType;
  typedef std::deque<PolychordPoints>                 PolychordContainer;
  typedef typename PolychordContainer::iterator       PolychordContainerIterator;
  typedef typename PolychordContainer::const_iterator PolychordContainerConstIterator;
  typedef std::map<FlowTopology, PolychordContainer>  PolychordMap;
  typedef typename PolychordMap::iterator             PolychordMapIterator;
  typedef typename PolychordMap::const_iterator       PolychordMapConstIterator;

  /**
   * @brief FindPolygonCorners finds all the corners of a (disk-like) mesh.
   * @note The mesh must be homeomorphic to a disk, but no check is made.
   * @param mesh The input mesh.
   * @param corners Vector of found corners to return.
   * @param convexOnly true if convex only corners must be found, false otherwise.
   * @param collectRegular true if all border vertices must be collected, false otherwise.
   * @param startPos The starting position (if not on boundary, the first border edge is searched).
   * @param endPos The ending position (not necessary, default null).
   * @return The number of corners found, -1 if convexOnly==true and a concave corner has been found (and pushed back).
   */
  static num_type FindPolygonCorners(PolyMeshType &mesh,
                                     std::deque<PosType> &corners,
                                     const bool convexOnly,
                                     const bool collectRegular,
                                     const vcg::face::Pos<FaceType> &startPos,
                                     vcg::face::Pos<FaceType> endPos = PosType()) {
    PosType tmpPos, runPos = startPos;

    // empty the list of corners
    corners.clear();

    // check input
    if (startPos.IsNull()) {
      std::cout << "Null start Pos. Impossible to find corners." << std::endl;
      return -1;
    }

    // find a (starting) polygon border edge
    while (!runPos.IsBorder() && vcg::tri::IsMarked(mesh, runPos.FFlip())) {
      // not on a polygon border edge, so go on next edge in _current face (counterclockwise)
      runPos.FlipV();
      runPos.FlipE();
      // if not on border yet, go on the adjacent face
      if (!runPos.IsBorder() && vcg::tri::IsMarked(mesh, runPos.FFlip()) && runPos != startPos) {
        runPos.FlipF();
        runPos.FlipE();
      }
      // if returned to the starting position, exit (should try another way)
      if (runPos == startPos)
        break;
    }
    // if it's a loop, try another way
    if (!runPos.IsBorder() && vcg::tri::IsMarked(mesh, runPos.FFlip())) {
      // not on a polygon border edge, so go on next edge in _current face (counterclockwise)
      runPos.FlipV();
      runPos.FlipE();
      tmpPos = runPos;
      // find a (starting) polygon border edge (on another way)
      while (!runPos.IsBorder() && vcg::tri::IsMarked(mesh, runPos.FFlip())) {
        // not on a polygon border edge, so go on next edge in _current face (counterclockwise)
        runPos.FlipV();
        runPos.FlipE();
        // if not on border yet, go on the adjacent face
        if (!runPos.IsBorder() && vcg::tri::IsMarked(mesh, runPos.FFlip()) && runPos != tmpPos) {
          runPos.FlipF();
          runPos.FlipE();
        }
        // if returned to the starting position, exit (no way)
        if (runPos == tmpPos)
          break;
      }
    }

    if (!runPos.IsBorder() && vcg::tri::IsMarked(mesh, runPos.FFlip())) {
      runPos.SetNull();
      for (size_t f = 0; f < mesh.face.size(); ++f)
        if (vcg::tri::IsMarked(mesh, &mesh.face[f])) {
          int v = 0;
          for (; v < mesh.face[f].VN(); ++v)
            if (vcg::face::IsBorder(mesh.face[f], v) || !vcg::tri::IsMarked(mesh, mesh.face[f].cFFp(v))) {
              runPos.Set(&mesh.face[f], v, mesh.face[f].V(v));
              break;
            }
          if (v < mesh.face[f].VN())
            break;
        }
      // if no border is found, this isn't a disk-like polygon, so exit
      if (runPos.IsNull()) {
        std::cout << "No boundary for this sub-mesh. No corners found." << std::endl;
        assert(false);
        return -1;
      }
    }

    // set end corner
    PosType initPos = runPos;
    if (endPos.IsNull() ||
        (!endPos.IsBorder() && vcg::tri::IsMarked(mesh, endPos.FFlip())))
      endPos = runPos;
    // navigate on the boundary and find the corners
    do {
      if (collectRegular)
        corners.push_back(runPos);
      // go to next vertex
      runPos.FlipV();
      runPos.FlipE();
      // if it has to turn to the let, it's a convex corner
      if (runPos.IsBorder() || !vcg::tri::IsMarked(mesh, runPos.FFlip())) {
        if (!collectRegular)
          corners.push_back(runPos);
      } else {
        // else go straight
        runPos.FlipF();
        runPos.FlipE();
        // if it has to turn right, it's a concave corner
        if (!runPos.IsBorder() && vcg::tri::IsMarked(mesh, runPos.FFlip())) {
          // go to the border edge
          while (!runPos.IsBorder() && vcg::tri::IsMarked(mesh, runPos.FFlip())) {
            runPos.FlipF();
            runPos.FlipE();
          }
          // push this concave corner
          if (!collectRegular)
            corners.push_back(runPos);
          if (convexOnly)
            return -1;
        }
      }
    } while (runPos != endPos && runPos != initPos);

    // if no corner's been found (polygon with only one side), use endPos as starting/ending position of the side
    if (corners.empty()) {
      corners.push_back(runPos);
      return 0;
    }
    return corners.size();
  }

  /**
   * @brief CountConcaveCorners counts how many concave veritices there are among those in corners.
   * @param mesh The input mesh.
   * @param corners The corners of the patch.
   * @return The number of concave corners.
   */
  static n_corners_type CountConcaveCorners(const PolyMeshType &mesh, const std::deque<PosType> &corners) {
    n_corners_type nConcaveCorners = 0;
    PosType cornerPos;
    for (size_t c = 0; c < corners.size(); ++c) {
      cornerPos = corners[c];
      cornerPos.FlipE();
      if (!cornerPos.IsBorder() && vcg::tri::IsMarked(mesh, cornerPos.FFlip()))
        ++nConcaveCorners;
    }
    return nConcaveCorners;
  }

  /**
   * @brief FindPolychordsWithTopology lists all the polychords with the given topology, that is,
   * those starting at flow.startSide side and, eventually, ending at flow.endSide side of the
   * convex patch whose sides are defined by corners, or loops, and having flow.nSelfIntersections.
   * @param mesh The mesh representing the patch.
   * @param corners The start edge of each side of the patch.
   * @param flow The topology of the flow: start and end side, number of self intersections.
   * @param polychords The output polychords.
   */
  static void FindPolychordsWithTopology(PolyMeshType &mesh,
                                         const std::deque<PosType> &corners,
                                         const FlowTopology &flow,
                                         std::deque<PosType> &polychords) {
    polychords.clear();

    std::list<FacePointer> selected;
    PosType runPos, startPos;

    if (flow.startSide >= corners.size()) {
      vcg::tri::PolychordCollapse<PolyMeshType>::FindPolychords(mesh, polychords, true);
      return;
    }

    startPos = corners[flow.startSide];
    do {
      if (!startPos.IsBorder()) {
        polychords.clear();
        return;
      }

      if (!startPos.F()->IsS()) {
        // from the current edge, go to the other side
        runPos = startPos;
        runPos.F()->SetS();
        selected.push_back(runPos.F());
        do {
          runPos.FlipE();
          runPos.FlipV();
          runPos.FlipE();
          runPos.FlipF();
        } while (!runPos.IsBorder());
        runPos.F()->SetS();
        selected.push_back(runPos.F());
        if (flow.endSide < corners.size()) {
          // go clockwisily along the border back to this side's first corner
          runPos.FlipV();
          runPos.FlipE();
          while (!runPos.IsBorder()) {
            runPos.FlipF();
            runPos.FlipE();
            if (!runPos.IsBorder()) {   // concave corner
              runPos.FlipE();
              runPos.FlipF();
              break;
            } else {
              runPos.FlipV();
              runPos.FlipE();
            }
          }
          runPos.FlipE();
          // check end side
          if (runPos == corners[flow.endSide])
            polychords.push_back(startPos);
        } else {
          polychords.push_back(startPos);
        }
      }

      startPos.FlipV();
      startPos.FlipE();
      while (!startPos.IsBorder()) {
        startPos.FlipF();
        startPos.FlipE();
      }
    } while (startPos != corners[(flow.startSide + 1) % corners.size()]);

    for (typename std::list<FacePointer>::iterator it = selected.begin(); it != selected.end(); ++it)
      (*it)->ClearS();
  }

  /**
   * @brief PolychordToFaceIndexes gives the sequence of indices of the faces of the polychord from startPos to the end.
   * @param mesh The mesh the polychord is owned by.
   * @param startPos The starting Pos of the polychord.
   * @param indexes The output indices.
   */
  static void PolychordToFaceIndexes(const PolyMeshType &mesh, const PosType &startPos, std::deque<size_t> &indexes) {
    indexes.clear();
    if (startPos.IsNull())
      return;
    PosType runPos = startPos;
    do {
      indexes.push_back(vcg::tri::Index(mesh, runPos.F()));
      runPos.FlipE();
      runPos.FlipV();
      runPos.FlipE();
      runPos.FlipF();
    } while (!runPos.IsBorder() && runPos != startPos);
  }
};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // POLYCHORD_SUPPORT_H
