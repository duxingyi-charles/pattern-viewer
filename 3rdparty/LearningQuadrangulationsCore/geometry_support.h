#ifndef GEOMETRY_SUPPORT_H
#define GEOMETRY_SUPPORT_H

#include <vcg/complex/algorithms/implicit_smooth.h>
#include "patch.h"
#include "uv_grid.h"

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief The PatchGeometrizator class is responsible of assigning a geometry to a Patch.
 * Inherit this and override PatchGeometrizator::geometrizePatch() to assign positions to
 * its vertices.
 */
template < typename PolyMeshType >
class PatchGeometrizator {
public:
  typedef typename PolychordSupport<PolyMeshType>::FlowMap    FlowMap;

  /**
   * @brief PatchGeometrizator Default constructor.
   */
  PatchGeometrizator() { }

  /**
   * @brief ~PatchGeometrizator Default destructor.
   */
  virtual ~PatchGeometrizator() { }

  /**
   * @brief geometrizePatch assigns vertex positions to the Patch patch.
   * @param patch The Patch to geometrize.
   * @param flows The flows in some space.
   * @return true if flows and polychords successfully match, false otherwise.
   */
  virtual bool geometrizePatch(Patch<PolyMeshType> &patch, const FlowMap &flows) const = 0;
};



/************************************* A simple example of a 2D Geometrizator ***********************************/
/**
 * @brief The PatchFlattener class provides methods for shaping a given patch as a planar polygon.
 */
template < typename PolyMeshType >
class PatchFlattener : public PatchGeometrizator<PolyMeshType> {
public:
  typedef typename PolyMeshType::VertexType               VertexType;
  typedef typename PolyMeshType::VertexPointer            VertexPointer;
  typedef typename PolyMeshType::FaceType                 FaceType;
  typedef typename PolyMeshType::FacePointer              FacePointer;
  typedef pl::Patch<PolyMeshType>                         Patch;
  typedef typename Patch::PosType                         PosType;
  typedef pl::PolychordSupport<PolyMeshType>              PolychordSupport;
  typedef typename PolychordSupport::ScalarType           ScalarType;
  typedef typename PolychordSupport::ScalarContainer      ScalarContainer;
  typedef typename PolychordSupport::PointType            PointType;
  typedef typename PolychordSupport::PointContainer       PointContainer;
  typedef typename PolychordSupport::FlowPoints           FlowPoints;
  typedef typename PolychordSupport::FlowContainer        FlowContainer;
  typedef typename PolychordSupport::FlowMap              FlowMap;
  typedef typename PolychordSupport::FlowMapIterator      FlowMapIterator;
  typedef typename PolychordSupport::FlowMapConstIterator FlowMapConstIterator;
  typedef typename PolychordSupport::PolychordPoints      PolychordPoints;
  typedef typename PolychordSupport::PolychordContainer   PolychordContainer;
  typedef typename PolychordSupport::PolychordMap         PolychordMap;
  typedef typename PolychordSupport::PolychordMapIterator PolychordMapIterator;

  static const unsigned int RADIUS = 300;   ///< Default radius.

  /**
   * @brief PatchFlattener Default constructor.
   */
  PatchFlattener() : _radius(RADIUS), _constrainedSmooth(true), _lapWeight(0.8), _qualityThreshold(1)/*, _lambda(0.8)*/ { }

  /**
   * @brief PatchFlattener Copy constructor.
   * @param pf The PatchFlattener to copy.
   */
  PatchFlattener(const PatchFlattener &pf) : _radius(pf._radius), _constrainedSmooth(pf._constrainedSmooth), _lapWeight(pf._lapWeight),
                                             _qualityThreshold(pf._qualityThreshold)/*, _lambda(pf._lambda)*/ { }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief PatchFlattener Move constructor.
   * @param pf The PatchFlattener to move.
   */
  PatchFlattener(PatchFlattener&& pf) : _radius(pf._radius), _constrainedSmooth(pf._constrainedSmooth), _lapWeight(pf._lapWeight),
                                        _qualityThreshold(pf._qualityThreshold)/*, _lambda(pf._lambda)*/ { }
#endif

  /**
   * @brief operator = Copy assignment operator.
   * @param pf The PatchFlattener to copy.
   * @return A reference to this.
   */
  PatchFlattener<PolyMeshType> & operator = (const PatchFlattener &pf) {
    if (this != &pf) {
      _radius = pf._radius;
      _constrainedSmooth = pf._constrainedSmooth;
      _lapWeight = pf._lapWeight;
      _qualityThreshold = pf._qualityThreshold;
    }
    return *this;
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief operator = Move assignment operator.
   * @param pf The PatchFlattener to move.
   * @return A reference to this.
   */
  PatchFlattener & operator = (PatchFlattener&& pf) {
    if (this != &pf) {
      _radius = pf._radius;
      _constrainedSmooth = pf._constrainedSmooth;
      _lapWeight = pf._lapWeight;
      _qualityThreshold = pf._qualityThreshold;
      pf._radius = RADIUS;
      pf._constrainedSmooth = false;
      pf._lapWeight = ScalarType(0.8);
      pf._qualityThreshold = ScalarType(1);
    }
    return *this;
  }
#endif

  /**
   * @brief geometrizePatch assigns vertex positions to the Patch patch in order to
   * lie on a planar polygon defined by the patch itself.
   * @param patch The Patch to geometrize.
   * @param flows The flows in planar space.
   * @return true if flows and polychords successfully match, false otherwise.
   */
  bool geometrizePatch(Patch &patch, const FlowMap &flows) const {
    // compute initial geometry
    PatchFlatteningBySmoothing(patch.getMesh(), patch.getStartCorner(), patch.getBoundaryLenVec(), _radius);
    // compute quad quality
    patch.setQuadsQuality(ComputeQuadsQuality(patch.getMesh()));
    // filter by quality
    if (patch.getQuadsQuality() > _qualityThreshold)
      return false;
    // if no flows, ok with this, return
    if (flows.empty())
      return true;
    // otherwise, compute flow-polychord matches, and set the new flow distance
    if (ComputeFlowsPolychords(patch, flows)) {
      // if ok, update the geometry according to the matches
      if (_constrainedSmooth) {
        SmoothWithFlowConstraints(patch, flows, _lapWeight);
        // update polychords
        UpdateFlowPolychords(patch, flows);
        // compute quad quality
        patch.setQuadsQuality(ComputeQuadsQuality(patch.getMesh()));
        // filter by quality
        if (patch.getQuadsQuality() > _qualityThreshold)
          return false;
      }
      // then return
      return true;
    }
    return false;
  }

  /**
   * @brief setRadius See _radius.
   * @param radius
   */
  void setRadius(const ScalarType radius) {
    _radius = radius;
  }

  /**
   * @brief getRadius See _radius.
   * @return
   */
  ScalarType getRadius() const {
    return _radius;
  }

  /**
   * @brief setSmoothWithConstraints See _constrainedSmooth.
   * @param smooth
   */
  void setSmoothWithConstraints(const bool smooth) {
    _constrainedSmooth = smooth;
  }

  /**
   * @brief getSmoothWithConstraints See _constrainedSmooth.
   * @return
   */
  bool getSmoothWithConstraints() const {
    return _constrainedSmooth;
  }

  /**
   * @brief setLaplacianWeight See _lapWeight.
   * @param weight
   */
  void setLaplacianWeight(const ScalarType weight = ScalarType(0.8)) {
    assert(weight > ScalarType(0) && weight < ScalarType(1));
    _lapWeight = weight;
  }

  /**
   * @brief getLaplacianWeight See _lapWeight.
   * @return
   */
  ScalarType getLaplacianWeight() const {
    return _lapWeight;
  }

  /**
   * @brief getQualityThreshold See _qualityThreshold.
   * @return
   */
  ScalarType getQualityThreshold() const {
    return _qualityThreshold;
  }

  /**
   * @brief setQualityThreshold See _qualityThreshold.
   * @param threshold
   */
  void setQualityThreshold(const ScalarType threshold) {
    _qualityThreshold = threshold;
    if (_qualityThreshold < ScalarType(0))
      _qualityThreshold = ScalarType(0);
    else if (_qualityThreshold > ScalarType(1))
      _qualityThreshold = ScalarType(1);
  }

  /**
   * @brief ComputeQuadsQuality computes the quad quality of the whole mesh.
   * @param mesh The input mesh.
   * @return A (normalized) value of quality, with the lower the better.
   */
  static ScalarType ComputeQuadsQuality(PolyMeshType &mesh) {
    ScalarType finalAvg = ScalarType(0), finalMax = ScalarType(0);
    ScalarType avgDev, maxDev;
    for (size_t f = 0; f < mesh.face.size(); ++f)
      if (!mesh.face[f].IsD()) {
        vcg::PolyAngleDeviation(mesh.face[f], avgDev, maxDev);
        mesh.face[f].Q() = ScalarType(1) - maxDev;
        finalAvg += avgDev;
        if (maxDev > finalMax)
          finalMax = maxDev;
      }
    finalAvg /= mesh.FN();
    return finalMax;
  }

  /**
   * @brief UpdateFlowPolychords updates polychord points and distances.
   * @note ComputeFlowsPolychords must be called before.
   * @param patch The input Patch whose polychords must be updated.
   * @param flows The input FlowMap the polychords are associated to.
   */
  static void UpdateFlowPolychords(Patch &patch, const FlowMap &flows) {
    ScalarType distance, totalDistance = ScalarType(0);
    assert(patch.getPolychordMap().size() == flows.size());
    FlowMapConstIterator flowsIt;
    PolychordMapIterator polychordsIt;
    for (flowsIt = flows.begin(), polychordsIt = patch.getPolychordMap().begin();
         flowsIt != flows.end() && polychordsIt != patch.getPolychordMap().end();
         ++flowsIt, ++polychordsIt) {
      const FlowContainer &flowContainer = flowsIt->second;
      PolychordContainer &polychordContainer = polychordsIt->second;
      assert(flowContainer.size() == polychordContainer.size());
      for (size_t i = 0; i < flowContainer.size(); ++i) {
        const FlowPoints &flowPoints = flowContainer[i];
        PolychordPoints &polychordPoints = polychordContainer[i];
        polychordPoints.polychordToPoints();
        polychordPoints.computeSamples(flowPoints.samplePoints.size());
        distance = ComputePointsToPointsDistance(flowPoints.samplePoints.begin(), flowPoints.samplePoints.end(),
                                                 polychordPoints.samplePoints.begin(), polychordPoints.samplePoints.end(),
                                                 &polychordPoints.distances, false);
        // add a stretch factor
        if (flowPoints.arcLengths.back() >= polychordPoints.arcLengths.back()) {
          if (polychordPoints.arcLengths.back() == ScalarType(0))
            distance = std::numeric_limits<ScalarType>::max();
          else
            distance *= flowPoints.arcLengths.back() / polychordPoints.arcLengths.back();
        } else {
          if (flowPoints.arcLengths.back() == ScalarType(0))
            distance = std::numeric_limits<ScalarType>::max();
          else
            distance *= polychordPoints.arcLengths.back() / flowPoints.arcLengths.back();
        }
        totalDistance += distance;
      }
    }

    patch.setFlowDistance(totalDistance);
  }

  /**
   * @brief SmoothWithFlowConstraints smoothes the mesh of patch with the constraints of flows.
   * @param patch The input Patch.
   * @param flows The input FlowMap constraints.
   * @param lapWeight Weight for the laplacian matrix entry. It must be in [0, 1].
   * @param flowFaceWeight The weight for each face barycenter of the flows.
   */
  static void SmoothWithFlowConstraints(Patch &patch, const FlowMap &flows,
                                        const ScalarType lapWeight = ScalarType(0.8)) {
    // check input
    if (patch.getMesh().IsEmpty() || patch.getStartCorner().IsNull()) {
      std::cout << "The patch is empty. Not smoothed with flow constraints." << std::endl;
      return;
    }

    FlowMapConstIterator flowsIt;
    PolychordMapIterator polychordsIt;
    PosType runPos;
    size_t bInd = 0, sInd = 0, fInd = 0;
    ScalarType step = ScalarType(0);
    std::vector<ScalarType> barycentricCoords(4, ScalarType(0.25));
    std::list<typename vcg::ImplicitSmoother<PolyMeshType>::FaceConstraint> faceConstraints;
    typename std::list<typename vcg::ImplicitSmoother<PolyMeshType>::FaceConstraint>::iterator faceConstraintsIt;
    typename vcg::ImplicitSmoother<PolyMeshType>::Parameter parameter;
    parameter.useMassMatrix = false;
    parameter.useCotWeight = false;
    parameter.lapWeight = ScalarType(1) - lapWeight;
    assert(parameter.lapWeight > ScalarType(0) && parameter.lapWeight < ScalarType(1));
    parameter.fixBorder = true;

    // for each flow topology
    for (polychordsIt = patch.getPolychordMap().begin(); polychordsIt != patch.getPolychordMap().end(); ++polychordsIt) {
      flowsIt = flows.find(polychordsIt->first);
      assert(flowsIt != flows.end());
      const PolychordContainer &polychordContainer = polychordsIt->second;
      const FlowContainer &flowContainer = flowsIt->second;
      assert(polychordContainer.size() == flowContainer.size());

      // for each polychord
      for (size_t p = 0; p < polychordContainer.size(); ++p) {
        const PolychordPoints &polychordPoints = polychordContainer[p];
        const FlowPoints &flowPoints = flowContainer[p];
        assert(polychordPoints.samplePoints.size() == flowPoints.samplePoints.size());

        // get the step length of samples
        step = polychordPoints.arcLengths.back() / ScalarType(polychordPoints.samplePoints.size() - 1);
        assert(step > ScalarType(0));

        // for each face
        runPos = polychordPoints.startPos;
        bInd = 1;
        do {
          if (!runPos.F()->IsV()) {
            fInd = vcg::tri::Index(patch.getMesh(), runPos.F());
            sInd = std::floor(polychordPoints.arcLengths.at(bInd) / step + ScalarType(0.5));
            if (!ComputeQuadFaceBarycentricCoordinates(*runPos.F(), polychordPoints.samplePoints.at(sInd), barycentricCoords))
              barycentricCoords.assign(4, ScalarType(0.25));
            faceConstraints.push_back(typename vcg::ImplicitSmoother<PolyMeshType>::FaceConstraint(fInd, barycentricCoords,
                                                                                                   flowPoints.samplePoints.at(sInd)));
            runPos.F()->SetV();
          }
          bInd += 2;
          runPos.FlipE();
          runPos.FlipV();
          runPos.FlipE();
          runPos.FlipF();
        } while (!runPos.IsBorder() && runPos != polychordPoints.startPos);

        // clear flag
        runPos = polychordPoints.startPos;
        do {
          runPos.F()->ClearV();
          runPos.FlipE();
          runPos.FlipV();
          runPos.FlipE();
          runPos.FlipF();
        } while (!runPos.IsBorder() && runPos != polychordPoints.startPos);
      }
    }

    // add constraints
    parameter.ConstrainedF.resize(faceConstraints.size());
    for (fInd = 0, faceConstraintsIt = faceConstraints.begin();
         fInd < parameter.ConstrainedF.size() && faceConstraintsIt != faceConstraints.end();
         ++fInd, ++faceConstraintsIt)
      parameter.ConstrainedF[fInd] = *faceConstraintsIt;
    faceConstraints.clear();

    // smooth
    vcg::ImplicitSmoother<PolyMeshType>::Compute(patch.getMesh(), parameter);

    vcg::tri::UpdateBounding<PolyMeshType>::Box(patch.getMesh());
    vcg::tri::UpdateNormal<PolyMeshType>::PerPolygonalFaceNormalized(patch.getMesh());
  }

  /**
   * @brief ComputeQuadFaceBarycentricCoordinates finds the four barycentrix coordinates of a point wrt a quad face.
   * @param face The input quad face.
   * @param point The input point.
   * @param coords The ouput four barycentric coordinates of point.
   * @return true if the point is inside face, false otherwise.
   */
  static bool ComputeQuadFaceBarycentricCoordinates(const FaceType &face, const PointType &point, std::vector<ScalarType> &coords) {
    coords.resize(4);
    std::vector<PointType> points(4);
    PointType p = point, m(0, 0, 0);
    for (size_t i = 0; i < 4; ++i) {
      points[i] = face.cP(i);
      m += points[i];
    }
    m /= 4;
    for (size_t i = 0; i < 4; ++i)
      points[i] -= m;
    vcg::Matrix33<ScalarType> rot = vcg::RotationMatrix<ScalarType>(face.cN(), PointType(0, 0, 1));
    for (size_t i = 0; i < 4; ++i)
      points[i] = rot * points[i];
    p = rot * p;
    vcg::Point2<ScalarType> v00, v01, v10, v11, x;
    vcg::Point2<ScalarType> aa, bb, cc, dd;
    ScalarType a, b, c, l, u, l1, l2, u1, u2;
    x[0] = p[0];
    x[1] = p[1];
    v00[0] = points[0][0];
    v00[1] = points[0][1];
    v10[0] = points[1][0];
    v10[1] = points[1][1];
    v01[0] = points[3][0];
    v01[1] = points[3][1];
    v11[0] = points[2][0];
    v11[1] = points[2][1];
    aa = v00 - x;
    bb = v10 - v00;
    cc = v01 - v00;
    dd = v00 - v01 - v10 + v11;
    // solve for l
    l = l1 = l2 = -ScalarType(1);
    a = bb ^ dd;
    b = (bb ^ cc) + (aa ^ dd);
    c = aa ^ cc;
    if (a != ScalarType(0) && (b * b - ScalarType(4) * a * c) >= ScalarType(0)) {
      l1 = (-b + std::sqrt(b * b - ScalarType(4) * a * c)) / (ScalarType(2) * a);
      l2 = (-b - std::sqrt(b * b - ScalarType(4) * a * c)) / (ScalarType(2) * a);
    }
    if (l1 >= ScalarType(0) && l1 <= ScalarType(1))
      l = l1;
    else if (l2 >= ScalarType(0) && l2 <= ScalarType(1))
      l = l2;
    // solve for u
    u = u1 = u2 = -ScalarType(1);
    a = cc ^ dd;
    b = (cc ^ bb) + (aa ^ dd);
    c = aa ^ bb;
    if (a != ScalarType(0) && b * b - ScalarType(4) * a * c >= ScalarType(0)) {
      u1 = (-b + std::sqrt(b * b - ScalarType(4) * a * c)) / (ScalarType(2) * a);
      u2 = (-b - std::sqrt(b * b - ScalarType(4) * a * c)) / (ScalarType(2) * a);
    }
    if (u1 >= ScalarType(0) && u1 <= ScalarType(1))
      u = u1;
    else if (u2 >= ScalarType(0) && u2 <= ScalarType(1))
      u = u2;
    // check
    if (l >= ScalarType(0) && l <= ScalarType(1) && u >= ScalarType(0) && u <= ScalarType(1)) {
      coords[0] = (ScalarType(1) - l) * (ScalarType(1) - u);
      coords[1] = l * (ScalarType(1) - u);
      coords[2] = l * u;
      coords[3] = (ScalarType(1) - l) * u;
      return true;
    }
    coords.assign(4, ScalarType(0.25));
    return false;
  }

  /**
   * @brief ComputeFlowsPolychords computes the flow-polychord best matches and the distance between them.
   * @param patch The input Patch.
   * @param flows The input flows map.
   * @return true if all flows are matched with polychords, false otherwise.
   */
  static bool ComputeFlowsPolychords(Patch &patch, const FlowMap &flows) {
    std::deque<PosType> polychordVec;
    std::deque<PosType> corners;
    PosType startPos;
    n_corners_type nc = 0;
    ScalarType dist = ScalarType(0);
    FlowMapConstIterator flowsIt;
    PolychordMapIterator polychordsIt;

    // init distance
    ScalarType distance = ScalarType(0);

    // check input
    if (patch.getMesh().IsEmpty() || patch.getStartCorner().IsNull()) {
      std::cout << "The patch is empty. No polychord associated with flows." << std::endl;
      return false;
    }

    n_corners_type nSides = patch.getNumberOfCorners();
    if (nSides == 0)
      nSides = 1;

    // get polychord map
    PolychordMap &polychords = patch.getPolychordMap();
    polychords.clear();

    // compute corners
    startPos = patch.getStartCorner();
    startPos.FlipE();
    while (!startPos.IsBorder()) {
      startPos.FlipF();
      startPos.FlipE();
    }
    startPos.FlipV();
    nc = PolychordSupport::FindPolygonCorners(patch.getMesh(), corners, false, false, startPos);
    assert(nc == patch.getNumberOfCorners());

    // for each flow topology
    for (flowsIt = flows.begin(); flowsIt != flows.end(); ++flowsIt) {
      const FlowTopology &flowTopology = flowsIt->first;
      const FlowContainer &flowContainer = flowsIt->second;

      // find all polychords with the current flow topology
      PolychordSupport::FindPolychordsWithTopology(patch.getMesh(), corners, flowTopology, polychordVec);
      if (polychordVec.size() < flowContainer.size()) {
        std::cout << "There are not as many polychords as flows. Impossible to map." << std::endl;
        polychords.clear();
        distance = std::numeric_limits<ScalarType>::max();
        patch.setFlowDistance(distance);
        return false;
      }

      // insert new sequence of polychords for the current flow topology
      polychordsIt = polychords.insert(std::make_pair(flowTopology, PolychordContainer())).first;
      assert(polychordsIt != polychords.end());
      PolychordContainer &polychordContainer = polychordsIt->second;

      // compute distances between flows and polychords
      dist = ComputeFlowsPolychordsDistance(nSides, flowTopology, flowContainer, polychordContainer, polychordVec);
      assert(dist != std::numeric_limits<ScalarType>::max() && polychordContainer.size() == flowContainer.size());

      // update distance
      distance += dist;
    }

    patch.setFlowDistance(distance);
    return true;
  }

  /**
   * @brief ComputeFlowsPolychordsDistance computes the best flows-polychords matching and the total distance.
   * @param nSides The number of sides.
   * @param flowTopology The topology of the flow. It is needed to distinguish between loops and border-to-border.
   * @param flowContainer The flows.
   * @param polychordContainer The output polychords.
   * @param polychordVec The input polychords' start Pos.
   * @return The total distance.
   */
  static ScalarType ComputeFlowsPolychordsDistance(const n_corners_type nSides,
                                                   const FlowTopology &flowTopology,
                                                   const FlowContainer &flowContainer,
                                                   PolychordContainer &polychordContainer,
                                                   const std::deque<PosType> &polychordVec) {
    polychordContainer.clear();

    if (polychordVec.size() < flowContainer.size() || flowContainer.empty())
      return std::numeric_limits<ScalarType>::max();

    polychordContainer.resize(flowContainer.size());

    // total distance to return
    ScalarType totalDistance = ScalarType(0);

    // build a matrix of distance values per-polychord and per-flow
    std::deque< std::deque<ScalarType> > distances(polychordVec.size(),
                                                   std::deque<ScalarType>(flowContainer.size(),
                                                                          std::numeric_limits<ScalarType>::max()));

    // some variables for the best maching
    std::list<size_t> rows;
    for (size_t i = 0; i < polychordVec.size(); ++i)
      rows.push_back(i);
    std::list<size_t> cols;
    for (size_t i = 0; i < flowContainer.size(); ++i)
      cols.push_back(i);
    ScalarType bestMatchDistance = std::numeric_limits<ScalarType>::max();
    std::list<size_t>::iterator bestRow, bestCol, rowsIt, colsIt;

    // loops or same side U's
    if (flowTopology.startSide >= nSides || flowTopology.endSide >= nSides || flowTopology.startSide == flowTopology.endSide) {
      PosType runPos;
      PointType polychordCenter(0, 0, 0), point(0, 0, 0), point2(0, 0, 0);
      size_t nEdges = 0;
      ScalarType value1, value2, dist, minDistance = std::numeric_limits<ScalarType>::max();

      // compute flow centers
      PointContainer flowCenters(flowContainer.size());
      for (size_t j = 0; j < flowCenters.size(); ++j)
        if (flowTopology.startSide >= nSides || flowTopology.endSide >= nSides) {
          flowCenters[j].SetZero();
          for (size_t k = 0; k < flowContainer[j].points.size() - 1; ++k)
            flowCenters[j] += flowContainer[j].points[k];
          flowCenters[j] /= flowContainer[j].points.size() - 1;
        } else {
          flowCenters[j] = (flowContainer[j].points.front() + flowContainer[j].points.back()) / ScalarType(2);
        }

      // build a matrix of PolychordPoints per-polychord and per-flow
      std::deque<PolychordContainer> polychords(polychordVec.size(), PolychordContainer(flowContainer.size()));

      // compute PolychordPoints and distances
      for (size_t i = 0; i < polychordVec.size(); ++i)
        for (size_t j = 0; j < flowContainer.size(); ++j) {
          const FlowPoints &flowPoints = flowContainer[j];
          PolychordPoints &polychordPoints = polychords[i][j];

          if (flowTopology.startSide >= nSides || flowTopology.endSide >= nSides) {
            // compute polychord center
            polychordCenter.SetZero();
            nEdges = 0;
            runPos = polychordVec[i];
            do {
              polychordCenter += (runPos.V()->P() + runPos.VFlip()->P()) / 2.0;
              ++nEdges;
              runPos.FlipE();
              runPos.FlipV();
              runPos.FlipE();
              runPos.FlipF();
            } while(runPos != polychordVec[i]);
            polychordCenter /= nEdges;

            // find start pos
            minDistance = std::numeric_limits<ScalarType>::max();
            runPos = polychordVec[i];
            do {
              point = (runPos.V()->P() + runPos.VFlip()->P()) / ScalarType(2);
              point -= polychordCenter - flowCenters[j];
              dist = (point - flowPoints.points.front()).Norm();
              if (dist < minDistance) {
                minDistance = dist;
                polychordPoints.startPos = runPos;
              }
              runPos.FlipE();
              runPos.FlipV();
              runPos.FlipE();
              runPos.FlipF();
            } while(runPos != polychordVec[i]);

            // compute points and arc lengths
            polychordPoints.polychordToPoints();

            // move to flow barycenter
            for (size_t k = 0; k < polychordPoints.points.size(); ++k)
              polychordPoints.points[k] -= polychordCenter - flowCenters[j];

            // compute samples
            polychordPoints.computeSamples(flowPoints.samplePoints.size());

            // compute distance to decide polychord orientation
            value1 = ComputePointsToPointsDistance(flowPoints.samplePoints.begin(), flowPoints.samplePoints.end(),
                                                   polychordPoints.samplePoints.begin(), polychordPoints.samplePoints.end());
            value2 = ComputePointsToPointsDistance(flowPoints.samplePoints.begin(), flowPoints.samplePoints.end(),
                                                   polychordPoints.samplePoints.rbegin(), polychordPoints.samplePoints.rend());
            if (value2 < value1) {
              polychordPoints.startPos.FlipF();
              polychordPoints.invertPoints();
              polychordPoints.updateArcLengths();
            }

            // move back to polychord barycenter
            for (size_t k = 0; k < polychordPoints.points.size(); ++k)
              polychordPoints.points[k] += polychordCenter - flowCenters[j];

          } else {
            // compute polychord center
            polychordPoints.startPos = runPos = polychordVec[i];
            polychordCenter = (runPos.V()->P() + runPos.VFlip()->P()) / ScalarType(2);
            do {
              runPos.FlipE();
              runPos.FlipV();
              runPos.FlipE();
              runPos.FlipF();
            } while(!runPos.IsBorder());
            polychordCenter = (polychordCenter + (runPos.V()->P() + runPos.VFlip()->P()) / ScalarType(2)) / ScalarType(2);

            // find best starting pos
            point  = (polychordPoints.startPos.V()->P() + polychordPoints.startPos.VFlip()->P()) / ScalarType(2);
            point  -= polychordCenter - flowCenters[j];
            point2 = (runPos.V()->P() + runPos.VFlip()->P()) / ScalarType(2);
            point2 -= polychordCenter - flowCenters[j];
            if ((point2 - flowPoints.points.front()).Norm() < (point - flowPoints.points.front()).Norm())
              polychordPoints.startPos = runPos;

            // compute polychord points
            polychordPoints.polychordToPoints();
          }

          // re-compute samples
          polychordPoints.computeSamples(flowPoints.samplePoints.size());

          // compute real distance
          distances[i][j] = ComputePointsToPointsDistance(flowPoints.samplePoints.begin(), flowPoints.samplePoints.end(),
                                                          polychordPoints.samplePoints.begin(), polychordPoints.samplePoints.end(),
                                                          &polychordPoints.distances, false);

          // add a stretch factor
          if (flowPoints.arcLengths.back() >= polychordPoints.arcLengths.back()) {
            if (polychordPoints.arcLengths.back() == ScalarType(0))
              distances[i][j] = std::numeric_limits<ScalarType>::max();
            else
              distances[i][j] *= flowPoints.arcLengths.back() / polychordPoints.arcLengths.back();
          } else {
            if (flowPoints.arcLengths.back() == ScalarType(0))
              distances[i][j] = std::numeric_limits<ScalarType>::max();
            else
              distances[i][j] *= polychordPoints.arcLengths.back() / flowPoints.arcLengths.back();
          }
        }

      // find the best matches
      while (!cols.empty()) {
        // find current minimum
        bestMatchDistance = std::numeric_limits<ScalarType>::max();
        for (rowsIt = rows.begin(); rowsIt != rows.end(); ++rowsIt)
          for (colsIt = cols.begin(); colsIt != cols.end(); ++colsIt)
            if (distances[*rowsIt][*colsIt] < bestMatchDistance) {
              bestMatchDistance = distances[*rowsIt][*colsIt];
              bestRow = rowsIt;
              bestCol = colsIt;
            }
        // add current match distance to total
        totalDistance += bestMatchDistance;
        // push output current best match
#if __cplusplus >= 201103L || _MSC_VER >= 1800
        polychordContainer[*bestCol] = std::move(polychords[*bestRow][*bestCol]);
#else
        polychordContainer[*bestCol] = polychords[*bestRow][*bestCol];
#endif
        // delete the current best row and column
        rows.erase(bestRow);
        cols.erase(bestCol);
      }

    } else {      // border to border
      PolychordContainer polychords(polychordVec.size());

      // compute PolychordPoints and distances
      for (size_t i = 0; i < polychords.size(); ++i) {
        PolychordPoints &polychordPoints = polychords[i];
        polychordPoints.startPos = polychordVec[i];
        polychordPoints.polychordToPoints();
        polychordPoints.computeSamples(flowContainer.front().samplePoints.size());
      }

      // for each polychord and for each flow
      for (size_t i = 0; i < polychords.size(); ++i) {
        PolychordPoints &polychordPoints = polychords[i];
        for (size_t j = 0; j < flowContainer.size(); ++j) {
          const FlowPoints &flowPoints = flowContainer[j];

          // compute real distance
          distances[i][j] = ComputePointsToPointsDistance(flowPoints.samplePoints.begin(), flowPoints.samplePoints.end(),
                                                          polychordPoints.samplePoints.begin(), polychordPoints.samplePoints.end(),
                                                          &polychordPoints.distances, false);

          // add a stretch factor
          if (flowPoints.arcLengths.back() >= polychordPoints.arcLengths.back()) {
            if (polychordPoints.arcLengths.back() == ScalarType(0))
              distances[i][j] = std::numeric_limits<ScalarType>::max();
            else
              distances[i][j] *= flowPoints.arcLengths.back() / polychordPoints.arcLengths.back();
          } else {
            if (flowPoints.arcLengths.back() == ScalarType(0))
              distances[i][j] = std::numeric_limits<ScalarType>::max();
            else
              distances[i][j] *= polychordPoints.arcLengths.back() / flowPoints.arcLengths.back();
          }
        }
      }

      // find the best matches
      while (!cols.empty()) {
        // find current minimum
        bestMatchDistance = std::numeric_limits<ScalarType>::max();
        for (rowsIt = rows.begin(); rowsIt != rows.end(); ++rowsIt)
          for (colsIt = cols.begin(); colsIt != cols.end(); ++colsIt)
            if (distances[*rowsIt][*colsIt] < bestMatchDistance) {
              bestMatchDistance = distances[*rowsIt][*colsIt];
              bestRow = rowsIt;
              bestCol = colsIt;
            }
        // add current match distance to total
        totalDistance += bestMatchDistance;
        // push output current best match
#if __cplusplus >= 201103L || _MSC_VER >= 1800
        polychordContainer[*bestCol] = std::move(polychords[*bestRow]);
#else
        polychordContainer[*bestCol] = polychords[*bestRow];
#endif
        // delete the current best row and column
        rows.erase(bestRow);
        cols.erase(bestCol);
      }
    }

    // finally return the total distance
    return totalDistance;
  }

  /**
   * @brief computePointsToPointsDistance computes the distance between a sequence of point pairs.
   * @param start1 The first point of the first sequence.
   * @param end1 The end point of the first sequence.
   * @param start2 The first point of the second sequence.
   * @param end2 The end point of the second sequence.
   * @return The total distance.
   */
  template < typename Iterator1, typename Iterator2 >
  static ScalarType ComputePointsToPointsDistance(Iterator1 start1, Iterator1 end1, Iterator2 start2, Iterator2 end2,
                                                  ScalarContainer *values = 0, const bool geodesicApprox = false) {
    if (values)
      values->clear();
    ScalarType distance = ScalarType(0), value = ScalarType(0);
    PointType dir(0, 0, 0), dir1(0, 0, 0), dir2(0, 0, 0);
    PointType last1(0, 0, 0), last2(0, 0, 0);
    if (start1 != end1 && start2 != end2) {
      dir = *start2 - *start1;
      last1 = *start1;
      last2 = *start2;
      ++start1;
      ++start2;
      if (start1 != end1 && start2 != end2) {
        dir1 = (*start1 - last1).Normalize();
        dir2 = (*start2 - last2).Normalize();
        if (geodesicApprox)
          value = vcg::ApproximateGeodesicDistance(last1, dir1, last2, dir2);
        else
          value = dir.Norm();
      } else {
        value = dir.Norm();
      }
      if (values)
        values->push_back(value);
      distance += value * value;
    }
    for (; start1 != end1 && start2 != end2; last1 = *start1, last2 = *start2, ++start1, ++start2) {
      dir = *start2 - *start1;
      dir1 = (*start1 - last1).Normalize();
      dir2 = (*start2 - last2).Normalize();
      if (geodesicApprox)
        value = vcg::ApproximateGeodesicDistance(*start1, dir1, *start2, dir2);
      else
        value = dir.Norm();
      if (values)
        values->push_back(value);
      distance += value * value;
    }
    distance = std::sqrt(distance);
    return distance;
  }

  /**
   * @brief PatchFlatteningBySmoothing assigns vertex positions to mesh s.t. boundaries lie on a polygon.
   * Internal vertex positions are then found by smoothing.
   * @param mesh The input mesh.
   * @param startPos The starting corner.
   * @param radius The radius of the circle for the boundary vertices.
   */
  static void PatchFlatteningBySmoothing(PolyMeshType &mesh,
                                         const PosType &startPos,
                                         const std::vector<num_type> &sidesLengthVec,
                                         const ScalarType radius) {
    // check input
    vcg::tri::RequirePerFaceFlags(mesh);
    vcg::tri::RequirePerVertexFlags(mesh);
    if (mesh.IsEmpty())
      return;
    if (startPos.IsNull())
      return;
    if (radius < 1)
      return;

    // reset coords
    for (size_t v = 0; v < mesh.vert.size(); v++)
      mesh.vert[v].P().SetZero();

    // boundary on polygon
    BoundaryToPolygon(startPos, sidesLengthVec, radius);

    // set smoother parameters
    typename vcg::ImplicitSmoother<PolyMeshType>::Parameter parameters;
    parameters.useMassMatrix = false;
    parameters.useCotWeight = false;
    parameters.fixBorder = true;
    parameters.lapWeight = ScalarType(1);

    // smooth
    vcg::ImplicitSmoother<PolyMeshType>::Compute(mesh, parameters);

    vcg::tri::UpdateBounding<PolyMeshType>::Box(mesh);
    vcg::tri::UpdateNormal<PolyMeshType>::PerPolygonalFaceNormalized(mesh);
  }

  /**
   * @brief BoundaryToPolygon assigns boundary vertex positions to lie on a polygon.
   * That is, corners are assigned positions on a circle, vertices between two corners
   * are assigned positions on the straight line connecting them.
   * @param startPos The starting corner.
   * @param sidesLengthVec A vector where each element is the length (number of edges) of the corresponding side.
   * @param radius The radius of the circle.
   */
  static void BoundaryToPolygon(const PosType &startPos,
                                const std::vector<num_type> &sidesLengthVec,
                                const ScalarType radius) {
    // check input
    if (startPos.IsNull())
      return;
    if (sidesLengthVec.size() < 1)
      return;
    if (sidesLengthVec.size() < 3) {
      BoundaryToCircle(startPos, sidesLengthVec, radius);
      return;
    }

    PosType runPos = startPos;
    if (!runPos.IsBorder())
      return;

    num_type perimeter = std::accumulate(sidesLengthVec.begin(), sidesLengthVec.end(), 0);
    ScalarType stepAngle = 2 * M_PI / ScalarType(perimeter);
    ScalarType angle = M_PI + (M_PI - stepAngle * sidesLengthVec[0]) / ScalarType(2);
    PointType currentCorner(radius * cos(angle), radius * sin(angle), 0);
    angle += stepAngle * sidesLengthVec[0];
    PointType nextCorner(radius * cos(angle), radius * sin(angle), 0);
    size_t side = 0;
    num_type edge = 0;
    do {
      runPos.V()->P() = currentCorner + (nextCorner - currentCorner) * edge / ScalarType(sidesLengthVec[side]);
      runPos.FlipV();
      runPos.FlipE();
        ///dxy add
        if (edge == sidesLengthVec[side] - 1) { //last edge on this side
            while (!runPos.IsBorder()) {
                runPos.FlipF();
                runPos.FlipE();
            }
            currentCorner = nextCorner;
            side += 1;
            angle += stepAngle * sidesLengthVec[side];
            nextCorner = PointType(radius * cos(angle), radius * sin(angle), 0);
            edge = 0;
        }
        else {
            while (!runPos.IsBorder()) {
                runPos.FlipF();
                runPos.FlipE();
            }
            edge += 1;
        }
        
        ///dxy add end
        
        ///dxy delete
//      if (!runPos.IsBorder()) {
//        runPos.FlipF();
//        runPos.FlipE();
//        if (!runPos.IsBorder()) {   // concave corner
//          while (!runPos.IsBorder()) {
//            runPos.FlipF();
//            runPos.FlipE();
//          }
//          currentCorner = nextCorner;
//          side = (side + 1) % sidesLengthVec.size();
//          angle += stepAngle * sidesLengthVec[side];
//          nextCorner = PointType(radius * cos(angle), radius * sin(angle), 0);
//          edge = 0;
//        } else
//          ++edge;
//      } else {    // convex corner
//        currentCorner = nextCorner;
//        side = (side + 1) % sidesLengthVec.size();
//        angle += stepAngle * sidesLengthVec[side];
//        nextCorner = PointType(radius * cos(angle), radius * sin(angle), 0);
//        edge = 0;
//      }
        ///dxy delete end
    } while (runPos != startPos);
  }

  /**
   * @brief BoundaryToCircle assigns boundary vertex positions lying on a circle.
   * @param startPos The starting corner.
   * @param sidesLengthVec A vector where each element is the length (number of edges) of the corresponding side.
   * @param radius The radius of the circle.
   */
  static void BoundaryToCircle(const PosType &startPos,
                               const std::vector<num_type> &sidesLengthVec,
                               const ScalarType radius) {
    // check input
    if (startPos.IsNull())
      return;
    if (sidesLengthVec.size() < 1)
      return;

    ScalarType angle = ScalarType(M_PI);
    num_type perimeter = std::accumulate(sidesLengthVec.begin(), sidesLengthVec.end(), 0);
    ScalarType stepAngle = 2 * M_PI / ScalarType(perimeter);
    if (sidesLengthVec.size() > 2)
      angle += (2 * M_PI - (angle + stepAngle * sidesLengthVec[0])) / ScalarType(2);

    PosType runPos = startPos;
    if (!runPos.IsBorder())
      return;
    do {
      runPos.V()->P() = PointType(radius * cos(angle), radius * sin(angle), ScalarType(0));
      angle += stepAngle;
      runPos.FlipV();
      runPos.FlipE();
      while (!runPos.IsBorder()) {
        runPos.FlipF();
        runPos.FlipE();
      }
    } while (runPos != startPos);
  }

  /**
   * @brief BoundaryToPolygon assigns boundary vertex positions to lie on a polygon.
   * That is, corners are assigned positions on a circle, vertices between two corners
   * are assigned positions on the straight line connecting them.
   * @param sidesLengthVec A vector where each element is the length (number of edges) of the corresponding side.
   * @param sides The output corner points.
   * @param radius The radius of the circle.
   */
  static void BoundaryToPolygon(const std::vector<num_type> &sidesLengthVec,
                                std::vector<PointContainer> &sides,
                                ScalarType radius) {
    sides.clear();

    if (sidesLengthVec.size() < 1)
      return;
    if (sidesLengthVec.size() < 3) {
      BoundaryToCircle(sidesLengthVec, sides, radius);
      return;
    }

    sides.resize(sidesLengthVec.size());

    num_type perimeter = std::accumulate(sidesLengthVec.begin(), sidesLengthVec.end(), 0);
    ScalarType stepAngle = 2 * M_PI / ScalarType(perimeter);
    ScalarType angle = M_PI + (M_PI - stepAngle * sidesLengthVec[0]) / ScalarType(2);
    PointType currentCorner(radius * cos(angle), radius * sin(angle), 0);
    angle += stepAngle * sidesLengthVec[0];
    PointType nextCorner(radius * cos(angle), radius * sin(angle), 0);

    for (size_t i = 0; i < sidesLengthVec.size(); ++i) {
      for (num_type j = 0; j < sidesLengthVec[i]; ++j)
        sides[i].push_back(currentCorner + (nextCorner - currentCorner) * j / ScalarType(sidesLengthVec[i]));
      sides[i].push_back(nextCorner);
      currentCorner = nextCorner;
      angle += stepAngle * sidesLengthVec[(i + 1) % sidesLengthVec.size()];
      nextCorner = PointType(radius * cos(angle), radius * sin(angle), 0);
    }
  }

  /**
   * @brief BoundaryToCircle assigns boundary vertex positions to lie on a circle.
   * @param sidesLengthVec A vector where each element is the length (number of edges) of the corresponding side.
   * @param sides The output corner points.
   * @param radius The radius of the circle.
   */
  static void BoundaryToCircle(const std::vector<num_type> &sidesLengthVec,
                               std::vector<PointContainer> &sides,
                               ScalarType radius) {
    sides.clear();

    if (sidesLengthVec.size() < 1)
      return;

    sides.resize(sidesLengthVec.size());

    ScalarType angle = ScalarType(M_PI);
    num_type perimeter = std::accumulate(sidesLengthVec.begin(), sidesLengthVec.end(), 0);
    ScalarType stepAngle = 2 * M_PI / ScalarType(perimeter);
    if (sidesLengthVec.size() > 2)
      angle += (2 * M_PI - (angle + stepAngle * sidesLengthVec[0])) / ScalarType(2);
    for (size_t i = 0; i < sidesLengthVec.size(); ++i) {
      for (num_type j = 0; j < sidesLengthVec[i]; ++j) {
        sides[i].push_back(PointType(radius * cos(angle), radius * sin(angle), ScalarType(0)));
        angle += stepAngle;
      }
      sides[i].push_back(PointType(radius * cos(angle), radius * sin(angle), ScalarType(0)));
    }
  }

protected:
  ScalarType      _radius;                ///< The radius of the circle circumscribing the polygon.
  bool            _constrainedSmooth;     ///< Smooth with flows as constraints, or not.
  ScalarType      _lapWeight;             ///< Weight for the smoother.
  ScalarType      _qualityThreshold;      ///< Threshold for the quad quality: patches with quality above this are discarded.
};




/**
 * @brief The PatchRetopologizer class provides methods for shaping a given patch in a way such that
 * it retopologizes a (part of a) given triangle mesh.
 */
template < typename PolyMeshType, typename TriMeshType >
class PatchRetopologizer : public PatchFlattener<PolyMeshType> {
public:
  typedef pl::PatchFlattener<PolyMeshType>                PatchFlattener;
  typedef typename PatchFlattener::VertexType             PVertexType;
  typedef typename PatchFlattener::VertexPointer          PVertexPointer;
  typedef typename PatchFlattener::FaceType               PFaceType;
  typedef typename PatchFlattener::FacePointer            PFacePointer;
  typedef typename PatchFlattener::Patch                  Patch;
  typedef typename PatchFlattener::PosType                PosType;
  typedef typename PatchFlattener::PolychordSupport       PolychordSupport;
  typedef typename PatchFlattener::ScalarType             ScalarType;
  typedef typename PatchFlattener::ScalarContainer        ScalarContainer;
  typedef typename PatchFlattener::PointType              PointType;
  typedef typename PatchFlattener::PointContainer         PointContainer;
  typedef typename PatchFlattener::FlowPoints             FlowPoints;
  typedef typename PatchFlattener::FlowContainer          FlowContainer;
  typedef typename PatchFlattener::FlowMap                FlowMap;
  typedef typename PatchFlattener::FlowMapIterator        FlowMapIterator;
  typedef typename PatchFlattener::FlowMapConstIterator   FlowMapConstIterator;
  typedef typename PatchFlattener::PolychordPoints        PolychordPoints;
  typedef typename PatchFlattener::PolychordContainer     PolychordContainer;
  typedef typename PatchFlattener::PolychordMap           PolychordMap;
  typedef typename PatchFlattener::PolychordMapIterator   PolychordMapIterator;
  typedef typename TriMeshType::VertexType                TVertexType;
  typedef typename TriMeshType::VertexPointer             TVertexPointer;
  typedef typename TriMeshType::FaceType                  TFaceType;
  typedef typename TriMeshType::FacePointer               TFacePointer;
  typedef typename TriMeshType::ScalarType                TScalarType;
  typedef typename TriMeshType::CoordType                 TPointType;
  typedef typename TVertexType::TexCoordType::PointType   TUVPointType;

  using PatchFlattener::ComputeQuadsQuality;
  using PatchFlattener::ComputeFlowsPolychords;
  using PatchFlattener::SmoothWithFlowConstraints;
  using PatchFlattener::UpdateFlowPolychords;
  /**
   * @brief PatchRetopologizer Default constructor;
   */
  PatchRetopologizer() : PatchFlattener(), _tMesh(0), _nCorners(std::numeric_limits<n_corners_type>::max()) { }

  /**
   * @brief PatchRetopologizer Copy constructor.
   * @param pr The PatchRetopologizer to copy.
   */
  PatchRetopologizer(const PatchRetopologizer &pr) : PatchFlattener(pr),
                                                     _tMesh(0),
                                                     _boundary(pr._boundary),
                                                     _boundaryUV(pr._boundaryUV),
                                                     _flowMap(pr._flowMap),
                                                     _flowMapUV(pr._flowMapUV),
                                                     _nCorners(pr._nCorners),
                                                     _boundaryLen(pr._boundaryLen) {
    if (pr._tMesh) {
      _tMesh = new TriMeshType;
      pr.getTriMesh(*_tMesh);
    }
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief PatchRetopologizer Move constructor.
   * @param pr The PatchRetopologizer to move.
   */
  PatchRetopologizer(PatchRetopologizer&& pr) : PatchFlattener(std::move(pr)),
                                                _tMesh(pr._tMesh),
                                                _boundary(std::move(pr._boundary)),
                                                _boundaryUV(std::move(pr._boundaryUV)),
                                                _flowMap(std::move(pr._flowMap)),
                                                _flowMapUV(std::move(pr._flowMapUV)),
                                                _nCorners(pr._nCorners),
                                                _boundaryLen(std::move(pr._boundaryLen)) {
    if (pr._tMesh)
      pr._tMesh = 0;
  }
#endif

  /**
   * @brief ~PatchRetopologizer Destructor.
   */
  ~PatchRetopologizer() {
    if (_tMesh)
      delete _tMesh;
  }

  /**
   * @brief operator = Copy assignment operator.
   * @param pr The PatchRetopologizer to copy.
   * @return A reference to this.
   */
  PatchRetopologizer & operator = (const PatchRetopologizer<PolyMeshType,TriMeshType> &pr) {
    if (this != &pr) {
      PatchFlattener::operator =(pr);
      if (pr._tMesh) {
        if (!_tMesh)
          _tMesh = new TriMeshType;
        pr.getTriMesh(*_tMesh);
      } else {
        if (_tMesh)
          delete _tMesh;
        _tMesh = 0;
      }
      _boundary = pr._boundary;
      _boundaryUV = pr._boundaryUV;
      _flowMap = pr._flowMap;
      _flowMapUV = pr._flowMapUV;
      _nCorners = pr._nCorners;
      _boundaryLen = pr._boundaryLen;
    }
    return *this;
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief operator = Move assignment operator.
   * @param pr The PatchRetopologizer to move.
   * @return A reference to this.
   */
  PatchRetopologizer & operator = (PatchRetopologizer<PolyMeshType,TriMeshType>&& pr) {
    if (this != &pr) {
      PatchFlattener::operator =(std::move(pr));
      if (_tMesh)
        delete _tMesh;
      _tMesh = pr._tMesh;
      _boundary = std::move(pr._boundary);
      _boundaryUV = std::move(pr._boundaryUV);
      _flowMap = std::move(pr._flowMap);
      _flowMapUV = std::move(pr._flowMapUV);
      _nCorners = pr._nCorners;
      _boundaryLen = std::move(pr._boundaryLen);
      pr._tMesh = 0;
      pr._nCorners = std::numeric_limits<n_corners_type>::max();
    }
    return *this;
  }
#endif

  /**
   * @brief getTriMesh gets a copy of _tMesh.
   * @param tMesh
   */
  void getTriMesh(TriMeshType &tMesh) const {
    tMesh.Clear();
    if (!_tMesh)
      return;
    vcg::tri::Append<TriMeshType,TriMeshType>::MeshCopy(tMesh, const_cast<TriMeshType &>(*_tMesh), false, true);
    vcg::tri::UpdateTopology<TriMeshType>::FaceFace(tMesh);
  }

  /**
   * @brief getTriMesh gets access to _tMesh.
   * @return
   */
  TriMeshType & getTriMesh() {
    if (!_tMesh)
      _tMesh = new TriMeshType;
    return *_tMesh;
  }

  /**
   * @brief getBoundary gets a copy of _boundary.
   * @param boundary
   */
  void getBoundary(std::vector<PointContainer> &boundary) const {
    boundary = _boundary;
  }

  /**
   * @brief getBoundary gets access to _boundary.
   * @return
   */
  std::vector<PointContainer> & getBoundary() {
    return _boundary;
  }

  /**
   * @brief getBoundaryUV gets a copy of _boundaryUV.
   * @param boundaryUV
   */
  void getBoundaryUV(std::vector<PointContainer> &boundaryUV) const {
    boundaryUV = _boundaryUV;
  }

  /**
   * @brief getBoundaryUV gets access to _boundaryUV.
   * @return
   */
  std::vector<PointContainer> & getBoundaryUV() {
    return _boundaryUV;
  }

  /**
   * @brief getFlowMap gets a copy of _flowMap.
   * @param flowMap
   */
  void getFlowMap(FlowMap &flowMap) const {
    flowMap = _flowMap;
  }

  /**
   * @brief getFlowMap gets access to _flowMap.
   * @return
   */
  FlowMap & getFlowMap() {
    return _flowMap;
  }

  /**
   * @brief getFlowMapUV gets a copy of _flowMapUV.
   * @param flowMapUV
   */
  void getFlowMapUV(FlowMap &flowMapUV) const {
    flowMapUV = _flowMapUV;
  }

  /**
   * @brief getFlowMapUV gets access to _flowMapUV.
   * @return
   */
  FlowMap & getFlowMapUV() {
    return _flowMapUV;
  }

  /**
   * @brief getNumberOfCorners gets a copy of _nCorners.
   * @return
   */
  n_corners_type getNumberOfCorners() const {
    return _nCorners;
  }

  /**
   * @brief getNumberOfCorners gets access to _nCorners.
   * @return
   */
  n_corners_type & getNumberOfCorners() {
    return _nCorners;
  }

  /**
   * @brief getBoundaryLen gets a copy of _boundaryLen.
   * @param boundaryLen
   */
  void getBoundaryLen(std::vector<num_type> &boundaryLen) const {
    boundaryLen = _boundaryLen;
  }

  /**
   * @brief getBoundaryLen gets access to _boundaryLen.
   * @return
   */
  std::vector<num_type> & getBoundaryLen() {
    return _boundaryLen;
  }

  /**
   * @brief geometrizePatch assigns vertex positions to the Patch patch in order to
   * lie on a 3D polygon defined by the patch itself.
   * @param patch The Patch to geometrize.
   * @param flows The flows in parameter space.
   * @return true if flows and polychords successfully match, false otherwise.
   */
  bool geometrizePatch(Patch &patch, const FlowMap &flows) const {
    // compute initial geometry
    PatchSmoothingFromBorder(patch.getMesh(), patch.getStartCorner(), _boundaryUV);
    // if no flows, ok with this, return
    if (flows.empty()) {
      // from UV parametrization to 3D coordinates
      if (_tMesh)
        PatchSmoothingFromUVto3D(patch.getMesh(), *_tMesh);
      // compute quad quality
      patch.setQuadsQuality(ComputeQuadsQuality(patch.getMesh()));
      // filter by quality
      if (patch.getQuadsQuality() > _qualityThreshold)
        return false;
      return true;
    }
    // otherwise, compute flow-polychord matches, and set the new flow distance
    if (ComputeFlowsPolychords(patch, flows)) {
      // if ok, update the geometry according to the matches
      if (_constrainedSmooth) {
        SmoothWithFlowConstraints(patch, flows, _lapWeight);
        // from UV parametrization to 3D coordinates
        if (_tMesh)
          PatchSmoothingFromUVto3D(patch.getMesh(), *_tMesh);
        // update polychords
        UpdateFlowPolychords(patch, _flowMap);
      } else {
        // from UV parametrization to 3D coordinates
        if (_tMesh)
          PatchSmoothingFromUVto3D(patch.getMesh(), *_tMesh);
      }
      // compute quad quality
      patch.setQuadsQuality(ComputeQuadsQuality(patch.getMesh()));
      // filter by quality
      if (patch.getQuadsQuality() > _qualityThreshold)
        return false;
      // then return
      return true;
    }
    return false;
  }

  /**
   * @brief PatchSmoothingFromUVto3D transforms a PolyMeshType UV coordinates into 3D coordinates, using
   * a TriMeshType in parameter space.
   * @param pMesh The input/output PolyMeshType.
   * @param tMesh The input TriMeshType in parameter space.
   */
  static void PatchSmoothingFromUVto3D(PolyMeshType &pMesh, const TriMeshType &tMesh) {
    TUVPointType uvPoint;
    TPointType point;
    TPointType normal;
    TFacePointer tFaceP;
    TPointType baryC;
    bool found = false;
    UVGrid<TriMeshType> uvGrid;
    uvGrid.Init(const_cast<TriMeshType &>(tMesh));
    for (size_t i = 0; i < pMesh.vert.size(); ++i)
      if (!pMesh.vert[i].IsD()) {
        uvPoint[0] = pMesh.vert[i].P()[0];
        uvPoint[1] = pMesh.vert[i].P()[1];
        found = uvGrid.CoordinatesPointUnique(uvPoint, point, normal);
        if (!found) {
          found = uvGrid.getClosest(uvPoint, tFaceP, baryC);
          assert(found);
          point = uvGrid.InterpolateRPos(tFaceP, baryC);
          normal = uvGrid.InterpolateRNormal(tFaceP, baryC);
        }
        pMesh.vert[i].P().Import(point);
        pMesh.vert[i].N().Import(normal);
      }
    vcg::tri::UpdateBounding<PolyMeshType>::Box(pMesh);
    vcg::tri::UpdateNormal<PolyMeshType>::PerPolygonalFaceNormalized(pMesh);
  }

  /**
   * @brief PatchSmoothingFromBorder assigns vertex positions with boundary and then smoothes the rest.
   * @param mesh The input mesh.
   * @param startPos The starting corner.
   * @param boundary The boundary positions to fix to.
   */
  static void PatchSmoothingFromBorder(PolyMeshType &mesh,
                                       const PosType &startPos,
                                       const std::vector<PointContainer> &boundary) {
    // check input
    if (mesh.IsEmpty())
      return;
    if (startPos.IsNull())
      return;
    if (boundary.size() < 1)
      return;

    // reset coords
    for (size_t v = 0; v < mesh.vert.size(); v++)
      mesh.vert[v].P().SetZero();

    // boundary on polygon
    BoundaryToBoundary(startPos, boundary);

    // set smoother parameters
    typename vcg::ImplicitSmoother<PolyMeshType>::Parameter parameters;
    parameters.useMassMatrix = false;
    parameters.useCotWeight = false;
    parameters.fixBorder = true;
    parameters.lapWeight = ScalarType(1);

    // smooth
    vcg::ImplicitSmoother<PolyMeshType>::Compute(mesh, parameters);

    vcg::tri::UpdateBounding<PolyMeshType>::Box(mesh);
    vcg::tri::UpdateNormal<PolyMeshType>::PerPolygonalFaceNormalized(mesh);
  }

  /**
   * @brief BoundaryToBoundary assigns border vertex positions to those of boundary.
   * @param startPos The starting corner.
   * @param boundary The boundary vertex positions.
   */
  static void BoundaryToBoundary(const PosType &startPos,
                                 const std::vector<PointContainer> &boundary) {
    // check input
    if (startPos.IsNull())
      return;
    if (boundary.size() < 1)
      return;

    PosType runPos = startPos;
    if (!runPos.IsBorder())
      return;

    size_t side = 0;
    num_type edge = 0;
    do {
      runPos.V()->P() = boundary[side][edge];
      runPos.FlipV();
      runPos.FlipE();
      if (!runPos.IsBorder()) {
        runPos.FlipF();
        runPos.FlipE();
        if (!runPos.IsBorder()) {   // concave corner
          while (!runPos.IsBorder()) {
            runPos.FlipF();
            runPos.FlipE();
          }
          side = (side + 1) % boundary.size();
          edge = 0;
        } else
          ++edge;
      } else {    // convex corner
        side = (side + 1) % boundary.size();
        edge = 0;
      }
    } while (runPos != startPos);
  }

  /**
   * @brief PointsToFlow converts points into FlowPoints to put into flowMap.
   * @param points The points of the flow.
   * @param isLoop true if it is a loop flow, false if it is border-to-border.
   * @param boundary The boundary points used to find the start and end sides of the flow.
   * @param flowMap The FlowMap where to store the FlowPoints.
   * @param nSamples Number of equally-spaced points the flow must be sampled with.
   */
  static void PointsToFlow(const PointContainer &points, const bool isLoop, const std::vector<PointContainer> &boundary,
                           FlowMap &flowMap, const count_type nSamples) {
    if (points.empty()) {
      std::cout << "No points to convert into flow." << std::endl;
      return;
    }
    if (nSamples < 2) {
      std::cout << "The number of samples must be greater than 1. Points not converted into flow" << std::endl;
      return;
    }

    ScalarType dist, bestDist;
    PointType tmpPoint, closestPoint;
    FlowTopology flowTopology;
    FlowPoints flowPoints;
    flowPoints.points = points;
    if (isLoop) {
      if (flowPoints.points.front() != flowPoints.points.back())
        flowPoints.points.push_back(flowPoints.points.front());
    } else {
      // find start side
      bestDist = std::numeric_limits<ScalarType>::max();
      for (n_corners_type s = 0; s < boundary.size(); ++s) {
        for (size_t p = 0; p < boundary[s].size() - 1; ++p) {
          vcg::SegmentPointDistance(vcg::Segment3<ScalarType>(boundary[s][p], boundary[s][p+1]),
                                    flowPoints.points.front(), tmpPoint, dist);
          if (dist < bestDist) {
            bestDist = dist;
            closestPoint = tmpPoint;
            flowTopology.startSide = s;
            flowPoints.startEdge = p;
          }
        }
        if (boundary[s].back() != boundary[(s + 1) % boundary.size()].front()) {
          vcg::SegmentPointDistance(vcg::Segment3<ScalarType>(boundary[s].back(), boundary[(s + 1) % boundary.size()].front()),
                                    flowPoints.points.front(), tmpPoint, dist);
          if (dist < bestDist) {
            bestDist = dist;
            closestPoint = tmpPoint;
            flowTopology.startSide = s;
            flowPoints.startEdge = boundary[s].size();
          }
        }
      }
      if (flowPoints.points.front() != closestPoint)
        flowPoints.points.push_front(closestPoint);
      // find end side
      bestDist = std::numeric_limits<ScalarType>::max();
      for (n_corners_type s = 0; s < boundary.size(); ++s) {
        for (size_t p = 0; p < boundary[s].size() - 1; ++p) {
          vcg::SegmentPointDistance(vcg::Segment3<ScalarType>(boundary[s][p], boundary[s][p+1]),
                                    flowPoints.points.back(), tmpPoint, dist);
          if (dist < bestDist) {
            bestDist = dist;
            closestPoint = tmpPoint;
            flowTopology.endSide = s;
            flowPoints.endEdge = p;
          }
        }
        if (boundary[s].back() != boundary[(s + 1) % boundary.size()].front()) {
          vcg::SegmentPointDistance(vcg::Segment3<ScalarType>(boundary[s].back(), boundary[(s + 1) % boundary.size()].front()),
                                    flowPoints.points.back(), tmpPoint, dist);
          if (dist < bestDist) {
            bestDist = dist;
            closestPoint = tmpPoint;
            flowTopology.endSide = s;
            flowPoints.endEdge = boundary[s].size();
          }
        }
      }
      if (flowPoints.points.back() != closestPoint)
        flowPoints.points.push_back(closestPoint);
      if (flowTopology.endSide < flowTopology.startSide) {
        flowTopology.swap();
        flowPoints.invertPoints();
      } else if (flowTopology.endSide == flowTopology.startSide && flowPoints.endEdge < flowPoints.startEdge) {
        flowPoints.invertPoints();
      }
    }
    flowPoints.updateArcLengths();
    flowPoints.computeSamples(nSamples);
    // put the FlowPoints into the FlowMap
    FlowMapIterator flowIt = flowMap.find(flowTopology);
    if (flowIt == flowMap.end())
      flowMap.insert(std::make_pair(flowTopology, FlowContainer(1, flowPoints)));
    else
      flowIt->second.push_back(flowPoints);
  }

private:
  using PatchFlattener::_constrainedSmooth;
  using PatchFlattener::_lapWeight;
  using PatchFlattener::_qualityThreshold;
  TriMeshType                              *_tMesh;           ///< The TriMeshType to retopologize.
  std::vector<PointContainer>               _boundary;        ///< The boundary vertex positions in 3D space.
  std::vector<PointContainer>               _boundaryUV;      ///< The boundary vertex positions in parameter space.
  FlowMap                                   _flowMap;         ///< The flows in 3D space.
  FlowMap                                   _flowMapUV;       ///< The flows in parameter space.
  n_corners_type                            _nCorners;        ///< The number of corners of the query.
  std::vector<num_type>                     _boundaryLen;     ///< The boundary subdivision of the query.
};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // GEOMETRY_SUPPORT_H
