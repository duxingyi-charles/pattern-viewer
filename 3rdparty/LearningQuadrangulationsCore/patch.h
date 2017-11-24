#ifndef PATCH_H
#define PATCH_H

#include <string>
#include <iomanip>
#include <sstream>
#include <queue>
#include <numeric>
#if __cplusplus >= 201103L || _MSC_VER >= 1800
#include <functional>
#endif
#include "polychord_support.h"

namespace vcg {
    namespace tri {
        
        // namespace pl: patch learning
        namespace pl {
            
            /**
             * @brief The Patch class keeps information about a polygonal patch.
             *
             * A Patch is a quad (sub)mesh representing a polygon, in a topological meaning.
             *
             * A single quad face is the most simple polygonal patch.
             * Any quad (sub)mesh is a n-sided polygon if all vertices on the boundary
             * have valence 3 (regular) or not 3 (corner).
             * Note that the valence of a boundary vertex is based on the edges belonging to the patch
             * (i.e. internal or boundary edges).
             *
             * Usage:
             * - Construct a patch from a mesh
             * \code
             * PolyMeshType mesh;
             * std::vector< vcg::face::Pos<PolyMeshType::FaceType> > corners;
             * n_corners_type nCorners;
             * num_type nFaces, nVertices;
             * std::vector< std::pair<valence_type,num_type> > singularityVec;
             * ... // define the mesh and other data here
             * Patch<PolyMeshType> patch;
             * patch.setPatch(mesh, corners, nCorners, nFaces, nVertices, singularityVec);
             * \endcode
             * - Construct a patch from known data
             * Just use the \code Patch<PolyMeshType>::setSomething(something); \endcode functions.
             * - Read data from a patch
             * Just use the \code Patch<PolyMeshType>::getSomething(); \endcode functions. If you want to
             * reconstruct a mesh use \code Patch<PolyMeshType>::reconstructPatch(); \endcode .
             */
            template < typename PolyMeshType >
            class Patch {
            public:
                // data types and members (public):
                typedef typename PolyMeshType::ScalarType                         ScalarType;
                typedef typename PolyMeshType::VertexType                         VertexType;
                typedef typename PolyMeshType::VertexPointer                      VertexPointer;
                typedef typename PolyMeshType::FaceType                           FaceType;
                typedef typename PolyMeshType::FacePointer                        FacePointer;
                typedef std::vector<var_label_type>                               BoundaryIDsOnSideType;
                typedef std::vector<BoundaryIDsOnSideType>                        BoundaryIDsVecType;
                typedef std::pair<n_corners_type, n_corners_type>                 BoundaryIDHeadTailType;
                typedef std::vector<BoundaryIDHeadTailType>                       BoundaryIDsEndsVecType;
                typedef pl::PolychordSupport<PolyMeshType>                        PolychordSupport;
                typedef typename PolychordSupport::PosType                        PosType;
                typedef typename PolychordSupport::PolychordMap                   PolychordMap;
                typedef typename PolychordSupport::PolychordMapIterator           PolychordMapIterator;
                typedef typename PolychordSupport::PolychordMapConstIterator      PolychordMapConstIterator;
                typedef typename MeshSupport<PolyMeshType>::HalfEdge              HalfEdge;
#if __cplusplus >= 201103L || _MSC_VER >= 1800
                typedef typename MeshSupport<PolyMeshType>::HalfEdgeHash          HalfEdgeHash;
#endif
                
                // data methods (public):
                /**
                 * @brief Patch Default constructor.
                 */
                Patch() : _nSingularities(0),
                _perimeter(0),
                _nCorners(0),
                _nFaces(0),
                _nVertices(0),
                _nRotations(0),
                _nOccurrences(0),
                _nConcaveCorners(0),
                _patchMesh(0),
                _startPos(0),
                _flowDistance(std::numeric_limits<ScalarType>::max()),
                _quadsQuality(std::numeric_limits<ScalarType>::max()) { }
                
                /**
                 * @brief Patch Copy constructor.
                 * @param p
                 */
                Patch(const Patch &p) : _topologyStr(p._topologyStr),
                _expandedTopologyStr(p._expandedTopologyStr),
                _boundaryIDsVec(p._boundaryIDsVec),
                _boundaryIDsEndsVec(p._boundaryIDsEndsVec),
                _boundaryVarVec(p._boundaryVarVec),
                _boundaryLenVec(p._boundaryLenVec),
                _singularityStr(p._singularityStr),
                _nSingularities(p._nSingularities),
                _perimeter(p._perimeter),
                _nCorners(p._nCorners),
                _nFaces(p._nFaces),
                _nVertices(p._nVertices),
                _nRotations(p._nRotations),
                _nOccurrences(p._nOccurrences),
                _nConcaveCorners(p._nConcaveCorners),
                _patchMesh(0),
                _startPos(0),
                _polychordMap(p._polychordMap),
                _flowDistance(p._flowDistance),
                _quadsQuality(p._quadsQuality) {
                    _patchMesh = new PolyMeshType;
                    _startPos = new vcg::face::Pos<FaceType>;
                    p.getMeshAndStartCorner(*_patchMesh, *_startPos);
                    p.getPolychordMap(*_patchMesh, _polychordMap);
                }
                
#if __cplusplus >= 201103L || _MSC_VER >= 1800
                /**
                 * @brief Patch Move constructor.
                 * @param p The Patch to move and reset.
                 */
                Patch(Patch&& p) : _topologyStr(std::move(p._topologyStr)),
                _expandedTopologyStr(std::move(p._expandedTopologyStr)),
                _boundaryIDsVec(std::move(p._boundaryIDsVec)),
                _boundaryIDsEndsVec(std::move(p._boundaryIDsEndsVec)),
                _boundaryVarVec(std::move(p._boundaryVarVec)),
                _boundaryLenVec(std::move(p._boundaryLenVec)),
                _singularityStr(std::move(p._singularityStr)),
                _nSingularities(p._nSingularities),
                _perimeter(p._perimeter),
                _nCorners(p._nCorners),
                _nFaces(p._nFaces),
                _nVertices(p._nVertices),
                _nRotations(p._nRotations),
                _nOccurrences(p._nOccurrences),
                _nConcaveCorners(p._nConcaveCorners),
                _patchMesh(p._patchMesh),
                _startPos(p._startPos),
                _polychordMap(std::move(p._polychordMap)),
                _flowDistance(p._flowDistance),
                _quadsQuality(p._quadsQuality) {
                    p._patchMesh = 0;
                    p._startPos = 0;
                    p.reset();
                }
#endif
                
                /**
                 * @brief ~Patch Destructor.
                 */
                ~Patch() {
                    reset();
                }
                
                /**
                 * @brief operator = Copy assignment operator.
                 * @param p The Patch to copy.
                 * @return A reference to this.
                 */
                Patch & operator = (const Patch &p) {
                    if (&p != this) {
                        _topologyStr = p._topologyStr;
                        _expandedTopologyStr = p._expandedTopologyStr;
                        _boundaryIDsVec = p._boundaryIDsVec;
                        _boundaryIDsEndsVec = p._boundaryIDsEndsVec;
                        _boundaryVarVec = p._boundaryVarVec;
                        _boundaryLenVec = p._boundaryLenVec;
                        _singularityStr = p._singularityStr;
                        _nSingularities = p._nSingularities;
                        _perimeter = p._perimeter;
                        _nCorners = p._nCorners;
                        _nFaces = p._nFaces;
                        _nVertices = p._nVertices;
                        _nRotations = p._nRotations;
                        _nOccurrences = p._nOccurrences;
                        _nConcaveCorners = p._nConcaveCorners;
                        if (_patchMesh == 0)
                            _patchMesh = new PolyMeshType;
                        if (_startPos == 0)
                            _startPos = new PosType;
                        p.getMeshAndStartCorner(*_patchMesh, *_startPos);
                        p.getPolychordMap(*_patchMesh, _polychordMap);
                        _flowDistance = p._flowDistance;
                        _quadsQuality = p._quadsQuality;
                    }
                    return *this;
                }
                
#if __cplusplus >= 201103L || _MSC_VER >= 1800
                /**
                 * @brief operator = Move assignment operator.
                 * @param p The Patch to move and reset.
                 * @return A reference to this.
                 */
                Patch & operator = (Patch&& p) {
                    if (&p != this) {
                        _topologyStr = std::move(p._topologyStr);
                        _expandedTopologyStr = std::move(p._expandedTopologyStr);
                        _boundaryIDsVec = std::move(p._boundaryIDsVec);
                        _boundaryIDsEndsVec = std::move(p._boundaryIDsEndsVec);
                        _boundaryVarVec = std::move(p._boundaryVarVec);
                        _boundaryLenVec = std::move(p._boundaryLenVec);
                        _singularityStr = std::move(p._singularityStr);
                        _nSingularities = p._nSingularities;
                        _perimeter = p._perimeter;
                        _nCorners = p._nCorners;
                        _nFaces = p._nFaces;
                        _nVertices = p._nVertices;
                        _nRotations = p._nRotations;
                        _nOccurrences = p._nOccurrences;
                        _nConcaveCorners = p._nConcaveCorners;
                        _patchMesh = p._patchMesh;
                        _startPos = p._startPos;
                        _polychordMap = std::move(p._polychordMap);
                        _flowDistance = p._flowDistance;
                        _quadsQuality = p._quadsQuality;
                        p._patchMesh = 0;
                        p._startPos = 0;
                        p.reset();
                    }
                    return *this;
                }
#endif
                
                /**
                 * @brief operator < determines the relational order between this patch and p.
                 * @param p The patch to compare.
                 * @return true if this patch < p, false otherwise.
                 */
                bool operator < (const Patch &p) const {
                    if (_expandedTopologyStr == p._expandedTopologyStr)
                        return false;
                    if (_nSingularities != p._nSingularities)
                        return _nSingularities < p._nSingularities;
                    if (_nOccurrences != p._nOccurrences)
                        return _nOccurrences > p._nOccurrences;
                    if (_nFaces != p._nFaces)
                        return _nFaces < p._nFaces;
                    if (_nVertices != p._nVertices)
                        return _nVertices < p._nVertices;
                    if (_nCorners != p._nCorners)
                        return _nCorners < p._nCorners;
                    if (_perimeter != p._perimeter)
                        return _perimeter < p._perimeter;
                    if (_boundaryLenVec != p._boundaryLenVec)
                        return _boundaryLenVec < p._boundaryLenVec;
                    if (_topologyStr != p._topologyStr)
                        return _topologyStr < p._topologyStr;
                    if (_boundaryIDsVec != p._boundaryIDsVec)
                        return _boundaryIDsVec < p._boundaryIDsVec;
                    if (_nRotations != p._nRotations)
                        return _nRotations < p._nRotations;
                    return _boundaryVarVec < p._boundaryVarVec;
                }
                
                /**
                 * @brief operator > determines the relational order between this patch and p.
                 * @param p The patch to compare.
                 * @return true if this patch > p, false otherwise.
                 */
                bool operator > (const Patch &p) const {
                    return p < *this;
                }
                
                /**
                 * @brief operator <= determines the relational order between this patch and p.
                 * @param p The patch to compare.
                 * @return true if this patch <= p, false otherwise.
                 */
                bool operator <= (const Patch &p) const {
                    return !(p < *this);
                }
                
                /**
                 * @brief operator >= determines the relational order between this patch and p.
                 * @param p The patch to compare.
                 * @return true if this patch >= p, false otherwise.
                 */
                bool operator >= (const Patch &p) const {
                    return !(*this < p);
                }
                
                /**
                 * @brief operator == determines the relational order between this patch and p.
                 * @param p The patch to compare.
                 * @return true if this patch == p, false otherwise.
                 */
                bool operator == (const Patch &p) const {
                    return !(*this < p || p < *this);
                }
                
                /**
                 * @brief operator != determines the relational order between this patch and p.
                 * @param p The patch to compare.
                 * @return true if this patch == p, false otherwise.
                 */
                bool operator != (const Patch &p) const {
                    return *this < p || p < *this;
                }
                
                // access methods:
                /**
                 * @brief isNull says whether this patch is empty or not.
                 * @return true if it is empty, false otherwise.
                 */
                bool isNull() const {
                    return _topologyStr.empty()
                    || _boundaryIDsVec.empty()
                    || _boundaryLenVec.empty()
                    || _boundaryVarVec.empty()
                    || _perimeter == 0
                    || _nFaces == 0
                    || _nVertices == 0;
                }
                
                /**
                 * @brief getTopologyStr See _topologyStr.
                 * @return
                 */
                const std::string & getTopologyStr() const {
                    return _topologyStr;
                }
                
                /**
                 * @brief getExpandedTopologyStr See _expandedTopologyStr.
                 * @return
                 */
                const std::string & getExpandedTopologyStr() const {
                    return _expandedTopologyStr;
                }
                
                /**
                 * @brief getBoundaryIDsVec See _boundaryIDsVec.
                 * @return
                 */
                const BoundaryIDsVecType & getBoundaryIDsVec() const {
                    return _boundaryIDsVec;
                }
                
                /**
                 * @brief getBoundaryIDsStr See _boundaryIDsVec.
                 * @param boundaryIDsStr
                 */
                void getBoundaryIDsStr(std::string &boundaryIDsStr) const {
                    std::stringstream stream;
                    for (size_t side = 0; side < _boundaryIDsVec.size(); ++side) {
                        for (size_t edge = 0; edge < _boundaryIDsVec[side].size(); ++edge) {
                            // write it into the string
                            stream << VARLABEL;
                            stream << std::setw(2 * sizeof(var_label_type)) << std::setfill('0') << std::hex << std::noshowbase;
                            stream << (num_type)_boundaryIDsVec[side][edge];
                        }
                        if (side < _boundaryIDsVec.size() - 1)
                            stream << SEPARATOR;
                    }
                    boundaryIDsStr = stream.str();
                }
                
                /**
                 * @brief getBoundaryVarVec See _boundaryVarVec.
                 * @return
                 */
                const std::vector<num_type> & getBoundaryVarVec() const {
                    return _boundaryVarVec;
                }
                
                /**
                 * @brief getBoundaryVarStr See _boundaryVarVec.
                 * @param boundaryVarStr
                 */
                void getBoundaryVarStr(std::string &boundaryVarStr) const {
                    std::stringstream stream;
                    for (size_t i = 0; i < _boundaryVarVec.size(); ++i) {
                        stream << VARLABEL;
                        stream << std::setw(2 * sizeof(var_label_type)) << std::setfill('0') << std::hex << std::noshowbase;
                        stream << i;
                        stream << VARSEPARATOR;
                        stream << std::setw(2 * sizeof(num_type)) << std::setfill('0') << std::hex << std::noshowbase;
                        stream << _boundaryVarVec[i];
                        if (i < _boundaryVarVec.size() - 1)
                            stream << SEPARATOR;
                    }
                    boundaryVarStr = stream.str();
                }
                
                /**
                 * @brief getBoundaryLenVec See _boundaryLenVec.
                 * @return
                 */
                const std::vector<num_type> & getBoundaryLenVec() const {
                    return _boundaryLenVec;
                }
                
                /**
                 * @brief getBoundaryLenStr See _boundaryLenVec.
                 * @param boundaryLenStr
                 */
                void getBoundaryLenStr(std::string &boundaryLenStr) const {
                    std::stringstream stream;
                    for (size_t i = 0; i < _boundaryLenVec.size(); ++i) {
                        stream << std::setw(2 * sizeof(num_type)) << std::setfill('0') << std::hex << std::noshowbase;
                        stream << _boundaryLenVec[i];
                    }
                    boundaryLenStr = stream.str();
                }
                
                /**
                 * @brief getBoundaryIDsEnds See _boundaryIDsEndsVec.
                 * @return
                 */
                const BoundaryIDsEndsVecType & getBoundaryIDsEnds() const {
                    return _boundaryIDsEndsVec;
                }
                
                /**
                 * @brief getBoundaryIDsEndsStr See _boundaryIDsEndsVec.
                 * @param boundaryIDsEndsStr
                 */
                void getBoundaryIDsEndsStr(std::string &boundaryIDsEndsStr) const {
                    std::stringstream stream;
                    for (size_t i = 0; i < _boundaryIDsEndsVec.size(); ++i) {
                        stream << VARLABEL;
                        stream << std::setw(2 * sizeof(var_label_type)) << std::setfill('0') << std::hex << std::noshowbase;
                        stream << i;
                        stream << VARSEPARATOR;
                        stream << std::setw(2 * sizeof(num_type)) << std::setfill('0') << std::hex << std::noshowbase;
                        stream << _boundaryIDsEndsVec[i].first;
                        stream << VARSEPARATOR;
                        stream << std::setw(2 * sizeof(num_type)) << std::setfill('0') << std::hex << std::noshowbase;
                        stream << _boundaryIDsEndsVec[i].second;
                        if (i < _boundaryIDsEndsVec.size() - 1)
                            stream << SEPARATOR;
                    }
                    boundaryIDsEndsStr = stream.str();
                }
                
                ///dxy add
                std::string getBoundaryTopoStr() {
                    return _boundaryTopoStr;
                }
                ///dxy add end
                
                /**
                 * @brief getSingularityStr See _singularityStr.
                 * @return
                 */
                const std::string & getSingularityStr() const {
                    return _singularityStr;
                }
                
                /**
                 * @brief getSingularityVec extracts a vector of pairs <valence_V,number_of_vertices_with_valence_V> from _singularityStr.
                 * @param singularityVec
                 */
                void getSingularityVec(std::vector< std::pair<valence_type,num_type> > &singularityVec) const {
                    std::stringstream stream;
                    singularityVec.clear();
                    singularityVec.resize((_singularityStr.length() + 1/*last SEPARATOR is absent*/) / (1/*VALENCELABEL*/ + 2 * sizeof(valence_type) + 1/*VALENCESEP*/ + 2 * sizeof(num_type)/*num*/ + 1/*SEPARATOR*/));
                    for (size_t i = 1/*VALENCELABEL*/; i < _singularityStr.size(); ) {
                        stream.str(_singularityStr.substr(i, 2 * sizeof(valence_type)) + std::string(" "));
                        stream >> std::hex >> singularityVec[i / (1/*VALENCELABEL*/ + 2 * sizeof(valence_type)/*valence*/ + 1/*VALENCESEP*/ + 2 * sizeof(num_type)/*num*/ + 1/*SEPARATOR*/)].first;
                        i += 2 * sizeof(valence_type)/*valence*/ + 1/*VALENCESEP*/;
                        stream.str(_singularityStr.substr(i, 2 * sizeof(num_type)) + std::string(" "));
                        stream >> std::hex >> singularityVec[i / (1/*VALENCELABEL*/ + 2 * sizeof(valence_type)/*valence*/ + 1/*VALENCESEP*/ + 2 * sizeof(num_type)/*num*/ + 1/*SEPARATOR*/)].second;
                        i += 2 * sizeof(num_type)/*num*/ + 1/*SEPARATOR*/ + 1/*VALENCELABEL*/;
                    }
                }
                
                /**
                 * @brief getNumberOfSingularities See _nSingularities.
                 * @return
                 */
                num_type getNumberOfSingularities() const {
                    return _nSingularities;
                }
                
                /**
                 * @brief getPerimeter See _perimeter.
                 * @return
                 */
                num_type getPerimeter() const {
                    return _perimeter;
                }
                
                /**
                 * @brief getNumberOfCorners See _nCorners.
                 * @return
                 */
                n_corners_type getNumberOfCorners() const {
                    return _nCorners;
                }
                
                /**
                 * @brief getNumberOfFaces See _nFaces.
                 * @return
                 */
                num_type getNumberOfFaces() const {
                    return _nFaces;
                }
                
                /**
                 * @brief getNumberOfVertices See _nVertices.
                 * @return
                 */
                num_type getNumberOfVertices() const {
                    return _nVertices;
                }
                
                /**
                 * @brief getNumberOfRotations See _nRotations.
                 * @return
                 */
                n_corners_type getNumberOfRotations() const {
                    return _nRotations;
                }
                
                /**
                 * @brief getNumberOfOccurrences See _nOccurrences.
                 * @return
                 */
                count_type getNumberOfOccurrences() const {
                    return _nOccurrences;
                }
                
                /**
                 * @brief getNumberOfConcaveCorners See _nConcaveCorners.
                 * @return
                 */
                n_corners_type getNumberOfConcaveCorners() const {
                    return _nConcaveCorners;
                }
                
                /**
                 * @brief getMeshAndStartCorner copies _patchMesh and _startPos into m and sp.
                 * @param m The new mesh.
                 * @param sp The new start Pos.
                 */
                void getMeshAndStartCorner(PolyMeshType &m, vcg::face::Pos<FaceType> &sp) const {
                    m.Clear();
                    sp.SetNull();
                    if (_patchMesh && !_patchMesh->IsEmpty() && _startPos && !_startPos->IsNull()) {
                        vcg::tri::Append<PolyMeshType, PolyMeshType>::MeshCopy(m, const_cast<PolyMeshType &>(*_patchMesh), false, true);
                        vcg::tri::UpdateTopology<PolyMeshType>::FaceFace(m);
                        size_t fInd = vcg::tri::Index(*_patchMesh, _startPos->F());
                        size_t vInd = vcg::tri::Index(*_patchMesh, _startPos->V());
                        sp.Set(&m.face[fInd], _startPos->E(), &m.vert[vInd]);
                    }
                }
                
                /**
                 * @brief getMesh accesses _patchMesh.
                 * @return A reference to the mesh.
                 */
                PolyMeshType & getMesh() {
                    if (_patchMesh == 0 || _startPos == 0) {
                        reconstructPatch();
                    }
                    return *_patchMesh;
                }
                
                /**
                 * @brief getStartCorner accesses _startPos.
                 * @return A reference to the start corner.
                 */
                vcg::face::Pos<FaceType> & getStartCorner() {
                    if (_patchMesh == 0 || _startPos == 0)
                        reconstructPatch();
                    return *_startPos;
                }
                
                /**
                 * @brief getPolychordMap See _polychordMap.
                 * @param mesh The mesh (copy of this) used to convert pointers.
                 * @param polychordMap The output PolychordMap.
                 */
                void getPolychordMap(PolyMeshType &mesh, PolychordMap &polychordMap) const {
                    polychordMap.clear();
                    if (mesh.IsEmpty())
                        return;
                    polychordMap = _polychordMap;
                    size_t fInd = 0, vInd = 0;
                    for (PolychordMapIterator it = polychordMap.begin(); it != polychordMap.end(); ++it)
                        for (size_t polyInd = 0; polyInd < it->second.size(); ++polyInd) {
                            fInd = vcg::tri::Index(*_patchMesh, it->second[polyInd].startPos.F());
                            vInd = vcg::tri::Index(*_patchMesh, it->second[polyInd].startPos.V());
                            it->second[polyInd].startPos.Set(&mesh.face[fInd], it->second[polyInd].startPos.E(), &mesh.vert[vInd]);
                        }
                }
                
                /**
                 * @brief getPolychordMap See _polychordMap.
                 * @return
                 */
                PolychordMap & getPolychordMap() {
                    return _polychordMap;
                }
                
                /**
                 * @brief getFlowDistance See _flowDistance.
                 * @return
                 */
                ScalarType getFlowDistance() const {
                    return _flowDistance;
                }
                
                /**
                 * @brief getQuadsQuality See _quadsQuality.
                 * @return
                 */
                ScalarType getQuadsQuality() const {
                    return _quadsQuality;
                }
                
                /**
                 * @brief reconstructTopologyFromString translates a string of (hexadecimal) numbers into a vector of adjacency.
                 *
                 * Assuming topology to be a string of hexadecimal numbers, every i-th subsequence of 4 numbers represents the indices
                 * of the 4 faces adjacent to the i-th face. Each number is stored without the prefix '0x', and corresponds
                 * to a not negative number in the range [0,<max of num_type>], with <max of num_type> meaning a border edge.
                 *
                 * The resulting vector has a cell for each face. The i-th cell is a vector of 4 elements. The j-th element (0 <= j <= 3)
                 * is the index of the face adjacent to the j-th edge.
                 *
                 * @param out The output vector of adjacency.
                 */
                void reconstructTopologyFromString(std::vector< std::vector<num_type> > &out) const {
                    std::stringstream stream;
                    std::string str;
                    num_type lastSeen = 0;
                    out.assign(_topologyStr.size() / (3 * 2 * sizeof(num_type)), std::vector<num_type>(4, -1));
                    for (size_t i = 0; i < out.size(); ++i) {
                        // out[i][0] is already known
                        for (size_t j = 1; j < 4; ++j) {
                            str = _topologyStr.substr(i * 3 * 2 * sizeof(num_type) + (j-1) * 2 * sizeof(num_type), 2 * sizeof(num_type));
                            if (str != std::string(2 * sizeof(num_type), BORDER)) {
                                stream.str(str + std::string(" "));
                                stream >> std::hex >> out[i][j];
                                if (out[i][j] > lastSeen) {
                                    lastSeen = out[i][j];
                                    out[lastSeen][0] = i;
                                }
                            } else {
                                out[i][j] = -1;
                            }
                        }
                    }
                }
                
                /**
                 * @brief reconstructPatch creates a vcg PolyMeshType mesh from this patch's information, also taking into account
                 * _boundaryVarStr values to split polychords and to have the full patch.
                 * @param patchMesh The output vcg PolyMeshType mesh.
                 * @param startPos The output starting corner.
                 */
                void reconstructPatch(PolyMeshType &patchMesh, PosType &startPos) const {
                    ///dxy test
//                    std::cout << "reconstructPatch() begin..." << std::endl;
                    ///
                    
                    vcg::tri::RequireFFAdjacency(patchMesh);
                    std::vector< std::vector<num_type> > ffAdj;
                    std::vector<bool> hasVertexesAssigned;
                    num_type lastAssignedVertex = 0;
                    num_type currentFace, currentEdge, tmp, ffi;
                    patchMesh.Clear();
                    startPos.SetNull();
                    
                    if (isNull())
                        return;
                    
                    // allocate vertexes and faces
                    vcg::tri::Allocator<PolyMeshType>::AddVertices(patchMesh, _nVertices);
                    if ((num_type)patchMesh.VN() != _nVertices) {
                        std::cout << "Error while allocating vertices. Patch not reconstructed." << std::endl;
                        patchMesh.Clear();
                        return;
                    }
                    vcg::tri::Allocator<PolyMeshType>::AddFaces(patchMesh, _nFaces);
                    if ((num_type)patchMesh.FN() != _nFaces) {
                        std::cout << "Error while allocating faces. Patch not reconstructed." << std::endl;
                        patchMesh.Clear();
                        return;
                    }
                    
                    // decode the topology
                    reconstructTopologyFromString(ffAdj);
                    hasVertexesAssigned.assign(_nFaces, false);
                    
                    // apply the decoded topology to the mesh
                    for (num_type i = 0; i < _nFaces; ++i) {
                        // allocate memory for adjacency
                        patchMesh.face[i].Alloc(4);
                        // set the adjacency
                        for (num_type j = 0; j < 4; ++j) {
                            if (ffAdj[i][j] >= 0) {
                                // set face pointer
                                patchMesh.face[i].FFp(j) = &patchMesh.face[ffAdj[i][j]];
                                // set face corner index:
                                // 1) find a common edge
                                ffi = 0;
                                while (ffi < 4 && ffAdj[ffAdj[i][j]][ffi] != i)
                                    ++ffi;
                                assert(ffi < 4);
                                // 2) find LAST common edge
                                while (ffAdj[ffAdj[i][j]][(ffi + 1) % 4] == i)
                                    ffi = (ffi + 1) % 4;
                                // 3) go back a number of common edges equal to how man previous common edge has been already processed
                                for (num_type k = (j + 4 - 1) % 4; k < 4 && ffAdj[i][k] == ffAdj[i][j]; k = (k + 4 - 1) % 4)
                                    ffi = (ffi + 4 - 1) % 4;
                                // 4) store this index
                                patchMesh.face[i].FFi(j) = ffi;
                            } else {
                                patchMesh.face[i].FFp(j) = &patchMesh.face[i];
                                patchMesh.face[i].FFi(j) = j;
                            }
                        }
                    }
                    // set face-vertex adjacency
                    for (num_type i = 0; i < _nFaces; ++i) {
                        // set the adjacency
                        for (num_type j = 0; j < 4; ++j) {
                            // check if an adjacent face has already assigned its vertexes
                            currentFace = (num_type)vcg::tri::Index(patchMesh, patchMesh.face[i].FFp(j));
                            currentEdge = (num_type)patchMesh.face[currentFace].Next(patchMesh.face[i].FFi(j));
                            while (!hasVertexesAssigned[currentFace] && currentFace != i
                                   && !vcg::face::IsBorder(patchMesh.face[currentFace], currentEdge)) {
                                tmp = currentFace;
                                currentFace = (num_type)vcg::tri::Index(patchMesh, patchMesh.face[currentFace].FFp(currentEdge));
                                currentEdge = (num_type)patchMesh.face[currentFace].Next(patchMesh.face[tmp].FFi(currentEdge));
                            }
                            if (!hasVertexesAssigned[currentFace]) {
                                currentFace = (num_type)vcg::tri::Index(patchMesh, patchMesh.face[i].FFp(patchMesh.face[i].Prev(j)));
                                currentEdge = (num_type)patchMesh.face[i].FFi(patchMesh.face[i].Prev(j));
                                while (!hasVertexesAssigned[currentFace] && currentFace != i
                                       && !vcg::face::IsBorder(patchMesh.face[currentFace], patchMesh.face[currentFace].Prev(currentEdge))) {
                                    tmp = currentFace;
                                    currentFace = (num_type)vcg::tri::Index(patchMesh, patchMesh.face[currentFace].FFp(patchMesh.face[currentFace].Prev(currentEdge)));
                                    currentEdge = (num_type)patchMesh.face[tmp].FFi(patchMesh.face[tmp].Prev(currentEdge));
                                }
                            }
                            // if the corresponding vertex has already been assigned
                            if (hasVertexesAssigned[currentFace])
                                // then assign the same vertex
                                // even though currentFace==i, it can't be catched here because hasVertexesAssigned[i]==false
                                patchMesh.face[i].V(j) = patchMesh.face[currentFace].V(currentEdge);
                            else
                                patchMesh.face[i].V(j) = &patchMesh.vert[lastAssignedVertex++];
                        }
                        // now this face has all vertexes assigned
                        hasVertexesAssigned[i] = true;
                    }
                    
                    // find the polychords' starting positions
                    startPos.Set(&patchMesh.face[0], 0, patchMesh.face[0].V(0));
                    
                    // extract variables values
                    PosType runPos = startPos;
                    // take care of updating pointers while reallocating the structures
                    std::vector<FacePointer *> facesToUpdate(1);
                    std::vector<VertexPointer *> verticesToUpdate(1);
                    facesToUpdate[0] = &startPos.F();
                    verticesToUpdate[0] = &startPos.V();
                    // run on the border and split the polychords when needed
                    num_type lastVar = -1;
                    size_t side = 0;
                    num_type edge = 0;
                    
                    
                    
                    ///get corners[_nRotations] in the new mesh
                    ///dxy add
                    //      if (!_boundaryTopoStr.empty()) {
                    //          // 1. move startPos to border
                    //          bool flipped = false;
                    //          runPos = startPos;
                    //          if (!runPos.IsBorder()) {
                    //              while (!runPos.IsBorder()) {
                    //                  runPos.FlipE();
                    //                  runPos.FlipV();
                    //                  runPos.FlipE();
                    //                  runPos.FlipF();
                    //                  if (runPos == startPos) {
                    //                      if (flipped) {
                    //                          std::cout << "reconstructPatch Error: can't find border" << std::endl;
                    //                          return;
                    //                      }
                    //                      else {
                    //                          runPos.FlipE();
                    //                          flipped = true;
                    //                      }
                    //                  }
                    //              }
                    //          }
                    //          startPos = runPos;  //on border
                    //          // 2. build boundary TopoStr from startPos
                    //          std::string btstr;
                    //          buildBoundaryTopoStr(patchMesh, startPos, btstr);
                    //          // 3. find the shift to match _boundaryTopoStr
                    //          int strlen = _boundaryTopoStr.size();
                    //          int shift = 0;
                    //          bool reversed = true;
                    //          for (; shift < strlen; ++shift) {
                    //              std::string ts = _boundaryTopoStr.substr(shift, strlen-shift) + _boundaryTopoStr.substr(0, shift);
                    //              if (ts == btstr) {
                    //                  reversed = false;
                    //                  break;
                    //              }
                    //          }
                    //          if (reversed) {
                    //              //reverse startPos
                    //              startPos.FlipE();
                    //              while(!startPos.IsBorder()) {
                    //                  startPos.FlipF();
                    //                  startPos.FlipE();
                    //              }
                    //              //
                    //              bool match = false;
                    //              buildBoundaryTopoStr(patchMesh, startPos, btstr);
                    //              shift = 0;
                    //              for (; shift < strlen; ++shift) {
                    //                  std::string ts = _boundaryTopoStr.substr(shift, strlen-shift) + _boundaryTopoStr.substr(0, shift);
                    //                  if (ts == btstr) {
                    //                      match = true;
                    //                      break;
                    //                  }
                    //              }
                    //              if (!match) {
                    //                  std::cout << "reconstructPatch Error: can't match _boundaryTopoStr" << std::endl;
                    //                  return;
                    //              }
                    //          }
                    //          // 4. shift the startPos to proper pos
                    //          while(shift > 0) {
                    //              startPos.FlipE();
                    //              while (!startPos.IsBorder()) {
                    //                  startPos.FlipF();
                    //                  startPos.FlipE();
                    //              }
                    //              startPos.FlipV();
                    //              shift--;
                    //          }
                    //          // 5. test: check if topostr match
                    //          buildBoundaryTopoStr(patchMesh, startPos, btstr);
                    //          if (btstr != _boundaryTopoStr) {
                    //              std::cout << "reconstructPatch Error: wrong, check your code" << std::endl;
                    //              return;
                    //          }
                    //          // 6. set runPos
                    //          runPos = startPos;
                    //      }
                    ///dxy add end
                    
                    
                    ///dxy add
                    /* calc expanded boundaryLenVec */
                    auto expandLenVec = _boundaryLenVec;
//                    std::cout << "Patch sides num: " << expandLenVec.size() << std::endl;
                    //      for(int i=0; i<_boundaryVarVec.size(); ++i) {
                    //          expandLenVec[_boundaryIDsEndsVec[i].first]  += _boundaryVarVec[i] - 1;
                    //          expandLenVec[_boundaryIDsEndsVec[i].second] += _boundaryVarVec[i] - 1;
                    //      }
                    
                    ///dxy add end
                    
                    
                    ///dxy test
//                    std::cout << "do{...} begin :" << std::endl;
                    ///
                    
                    do {
                        ///dxy test
//                        std::cout << "side: " << side << ", edge: " << edge << std::endl;
                        ///
                        // split the polychord, if needed
                        if (static_cast<num_type>(_boundaryIDsVec[side][edge]) > lastVar) {
                            lastVar = _boundaryIDsVec[side][edge];
                            if (_boundaryVarVec[lastVar] > 1) {
                                ///dxy test
                                std::cout << "going to split polychord..." << std::endl;
                                //              std::cout << "side: " << side << ", edge: " << edge << std::endl;
                                std::cout << "lastVar: " << lastVar << std::endl;
                                std::cout << "_boundaryVarVec[lastVar] = " << _boundaryVarVec[lastVar] << std::endl;
                                ///
                                vcg::tri::PolychordCollapse<PolyMeshType>::SplitPolychord(patchMesh, runPos, _boundaryVarVec[lastVar],
                                                                                          facesToUpdate, verticesToUpdate);
                            }
                        }
                        // jump to the last inserted polychord
                        for (num_type i = 1; i < _boundaryVarVec[_boundaryIDsVec[side][edge]]; ++i) {
                            ///dxy test
//                            if (i == 1) { //first iter
//                                std::cout << "jump to the last inserted polychord..." << std::endl;
//                                std::cout << "_boundaryVarVec[_boudaryIDsVec[" << side << "][" << edge << "]] = ";
//                                std::cout << _boundaryVarVec[_boundaryIDsVec[side][edge]] << std::endl;
//                            }
                            ///
                            runPos.FlipV();
                            runPos.FlipE();
                            runPos.FlipF();
                            runPos.FlipE();
                            
                            ///dxy test
//                            if (!runPos.IsBorder()) {
//                                std::cout << "Warning: runPos should be on border ! ! !" << std::endl;
//                            }
                            
                            ///dxy test
//                            if (i == _boundaryVarVec[_boundaryIDsVec[side][edge]] - 1) {
//                                std::cout << "jump end" << std::endl;
//                            }
                        }
                        
                        ///TODO: run around border based on boundaryLenVec
                        
                        // go to the next polychord
                        runPos.FlipV();
                        runPos.FlipE();
                        ///dxy add
                        //        if (!_boundaryTopoStr.empty()) {
                        //        if(_bts) {
                        if (edge == expandLenVec[side] - 1) { //last edge on this side
                            while (!runPos.IsBorder()) {
                                runPos.FlipF();
                                runPos.FlipE();
                            }
                            side += 1;
                            edge = 0;
                        }
                        else {
                            while (!runPos.IsBorder()) {
                                runPos.FlipF();
                                runPos.FlipE();
                            }
                            edge += 1;
                        }
                        //continue;
                        //        }
                        ///dxy add end
                        
                        ///dxy delete
                        //      if (!runPos.IsBorder()) {
                        //        runPos.FlipF();
                        //        runPos.FlipE();
                        //        ++edge;
                        //        if (!runPos.IsBorder()) {     // concave corner
                        //          while (!runPos.IsBorder()) {
                        //            runPos.FlipF();
                        //            runPos.FlipE();
                        //          }
                        //          ++side;
                        //          edge = 0;
                        //        }
                        //      } else {
                        //        ++side;
                        //        edge = 0;
                        //      }
                        ///dxy delete end
                    } while (runPos != startPos);
                    
                    ///dxy delete
                    //      if (_boundaryTopoStr.empty()) {
                    //      ///
                    //          // update the starting corner
                    //          startPos.FlipE();
                    //          while (!startPos.IsBorder() && _nCorners > 0) {
                    //              startPos.FlipF();
                    //              startPos.FlipE();
                    //              if (!startPos.IsBorder()) {     // concave corner
                    //                  startPos.FlipE();
                    //                  while (!startPos.IsBorder()) {
                    //                      startPos.FlipF();
                    //                      startPos.FlipE();
                    //                  }
                    //                  startPos.FlipE();
                    //                  break;
                    //              } else {
                    //                  startPos.FlipV();
                    //                  startPos.FlipE();
                    //              }
                    //          }
                    //          startPos.FlipE();
                    //          // rotate to match the boundary contraints
                    //          ///dxy test
                    //          if (_nRotations > 0) {
                    //              std::cout << "_nRot = " << _nRotations << std::endl;
                    //          }
                    //          ///
                    //          for (n_corners_type i = 0; i < _nRotations; ++i) {
                    //              startPos.FlipV();
                    //              startPos.FlipE();
                    //              while (!startPos.IsBorder()) {
                    //                  startPos.FlipF();
                    //                  startPos.FlipE();
                    //                  if (!startPos.IsBorder()) {   // concave corner
                    //                      while (!startPos.IsBorder()) {
                    //                          startPos.FlipF();
                    //                          startPos.FlipE();
                    //                      }
                    //                      break;
                    //                  } else {
                    //                      startPos.FlipV();
                    //                      startPos.FlipE();
                    //                  }
                    //              }
                    //          }
                    //      }
                    ///dxy delete end
                    
                    // update boundary vertices flags
                    if (vcg::tri::HasPerVertexFlags(patchMesh)) {
                        runPos = startPos;
                        do {
                            runPos.V()->SetB();
                            runPos.FlipV();
                            runPos.FlipE();
                            while (!runPos.IsBorder()) {
                                runPos.FlipF();
                                runPos.FlipE();
                            }
                        } while (runPos != startPos);
                    }
                    
                    // initialize mark
                    vcg::tri::InitFaceIMark(patchMesh);
                    vcg::tri::InitVertexIMark(patchMesh);
                    vcg::tri::IMark(patchMesh) = 0;
                    
                    ///dxy test
//                    std::cout << "reconstructPatch end" << std::endl << std::endl;
                    ///
                    
                }
                
                // set methods:
                /**
                 * @brief reset deletes the previous convex patch's information.
                 */
                void reset() {
                    _topologyStr.clear();
                    _expandedTopologyStr.clear();
                    _boundaryIDsVec.clear();
                    _boundaryIDsEndsVec.clear();
                    _boundaryVarVec.clear();
                    _boundaryLenVec.clear();
                    _singularityStr.clear();
                    _nSingularities = 0;
                    _perimeter = 0;
                    _nCorners = 0;
                    _nFaces = 0;
                    _nVertices = 0;
                    _nRotations = 0;
                    _nOccurrences = 0;
                    _nConcaveCorners = 0;
                    if (_patchMesh != 0)
                        delete _patchMesh;
                    _patchMesh = 0;
                    if (_startPos != 0)
                        delete _startPos;
                    _startPos = 0;
                    _polychordMap.clear();
                    _flowDistance = std::numeric_limits<ScalarType>::max();
                    _quadsQuality = std::numeric_limits<ScalarType>::max();
                    ///dxy add
                    _boundaryTopoStr.clear();
                    ///
                }
                
                /**
                 * @brief setPatch stores all the information about the patch (sub)mesh bounded by the nCorners corners.
                 * @param mesh The vcg PolyMeshType mesh.
                 * @param corners A vector of "face-edge-vertex" positions identifying the corners (or a starting point, if nCorners==0).
                 * @param nCorners The number of corners: 0 and 1 corner = 1 side, n corners = n sides.
                 * @param nFaces The number of internal faces of the patch (sub)mesh.
                 * @param nVertices The total number of vertices (both internal and on boundaries).
                 * @param singularityVec A vector of pairs <valence_V,number_of_vertices_with_valence_V> (see _singularityStr).
                 */
#include <iostream>
                void setPatch(PolyMeshType &mesh,
                              const std::deque<PosType> &corners,
                              const n_corners_type nCorners,
                              const num_type nFaces,
                              const num_type nVertices,
                              const std::vector< std::pair<valence_type,num_type> > &singularityVec) {
                    // delete the previous convex patch
                    reset();
                    
                    if (corners.size() < 1 || (nCorners != (n_corners_type)corners.size() && nCorners == 0 && corners.size() > 1))
                        return;
                    
                    // clear flags
                    vcg::tri::UpdateFlags<PolyMeshType>::FaceClearV(mesh);
                    
                    // set number of corners, faces and vertices
                    _nCorners = nCorners;
                    _nFaces = nFaces;
                    _nVertices = nVertices;
                    if (_nCorners > 0)
                        _nConcaveCorners = PolychordSupport::CountConcaveCorners(mesh, corners);
                    else
                        _nConcaveCorners = 0;
                    
                    PosType runPos, startPos;
                    num_type length = 0;
                    n_corners_type b = 0;;
                    std::string str;
                    
                    // compute keys and boundaries
                    std::vector<num_type> boundaryLenVec(corners.size(), 0);
                    startPos = corners.front();
                    BuildTopologyStr(mesh, startPos, _topologyStr);
                    for (n_corners_type i = 0; i < (n_corners_type)corners.size(); ++i) {
                        length = 0;
                        runPos = corners[i];
                        do {
                            if (_nCorners == 0) {
                                BuildTopologyStr(mesh, runPos, str);
                                if (str < _topologyStr) {
                                    _topologyStr = str;
                                    startPos = runPos;
                                }
                            }
                            runPos.FlipV();
                            runPos.FlipE();
                            ++length;
                            while (!runPos.IsBorder()) {
                                runPos.FlipF();
                                runPos.FlipE();
                            }
                        } while (runPos != corners[(i + 1) % corners.size()]);
                        boundaryLenVec[i] = length;
                        _perimeter += length;
                        if (_nCorners > 0) {
                            BuildTopologyStr(mesh, corners[i], str);
                            if (str < _topologyStr) {
                                _topologyStr = str;
                                b = i;
                                startPos = corners[i];
                            }
                        }
                    }
                    
                    // build boundary Length vector
                    _boundaryLenVec.resize(boundaryLenVec.size());
                    for (size_t i = 0; i < boundaryLenVec.size(); ++i)
                        _boundaryLenVec[i] = boundaryLenVec[(i + b) % boundaryLenVec.size()];
                    
                    ///dxy add: build boundary topoStr
                    buildBoundaryTopoStr(mesh, startPos, _boundaryTopoStr);
                    ///
                    
                    // build the boundary IDs vector and ends
                    buildBoundaryIDsFromMesh(mesh, startPos);  ///this has been dxy modified
                    
                    // build the boundary IDs values
                    _boundaryVarVec.assign(std::accumulate(_boundaryLenVec.begin(), _boundaryLenVec.end(), 0) / 2, 1);
                    
                    // build the string representing the singular vertices count
                    buildSingularityStr(singularityVec);
                    
                    // count the total number of singularities
                    setNumberOfSingularities(singularityVec);
                    
                    // set the number of rotations w.r.t. the first corner
                    _nRotations = b;  // b
                    
                    // check string representations
                    if (_topologyStr.empty())
                        reset();
                }
                
                /**
                 * @brief setTopologyStr See _topologyStr.
                 * @param topologyStr
                 */
                void setTopologyStr(const std::string &topologyStr) {
                    _topologyStr = topologyStr;
                    _expandedTopologyStr.clear();
                }
                
                /**
                 * @brief setExpandedTopologyStr See _expandedTopologyStr.
                 * @param topologyStr
                 */
                void setExpandedTopologyStr(const std::string &expandedTopologyStr) {
                    _expandedTopologyStr = expandedTopologyStr;
                }
                
                /**
                 * @brief buildExpandedTopologyStr builds and sets the _expandedTopologyStr of this Patch.
                 * @note Call this after all ther setSomething() methods.
                 */
                void buildExpandedTopologyStr() {
                    // reconstruct the mesh
                    reconstructPatch();
                    // build the topology string key
                    BuildTopologyStr(*_patchMesh, *_startPos, _expandedTopologyStr);
                }
                
                /**
                 * @brief setBoundaryFromStr sets _boundaryIDsVec, _boundaryIDsEndsVec, _boundaryVarVec, _boundaryLenVec
                 * and _perimeter from the string representing _boundaryIDsVec.
                 * @note Ensure to set _nCorners before calling this method.
                 * @param boundaryIDsStr
                 */
                void setBoundaryFromStr(const std::string &boundaryIDsStr) {
                    buildBoundaryLenVecFromIDsStr(boundaryIDsStr);
                    _perimeter = std::accumulate(_boundaryLenVec.begin(), _boundaryLenVec.end(), 0);
                    _boundaryVarVec.assign(_perimeter / 2, 1);
                    buildBoundaryIDsAndEndsVecsFromIDsStr(boundaryIDsStr);
                    _boundaryVarVec.assign(_boundaryIDsEndsVec.size(), 1);
                    std::vector<num_type> bLen = _boundaryLenVec;
                    for (n_corners_type side = 0; side < bLen.size(); ++side)
                        bLen[side] = _boundaryLenVec[(side + _nRotations) % bLen.size()];
                    _boundaryLenVec = bLen;
                    _expandedTopologyStr.clear();
                }
                
                /**
                 * @brief setBoundaryVarVec See _boundaryVarVec.
                 * @param boundaryVarVec
                 */
                void setBoundaryVarVec(const std::vector<num_type> &boundaryVarVec) {
                    _boundaryVarVec = boundaryVarVec;
                    _expandedTopologyStr.clear();
                }
                
                /**
                 * @brief setBoundaryLenVec See _boundaryLenVec.
                 * @param boundaryLenVec
                 */
                void setBoundaryLenVec(const std::vector<num_type> &boundaryLenVec) {
                    _boundaryLenVec = boundaryLenVec;
                    _expandedTopologyStr.clear();
                }
                
                /**
                 * @brief setSingularityStr See _singularityStr.
                 * @param singularityStr
                 */
                void setSingularityStr(const std::string &singularityStr) {
                    _singularityStr = singularityStr;
                    _expandedTopologyStr.clear();
                }
                
                /**
                 * @brief setNumberOfSingularities See _nSingularities.
                 * @param nSingularities
                 */
                void setNumberOfSingularities(const num_type nSingularities) {
                    _nSingularities = nSingularities;
                    _expandedTopologyStr.clear();
                }
                
                /**
                 * @brief setNumberOfSingularities sums the second component of each singukarityVec element into _nSingularities.
                 * @param singularityVec is a vector of pairs <valence_type, number-of-occurrences>.
                 */
                void setNumberOfSingularities(const std::vector< std::pair<valence_type,num_type> > &singularityVec) {
                    _nSingularities = 0;
                    for (size_t i = 0; i < singularityVec.size(); ++i)
                        _nSingularities += singularityVec[i].second;
                    _expandedTopologyStr.clear();
                }
                
                /**
                 * @brief setPerimeter See _perimeter.
                 * @param perimeter
                 */
                void setPerimeter(const num_type perimeter) {
                    _perimeter = perimeter;
                    _expandedTopologyStr.clear();
                }
                
                /**
                 * @brief setNumberOfCorners See _nCorners.
                 * @param nCorners
                 */
                void setNumberOfCorners(const n_corners_type nCorners) {
                    _nCorners = nCorners;
                    _expandedTopologyStr.clear();
                }
                
                /**
                 * @brief setNumberOfFaces See _nFaces.
                 * @param nFaces
                 */
                void setNumberOfFaces(const num_type nFaces) {
                    _nFaces = nFaces;
                    _expandedTopologyStr.clear();
                }
                
                /**
                 * @brief setNumberOfVertices See _nVertices.
                 * @param nVertices
                 */
                void setNumberOfVertices(const num_type nVertices) {
                    _nVertices = nVertices;
                    _expandedTopologyStr.clear();
                }
                
                /**
                 * @brief setNumberOfRotations See _nRotations.
                 * @param nRotations
                 */
                void setNumberOfRotations(const n_corners_type nRotations) {
                    _nRotations = nRotations;
                    _expandedTopologyStr.clear();
                }
                
                /**
                 * @brief setNumberOfOccurrences See _nOccurrences.
                 * @param nOccurences
                 */
                void setNumberOfOccurrences(const count_type nOccurrences) {
                    _nOccurrences = nOccurrences;
                }
                
                /**
                 * @brief incrementNumberOfOccurrences increments _nOccurrences by one.
                 */
                void incrementNumberOfOccurrences() {
                    ++_nOccurrences;
                }
                
                /**
                 * @brief setNumberOfConcaveCorners See _nConcaveCorners.
                 * @param nConcaveCorners
                 */
                void setNumberOfConcaveCorners(const n_corners_type nConcaveCorners) {
                    _nConcaveCorners = nConcaveCorners;
                }
                
                /**
                 * @brief reconstructConvexPatch creates a vcg PolyMeshType mesh from this convex patch's information,
                 * also taking into account _boundaryVarStr values to split polychords and to have the full convex patch.
                 * The mesh and the starting Pos are stored internally.
                 */
                void reconstructPatch() {
                    if (_patchMesh == 0)
                        _patchMesh = new PolyMeshType;
                    if (_startPos == 0)
                        _startPos = new PosType;
                    reconstructPatch(*_patchMesh, *_startPos);
                }
                
                /**
                 * @brief setFlowDistance See _flowDistance.
                 * @param flowDistance
                 */
                void setFlowDistance(const ScalarType flowDistance) {
                    _flowDistance = flowDistance;
                }
                
                /**
                 * @brief setQuadsQuality See _quadsQuality.
                 * @param quadsQuality
                 */
                void setQuadsQuality(const ScalarType quadsQuality) {
                    _quadsQuality = quadsQuality;
                }
                
                // useful static methods:
                /**
                 * @brief BuildTopologyStr constructs a string representing the topology in a unique form, starting from a given vertex.
                 * @param mesh The mesh into which analyze the (marked) patch.
                 * @param startCornerPos The vertex-edge-face the unique topology starts from.
                 * @param out The resulting string.
                 */
                static void BuildTopologyStr(const PolyMeshType &mesh, const PosType &startCornerPos, std::string &out) {
                    // reset the output string
                    out.clear();
                    
                    // check parameters
                    if (startCornerPos.IsNull())
                        return;
                    vcg::tri::RequirePerFaceMark(mesh);
                    if (!vcg::tri::IsMarked(mesh, startCornerPos.F()))
                        return;
                    vcg::tri::RequirePerFaceFlags(mesh);
                    
                    std::queue<PosType> facesQueue, faceSequenceQueue;
                    std::vector<num_type> faces_id_Map(mesh.face.size());
                    size_t ind = 0;
                    PosType currentPos, tmpPos, newPos;
                    
                    // visit the starting face
                    startCornerPos.F()->SetV();
                    // initial ID
                    num_type currentID = 0;
                    // map face with ID
                    ind = vcg::tri::Index(mesh, startCornerPos.F());
                    faces_id_Map[ind] = currentID;
                    // push the starting pos on queue
                    facesQueue.push(startCornerPos);
                    faceSequenceQueue.push(startCornerPos);
                    // visit the faces
                    while (!facesQueue.empty()) {
                        // get the top face
                        currentPos = facesQueue.front();
                        facesQueue.pop();
                        
                        // run over the edges
                        tmpPos = currentPos;
                        do {
                            // if find a new unvisited face
                            if (!tmpPos.IsBorder() && vcg::tri::IsMarked(mesh, tmpPos.FFlip()) && !tmpPos.FFlip()->IsV()) {
                                // build a new pos
                                newPos = tmpPos;
                                newPos.FlipF();
                                newPos.FlipV();
                                // visit the new face
                                newPos.F()->SetV();
                                // increment the ID
                                if (currentID == std::numeric_limits<num_type>::max()) {
                                    std::cout << "Error. Patch too much big. " << sizeof(num_type) << " byte"
                                    << (sizeof(num_type) > 1 ? "s are" : " is")
                                    << " not enough to represent more than "
                                    << std::numeric_limits<num_type>::max()
                                    << " faces. Topology string not built." << std::endl;
                                    return;
                                }
                                ++currentID;
                                // map the new face with the new ID
                                ind = vcg::tri::Index(mesh, newPos.F());
                                faces_id_Map[ind] = currentID;
                                // push the new face on queue
                                facesQueue.push(newPos);
                                faceSequenceQueue.push(newPos);
                            }
                            // go to the next edge
                            tmpPos.FlipV();
                            tmpPos.FlipE();
                        } while (tmpPos != currentPos);
                    }
                    
                    // now list the faces and its adjacent ones in the given order
                    std::stringstream stream;
                    int quad_egde = 0;
                    // visit the faces
                    while (!faceSequenceQueue.empty()) {
                        // get the top face
                        currentPos = faceSequenceQueue.front();
                        faceSequenceQueue.pop();
                        
                        // reset the new face
                        currentPos.F()->ClearV();
                        
                        // run over the edges
                        tmpPos = currentPos;
                        do {
                            if (quad_egde != 0) {
                                // if find a new unvisited face
                                if (!tmpPos.IsBorder() && vcg::tri::IsMarked(mesh, tmpPos.FFlip())) {   // push the ID into the string
                                    stream << std::setw(2 * sizeof(num_type)) << std::setfill('0');
                                    stream << std::hex << std::noshowbase;
                                    ind = vcg::tri::Index(mesh, tmpPos.FFlip());
                                    stream << faces_id_Map[ind];
                                } else {
                                    // else push a special value representing a border into the string
                                    stream << std::string(2 * sizeof(num_type), BORDER);
                                }
                            }
                            quad_egde = (quad_egde + 1) % 4;
                            // go to the next edge
                            tmpPos.FlipV();
                            tmpPos.FlipE();
                        } while (tmpPos != currentPos);
                    }
                    
                    // finally, store the new string
                    out = stream.str();
                    
                }
                
                
                
            private:
                // private methods:
                /**
                 * @brief buildBoundaryIDsFromMesh stores the boundary IDs by walking on the border of the mesh.
                 * @param mesh The input mesh.
                 * @param startPos The starting position.
                 * @note _boundaryLenVec must be already built.
                 */
                void buildBoundaryIDsFromMesh(const PolyMeshType &mesh, const PosType &startPos) {
                    // check parameters
                    if (mesh.IsEmpty()) {
                        std::cout << "Mesh is empty. Boundary IDs not built. Skipped." << std::endl;
                        return;
                    }
                    if (startPos.IsNull()) {
                        std::cout << "Wrong start. Boundary IDs not built. Skipped." << std::endl;
                        return;
                    }
                    if (_boundaryLenVec.empty()) {
                        std::cout << "No boundary sizes. Boundary IDs not built. Skipped." << std::endl;
                        return;
                    }
                    vcg::tri::RequirePerFaceMark(mesh);
                    if (!vcg::tri::IsMarked(mesh, startPos.F())) {
                        std::cout << "Wrong start. Boundary IDs not built. Skipped." << std::endl;
                        return;
                    }
                    
                    // some variables
#if __cplusplus >= 201103L || _MSC_VER >= 1800
                    std::unordered_map<HalfEdge, var_label_type, HalfEdgeHash> boundary_map;
                    typename std::unordered_map<HalfEdge, var_label_type, HalfEdgeHash>::iterator bIt;
#else
                    std::map<HalfEdge, var_label_type> boundary_map;
                    typename std::map<HalfEdge, var_label_type>::iterator bIt;
#endif
                    PosType borderPos = startPos, runPos;
                    n_corners_type side = 0;
                    var_label_type edge = 0;
                    var_label_type currentPolychordID = 0;
                    var_label_type newPolychordID = 0;
                    
                    // resize the boundary IDs vector
                    _boundaryIDsVec.resize(_boundaryLenVec.size());
                    for (size_t i = 0; i < _boundaryIDsVec.size(); ++i)
                        _boundaryIDsVec[i].resize(_boundaryLenVec[i]);
                    // resize the boundary IDs ends
                    _boundaryIDsEndsVec.resize(std::accumulate(_boundaryLenVec.begin(), _boundaryLenVec.end(), 0) / 2);
                    
                    // run on the boundary
                    do {
                        bIt = boundary_map.find(HalfEdge(vcg::tri::Index(mesh, borderPos.F()), borderPos.E()));
                        if (bIt == boundary_map.end()) {
                            ///dxy test
                            //          std::cout << "border not in b_map" << std::endl;
                            ///
                            // assign a new polychord ID
                            if (currentPolychordID == std::numeric_limits<var_label_type>::max()) {
                                std::cout << "Error. Symbols " << VARLABEL << "[0-" << std::numeric_limits<var_label_type>::max()
                                << "] are not enough to represent all the boundary IDs. Skipped." << std::endl;
                                _boundaryIDsVec.clear();
                                return;
                            }
                            currentPolychordID = newPolychordID;
                            ++newPolychordID;
                            boundary_map.insert(std::make_pair(HalfEdge(vcg::tri::Index(mesh, borderPos.F()), borderPos.E()), currentPolychordID));
                            
                            // run towards other end of this polychord
                            runPos = borderPos;
                            runPos.FlipV();
                            do {
                                runPos.FlipE();
                                runPos.FlipV();
                                runPos.FlipE();
                                if (!runPos.IsBorder() && vcg::tri::IsMarked(mesh, runPos.FFlip()))
                                    runPos.FlipF();
                            } while (!runPos.IsBorder() && vcg::tri::IsMarked(mesh, runPos.FFlip()));
                            // tag the other end
                            boundary_map.insert(std::make_pair(HalfEdge(vcg::tri::Index(mesh, runPos.F()), runPos.E()), currentPolychordID));
                            
                            // set ID
                            ///dxy test
                            //          std::cout << "side " << side << " edge " << edge << std::endl;
                            ///
                            _boundaryIDsVec[side][edge] = currentPolychordID;
                            // set ID's end
                            _boundaryIDsEndsVec[currentPolychordID].first = side;
                        } else {
                            ///test
                            //          std::cout << "border already in b_map" << std::endl;
                            ///
                            ///dxy test
                            //          std::cout << "side " << side << " edge " << edge << std::endl;
                            ///
                            // set ID
                            ///dxy test
                            //          std::cout << "bIt->second = " << bIt->second << std::endl;
                            ///
                            _boundaryIDsVec[side][edge] = bIt->second;
                            // set ID's end
                            _boundaryIDsEndsVec[bIt->second].second = side;
                        }
                        
                        // go to the next edge
                        borderPos.FlipV();
                        borderPos.FlipE();
                        ///dxy add
                        if (edge == _boundaryLenVec[side] - 1) { //last edge on this side
                            while (!borderPos.IsBorder()) {
                                borderPos.FlipF();
                                borderPos.FlipE();
                            }
                            side += 1;
                            edge = 0;
                        }
                        else {
                            while (!borderPos.IsBorder()) {
                                borderPos.FlipF();
                                borderPos.FlipE();
                            }
                            edge += 1;
                        }
                        ///dxy add end
                        ///dxy delete
                        //      if (borderPos != startPos && !borderPos.IsBorder() && vcg::tri::IsMarked(mesh, borderPos.FFlip())) {
                        //        borderPos.FlipF();
                        //        borderPos.FlipE();
                        //        if (!borderPos.IsBorder() && vcg::tri::IsMarked(mesh, borderPos.FFlip())) {   // concave corner
                        //          while (!borderPos.IsBorder() && vcg::tri::IsMarked(mesh, borderPos.FFlip())) {
                        //            borderPos.FlipF();
                        //            borderPos.FlipE();
                        //          }
                        //          ++side;
                        //          edge = 0;
                        //        } else
                        //          ++edge;
                        //      } else {
                        //        ++side;
                        //        edge = 0;
                        //      }
                        ///dxy delete end
                    } while (borderPos != startPos);
                }
                
                /**
                 * @brief buildBoundaryLenVecFromIDsStr extracts a vector of boundary lengths from a boundary ID string.
                 * @note _nCorners must be already set.
                 * @param boundaryIDs A tring representing the boundary iDs.
                 */
                void buildBoundaryLenVecFromIDsStr(const std::string &boundaryIDs) {
                    _boundaryLenVec.clear();
                    if (boundaryIDs.empty())
                        return;
                    
                    n_corners_type side = 0;
                    _boundaryLenVec.resize((_nCorners > 0 ? _nCorners : 1), 0);
                    for (size_t i = 0; i < boundaryIDs.length(); i += 1/*VARLABEL*/ + 2 * sizeof(var_label_type)/*variable index*/) {
                        while (boundaryIDs[i] == SEPARATOR) {
                            ++side;
                            ++i;
                        }
                        ++_boundaryLenVec[side];
                    }
                }
                
                /**
                 * @brief buildBoundaryIDsAndEndsVecsFromIDsStr converts the boundaryIDsStr into a vector of vectors of variables,
                 * a vector for each side, and stores the head and tail of each variable. See _boundaryIDsVec and _boundaryIDsEndsVec.
                 * @note Ensure that !_boundaryLenVec.empty() when this method is called.
                 */
                void buildBoundaryIDsAndEndsVecsFromIDsStr(const std::string &boundaryIDsStr) {
                    _boundaryIDsVec.clear();
                    _boundaryIDsEndsVec.clear();
                    if (_boundaryLenVec.empty())
                        return;
                    
                    n_corners_type side = 0;
                    var_label_type edge = 0;
                    std::stringstream stream;
                    num_type num = 0;
                    
                    // get sizes and variables values
                    _boundaryIDsVec.resize(_boundaryLenVec.size());
                    for (size_t i = 0; i < _boundaryLenVec.size(); ++i)
                        _boundaryIDsVec[i].resize(_boundaryLenVec[i]);
                    _boundaryIDsEndsVec.resize(_perimeter / 2, BoundaryIDHeadTailType(std::numeric_limits<n_corners_type>::max(),
                                                                                      std::numeric_limits<n_corners_type>::max()));
                    // boundary IDs string to vector
                    for (size_t i = 0; i < boundaryIDsStr.length(); i += 1/*VARLABEL*/ + 2 * sizeof(var_label_type)/*variable index*/, ++edge) {
                        while (boundaryIDsStr[i] == SEPARATOR) {
                            ++side;
                            ++i;
                            edge = 0;
                        }
                        stream.str(boundaryIDsStr.substr(i + 1, 2 * sizeof(var_label_type)) + ' ');
                        stream >> std::hex >> num;
                        _boundaryIDsVec[side][edge] = static_cast<var_label_type>(num);
                        if (_boundaryIDsEndsVec[num].first == std::numeric_limits<n_corners_type>::max())
                            _boundaryIDsEndsVec[num].first = side;
                        else {
                            assert(_boundaryIDsEndsVec[num].second == std::numeric_limits<n_corners_type>::max());
                            _boundaryIDsEndsVec[num].second = static_cast<n_corners_type>(side);
                        }
                    }
                }
                
                /**
                 * @brief buildSingularityStr constructs a string representing the singularities (see _singularityStr).
                 * @param singularityVec A vector of pairs <valence_V,number_of_vertices_with_valence_V>.
                 */
                void buildSingularityStr(const std::vector< std::pair<valence_type,num_type> > &singularityVec) {
                    std::stringstream stream;
                    for (size_t i = 0; i < singularityVec.size(); ++i) {
                        stream << VALENCELABEL;
                        stream << std::setw(2 * sizeof(valence_type)) << std::setfill('0') << std::hex << std::noshowbase;
                        stream << (num_type)singularityVec[i].first;
                        stream << VALENCESEP;
                        stream << std::setw(2 * sizeof(num_type)) << std::setfill('0') << std::hex << std::noshowbase;
                        stream << singularityVec[i].second;
                        if (i < singularityVec.size() - 1)
                            stream << SEPARATOR;
                    }
                    _singularityStr = stream.str();
                }
                
                
                ///dxy add
                static void buildBoundaryTopoStr(const PolyMeshType &mesh, const PosType startPos,
                                                 std::string &result)
                {
                    //check input
                    if (!startPos.IsBorder()) {
                        return;
                    }
                    
                    // run around the border (clockwise)
                    std::vector<std::string> topoStrVec;
                    std::string tStr;
                    PosType runPos = startPos;
                    do {
                        BuildTopologyStr(mesh, runPos, tStr);
                        topoStrVec.push_back(tStr);
                        runPos.FlipV();
                        runPos.FlipE();
                        while (!runPos.IsBorder()) {
                            runPos.FlipF();
                            runPos.FlipE();
                        }
                    } while (runPos != startPos);
                    
                    // topoStrVec to topoIdxVec
                    //1. sort topoStrVec
                    std::map<std::string, int> topoStr_count;
                    for (int i=0; i<topoStrVec.size(); ++i) {
                        topoStr_count[topoStrVec[i]]++;
                    }
                    std::map<std::string, int> topoStr_order;
                    int order = 0;
                    for (std::map<std::string, int>::iterator it=topoStr_count.begin(); it!=topoStr_count.end(); ++it) {
                        topoStr_order[it->first] = order;
                        order++;
                    }
                    //2. convert to topoIdxVec
                    std::vector<int> topoIdxVec;
                    for (int i=0; i<topoStrVec.size(); ++i) {
                        topoIdxVec.push_back(topoStr_order[topoStrVec[i]]);
                    }
                    
                    
                    // topoIdxVec to string
                    std::string idxStr;
                    for (int i=0; i<topoIdxVec.size(); ++i) {
                        idxStr.push_back('a' + topoIdxVec[i]);
                    }
                    
                    // return
                    result = idxStr;
                }
                
                ///dxy add end
                
                
                
                // data members (private):
                std::string                 _topologyStr;         ///< String representing the connectivity in a unique form.
                std::string                 _expandedTopologyStr; ///< String representing the connectivity of the expanded patch in a unique form.
                BoundaryIDsVecType          _boundaryIDsVec;      ///< Vector of vectors of boundary's edges IDs, one for each side.
                BoundaryIDsEndsVecType      _boundaryIDsEndsVec;  ///< Vector of pairs, one for each boundary's edges ID, indicating head and tail edges.
                std::vector<num_type>       _boundaryVarVec;      ///< Vector of boundary's ID's values.
                std::vector<num_type>       _boundaryLenVec;      ///< Vector of boundary' sides' length in counterclockwise order.
                std::string                 _singularityStr;      ///< String of pairs <SS:ss> where SS is the (hexadecimal) valence of the singularity and ss is the (hexadecimal) number of vertices with valence SS (with separator).
                num_type                    _nSingularities;      ///< Total number of singularities.
                num_type                    _perimeter;           ///< Total length of the boundary.
                n_corners_type              _nCorners;            ///< Number of corners: 0 corners and 1 corner = 1 side, n corners = n sides.
                num_type                    _nFaces;              ///< Number of quad faces.
                num_type                    _nVertices;           ///< Number of vertices.
                n_corners_type              _nRotations;          ///< Number of side rotations w.r.t. the unique form of the connectivity.
                count_type                  _nOccurrences;        ///< Number of occurrences of this patch found during the learning phase.
                n_corners_type              _nConcaveCorners;     ///< NUmber of concave corners.
                
                PolyMeshType               *_patchMesh;           ///< The mesh of this patch.
                PosType                    *_startPos;            ///< The starting corner.
                PolychordMap                _polychordMap;        ///< The map of this patch's polychords with flows.
                ScalarType                  _flowDistance;        ///< A value of distance of this patch wrt to the flows.
                ScalarType                  _quadsQuality;        ///< A value summing up the overall quad quality of the patch.
                
                ///dxy add : variables
                std::string                 _boundaryTopoStr;
                ///dxy add end
            };
            
            
            
            /**
             * @brief The PatchCompare class provides a way of changing the default Patch sort criteria.
             */
            template < typename PolyMeshType >
            class PatchCompare {
            public:
                typedef pl::Patch<PolyMeshType>     Patch;          ///< The type of objects to compare.
                typedef typename Patch::ScalarType  ScalarType;
                
                /**
                 * @brief The CompareCriteria enum defines some criteria of comparison.
                 */
                enum CompareCriteria {
                    PATCH_DEFAULT         = 0x00,                     ///< Compare using the criterion of Patch.
                    SINGULARITYTYPES      = 0x01,                     ///< Compare by singularity types and numbers (string).
                    SINGULARITYNUM        = 0x02,                     ///< Compare by number of singularities.
                    FLOW_DISTANCE         = 0x04,                     ///< Compare by flow distance.
                    QUADS_QUALITY         = 0x08                      ///< Compare by quality of quads.
                };
                
                /**
                 * @brief PatchCompare Default constructor.
                 */
                PatchCompare() : _lambda(-1) {
                    _criteria.push_back(FLOW_DISTANCE);
                    _criteria.push_back(SINGULARITYNUM);
                    _criteria.push_back(SINGULARITYTYPES);
                    _criteria.push_back(QUADS_QUALITY);
                }
                
                /**
                 * @brief PatchCompare Constructor with copy initializer.
                 * @param cs The sequence of criteria di copy.
                 */
                PatchCompare(const std::list<CompareCriteria> &cs) : _criteria(cs), _lambda(-1) { }
                
#if __cplusplus >= 201103L || _MSC_VER >= 1800
                /**
                 * @brief PatchCompare Constructor with move initializer.
                 * @param cs The sequence of criteria to move.
                 */
                PatchCompare(std::list<CompareCriteria>&& cs) : _criteria(cs), _lambda(-1) { }
#endif
                
                /**
                 * @brief PatchCompare Copy constructor.
                 * @param pc The PatchCompare to copy.
                 */
                PatchCompare(const PatchCompare &pc) : _criteria(pc._criteria), _lambda(pc._lambda) { }
                
#if __cplusplus >= 201103L || _MSC_VER >= 1800
                /**
                 * @brief PatchCompare Move constructor.
                 * @param pc The PatchCompare to move.
                 */
                PatchCompare(PatchCompare&& pc) : _criteria(std::move(pc._criteria)), _lambda(pc._lambda) {
                    pc._criteria.clear();
                    pc._criteria.push_back(FLOW_DISTANCE);
                    pc._criteria.push_back(SINGULARITYNUM);
                    pc._criteria.push_back(SINGULARITYTYPES);
                    pc._criteria.push_back(QUADS_QUALITY);
                    pc._lambda = ScalarType(-1);
                }
#endif
                
                /**
                 * @brief operator = Copy assignment.
                 * @param pc The PatchCompare to copy.
                 * @return A reference to this.
                 */
                PatchCompare & operator = (const PatchCompare &pc) {
                    if (&pc != this) {
                        _criteria = pc._criteria;
                        _lambda = pc._lambda;
                    }
                    return *this;
                }
                
#if __cplusplus >= 201103L || _MSC_VER >= 1800
                /**
                 * @brief operator = Move assignment.
                 * @param pc The PatchCompare to move.
                 * @return A reference to this.
                 */
                PatchCompare & operator = (PatchCompare&& pc) {
                    if (&pc != this) {
                        _criteria = std::move(pc._criteria);
                        _lambda = pc._lambda;
                        pc._criteria.clear();
                        pc._criteria.push_back(FLOW_DISTANCE);
                        pc._criteria.push_back(SINGULARITYNUM);
                        pc._criteria.push_back(SINGULARITYTYPES);
                        pc._criteria.push_back(QUADS_QUALITY);
                        pc._lambda = ScalarType(-1);
                    }
                    return *this;
                }
#endif
                
                /**
                 * @brief setCriteria sets a new sequence of criteria.
                 * @param criteria The list of CompareCriteria to copy.
                 */
                void setCriteria(const std::list<CompareCriteria> &criteria) {
                    _criteria = criteria;
                }
                
#if __cplusplus >= 201103L || _MSC_VER >= 1800
                /**
                 * @brief setCriteria sets a new sequence of criteria.
                 * @param criteria The list of CompareCriteria to move.
                 */
                void setCriteria(std::list<CompareCriteria>&& criteria) {
                    _criteria = criteria;
                }
#endif
                
                /**
                 * @brief getCriteria gives the current sequence of criteria.
                 * @return The current criteria.
                 */
                const std::list<CompareCriteria> & getCriteria() const {
                    return _criteria;
                }
                
                /**
                 * @brief setLambda sets _lambda.
                 * @param lambda
                 */
                void setLambda(const ScalarType lambda) {
                    _lambda = lambda;
                }
                
                /**
                 * @brief getLambda See _lambda.
                 * @return
                 */
                ScalarType getLambda() const {
                    return _lambda;
                }
                
                /**
                 * @brief operator () compares two Patches, using the current sequence of comparison criteria.
                 * @note If the current sequence of criteria is empty, the PATCH_DEFAULT is still used at the end.
                 * @param l The left Patch.
                 * @param r The right Patch.
                 * @return true if l < r, false otherwise.
                 */
                bool operator () (const Patch &l, const Patch &r) const {
                    std::vector< std::pair<valence_type,num_type> > lSings, rSings;
                    for (typename std::list<CompareCriteria>::const_iterator cIt = _criteria.begin(); cIt != _criteria.end(); ++cIt) {
                        switch (*cIt) {
                            case PATCH_DEFAULT:
                                if (l != r)
                                    return l < r;
                                break;
                            case SINGULARITYTYPES:
                                l.getSingularityVec(lSings);
                                r.getSingularityVec(rSings);
                                if (lSings != rSings)
                                    return lSings < rSings;
                                break;
                            case SINGULARITYNUM:
                                if (l.getNumberOfSingularities() != r.getNumberOfSingularities())
                                    return l.getNumberOfSingularities() < r.getNumberOfSingularities();
                                break;
                            case FLOW_DISTANCE:
                                if (l.getFlowDistance() != r.getFlowDistance())
                                    return l.getFlowDistance() < r.getFlowDistance();
                                break;
                            case QUADS_QUALITY:
                                if (l.getQuadsQuality() != r.getQuadsQuality())
                                    return l.getQuadsQuality() < r.getQuadsQuality();
                                break;
                            default:
                                break;
                        }
                    }
                    return l < r;
                }
                
            private:
                std::list<CompareCriteria>          _criteria;      ///< A sequence of criteria do use in the comparison.
                ScalarType                          _lambda;        ///< A weight for blending two different measures.
            };
            
        } // end namespace pl
        
    } // end namespace tri
} // end namespace vcg

#endif // PATCH_H
