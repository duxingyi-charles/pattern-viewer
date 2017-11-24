#ifndef MESH_H
#define MESH_H

#include <vector>
#include <vcg/complex/complex.h>

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

class TVertex;
class TFace;

struct TUsedTypes : public vcg::UsedTypes<
        vcg::Use<TVertex>::AsVertexType,
        vcg::Use<TFace>::AsFaceType
        > {};

class TVertex : public vcg::Vertex<
        TUsedTypes,
        vcg::vertex::Coord3d,
        vcg::vertex::Normal3d,
        vcg::vertex::TexCoord2d,
        vcg::vertex::BitFlags,
        vcg::vertex::Mark,
        vcg::vertex::Color4b
        > {};

class TFace : public vcg::Face<
        TUsedTypes,
        vcg::face::Normal3d,
        vcg::face::BitFlags,
        vcg::face::VertexRef,
        vcg::face::FFAdj,
        vcg::face::Color4b
        > {};

class TMesh : public vcg::tri::TriMesh<
        std::vector<TVertex>,
        std::vector<TFace>
        > {};



class PolyVertex;
class PolyFace;

struct PolyUsedTypes : public vcg::UsedTypes<
        vcg::Use<PolyVertex>::AsVertexType,
        vcg::Use<PolyFace>::AsFaceType
        > {};

class PolyVertex : public vcg::Vertex<
        PolyUsedTypes,
        vcg::vertex::Coord3d,
        vcg::vertex::Normal3d,
        vcg::vertex::BitFlags,
        vcg::vertex::Mark,
        vcg::vertex::Color4b
        > {};

class PolyFace : public vcg::Face<
        PolyUsedTypes,
        vcg::face::PolyInfo,
        vcg::face::Normal3d,
        vcg::face::BitFlags,
        vcg::face::Mark,
        vcg::face::PFVAdj,
        vcg::face::PFFAdj,
        vcg::face::Qualityd,
        vcg::face::Color4b
        > {};

class PolyMesh : public vcg::tri::TriMesh<
        std::vector<PolyVertex>,
        std::vector<PolyFace>
        > {};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // MESH_H
