#ifndef MY_POLYGON_POLYCHORD_COLLAPSE_H
#define MY_POLYGON_POLYCHORD_COLLAPSE_H

#include <vcg/complex/algorithms/polygon_polychord_collapse.h>

/// for test only
#include <iostream>
#include <wrap/io_trimesh/export.h>
#include <string>

///

namespace vcg {
namespace tri {
  
    template < typename PolyMeshType >
    class myPolychordCollapse : public PolychordCollapse <PolyMeshType> {
    public:

        // enum PC_ResultCode
        using typename PolychordCollapse <PolyMeshType> :: PC_ResultCode;
        using          PolychordCollapse <PolyMeshType> :: PC_SUCCESS;
        using          PolychordCollapse <PolyMeshType> :: PC_NOTMANIF;
        using          PolychordCollapse <PolyMeshType> :: PC_NOTQUAD;
        using          PolychordCollapse <PolyMeshType> :: PC_NOLINKCOND;
        using          PolychordCollapse <PolyMeshType> :: PC_SINGSIDEA;
        using          PolychordCollapse <PolyMeshType> :: PC_SINGSIDEB;
        using          PolychordCollapse <PolyMeshType> :: PC_SINGBOTH;
        using          PolychordCollapse <PolyMeshType> :: PC_SELFINTERSECT;
        using          PolychordCollapse <PolyMeshType> :: PC_NOMOREMANIF;
        using          PolychordCollapse <PolyMeshType> :: PC_VOID;
        using          PolychordCollapse <PolyMeshType> :: PC_OTHER;
        //
        using typename PolychordCollapse <PolyMeshType> :: CoordType;
        using typename PolychordCollapse <PolyMeshType> :: FaceType;
        using typename PolychordCollapse <PolyMeshType> :: VertexPointer;
        using typename PolychordCollapse <PolyMeshType> :: FacePointer;
        //
        using typename PolychordCollapse <PolyMeshType> :: PC_Chords;
        using typename PolychordCollapse <PolyMeshType> :: LinkConditions;
        //methods
        using PolychordCollapse <PolyMeshType> :: VisitPolychord;
        using PolychordCollapse <PolyMeshType> :: IsVertexAdjacentToAnyNonManifoldEdge;
        using PolychordCollapse <PolyMeshType> :: WillPolychordBeManifold;
        using PolychordCollapse <PolyMeshType> :: IsPolychordSelfIntersecting;
        using PolychordCollapse <PolyMeshType> :: MarkPolychords;
        
        
        
        static void printResultCode(PC_ResultCode c)
        {
            std::string cname;
            switch (c) {
                case PC_SUCCESS:
                    cname = "PC_SUCCESS";
                    break;
                case PC_NOTMANIF:
                    cname = "PC_NOTMANIF";
                    break;
                case PC_NOTQUAD:
                    cname = "PC_NOTQUAD";
                    break;
                case PC_NOLINKCOND:
                    cname = "PC_NOLINKCOND";
                    break;
                case PC_SINGSIDEA:
                    cname = "PC_SINGSIDEA";
                    break;
                case PC_SINGSIDEB:
                    cname = "PC_SINGSIDEB";
                    break;
                case PC_SINGBOTH:
                    cname = "PC_SINGBOTH";
                    break;
                case PC_NOMOREMANIF:
                    cname = "PC_NOMOREMANIF";
                    break;
                case PC_VOID:
                    cname = "PC_VOID/PC_OTHER"; //PC_VOID = PC_OTHER = 0x100
                    break;
                default:
                    cname = "ERR: no matching ResultCode";
                    break;
            }
            std::cout << cname << std::endl;
            
        }
        

        
        static void say()
        {
            std::cout << "Hi, I am a myPolychordCollapse object" << std::endl;
        }
        
        
        /// dxy new type
        enum V_Type {
            V_Border    = 0,
            V_Corner    = 1,
            V_Inner     = 2,
            V_Abnormal  = 3
        };
        
        enum V_ValenceType {
            V_Regular   = 0,
            V_Singular  = 1
        };
        
        typedef std::pair<V_Type, V_ValenceType> V_Property;
        
        static void printVProperty(V_Property v)
        {
            std::cout << "Type: ";
            switch (v.first) {
                case V_Border:
                    std::cout << "Border";
                    break;
                case V_Corner:
                    std::cout << "Corner";
                    break;
                case V_Inner:
                    std::cout << "Inner";
                    break;
                default:
                    std::cout << "Abnormal";
                    break;
            }
            std::cout << " Valence: ";
            switch (v.second) {
                case V_Regular:
                    std::cout << "Regular";
                    break;
                default:
                    std::cout << "Singular";
                    break;
            }
            std::cout << std::endl;
        }
        
        /// dxy new functions:
        // find property of vertex binding to pos
        static V_Property getProperty(vcg::face::Pos<FaceType> &pos)
        {
            V_Property vp;
            int valence = pos.NumberOfIncidentVertices();
            if (valence == 1) {
                vp.first = V_Abnormal;
                vp.second = V_Singular;
                return vp;
            }
            vcg::face::JumpingPos<FaceType> jmpPos;
            jmpPos.Set(pos.F(), pos.E(), pos.V());
            if (jmpPos.FindBorder()) {
                switch (valence) {
                    case 2:
                        vp.first = V_Corner;
                        vp.second = V_Regular;
                        break;
                    case 3:
                        vp.first = V_Border;
                        vp.second = V_Regular;
                        break;
                    default:  // >=4
                        vp.first = V_Corner;
                        vp.second = V_Singular;
                        break;
                }
            }
            else { //inner vertex
                vp.first = V_Inner;
                vp.second = V_Regular;
                if (valence != 4) {
                    vp.second = V_Singular;
                }
            }
            return vp;
        }
        
        ///dxy add: when cornerBit is enabled
        static V_Property getProperty(vcg::face::Pos<FaceType> &pos, int cornerBit)
        {
            //            std::cout << "get start" << std::endl;
            V_Property vp;
            int valence = pos.NumberOfIncidentVertices();
            if (valence == 1) {
                vp.first = V_Abnormal;
                vp.second = V_Singular;
                //                std::cout << "get end" << std::endl;
                return vp;
            }
            vcg::face::JumpingPos<FaceType> jmpPos;
            jmpPos.Set(pos.F(), pos.E(), pos.V());
            if (jmpPos.FindBorder()) {
                if (pos.V()->IsUserBit(cornerBit)) {
                    vp.first = V_Corner;
                    vp.second = V_Regular;
                    valence += 2;
                    if (valence != 4) {
                        vp.second = V_Singular;
                    }
                }
                else {
                    vp.first = V_Border;
                    vp.second = V_Regular;
                    valence += 1;
                    if (valence != 4) {
                        vp.second = V_Singular;
                    }
                }
            }
            else { //inner vertex
                vp.first = V_Inner;
                vp.second = V_Regular;
                if (valence != 4) {
                    vp.second = V_Singular;
                }
            }
            //            std::cout << "get end" << std::endl;
            return vp;
        }

        
        
        
        /** dxy @override
         * @brief CheckPolychordFindStartPosition checks if it's a collapsable polychord.
         * @param pos Input The starting position.
         * @param startPos Output the new starting position (in case of borders).
         * @param checkSing true if singularities on both sides are not allowed.
         * @return PC_SUCCESS if it's a collapsable polychord, otherwise the code for the cause (startPos is on it).
         */
        static PC_ResultCode CheckPolychordFindStartPosition (const vcg::face::Pos<FaceType> &pos,
                                                              vcg::face::Pos<FaceType> &startPos,
                                                              const bool checkSing = true) {
            ///test
//            std::cout << "---------------------------------------------" << std::endl;
//            std::cout << "using my CheckPolychordFindStartPosition(...)" << std::endl;
//            std::cout << "---------------------------------------------" << std::endl;
            int counter = 0;
            ///
            assert(!pos.IsNull());
//            int valence = 0;
            bool singSideA = false, singSideB = false;
            bool borderSideA = false, borderSideB = false;
            bool polyBorderFound = false;
            vcg::face::JumpingPos<FaceType> jmpPos;
            
            startPos = pos;
            do {
                ///test
                counter++;
//                std::cout << std::endl << "Cross No." << counter << std::endl;
                ///
                
                // check if it is a quad
                if (startPos.F()->VN() != 4)
                    return PC_NOTQUAD;
                // check manifoldness
                if (IsVertexAdjacentToAnyNonManifoldEdge(startPos))
                    return PC_NOTMANIF;
                startPos.FlipV();
                if (IsVertexAdjacentToAnyNonManifoldEdge(startPos))
                    return PC_NOTMANIF;
                startPos.FlipV();
                
                // check if side A is on border
                startPos.FlipE();
                if (startPos.IsBorder())
                    borderSideA = true;
                startPos.FlipE();
                // check if side B is on border
                startPos.FlipV();
                startPos.FlipE();
                if (startPos.IsBorder())
                    borderSideB = true;
                startPos.FlipE();
                startPos.FlipV();
                
                /// dxy add
                if (checkSing) {
//                    vcg::face::Pos<FaceType> endA = startPos;
//                    startPos.FlipV();
//                    vcg::face::Pos<FaceType> endB = startPos;
//                    startPos.FlipV();
                    
                    V_Property pA = getProperty(startPos);
//                    std::cout << "A: ";
//                    printVProperty(pA);
                    int valenceA = startPos.NumberOfIncidentVertices();
//                    std::cout << "valence A: " << valenceA << std::endl;
                    startPos.FlipV();
                    V_Property pB = getProperty(startPos);
//                    std::cout << "B: ";
//                    printVProperty(pB);
                    int valenceB = startPos.NumberOfIncidentVertices();
//                    std::cout << "valence B: " << valenceB << std::endl;
                    startPos.FlipV();
                    
                    //1. check if both ends of cross are corners
                    if (pA.first == V_Corner && pB.first == V_Corner) {
//                        std::cout << "Bad cross: two corner ends" << std::endl;
                        singSideA = true;
                        singSideB = true;
                    }
//                    //2. check if one end is sigCorner, the other is singularity
//                    if (pA.first == V_Corner && pA.second == V_Singular && pB.second == V_Singular) {
//                        std::cout << "Bad cross: A is sigCorner and B is singular" << std::endl;
//                    }
//                    if (pB.first == V_Corner && pB.second == V_Singular && pA.second == V_Singular) {
//                        std::cout << "Bad cross: A is singular and B is sigCorner" << std::endl;
//                    }
                    //3. check if both ends of cross are singularities
                    if (pA.second == V_Singular && pB.second == V_Singular) {
//                        std::cout << "Bad cross: two singular ends" << std::endl;
                        singSideA = true;
                        singSideB = true;
                    }
                    //4. check if 2-inner end exist
//                    if (pA.first == V_Inner && valenceA == 2) {
//                        std::cout << "Bad cross: A is inner 2-singularity" << std::endl;
//                        singSideA = true;
//                        singSideB = true;
//                    }
//                    if (pB.first == V_Inner && valenceB == 2) {
//                        std::cout << "Bad cross: B is inner 2-singularity" << std::endl;
//                        singSideA = true;
//                        singSideB = true;
//                    }
                    //4. check if 2-singularity end exist
                    if (pA.first != V_Corner && valenceA == 2) {
//                        std::cout << "Bad cross: A is 2-singularity" << std::endl;
                        singSideA = true;
                        singSideB = true;
                    }
                    if (pB.first != V_Corner && valenceB == 2) {
//                        std::cout << "Bad cross: B is 2-singularity" << std::endl;
                        singSideA = true;
                        singSideB = true;
                    }

                    
                    
                    if (pA.second == V_Singular) {
                        singSideA = true;
                    }
                    if (pB.second == V_Singular) {
                        singSideB = true;
                    }
                    
                    //
                }
                
                // if the first border has been reached, go on the other direction to find the other border
                if (startPos != pos && startPos.IsBorder() && !polyBorderFound) {
                    startPos = pos;
                    startPos.FlipF();
                    polyBorderFound = true;
                }
                
                // if the other border has been reached, return
                if (polyBorderFound && startPos.IsBorder())
                    break;
                
                // go to the next edge
                startPos.FlipE();
                startPos.FlipV();
                startPos.FlipE();
                // check manifoldness
                if (IsVertexAdjacentToAnyNonManifoldEdge(startPos))
                    return PC_NOTMANIF;
                startPos.FlipV();
                if (IsVertexAdjacentToAnyNonManifoldEdge(startPos))
                    return PC_NOTMANIF;
                startPos.FlipV();
                // go to the next face
                startPos.FlipF();
            } while (startPos != pos);
            
            
            // polychord with singularities on both sides can not collapse
            /// dxy delete
            if ((singSideA && singSideB) //||
                //(singSideA && borderSideB) ||
                //(singSideB && borderSideA))
                )
                return PC_SINGBOTH;
            ///
            
            // polychords that are rings and have borders on both sides can not collapse
            if (!polyBorderFound && borderSideA && borderSideB)
                return PC_SINGBOTH;
            
            // if there are singularities or borders on the side A, remember to keep coordinates on it
            if (singSideA || borderSideA)
                return PC_SINGSIDEA;
            // if there are singularities or borders on the side B, remember to keep coordinates on it
            if (singSideB || borderSideB)
                return PC_SINGSIDEB;
            
            return PC_SUCCESS;
        }
        
        
        ///dxy add: when User cornerBit is enabled
        static PC_ResultCode CheckPolychordFindStartPosition (const vcg::face::Pos<FaceType> &pos,
                                                              vcg::face::Pos<FaceType> &startPos,
                                                              int cornerBit,
                                                              const bool checkSing = true) {
            ///test
//            std::cout << "---------------------------------------------" << std::endl;
//            std::cout << "using my CheckPolychordFindStartPosition(...)" << std::endl;
//            std::cout << "---------------------------------------------" << std::endl;
//            int counter = 0;
            ///
            assert(!pos.IsNull());
//            int valence = 0;
            bool singSideA = false, singSideB = false;
            bool borderSideA = false, borderSideB = false;
            bool polyBorderFound = false;
            vcg::face::JumpingPos<FaceType> jmpPos;
            
            startPos = pos;
            ///dxy add
            int sideA_repeat = 0;
            int sideB_repeat = 0;
            ///
            do {
                ///test
//                counter++;
//                std::cout << std::endl << "Cross No." << counter << std::endl;
                ///
                
                // check if it is a quad
                if (startPos.F()->VN() != 4)
                    return PC_NOTQUAD;
                // check manifoldness
                if (IsVertexAdjacentToAnyNonManifoldEdge(startPos))
                    return PC_NOTMANIF;
                startPos.FlipV();
                if (IsVertexAdjacentToAnyNonManifoldEdge(startPos))
                    return PC_NOTMANIF;
                startPos.FlipV();
                
                // check if side A is on border
                startPos.FlipE();
                if (startPos.IsBorder())
                    borderSideA = true;
                startPos.FlipE();
                // check if side B is on border
                startPos.FlipV();
                startPos.FlipE();
                if (startPos.IsBorder())
                    borderSideB = true;
                startPos.FlipE();
                startPos.FlipV();
                
                /// dxy add
                if (checkSing) {
                    V_Property pA = getProperty(startPos, cornerBit);
//                    std::cout << "A: ";
//                    printVProperty(pA);
                    int valenceA = startPos.NumberOfIncidentVertices();
//                    std::cout << "valence A: " << valenceA << std::endl;
                    
                    startPos.FlipV();
                    V_Property pB = getProperty(startPos, cornerBit);
//                    std::cout << "B: ";
//                    printVProperty(pB);
                    int valenceB = startPos.NumberOfIncidentVertices();
//                    std::cout << "valence B: " << valenceB << std::endl;
                    startPos.FlipV();
                    
                    //1. check if both ends of cross are corners
                    if (pA.first == V_Corner && pB.first == V_Corner) {
//                        std::cout << "Bad cross: two corner ends" << std::endl;
                        singSideA = true;
                        singSideB = true;
                    }
                    //2. check if both ends of cross are singularities
                    if (pA.second == V_Singular && pB.second == V_Singular) {
//                        std::cout << "Bad cross: two singular ends" << std::endl;
                        singSideA = true;
                        singSideB = true;
                    }
                    //3. check if 2-inner end exist
                    if (pA.first == V_Inner && valenceA == 2) {
//                        std::cout << "Bad cross: A is inner 2-singularity" << std::endl;
                        singSideA = true;
                        singSideB = true;
                    }
                    if (pB.first == V_Inner && valenceB == 2) {
//                        std::cout << "Bad cross: B is inner 2-singularity" << std::endl;
                        singSideA = true;
                        singSideB = true;
                    }
                    //4. check if 2-singularity end exist
                    if (pA.first != V_Corner && valenceA == 2) {
//                        std::cout << "Bad cross: A is 2-singularity" << std::endl;
                        singSideA = true;
                        singSideB = true;
                    }
                    if (pB.first != V_Corner && valenceB == 2) {
//                        std::cout << "Bad cross: B is 2-singularity" << std::endl;
                        singSideA = true;
                        singSideB = true;
                    }
                    //5. check if both ends are border but cross is not on border
                    if (pA.first != V_Inner && pB.first != V_Inner && !startPos.IsBorder()) {
                        //                        std::cout << "Bad cross: A, B are on border but cross is not border" << std::endl;
                        singSideA = true;
                        singSideB = true;
                    }
                    //6. check side edges
                    vcg::face::Pos<FaceType> lastPosA = startPos;
                    lastPosA.FlipF();
                    lastPosA.FlipE();
                    lastPosA.FlipF();
                    vcg::face::Pos<FaceType> sidePosA = startPos;
                    sidePosA.FlipE();
                    sidePosA.FlipF();
                    if (sidePosA.E() != lastPosA.E() && sidePosA.F() == lastPosA.F()) {
                        sideA_repeat += 1;
                    }
                    else {
                        sideA_repeat = 1;
                    }
                    startPos.FlipV();
                    vcg::face::Pos<FaceType> lastPosB = startPos;
                    lastPosB.FlipF();
                    lastPosB.FlipE();
                    lastPosB.FlipF();
                    vcg::face::Pos<FaceType> sidePosB = startPos;
                    sidePosB.FlipE();
                    sidePosB.FlipF();
                    startPos.FlipV();
                    if (sidePosB.E() != lastPosB.E() && sidePosB.F() == lastPosB.F()) {
                        sideB_repeat += 1;
                    }
                    else {
                        sideB_repeat = 1;
                    }
                    if (sideA_repeat >= 3 || sideB_repeat >= 3) {
//                        std::cout << "Bad Side: 3 consecutive side edges come from a common face" << std::endl;
                        singSideA = true;
                        singSideB = true;
                    }
                    //7. mark singularity side
                    if (pA.second == V_Singular) {
                        singSideA = true;
                    }
                    if (pB.second == V_Singular) {
                        singSideB = true;
                    }

                } ///
                
                
                // if the first border has been reached, go on the other direction to find the other border
                if (startPos != pos && startPos.IsBorder() && !polyBorderFound) {
                    startPos = pos;
                    startPos.FlipF();
                    polyBorderFound = true;
                    ///dxy add
                    sideA_repeat = 0;
                    sideB_repeat = 0;
                    ///
                }
                
                // if the other border has been reached, return
                if (polyBorderFound && startPos.IsBorder())
                    break;
                
                // go to the next edge
                startPos.FlipE();
                startPos.FlipV();
                startPos.FlipE();
                // check manifoldness
                if (IsVertexAdjacentToAnyNonManifoldEdge(startPos))
                    return PC_NOTMANIF;
                startPos.FlipV();
                if (IsVertexAdjacentToAnyNonManifoldEdge(startPos))
                    return PC_NOTMANIF;
                startPos.FlipV();
                // go to the next face
                startPos.FlipF();
            } while (startPos != pos);
            
            
            // polychord with singularities on both sides can not collapse
            /// dxy delete
            if ((singSideA && singSideB) //||
                //(singSideA && borderSideB) ||
                //(singSideB && borderSideA))
                )
                return PC_SINGBOTH;
            ///
            
            // polychords that are rings and have borders on both sides can not collapse
            if (!polyBorderFound && borderSideA && borderSideB)
                return PC_SINGBOTH;
            
            ///dxy add
            // consider border first
            if (borderSideA) {
                return PC_SINGSIDEA;
            }
            if (borderSideB) {
                return PC_SINGSIDEB;
            }
            // then singularities
            if (singSideA) {
                return PC_SINGSIDEA;
            }
            if (singSideB) {
                return PC_SINGSIDEB;
            }
            ///dxy add end
            
            
            ///dxy delete
//            // if there are singularities or borders on the side A, remember to keep coordinates on it
//            if (singSideA || borderSideA)
//                return PC_SINGSIDEA;
//            // if there are singularities or borders on the side B, remember to keep coordinates on it
//            if (singSideB || borderSideB)
//                return PC_SINGSIDEB;
            ///dxy delete end
            
            return PC_SUCCESS;
        }

        
        ///
        
        
        /**
         * @brief CollapsePolychord performs all checks and then collapses the polychord.
         *
         * @warning This function deletes faces and vertices by calling
         * vcg::tri::Allocator<PolyMeshType>::DeleteFace() and
         * vcg::tri::Allocator<PolyMeshType>::DeleteVertex().
         * The object PC_Chords chords is used to track the polychords, and it has got
         * a size proportional to that of the mesh face container. If you actually
         * delete faces and vertices by calling vcg::tri::Allocator<PolyMeshType>::CompactFaceVector()
         * and vcg::tri::Allocator<PolyMeshType>::CompactVertexVector() after this function,
         * object PC_Chords chords then is not valid any more, so you MUST rearrange it
         * by calling PC_Chords.Reset(). For the same reason, you MUST rearrange LinkConditions linkConditions
         * by calling LinkConditions.Resize().
         * However, for efficiency, you SHOULD compact vertex and face containers at the end of all your
         * polychord collapsing operations, without having to rearrange chords and linkConditions.
         * The function CollapseAllPolychords() does this for you.
         *
         * @note Vertex flags, face flags, FF adjacency and FV adjacency are required. Not anything else.
         * Such components are automatically updated here. If the mesh has other components that may be
         * affected by this editing, you should update them later by yourself.
         *
         * @param mesh The polygonal mesh used for getting the face index and deleting the faces
         * (it SHOULD have the vcg::face::PolyInfo component).
         * @param pos Position of the polychord.
         * @param mark Mark for the current polychord.
         * @param chords Vector of chords.
         * @param linkConditions Link conditions checker.
         * @param checkSing true if singularities on both sides are not allowed.
         * @return A PC_ResultCode resulting from checks or PC_SUCCESS if the collapse has been performed.
         */
        static PC_ResultCode CollapsePolychord (PolyMeshType &mesh,
                                                const vcg::face::Pos<FaceType> &pos,
                                                const unsigned long mark,
                                                PC_Chords &chords,
                                                LinkConditions &linkConditions,
                                                const bool checkSing = true) {
            vcg::tri::RequireFFAdjacency(mesh);
            
            if (mesh.IsEmpty())
                return PC_VOID;
            
            if (pos.IsNull())
                return PC_VOID;
            
            vcg::face::Pos<FaceType> tempPos, startPos;
            
            // check if the sequence of facets is a polychord and find the starting coord
            PC_ResultCode resultCode = CheckPolychordFindStartPosition(pos, startPos, checkSing);
            // if not successful, visit the sequence for marking it and return
            if (resultCode != PC_SUCCESS && resultCode != PC_SINGSIDEA && resultCode != PC_SINGSIDEB) {
                // if not manifold, visit the entire polychord ending on the non-manifold edge
                if (resultCode == PC_NOTMANIF) {
                    tempPos = pos;
                    VisitPolychord(mesh, tempPos, chords, mark, resultCode);
                    if (tempPos.IsManifold() && !tempPos.IsBorder()) {
                        tempPos.FlipF();
                        VisitPolychord(mesh, tempPos, chords, mark, resultCode);
                    }
                    return resultCode;
                }
                // if not quad, visit all the polychords passing through this coord
                if (resultCode == PC_NOTQUAD) {
                    tempPos = startPos;
                    do {
                        if (!tempPos.IsBorder()) {
                            tempPos.FlipF();
                            VisitPolychord(mesh, tempPos, chords, mark, resultCode);
                            tempPos.FlipF();
                        }
                        tempPos.FlipV();
                        tempPos.FlipE();
                    } while (tempPos != startPos);
                    VisitPolychord(mesh, startPos, chords, mark, resultCode);
                    return resultCode;
                }
                VisitPolychord(mesh, startPos, chords, mark, resultCode);
                return resultCode;
            }
            // check if the link conditions are satisfied
            // if not satisfied, visit the sequence for marking it and return
            if (!linkConditions.CheckLinkConditions(mesh, startPos)) {
                VisitPolychord(mesh, startPos, chords, mark, PC_NOLINKCOND);
                return PC_NOLINKCOND;
            }
            // mark the polychord's chords
            MarkPolychords(mesh, startPos, chords, mark);
            // check if the polychord does not intersect itself
            // if it self-intersects, visit the polychord for marking it and return
            if (IsPolychordSelfIntersecting(mesh, startPos, chords, mark)) {
                VisitPolychord(mesh, startPos, chords, mark, PC_SELFINTERSECT);
                return PC_SELFINTERSECT;
            }
            // check if manifoldness remains
            // if it will loose manifoldness, visit the sequence for marking it and return
            if (!WillPolychordBeManifold(mesh, startPos, chords, mark)) {
                VisitPolychord(mesh, startPos, chords, mark, PC_NOMOREMANIF);
                return PC_NOMOREMANIF;
            }
            // at this point the polychord is collapsable, visit it for marking
            VisitPolychord(mesh, startPos, chords, mark, PC_SUCCESS);
            
            // now collapse
            CoordType point;
            //    int valenceA = 0, valenceB = 0;
            vcg::face::Pos<FaceType> runPos = startPos;
            vcg::face::JumpingPos<FaceType> tmpPos;
            //    bool onSideA = false, onSideB = false;
            vcg::face::Pos<FaceType> sideA, sideB;
            typedef std::queue<VertexPointer *> FacesVertex;
            typedef std::pair<VertexPointer, FacesVertex> FacesVertexPair;
            typedef std::queue<FacesVertexPair> FacesVertexPairQueue;
            FacesVertexPairQueue vQueue;
            typedef std::pair<FacePointer *, FacePointer> FFpPair;
            typedef std::pair<char *, char> FFiPair;
            typedef std::pair<FFpPair, FFiPair> FFPair;
            typedef std::queue<FFPair> FFQueue;
            FFQueue ffQueue;
            std::queue<VertexPointer> verticesToDeleteQueue;
            std::queue<FacePointer> facesToDeleteQueue;
            
            runPos = startPos;
            do {
                // compute new vertex
                point = (runPos.V()->P() + runPos.VFlip()->P()) / 2.f;
                if (checkSing) {
                    if (resultCode == PC_SINGSIDEA)
                        point = runPos.V()->P();
                    else if (resultCode == PC_SINGSIDEB)
                        point = runPos.VFlip()->P();
                }
                runPos.V()->P() = point;
                // list the vertex pointer of the faces on the other side to be updated
                vQueue.push(FacesVertexPair());
                vQueue.back().first = runPos.V();
                tmpPos.Set(runPos.F(), runPos.E(), runPos.V());
                tmpPos.FlipV();
                tmpPos.NextFE();    // go to next face
                while (tmpPos.F() != runPos.F()) {
                    if (tmpPos.F() != runPos.FFlip())
                        vQueue.back().second.push(&tmpPos.F()->V(tmpPos.VInd()));
                    tmpPos.NextFE();    // go to next face
                }
                
                
                // enqueue to delete the other vertex
                verticesToDeleteQueue.push(runPos.VFlip());
                
                // list the adjacencies
                sideA = runPos;
                sideA.FlipE();
                sideA.FlipF();
                sideB = runPos;
                sideB.FlipV();
                sideB.FlipE();
                sideB.FlipF();
                // first side
                if (!sideA.IsBorder()) {
                    ffQueue.push(FFPair(FFpPair(),FFiPair()));
                    ffQueue.back().first.first = &sideA.F()->FFp(sideA.E());
                    ffQueue.back().second.first = &sideA.F()->FFi(sideA.E());
                    if (!sideB.IsBorder()) {
                        ffQueue.back().first.second = sideB.F();
                        ffQueue.back().second.second = sideB.E();
                    } else {
                        ffQueue.back().first.second = sideA.F();
                        ffQueue.back().second.second = sideA.E();
                    }
                }
                // second side
                if (!sideB.IsBorder()) {
                    ffQueue.push(FFPair(FFpPair(),FFiPair()));
                    ffQueue.back().first.first = &sideB.F()->FFp(sideB.E());
                    ffQueue.back().second.first = &sideB.F()->FFi(sideB.E());
                    if (!sideA.IsBorder()) {
                        ffQueue.back().first.second = sideA.F();
                        ffQueue.back().second.second = sideA.E();
                    } else {
                        ffQueue.back().first.second = sideB.F();
                        ffQueue.back().second.second = sideB.E();
                    }
                }
                
                // enqueue to delete the face
                facesToDeleteQueue.push(runPos.F());
                
                // go on next edge/face
                runPos.FlipE();
                runPos.FlipV();
                runPos.FlipE();
                runPos.FlipF();
            } while (runPos != startPos && !runPos.IsBorder());
            assert(runPos == startPos || vcg::face::IsBorder(*startPos.F(),startPos.E()));
            if (runPos.IsBorder()) {
                // compute new vertex on the last (border) edge
                point = (runPos.V()->P() + runPos.VFlip()->P()) / 2.f;
                if (checkSing) {
                    if (resultCode == PC_SINGSIDEA)
                        point = runPos.V()->P();
                    else if (resultCode == PC_SINGSIDEB)
                        point = runPos.VFlip()->P();
                }
                runPos.V()->P() = point;
                // list the vertex pointer of the faces on the other side to be updated
                vQueue.push(FacesVertexPair());
                vQueue.back().first = runPos.V();
                tmpPos.Set(runPos.F(), runPos.E(), runPos.V());
                tmpPos.FlipV();
                tmpPos.NextFE();    // go to next face
                while (tmpPos.F() != runPos.F()) {
                    vQueue.back().second.push(&tmpPos.F()->V(tmpPos.VInd()));
                    tmpPos.NextFE();
                }
                
                // enqueue to delete the other vertex
                verticesToDeleteQueue.push(runPos.VFlip());
            }
            
            // update vertices
            while (!vQueue.empty()) {
                while (!vQueue.front().second.empty()) {
                    *vQueue.front().second.front() = vQueue.front().first;
                    vQueue.front().second.pop();
                }
                vQueue.pop();
            }
            
            // update adjacencies
            while (!ffQueue.empty()) {
                *ffQueue.front().first.first = ffQueue.front().first.second;
                *ffQueue.front().second.first = ffQueue.front().second.second;
                ffQueue.pop();
            }
            
            // delete faces
            while (!facesToDeleteQueue.empty()) {
                vcg::tri::Allocator<PolyMeshType>::DeleteFace(mesh, *facesToDeleteQueue.front());
                facesToDeleteQueue.pop();
            }
            
            // delete vertices
            while (!verticesToDeleteQueue.empty()) {
                vcg::tri::Allocator<PolyMeshType>::DeleteVertex(mesh, *verticesToDeleteQueue.front());
                verticesToDeleteQueue.pop();
            }
            
            return PC_SUCCESS;
        }
        
        
        /****************************************************************/
        ///dxy add: when cornerBit is enabled
        static PC_ResultCode CollapsePolychord (PolyMeshType &mesh,
                                                const vcg::face::Pos<FaceType> &pos,
                                                const unsigned long mark,
                                                PC_Chords &chords,
                                                LinkConditions &linkConditions,
                                                int cornerBit,
                                                const bool checkSing = true) {
            vcg::tri::RequireFFAdjacency(mesh);
            
            if (mesh.IsEmpty())
                return PC_VOID;
            
            if (pos.IsNull())
                return PC_VOID;
            
            vcg::face::Pos<FaceType> tempPos, startPos;
            
            // check if the sequence of facets is a polychord and find the starting coord
            PC_ResultCode resultCode = CheckPolychordFindStartPosition(pos, startPos, cornerBit, checkSing);
            // if not successful, visit the sequence for marking it and return
            if (resultCode != PC_SUCCESS && resultCode != PC_SINGSIDEA && resultCode != PC_SINGSIDEB) {
                // if not manifold, visit the entire polychord ending on the non-manifold edge
                if (resultCode == PC_NOTMANIF) {
                    tempPos = pos;
                    VisitPolychord(mesh, tempPos, chords, mark, resultCode);
                    if (tempPos.IsManifold() && !tempPos.IsBorder()) {
                        tempPos.FlipF();
                        VisitPolychord(mesh, tempPos, chords, mark, resultCode);
                    }
                    return resultCode;
                }
                // if not quad, visit all the polychords passing through this coord
                if (resultCode == PC_NOTQUAD) {
                    tempPos = startPos;
                    do {
                        if (!tempPos.IsBorder()) {
                            tempPos.FlipF();
                            VisitPolychord(mesh, tempPos, chords, mark, resultCode);
                            tempPos.FlipF();
                        }
                        tempPos.FlipV();
                        tempPos.FlipE();
                    } while (tempPos != startPos);
                    VisitPolychord(mesh, startPos, chords, mark, resultCode);
                    return resultCode;
                }
                VisitPolychord(mesh, startPos, chords, mark, resultCode);
                return resultCode;
            }
            // check if the link conditions are satisfied
            // if not satisfied, visit the sequence for marking it and return
            if (!linkConditions.CheckLinkConditions(mesh, startPos)) {
                VisitPolychord(mesh, startPos, chords, mark, PC_NOLINKCOND);
                return PC_NOLINKCOND;
            }
            // mark the polychord's chords
            MarkPolychords(mesh, startPos, chords, mark);
            // check if the polychord does not intersect itself
            // if it self-intersects, visit the polychord for marking it and return
            if (IsPolychordSelfIntersecting(mesh, startPos, chords, mark)) {
                VisitPolychord(mesh, startPos, chords, mark, PC_SELFINTERSECT);
                return PC_SELFINTERSECT;
            }
            // check if manifoldness remains
            // if it will loose manifoldness, visit the sequence for marking it and return
            if (!WillPolychordBeManifold(mesh, startPos, chords, mark)) {
                VisitPolychord(mesh, startPos, chords, mark, PC_NOMOREMANIF);
                return PC_NOMOREMANIF;
            }
            // at this point the polychord is collapsable, visit it for marking
            VisitPolychord(mesh, startPos, chords, mark, PC_SUCCESS);
            
            
            ///dxy mod: update cornerBit while collapse
            // now collapse
            CoordType point;
            //    int valenceA = 0, valenceB = 0;
            vcg::face::Pos<FaceType> runPos = startPos;
            vcg::face::JumpingPos<FaceType> tmpPos;
            //    bool onSideA = false, onSideB = false;
            vcg::face::Pos<FaceType> sideA, sideB;
            typedef std::queue<VertexPointer *> FacesVertex;
            typedef std::pair<VertexPointer, FacesVertex> FacesVertexPair;
            typedef std::queue<FacesVertexPair> FacesVertexPairQueue;
            FacesVertexPairQueue vQueue;
            typedef std::pair<FacePointer *, FacePointer> FFpPair;
            typedef std::pair<char *, char> FFiPair;
            typedef std::pair<FFpPair, FFiPair> FFPair;
            typedef std::queue<FFPair> FFQueue;
            FFQueue ffQueue;
            std::queue<VertexPointer> verticesToDeleteQueue;
            std::queue<FacePointer> facesToDeleteQueue;
            
            runPos = startPos;
            do {
                // compute new vertex
                point = (runPos.V()->P() + runPos.VFlip()->P()) / 2.f;
                if (checkSing) {
                    if (resultCode == PC_SINGSIDEA)
                        point = runPos.V()->P();
                    else if (resultCode == PC_SINGSIDEB)
                        point = runPos.VFlip()->P();
                }
                runPos.V()->P() = point;
                // list the vertex pointer of the faces on the other side to be updated
                vQueue.push(FacesVertexPair());
                vQueue.back().first = runPos.V();
                tmpPos.Set(runPos.F(), runPos.E(), runPos.V());
                tmpPos.FlipV();
                tmpPos.NextFE();    // go to next face
                while (tmpPos.F() != runPos.F()) {
                    if (tmpPos.F() != runPos.FFlip())
                        vQueue.back().second.push(&tmpPos.F()->V(tmpPos.VInd()));
                    tmpPos.NextFE();    // go to next face
                }
                
                ///dxy: update cornerBit flag
                if (runPos.VFlip()->IsUserBit(cornerBit)) {
                    runPos.V()->SetUserBit(cornerBit);
                }
                ///
                
                // enqueue to delete the other vertex
                verticesToDeleteQueue.push(runPos.VFlip());
                
                // list the adjacencies
                sideA = runPos;
                sideA.FlipE();
                sideA.FlipF();
                sideB = runPos;
                sideB.FlipV();
                sideB.FlipE();
                sideB.FlipF();
                // first side
                if (!sideA.IsBorder()) {
                    ffQueue.push(FFPair(FFpPair(),FFiPair()));
                    ffQueue.back().first.first = &sideA.F()->FFp(sideA.E());
                    ffQueue.back().second.first = &sideA.F()->FFi(sideA.E());
                    if (!sideB.IsBorder()) {
                        ffQueue.back().first.second = sideB.F();
                        ffQueue.back().second.second = sideB.E();
                    } else {
                        ffQueue.back().first.second = sideA.F();
                        ffQueue.back().second.second = sideA.E();
                    }
                }
                // second side
                if (!sideB.IsBorder()) {
                    ffQueue.push(FFPair(FFpPair(),FFiPair()));
                    ffQueue.back().first.first = &sideB.F()->FFp(sideB.E());
                    ffQueue.back().second.first = &sideB.F()->FFi(sideB.E());
                    if (!sideA.IsBorder()) {
                        ffQueue.back().first.second = sideA.F();
                        ffQueue.back().second.second = sideA.E();
                    } else {
                        ffQueue.back().first.second = sideB.F();
                        ffQueue.back().second.second = sideB.E();
                    }
                }
                
                // enqueue to delete the face
                facesToDeleteQueue.push(runPos.F());
                
                // go on next edge/face
                runPos.FlipE();
                runPos.FlipV();
                runPos.FlipE();
                runPos.FlipF();
            } while (runPos != startPos && !runPos.IsBorder());
            assert(runPos == startPos || vcg::face::IsBorder(*startPos.F(),startPos.E()));
            if (runPos.IsBorder()) {
                // compute new vertex on the last (border) edge
                point = (runPos.V()->P() + runPos.VFlip()->P()) / 2.f;
                if (checkSing) {
                    if (resultCode == PC_SINGSIDEA)
                        point = runPos.V()->P();
                    else if (resultCode == PC_SINGSIDEB)
                        point = runPos.VFlip()->P();
                }
                runPos.V()->P() = point;
                // list the vertex pointer of the faces on the other side to be updated
                vQueue.push(FacesVertexPair());
                vQueue.back().first = runPos.V();
                tmpPos.Set(runPos.F(), runPos.E(), runPos.V());
                tmpPos.FlipV();
                tmpPos.NextFE();    // go to next face
                while (tmpPos.F() != runPos.F()) {
                    vQueue.back().second.push(&tmpPos.F()->V(tmpPos.VInd()));
                    tmpPos.NextFE();
                }
                
                ///dxy: update cornerBit flag
                if (runPos.VFlip()->IsUserBit(cornerBit)) {
                    runPos.V()->SetUserBit(cornerBit);
                }
                ///
                
                // enqueue to delete the other vertex
                verticesToDeleteQueue.push(runPos.VFlip());
            }
            
            // update vertices
            while (!vQueue.empty()) {
                while (!vQueue.front().second.empty()) {
                    *vQueue.front().second.front() = vQueue.front().first;
                    vQueue.front().second.pop();
                }
                vQueue.pop();
            }
            
            // update adjacencies
            while (!ffQueue.empty()) {
                *ffQueue.front().first.first = ffQueue.front().first.second;
                *ffQueue.front().second.first = ffQueue.front().second.second;
                ffQueue.pop();
            }
            
            // delete faces
            while (!facesToDeleteQueue.empty()) {
                vcg::tri::Allocator<PolyMeshType>::DeleteFace(mesh, *facesToDeleteQueue.front());
                facesToDeleteQueue.pop();
            }
            
            // delete vertices
            while (!verticesToDeleteQueue.empty()) {
                vcg::tri::Allocator<PolyMeshType>::DeleteVertex(mesh, *verticesToDeleteQueue.front());
                verticesToDeleteQueue.pop();
            }
            
            return PC_SUCCESS;
        }
        ///

        
        /**
         * @brief CollapseAllPolychords finds and collapses all the polychords.
         * @param mesh The input polygonal mesh (it SHOULD have the vcg::face::PolyInfo component).
         * @param checkSing true if singularities on both sides of a polychord are not allowed.
         */
        static void CollapseAllPolychords (PolyMeshType &mesh, const bool checkSing = true) {
            vcg::tri::RequireFFAdjacency(mesh);
            
            if (mesh.IsEmpty())
                return;
            
            vcg::face::Pos<FaceType> pos;
            PC_ResultCode resultCode;
            std::pair<size_t, unsigned char> face_edge;
            // construct the link conditions checker
            LinkConditions linkConditions(mesh.vert.size());
            // construct the vector of chords
            PC_Chords chords(mesh);
            unsigned long mark = 0;
            
            ///test
            int nCollapse = 0;
            
            // iterate over all the chords
            while (!chords.End()) {
                // get the current coord
                chords.GetCurrent(face_edge);
                resultCode = chords[face_edge].q;
                assert(resultCode == PC_VOID);
                // construct a pos on the face and edge of the current coord
                pos.Set(&mesh.face[face_edge.first], face_edge.second, mesh.face[face_edge.first].V(face_edge.second));
                // (try to) collapse the polychord
                ///test
//                std::cout << nCollapse << std::endl;
                ///
                resultCode = CollapsePolychord(mesh, pos, mark, chords, linkConditions, checkSing);
                /// test: save mesh every time after a collapse
//                printResultCode(resultCode);
//                std::cout << std::endl;
//                std::string filename = "output_" + std::to_string(nCollapse) + ".obj";
//                vcg::tri::io::Exporter<PolyMeshType>::Save(mesh, filename.data());
                nCollapse += 1;
                ///
                // go to the next coord
                chords.Next();
                
                // increment the mark
                ++mark;
                if (mark == std::numeric_limits<unsigned long>::max()) {
                    chords.ResetMarks();
                    mark = 0;
                }
            }
        }
        
        
        
        /******************************************************/
        ///dxy add: when User cornerBit is enabled
        static void CollapseAllPolychords (PolyMeshType &mesh, int cornerBit, const bool checkSing = true) {
            vcg::tri::RequireFFAdjacency(mesh);
            
            if (mesh.IsEmpty())
                return;
            
            vcg::face::Pos<FaceType> pos;
            PC_ResultCode resultCode;
            std::pair<size_t, unsigned char> face_edge;
            // construct the link conditions checker
            LinkConditions linkConditions(mesh.vert.size());
            // construct the vector of chords
            PC_Chords chords(mesh);
            unsigned long mark = 0;
            
            ///test
            int nCollapse = 0;
            
            // iterate over all the chords
            while (!chords.End()) {
                // get the current coord
                chords.GetCurrent(face_edge);
                resultCode = chords[face_edge].q;
                assert(resultCode == PC_VOID);
                // construct a pos on the face and edge of the current coord
                pos.Set(&mesh.face[face_edge.first], face_edge.second, mesh.face[face_edge.first].V(face_edge.second));
                // (try to) collapse the polychord
                ///test
//                std::cout << nCollapse << std::endl;
                ///
                resultCode = CollapsePolychord(mesh, pos, mark, chords, linkConditions, cornerBit, checkSing);
                /// test: save mesh every time after a collapse
//                printResultCode(resultCode);
//                std::cout << std::endl;
//                std::string filename = "output_" + std::to_string(nCollapse) + ".obj";
//                vcg::tri::io::Exporter<PolyMeshType>::Save(mesh, filename.data());
                nCollapse += 1;
                ///
                // go to the next coord
                chords.Next();
                
                // increment the mark
                ++mark;
                if (mark == std::numeric_limits<unsigned long>::max()) {
                    chords.ResetMarks();
                    mark = 0;
                }
            }
        }
        ///
        
    };
    
    
    
    
}
}


#endif