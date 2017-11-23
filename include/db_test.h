#ifndef DB_TEST_H
#define DB_TEST_H

#include <string>



#include "common_header.h"


// mesh save
#include <wrap/io_trimesh/export.h>

#include <iostream>
#include <vector>




template <typename PolyMeshType>
class DB_Test {
public:
    
    static int cornerBit;
    
    typedef vcg::tri::pl::PatchFlattener<PolyMeshType> Flattener;
    
    
    
    
    
    template <typename MatrixType>
    static void print_Mat(MatrixType &m)
    {
        for (int i=0; i<m.size(); ++i) {
            auto &row = m[i];
            for (int j=0; j<row.size(); ++j) {
                std::cout << row[j] << "\t";
            }
            std::cout << std::endl;
        }
    }
    
    
    template <typename VectorType>
    static void print_Vec(VectorType &v)
    {
        for (int i=0; i<v.size(); ++i) {
            std::cout << v[i] << "\t";
        }
        std::cout << std::endl;
    }
    
    
    
    static void print_points(PolyMeshType &m)
    {
        std::cout << "x" << "\t\t" << "y" << std::endl;
        for (int i=0; i<m.vert.size(); ++i) {
            auto &p = m.vert[i].P();
            std::cout << p[0] << "\t\t" << p[1] << std::endl;
        }
        std::cout << std::endl;
    }
    
    
    static void print_corners(PolyMeshType &m)
    {
        if (cornerBit == -1) {
            std::cout << "No corner bit" << std::endl;
            return;
        }
        std::cout << "x" << "\t\t" << "y" << std::endl;
        for (int i=0; i<m.vert.size(); ++i) {
            if (m.vert[i].IsUserBit(cornerBit)) {
                auto &p = m.vert[i].P();
                std::cout << p[0] << "\t\t" << p[1] << std::endl;
            }
        }
        std::cout << std::endl;
    }
    
    static void print_corners(Patch<PolyMeshType> &p)
    {
        PolyMeshType &mesh = p.getMesh();
        auto &sidesLengthVec = p.getBoundaryLenVec();
        auto startPos = p.getStartCorner();  //not reference
        
        //check input
        if (startPos.IsNull()) {
            return;
        }
        if (sidesLengthVec.size() < 2) {
            return;
        }
        if (mesh.VN() != mesh.vert.size()) {  //some vertices deleted
            return;
        }
        if (!startPos.IsBorder()) {
            return;
        }
        
        auto runPos = startPos;
        std::cout << "x" << "\t\t" << "y" << std::endl;
        
        for(int i=0; i<sidesLengthVec.size(); ++i) {
            auto &coord = runPos.V()->P();
            std::cout << coord[0] << "\t\t" << coord[1] << std::endl;
            int sideLen = sidesLengthVec[i];
            while (sideLen > 0) {
                runPos.FlipV();
                runPos.FlipE();
                while (!runPos.IsBorder()) {
                    runPos.FlipF();
                    runPos.FlipE();
                }
                sideLen--;
            }
        }
        std::cout << std::endl;
    }
    
    
    
    template <typename MatrixType>
    static void calc_VVAdj(PolyMeshType &mesh, MatrixType &Adj)
    {
        int nVert = mesh.VN();
        Adj.resize(nVert);
        for (int i=0; i<nVert; ++i) {
            Adj[i].assign(nVert, 0);
        }

        int nFace = mesh.face.size();
        for (int i=0; i < nFace; ++i) {
            for (int j=0; j<4; ++j) {
                int v1 = std::distance(&mesh.vert[0], mesh.face[i].V(j));
                int v2 = std::distance(&mesh.vert[0], mesh.face[i].V((j+1)%4));
                Adj[v1][v2] = 1;
                Adj[v2][v1] = 1;
            }
        }
    }
    
    static void geometrize(Patch<PolyMeshType> &p, unsigned int radius)
    {
        auto &mesh = p.getMesh();
        auto varVec = p.getBoundaryVarVec();
        auto IDsEnds = p.getBoundaryIDsEnds();
        auto LenVec = p.getBoundaryLenVec();
        for(int i=0; i<varVec.size(); ++i) {
            LenVec[IDsEnds[i].first]  += varVec[i] - 1;
            LenVec[IDsEnds[i].second] += varVec[i] - 1;
        }
        Flattener::PatchFlatteningBySmoothing(mesh, p.getStartCorner(), LenVec, radius);

//        Flattener::PatchFlatteningBySmoothing(mesh, p.getStartCorner(), p.getBoundaryLenVec(), radius);
    }
    
    
    static void display(Patch<PolyMeshType> & p)
    {
        std::cout << "Number of faces: " << p.getNumberOfFaces() << std::endl;
        
        //    std::cout << "Topo String: " << p.getTopologyStr() << std::endl;
        
        
        std::cout << "Boundary IDs Vec: " << std::endl;
        auto &BID = p.getBoundaryIDsVec();
        print_Mat(BID);
        
        
        std::cout << "Boundary IDs Ends Vec: " << std::endl;
        auto &BIDE = p.getBoundaryIDsEnds();
        for (int i=0; i<BIDE.size(); ++i) {
            std::cout << "(" << BIDE[i].first << ", " << BIDE[i].second << ")" << std::endl;
        }
        
        std::cout << "Boundary Len Vec: " << std::endl;
        auto &BLen = p.getBoundaryLenVec();
        print_Vec(BLen);
        
        
        std::cout << "Boundary Var Vec: " << std::endl;
        auto &BVar = p.getBoundaryVarVec();
        print_Vec(BVar);
        
        std::cout << "Singularity Vec: " << std::endl;
        std::vector<std::pair<vcg::tri::pl::valence_type, vcg::tri::pl::num_type>> Sig;
        p.getSingularityVec(Sig);
        for (int i=0; i<Sig.size(); ++i) {
            std::cout << "(" << Sig[i].first << ", " << Sig[i].second << ")" << std::endl;
        }
        
        //Topology: Face-Face Adjacency
        std::cout << "Topo FFAdj: " << std::endl;
//        std::vector< std::vector<long> > FFAdj;
        std::vector< std::vector<vcg::tri::pl::num_type> > FFAdj;
        p.reconstructTopologyFromString(FFAdj);
        print_Mat(FFAdj);
        
        //Geometrize Patch
        auto &mesh = p.getMesh();
        unsigned int radius = 10;
//        Flattener::PatchFlatteningBySmoothing(mesh, p.getStartCorner(), p.getBoundaryLenVec(), radius);
        geometrize(p, radius);
        int vn = mesh.VN();
        std::cout << "number of vertices: " << mesh.VN() << std::endl;
        print_points(mesh);
        
        std::cout << "Corners: " << std::endl;
        initCornerBit(p);
        print_corners(mesh);
        ///test
        //print_corners(p);
        
        std::vector<std::vector<int>> VVAdj(vn);
        for (int i=0; i<VVAdj.size(); ++i) {
            VVAdj[i].resize(vn, 0);
        }
        calc_VVAdj(mesh, VVAdj);
        std::cout << "VVAdj: " << std::endl;
        print_Mat(VVAdj);
        
        
        
        
        std::cout << std::endl;
        
        
    }
    
    
    static void save_mesh(PolyMeshType & mesh, int nCorners, int index)
    {
        std::string filename = std::to_string(nCorners) + "_" + std::to_string(index) + ".obj";
        vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(mesh, filename.data(), 0);
    }
    
    
    static void save_mesh(PolyMeshType &mesh, std::string fname)
    {
        vcg::tri::io::ExporterOBJ<PolyMeshType>::Save(mesh, fname.data(), 0);
    }
    
    
    ///
    
    static void enable_corner_flag(PolyMeshType &mesh)
    {
        if (mesh.IsEmpty()) {
            return;
        }
        if (cornerBit == -1) {
            cornerBit = mesh.vert[0].NewBitFlag();
            ///test
            std::cout << "User Corner Bit: " << cornerBit << std::endl;
        }
    }
    
    static void disable_corner_flag(PolyMeshType &mesh)
    {
        if (mesh.IsEmpty()) {
            return;
        }
        if (cornerBit != -1) {
            //clear all corner bit
            for (int i=0; i<mesh.vert.size(); ++i) {
                mesh.vert[i].ClearUserBit(cornerBit);
            }
            mesh.vert[0].DeleteBitFlag(cornerBit);
            cornerBit = -1;
            ///test
            std::cout << "User Corner Bit deleted" << std::endl;
        }
    }
    
    
    static bool initCornerBit(Patch<PolyMeshType> &p)
    {
        PolyMeshType &mesh = p.getMesh();
        auto sidesLengthVec = p.getBoundaryLenVec();  //not reference
        auto startPos = p.getStartCorner();  //not reference
        
        //check input
        if (startPos.IsNull()) {
            return false;
        }
        if (sidesLengthVec.size() < 2) {
            return false;
        }
        if (mesh.VN() != mesh.vert.size()) {  //some vertices deleted
            return false;
        }
        if (!startPos.IsBorder()) {
            return false;
        }
        
        ///consider boundaryVarVec
        auto &varVec = p.getBoundaryVarVec();
        auto &IDsEnds = p.getBoundaryIDsEnds();
        for (int i=0; i<varVec.size(); ++i) {
            int end1 = IDsEnds[i].first;
            int end2 = IDsEnds[i].second;
            int wider = varVec[i] - 1;
            sidesLengthVec[end1] += wider;
            sidesLengthVec[end2] += wider;
        }
        
        enable_corner_flag(mesh);
        auto runPos = startPos;
        
        for(int i=0; i<sidesLengthVec.size(); ++i) {
            runPos.V()->SetUserBit(cornerBit);
            int sideLen = sidesLengthVec[i];
            while (sideLen > 0) {
                runPos.FlipV();
                runPos.FlipE();
                while (!runPos.IsBorder()) {
                    runPos.FlipF();
                    runPos.FlipE();
                }
                sideLen--;
            }
        }
        
        //                do {
        //                    runPos.FlipV();
        //                    runPos.FlipE();
        //                    if (!runPos.IsBorder()) {
        //                        runPos.FlipF();
        //                        runPos.FlipE();
        //                        if (!runPos.IsBorder()) {   // concave corner
        //                            while (!runPos.IsBorder()) {
        //                                runPos.FlipF();
        //                                runPos.FlipE();
        //                            }
        //                            runPos.V()->SetUserBit(cornerBit);
        //                        }
        //                    } else {    // convex corner
        //                        runPos.V()->SetUserBit(cornerBit);
        //                    }
        //                } while (runPos != startPos);
        
        return true;
    }
    
    
    static bool initCornerBit(PolyMeshType &mesh)
    {
        //init start corner
        vcg::face::Pos<typename PolyMeshType::FaceType> startPos;
        startPos.Set(&mesh.face[0], 0, mesh.face[0].V(0));
        
        
        //find a border
        vcg::face::Pos<typename PolyMeshType::FaceType> runPos = startPos;
        bool flipped = false;  //startPos has been flipped before?
        
        
        while (!runPos.IsBorder()) {
            //go to the next edge
            runPos.FlipE();
            runPos.FlipV();
            runPos.FlipE();
            runPos.FlipF();
            
            
            if (runPos == startPos) { //find a ring, change direction
                startPos.FlipE();
                if (flipped) {
                    return false;
                }
                else {
                    flipped = true;
                }
            }
        }
        
        startPos = runPos;
        // update the starting corner
        startPos.FlipE();
        while (!startPos.IsBorder()) {
            startPos.FlipF();
            startPos.FlipE();
            if (!startPos.IsBorder()) {     // concave corner
                startPos.FlipE();
                while (!startPos.IsBorder()) {
                    startPos.FlipF();
                    startPos.FlipE();
                }
                startPos.FlipE();
                break;
            } else {
                startPos.FlipV();
                startPos.FlipE();
            }
        }
        startPos.FlipE();
        
        
        //check input
        if (startPos.IsNull()) {
            return false;
        }
        //        if (sidesLengthVec.size() < 2) {
        //            return false;
        //        }
        if (mesh.VN() != mesh.vert.size()) {  //some vertices deleted
            return false;
        }
        if (!startPos.IsBorder()) {
            return false;
        }
        
        enable_corner_flag(mesh);
        runPos = startPos;
        do {
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
                    runPos.V()->SetUserBit(cornerBit);
                }
            } else {    // convex corner
                runPos.V()->SetUserBit(cornerBit);
            }
        } while (runPos != startPos);
        
        return true;
    }
    
    
    
    static bool collectCorners(PolyMeshType &mesh,
                               std::deque<vcg::face::Pos<typename PolyMeshType::FaceType>> &corners)
    {
        //clear
        corners.clear();
        
        //init startPos
        vcg::face::Pos<typename PolyMeshType::FaceType> startPos;
        startPos.Set(&mesh.face[0], 0, mesh.face[0].V(0));
        
        //find a border
        vcg::face::Pos<typename PolyMeshType::FaceType> runPos = startPos;
        bool flipped = false;  //startPos has been flipped before?
        
        
        while (!runPos.IsBorder()) {
            //go to the next edge
            runPos.FlipE();
            runPos.FlipV();
            runPos.FlipE();
            runPos.FlipF();
            
            if (runPos == startPos) { //find a ring, change direction
                startPos.FlipE();
                if (flipped) {
                    return false;
                }
                else {
                    flipped = true;
                }
            }
        }
        
        startPos = runPos;
        //run around the border to collect corners
        do {
            runPos.FlipV();
            runPos.FlipE();
            while (!runPos.IsBorder()) {
                runPos.FlipF();
                runPos.FlipE();
            }
            //            runPos.FlipV();
            if (runPos.V()->IsUserBit(cornerBit)) {
                corners.push_back(runPos);
            }
            
        } while(runPos != startPos);
        
        
        return true;
    }
    
    
    ///assume: cornerBit has been initialized
    static void calc_valenceVec(PolyMeshType &mesh,
                                std::vector<int> &vVec)
    {
        assert(mesh.face.size() == mesh.FN());
        assert(mesh.vert.size() == mesh.VN());
        
        vVec.clear();
        vVec.assign(mesh.VN(), -1);
        
        int nFace = mesh.FN();
        for (int i=0; i<nFace; ++i) {
            for (int j=0; j<4; ++j) {
                int v_id = std::distance(&mesh.vert[0], mesh.face[i].V(j));
                if (vVec[v_id] == -1) {
                    vcg::face::Pos<typename PolyMeshType::FaceType> pos(&mesh.face[i], j, mesh.face[i].V(j));
                    int valence = pos.NumberOfIncidentVertices();
                    if (pos.V()->IsUserBit(cornerBit)) { //pos.V is corner
                        valence += 2;
                    }
                    else {
                        vcg::face::JumpingPos<typename PolyMeshType::FaceType> jmpPos(pos.F(), pos.E(), pos.V());
                        if (jmpPos.FindBorder()) {
                            valence += 1;
                        }
                    }
                    vVec[v_id] = valence;
                }
            }
        }
    }
    
    
    void patchPrepare(PolyMeshType &mesh,
                      std::deque<vcg::face::Pos<typename PolyMeshType::FaceType>> &corners,
                      vcg::tri::pl::n_corners_type &nCorners)
    {
        //corners
        if(collectCorners(mesh, corners)){
            nCorners = corners.size();
        }
        else {
            std::cout << "Error in collectCorners()" << endl;
            nCorners = -1;
        }
    }
    
    
    
    /////////////////////////
    
    //check if each vertices on mesh border has an unique topoStr
    static bool IsTopostrUnique(PolyMeshType &mesh,
                                const std::deque<vcg::face::Pos<typename PolyMeshType::FaceType>> &corners)
    {
        std::map<std::string, int> topostr_count;
        std::string tmp_topostr;
        vcg::face::Pos<typename PolyMeshType::FaceType> startPos = corners.front();
        vcg::face::Pos<typename PolyMeshType::FaceType> runPos = startPos;
        do {
            runPos.FlipV();
            runPos.FlipE();
            while (!runPos.IsBorder()) {
                runPos.FlipF();
                runPos.FlipE();
            }
            Patch<PolyMeshType>::BuildTopologyStr(mesh, runPos, tmp_topostr);
            topostr_count[tmp_topostr]++;
        } while (runPos != startPos);
        
        //check uniqueness
        std::map<std::string, int>::iterator mIt = topostr_count.begin();
        for (; mIt != topostr_count.end(); ++mIt) {
            if (mIt->second > 1) {
                return false;
            }
        }
        
        return true;
    }
    
    
    //
    static void getBoundaryValenceStr(PolyMeshType &mesh,
                                      std::string &result, int& cnt)
    {
        //init start corner
        vcg::face::Pos<typename PolyMeshType::FaceType> startPos;
        startPos.Set(&mesh.face[0], 0, mesh.face[0].V(0));
        
        
        //find a border
        vcg::face::Pos<typename PolyMeshType::FaceType> runPos = startPos;
        bool flipped = false;  //startPos has been flipped before?
        
        
        while (!runPos.IsBorder()) {
            //go to the next edge
            runPos.FlipE();
            runPos.FlipV();
            runPos.FlipE();
            runPos.FlipF();
            
            if (runPos == startPos) { //find a ring, change direction
                startPos.FlipE();
                if (flipped) {
                    result = "err";
                    cnt = -1;
                    return;
                }
                else {
                    flipped = true;
                }
            }
        }
        
        startPos = runPos;
        
        // run around the border (clockwise?)
        std::vector<int> valenceVec;
        do {
            runPos.FlipV();
            runPos.FlipE();
            while (!runPos.IsBorder()) {
                runPos.FlipF();
                runPos.FlipE();
            }
            valenceVec.push_back(runPos.NumberOfIncidentVertices());
        } while (runPos != startPos);
        
        // valenceVec to string
        std::string valStr;
        for (int i=0; i<valenceVec.size(); ++i) {
            valStr.push_back('a' + valenceVec[i]);
        }
        
        // minimize dictionary order
        /// find min char
        char min = '~'; //ascii 126
        for (int i=0; i<valStr.size(); ++i) {
            if (valStr[i] < min) {
                min = valStr[i];
            }
        }
        
        /// locate all chars equal to min
        std::vector<int> minIndex;
        for (int i=0; i<valStr.size(); ++i) {
            if (valStr[i] == min) {
                minIndex.push_back(i);
            }
        }
        
        /// record all str begin with min char
        int len = valStr.size();
        std::vector<std::string> smallStrVec;
        std::string minStr = valStr.substr(minIndex[0], len-minIndex[0]) + valStr.substr(0, minIndex[0]);
        smallStrVec.push_back(minStr);
        for (int i=1; i<minIndex.size(); ++i) {
            std::string iStr = valStr.substr(minIndex[i], len-minIndex[i]) + valStr.substr(0, minIndex[i]);
            smallStrVec.push_back(iStr);
            if (iStr < minStr) {  //dictionary order?
                minStr = iStr;
            }
        }
        
        /// count min str
        int n_minStr = 0;
        for(int i=0; i<smallStrVec.size(); ++i) {
            if (smallStrVec[i] == minStr) {
                n_minStr++;
            }
        }
        
        // return
        result = minStr;
        cnt = n_minStr;
    }
    
    
    
    //
    static void getBoundaryTopoStr(PolyMeshType &mesh,
                                   std::string &result, int& cnt,
                                   int translation = 0,
                                   bool clockwise = true,
                                   bool flip_runPos = false)
    {
        //init start corner
        vcg::face::Pos<typename PolyMeshType::FaceType> startPos;
        startPos.Set(&mesh.face[0], 0, mesh.face[0].V(0));
        
        
        //find a border
        vcg::face::Pos<typename PolyMeshType::FaceType> runPos = startPos;
        bool flipped = false;  //startPos has been flipped before?
        
        
        while (!runPos.IsBorder()) {
            //go to the next edge
            runPos.FlipE();
            runPos.FlipV();
            runPos.FlipE();
            runPos.FlipF();
            
            if (runPos == startPos) { //find a ring, change direction
                startPos.FlipE();
                if (flipped) {
                    result = "err";
                    cnt = -1;
                    return;
                }
                else {
                    flipped = true;
                }
            }
        }
        
        startPos = runPos;
        
        /// test: perform startPos translation
        for (int i=0; i<translation; ++i) {
            startPos.FlipV();
            startPos.FlipE();
            while (!startPos.IsBorder()) {
                startPos.FlipF();
                startPos.FlipE();
            }
        }
        runPos = startPos;
        
        /// test: run around the border (counterclockwise)
        if (!clockwise) {
            startPos.FlipE();
            while (!startPos.IsBorder()) {
                startPos.FlipF();
                startPos.FlipE();
            }
            runPos = startPos;
        }
        
        // run around the border (clockwise?)
        std::vector<std::string> topoStrVec;
        do {
            runPos.FlipV();
            runPos.FlipE();
            while (!runPos.IsBorder()) {
                runPos.FlipF();
                runPos.FlipE();
            }
            std::string tStr;
            ///test
            if (!flip_runPos) {
                Patch<PolyMeshType>::BuildTopologyStr(mesh, runPos, tStr);
            }
            else {
                runPos.FlipV();
                Patch<PolyMeshType>::BuildTopologyStr(mesh, runPos, tStr);
                runPos.FlipV();
            }
            topoStrVec.push_back(tStr);
            
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
        
        // minimize dictionary order
        /// find min char
        char min = '~'; //ascii 126
        for (int i=0; i<idxStr.size(); ++i) {
            if (idxStr[i] < min) {
                min = idxStr[i];
            }
        }
        
        /// locate all chars equal to min
        std::vector<int> minIndex;
        for (int i=0; i<idxStr.size(); ++i) {
            if (idxStr[i] == min) {
                minIndex.push_back(i);
            }
        }
        
        /// record all str begin with min char
        int len = idxStr.size();
        std::vector<std::string> smallStrVec;
        std::string minStr = idxStr.substr(minIndex[0], len-minIndex[0]) + idxStr.substr(0, minIndex[0]);
        smallStrVec.push_back(minStr);
        for (int i=1; i<minIndex.size(); ++i) {
            std::string iStr = idxStr.substr(minIndex[i], len-minIndex[i]) + idxStr.substr(0, minIndex[i]);
            smallStrVec.push_back(iStr);
            if (iStr < minStr) {  //dictionary order?
                minStr = iStr;
            }
        }
        
        /// count min str
        int n_minStr = 0;
        for(int i=0; i<smallStrVec.size(); ++i) {
            if (smallStrVec[i] == minStr) {
                n_minStr++;
            }
        }
        
        // return
        result = minStr;
        cnt = n_minStr;
    }
    
    
    static void getBoundaryTopoStr(PolyMeshType &mesh,
                                   vcg::face::Pos<typename PolyMeshType::FaceType> startPos,
                                   std::string &result)
    {
        //check input
        if (!startPos.IsBorder()) {
            return;
        }
        
        // run around the border (clockwise)
        std::vector<std::string> topoStrVec;
        std::string tStr;
        vcg::face::Pos<typename PolyMeshType::FaceType> runPos = startPos;
        do {
            Patch<PolyMeshType>::BuildTopologyStr(mesh, runPos, tStr);
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
    
    
    
    //helper: if s1 and s2 equals after some translation
    static bool strEqual(std::string s1, std::string s2)
    {
        if (s1.size() != s2.size()) {
            return false;
        }
        int len = s1.size();
        std::string ts;
        for (int i=0; i<len; ++i) {
            ts = s1.substr(i, len-i) + s1.substr(0, i);
            if (ts == s2) {
                return true;
            }
        }
        return false;
    }
    
    
    
    
    
};


template <typename PolyMeshType>
int DB_Test<PolyMeshType>::cornerBit = -1;


#endif
