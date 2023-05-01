#pragma once

#include <map>
#include <queue>

#include "CommonMeshDefinitions.hh"

namespace OpenVolumeMesh{


class SubTriangleFace{

public:
    SubTriangleFace(const FaceVertices& vertices);

    SubTriangleFace(const std::initializer_list<VertexHandle> vertices);


    bool operator<(const SubTriangleFace& other) const;
    bool operator==(const SubTriangleFace& other) const;
    bool operator!=(const SubTriangleFace& other) const;

    /* \return VertexHandle(-1) if the index is out of bounds */
    VertexHandle operator[](const int& index) const;

    const FaceVertices& get_vertices() const;

    SubTriangleFace opposite_face() const;

    //replaces vertex 'source' with vertex 'target' in the SubTriangleFace.
    //does nothing if source is not present
    bool replace_vertex(const VertexHandle& source,
                        const VertexHandle& target);

    bool contains_edge(const VertexHandle& from_vertex,
                       const VertexHandle& to_vertex) const;

    bool contains_vertex(const VertexHandle& vertex) const;


private:
    std::vector<VertexHandle> vertices_;
};


using SubFaceSet = std::set<SubTriangleFace>;

//this is more of an esthetics alias
using EquatorialDisc = SubFaceSet;




class SubTriangleMap : public std::map<SubTriangleFace, std::set<SubTriangleFace>>{

public:

    using BaseMap = std::map<SubTriangleFace, std::set<SubTriangleFace>>;

    /* overrides the map::insert() method to insert the subfaces for both
     * face orientations at once
     * runs in O(n) with n being the number of subfaces */
    void insert(const SubTriangleFace& face,
                const std::vector<SubTriangleFace>& subfaces);


    /*void insert(const SubTriangleFace& face,
                const std::set<SubTriangleFace>& subfaces);*/

    void insert(const SubTriangleFace& face,
                const SubTriangleFace& subface);

    /* override map::erase() to erase both sides of the face */
    void erase(const SubTriangleFace& face);

    /* NOTE: based on set::find() so same complexity */
    bool contains(const SubTriangleFace& face) const;

#if 0
    bool merge_node_with_super_node(const SubTriangleFace& face,
                                    const SubTriangleFace& super_face);

#endif


    void extract_all_subvertices(const SubTriangleFace& face,
                                 std::set<VertexHandle>& equatorial_disc_vertices) const;

    void extract_all_subfaces(const SubTriangleFace& face,
                             EquatorialDisc& equatorial_disc) const;

    void extract_subvertices_intersection(const SubTriangleFace& first_face,
                                          const SubTriangleFace& second_face,
                                          std::set<VertexHandle>& subvertices_intersection) const;


    void replace_vertex_in_all_faces_containing_vertex(const VertexHandle& source,
                                                       const VertexHandle& target);


#if 0
    bool replace_vertex_in_map_and_update(const VertexHandle& source,
                                          const VertexHandle& target) const;

    SubTriangleMap extract_sub_map(const SubTriangleFace& face) const;
#endif
};



    std::ostream& operator<<(std::ostream& os, const SubTriangleFace& face);

}

