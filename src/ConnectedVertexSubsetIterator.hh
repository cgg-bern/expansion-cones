//
// Created by Valentin Nigolian on 29.10.21.
//

#ifndef OPENFLIPPER_CONNECTEDVERTEXSUBSETITERATOR_H
#define OPENFLIPPER_CONNECTEDVERTEXSUBSETITERATOR_H

#include <boost/fusion/include/hash.hpp>

#include "CommonMeshDefinitions.hh"
#include <map>
#include <queue>



namespace OpenVolumeMesh {


    class ConnectedVertexSubset : public std::vector<VertexHandle> {

    public:
        ConnectedVertexSubset(int k) : std::vector<VertexHandle>(k) {}

        ConnectedVertexSubset() : ConnectedVertexSubset(0) {}


        std::set<int> map_key() const {
            std::set<int> key;
            for (auto v: *this) {
                key.insert(v.idx());
            }
            return key;
        }

        void set_from_key(const std::set<int>& key){
            this->clear();
            for(auto v: key){
                this->push_back(VertexHandle(v));
            }
        }

    private:

    };



    /**
     * NOTE: Connected vertex subsets are called "clusters" internally because it's shorter */
    using VertexCluster = ConnectedVertexSubset;

    std::ostream &operator<<(std::ostream &stream, const ConnectedVertexSubset &cluster);
    std::ostream &operator<<(std::ostream &stream, const std::set<int> &key);
    //std::ostream& operator<<(std::ostream& stream, const std::vector<VertexHandle>& vertices);

    //*****************************************************************************  ITERATOR


    class ConnectedVertexSubsetRecursiveIterator;

#if 0//TODO: make wrapper to use without pointer to iterator k-1
    class ConnectedVertexSubsetIterator {

    public:


        ConnectedVertexSubsetIterator(int k,
                                      TetrahedralMesh &mesh,
                                      const VertexPropertyT<bool> &candidate_vertex_prop,
                                      bool negate_property = false);

        virtual ~ConnectedVertexSubsetIterator();

        ConnectedVertexSubset next();

        bool is_valid() const;

        ConnectedVertexSubsetRecursiveIterator* p_k_iterator_ = nullptr;
    };
#endif


    //*****************************************************************************  RECURSIVE ITERATOR
    /**
     * NOTE: this iterator works by considering that all connected k-clusters can be obtained by taking
     * all the (k-1)-clusters, and then iteratively adding each vertex that is a neighbor of any of the
     * vertices of the (k-1)-cluster.
     * Hence we can use this class recursively.
     * Then, the idea is that if you don't need to just find the k-clusters but also the (k-j)-clusters
     * (j={1,...,k-1)), then you should only need to enumerate cluster of each k only once and store
     * them in a list (or in this case, a map).
     * Hence you should always provide a pointer to an iterator of (k-1)-clusters so you don't have to
     * recompute them all again.
     *
     * If you really just need the k-clusters, then you can use the wrapper above (not done yet)
     * */
    class ConnectedVertexSubsetRecursiveIterator {

    public:



        ConnectedVertexSubsetRecursiveIterator(int k,
                                               TetrahedralMesh &mesh,
                                               const VertexPropertyT<bool> &candidate_vertex_prop,
                                               ConnectedVertexSubsetRecursiveIterator* p_done_sub_iterator_for_k_minus_one,
                                               bool negate_property = false);

        virtual ~ConnectedVertexSubsetRecursiveIterator();

        ConnectedVertexSubset next();

        bool is_valid() const;

    private:

        void update_current_cluster();

        //update the current cluster by using a sub-iterator to find the next sub-cluster
        void update_current_cluster_with_sub_iterator();

        //update the current cluster with a map look-up
        void update_current_cluster_with_exploration_map();

        void invalidate();

        void reactivate();

        const int k_;

        TetrahedralMesh &mesh_;

        //NOTE: copy so we can modify it
        VertexPropertyT<bool> candidate_vertex_prop_;

        bool negate_candidate_property_;


        std::queue<VertexHandle> neighbors_queue_;

        ConnectedVertexSubset current_cluster_;

        std::map<std::set<int>, bool> exploration_map_;

        bool is_valid_;

        ConnectedVertexSubsetRecursiveIterator* p_sub_iterator_;

        std::map<std::set<int>, bool>::iterator exploration_map_iterator_;

    };


}
#endif //OPENFLIPPER_CONNECTEDVERTEXSUBSETITERATOR_H
