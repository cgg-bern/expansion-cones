//
// Created by Valentin Nigolian on 29.10.21.
//

#include "ConnectedVertexSubsetIterator.hh"

namespace OpenVolumeMesh {


    std::ostream &operator<<(std::ostream &stream, const ConnectedVertexSubset &cluster) {
        for (auto v: cluster) {
            stream << v << " ";
        }
        return stream;
    }


    std::ostream &operator<<(std::ostream &stream, const std::set<int> &key) {
        for (auto v: key) {
            stream << v << " ";
        }
        return stream;
    }

    //*****************************************************************************  ITERATOR


#if 0
    ConnectedVertexSubsetIterator::ConnectedVertexSubsetIterator(int k,
                                                                 TetrahedralMesh &mesh,
                                                                 const VertexPropertyT<bool> &candidate_vertex_prop,
                                                                 bool negate_property){

    }

    ConnectedVertexSubsetIterator::~ConnectedVertexSubsetIterator(){

    }

    ConnectedVertexSubset ConnectedVertexSubsetIterator::next(){

    }

    bool ConnectedVertexSubsetIterator::is_valid() const{

    }
#endif


    //*****************************************************************************  RECURSIVE ITERATOR


    ConnectedVertexSubsetRecursiveIterator::ConnectedVertexSubsetRecursiveIterator(int k,
                                                                                   TetrahedralMesh &mesh,
                                                                                   const VertexPropertyT<bool> &candidate_vertex_prop,
                                                                                   ConnectedVertexSubsetRecursiveIterator* p_done_sub_iterator_for_k_minus_one,
                                                                                   bool negate_candidate_property)
            : k_(k), mesh_(mesh), candidate_vertex_prop_(candidate_vertex_prop),
              negate_candidate_property_(negate_candidate_property), is_valid_(true),
              p_sub_iterator_(p_done_sub_iterator_for_k_minus_one) {

        if (k < 1) {
            invalidate();
            std::cout<<" ERROR - "<<k<<"-sized Connected Vertex Subset"<<std::endl;
            return;
        }

        if(k > 1) {
            //std::cout<<" sub-iterator: "<<p_done_sub_iterator_for_k_minus_one<<std::endl;
            //std::cout<<" is valid: "<<p_done_sub_iterator_for_k_minus_one->is_valid()<<std::endl;
            if (p_done_sub_iterator_for_k_minus_one &&
                !p_done_sub_iterator_for_k_minus_one->is_valid()) {

                p_sub_iterator_->reactivate();
            } else {
                std::cout << " ERROR - nullptr sub-iterator for k-1" << std::endl;
                invalidate();
                return;
            }
        }

        int candidate_vertices_count(0);

        //copying the property
        //std::cout<<" candidate vertices (k="<<k_<<"): ";
        candidate_vertex_prop_ = mesh_.request_vertex_property<bool>();
        for (auto v: mesh_.vertices()) {
            candidate_vertex_prop_[v] = candidate_vertex_prop[v];
            candidate_vertices_count += (candidate_vertex_prop_[v] != negate_candidate_property_);
            /*if(candidate_vertex_prop_[v] != negate_candidate_property_) {
                std::cout << " " << v;
            }*/
        }
        //std::cout<<std::endl;


        if(candidate_vertices_count < k_) {
            std::cout << " ERROR - Requesting " << k_ << "-clusters but only " << candidate_vertices_count
                      << " candidate vertices. Invalidating." << std::endl;

            invalidate();
            return;
        }

        update_current_cluster();
        exploration_map_.insert({current_cluster_.map_key(), true});

        //std::cout<<" initial current cluster: "<<current_cluster_<<std::endl;

    }

    ConnectedVertexSubsetRecursiveIterator::~ConnectedVertexSubsetRecursiveIterator() {
        if (!p_sub_iterator_) {
            delete p_sub_iterator_;
        }
    }

    bool ConnectedVertexSubsetRecursiveIterator::is_valid() const {
        return is_valid_;
    }

    void ConnectedVertexSubsetRecursiveIterator::invalidate() {
        p_sub_iterator_ = nullptr;
        is_valid_ = false;
        current_cluster_.clear();

#warning TODO: figure out why they're invalidated multiple times'

        //std::cout<<" invalidated "<<this<<std::endl;

        /*std::cout<<" invalidated "<<this<<". map: "<<std::endl;

        for(auto pair = exploration_map_.begin(); pair != exploration_map_.end(); pair++){
            std::cout<<" - "<<pair->first<<std::endl;
        }*/
    }

    void ConnectedVertexSubsetRecursiveIterator::reactivate() {

        exploration_map_iterator_ = exploration_map_.begin();
        //std::cout<<" --> first key at reactivation: "<<exploration_map_iterator_->first<<std::endl;
        is_valid_ = true;

        //to mimick the constructor
        update_current_cluster();
    }


    ConnectedVertexSubset ConnectedVertexSubsetRecursiveIterator::next() {


        //std::cout << "----------------------------" << std::endl;
        //std::cout << " --------- gathering next cluster for " << k_ << "-cluster iterator " << this << std::endl;
        if (!is_valid()) {
            //std::cout << " --> no longer valid. done" << std::endl;
            return ConnectedVertexSubset(0);
        }


        auto current_cluster_copy = current_cluster_;

        this->update_current_cluster();

        if(k_>1 && p_sub_iterator_) {

            //NOTE: could be a do{}while; but I need the output for now
            while (exploration_map_.find(current_cluster_.map_key()) != exploration_map_.end() && is_valid()) {
                //std::cout << " ---> current cluster " << current_cluster_ << " already in map, looking for next one"<< std::endl;
                this->update_current_cluster();
            }
        }

        if(current_cluster_.empty()){
            //std::cout<<" ...done. Sub-iterator reached the end, terminating"<<std::endl;
            invalidate();
        }else {

            exploration_map_.insert({current_cluster_.map_key(), true});
            //std::cout << " done. Updated next cluster of iterator " << this << " to " << current_cluster_
            //          << " and inserted it in map" << std::endl;
            //std::cout << "----------------------------" << std::endl;
        }

        //std::cout<<" next "<<k_<<"-cluster is "<<current_cluster_<<std::endl;
        //std::cout << "----------------------------" << std::endl;

        return current_cluster_copy;

    }



    void ConnectedVertexSubsetRecursiveIterator::update_current_cluster(){

        //std::cout<<" - updating "<<k_<<"-cluster"<<std::endl;
        //base case
        if(k_ == 1){
            if(neighbors_queue_.empty()){
                //std::cout<<" -- empty neighbor queue"<<std::endl;
                if(current_cluster_.empty()) {

                    //std::cout<<"  --- 1-cluster queue: ";
                    for (auto v: mesh_.vertices()) {
                        if (candidate_vertex_prop_[v] != negate_candidate_property_) {
                            neighbors_queue_.push(v);
                            //std::cout<<" "<<v;
                        }
                    }
                    //std::cout<<std::endl;
                    //std::cout<<" --> "<<neighbors_queue_.size()<<" elements"<<std::endl;
                }else{
                    //std::cout<<" --> current cluster is not empty, invalidating"<<std::endl;
                    invalidate();
                    return;
                }
            }

            current_cluster_.clear();
            current_cluster_.push_back(neighbors_queue_.front());
            neighbors_queue_.pop();
            //std::cout<<" - updated 1-cluster to "<<current_cluster_<<std::endl;

        }else{
            if(p_sub_iterator_){
                update_current_cluster_with_sub_iterator();
            }else{
                update_current_cluster_with_exploration_map();
            }
        }
    }


    void ConnectedVertexSubsetRecursiveIterator::update_current_cluster_with_sub_iterator(){

        //std::cout<<" - updating with sub-iterator"<<std::endl;
        if(neighbors_queue_.empty()){

            current_cluster_ = p_sub_iterator_->next();
            //std::cout<<" - next sub-cluster: "<<current_cluster_<<std::endl;

            if(current_cluster_.empty()){
                //std::cout<<" --> reached end of sub-cluster. Deleting sub-cluster pointer and invalidating iterator."<<std::endl;
                invalidate();
                return;
            }

            auto added_prop = mesh_.request_vertex_property<bool>();

            for(auto v: current_cluster_) {
                added_prop[v] = true;
            }

            for(auto v: current_cluster_){
                for(auto out_he: mesh_.outgoing_halfedges(v)){
                    auto neighbor = mesh_.to_vertex_handle(out_he);
                    if(!added_prop[neighbor] &&
                       candidate_vertex_prop_[neighbor] != negate_candidate_property_){
                        neighbors_queue_.push(neighbor);
                    }
                }
            }


            //std::cout<<" - neighbors queue initial size = "<<neighbors_queue_.size()<<std::endl;

            //to get the cluster to the right size
            current_cluster_.push_back(VertexHandle(-1));
        }

        if(neighbors_queue_.empty()){
            //std::cout<<" ERROR - empty neighbors list"<<std::endl;
            invalidate();
            return;
        }


        //by now we should have the current cluster and the queue of neighbors to add ready.
        current_cluster_.back() = neighbors_queue_.front();
        neighbors_queue_.pop();
        //std::cout<<" - updated neighbors queue size = "<<neighbors_queue_.size()<<std::endl;

    }

    void ConnectedVertexSubsetRecursiveIterator::update_current_cluster_with_exploration_map(){

        //std::cout<<" - updating with exploration map"<<std::endl;
        if(exploration_map_iterator_ == exploration_map_.end()){
            invalidate();
            return;
        }

        current_cluster_.set_from_key(exploration_map_iterator_->first);
        //std::cout<<" - set current cluster from map as "<<current_cluster_<<std::endl;

        exploration_map_iterator_++;
    }

}
