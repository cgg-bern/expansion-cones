
#pragma once


#include <iostream>
#include <string>
#include <vector>

/**NOTE: I know there's a TON of JSON APIs out there but since my needs are very simple,
 * I prefer to write something than to add another dependency (even if it's just a header file)
 *
 * NOTE2: only works for writing numbers and arrays of numbers
 */

template<typename _STREAM>
class LightWeightJsonExporter{
public:


    LightWeightJsonExporter(_STREAM& stream) : stream_(stream){
        stream_ <<"{"<<std::endl;
    }


    void write(const std::string& field_name, const char* data){
        write(field_name, std::string(data));
    }

    void write(const std::string& field_name, const std::string& data){
        if(!first_write_){
            stream_ <<","<<std::endl;
        }
        first_write_ = false;
        stream_ << "\""<<field_name<<"\": \""<<data<<"\"";
    }

    template<typename _DATA>
    void write(const std::string& field_name, const _DATA& data){
        if(!first_write_){
            stream_ <<","<<std::endl;
        }
        first_write_ = false;
        stream_ << "\""<<field_name<<"\": "<<data;
    }

    template<typename _DATA>
    void write_vector(const std::string& field_name, const std::vector<_DATA>& data_vector){
        if(data_vector.empty()){
            return;
        }
        if(!first_write_){
            stream_ <<","<<std::endl;
        }
        first_write_ = false;

        stream_ << "\""<<field_name<<"\": [ ";
        for(int i(0); i<((int)data_vector.size())-1; i++){
            stream_ << data_vector[i]<<", ";
        }
        stream_<<data_vector.back()<<" ]";
    }


    void close(){
        stream_ << std::endl <<"}"<<std::endl;
    }

    _STREAM& stream(){ return stream_; }


private:
    _STREAM& stream_;
    bool first_write_ = true;

};
