// Derived Classes


#include <boost/serialization/base_object.hpp>
#include <string>
#include <fstream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

class bus_stop_corner : public bus_stop
{
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        // serialize base class information
        ar & boost::serialization::base_object<bus_stop>(*this);
        ar & street1;
        ar & street2;
    }
    std::string street1;
    std::string street2;
    virtual std::string description() const
    {
        return street1 + " and " + street2;
    }
public:
    bus_stop_corner(){}
    bus_stop_corner(const gps_position & lat_, const gps_position & long_,
        const std::string & s1_, const std::string & s2_
    ) :
        bus_stop(lat_, long_), street1(s1_), street2(s2_)
    {}
};

int main() {
    // create and open a character archive for output
    std::ofstream ofs("filename");

    // create class instance
    const gps_position g(35, 59, 24.567f);
    const gps_position g2(35, 59, 24.567f);
    const bus_stop_corner bsc(g, g2, "3th street", "5th street");

    // save data to archive
    {
        boost::archive::text_oarchive oa(ofs);
        // write class instance to archive
        oa << g;
        // archive and stream closed when destructors are called
    }

    // ... some time later restore the class instance to its orginal state
    gps_position newg;
    {
        // create and open an archive for input
        std::ifstream ifs("filename");
        boost::archive::text_iarchive ia(ifs);
        // read class state from archive
        ia >> newg;
        // archive and stream closed when destructors are called
    }
    return 0;
}