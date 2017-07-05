//Non Intrusive Version

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

class gps_position
{
public:
    int degrees;
    int minutes;
    float seconds;
    gps_position(){};
    gps_position(int d, int m, float s) :
        degrees(d), minutes(m), seconds(s)
    {}
};

namespace boost {
    namespace serialization {

    template<class Archive>
    void serialize(Archive & ar, gps_position & g, const unsigned int version)
    {
        ar & g.degrees;
        ar & g.minutes;
        ar & g.seconds;
    }

    } // namespace serialization
} // namespace boost