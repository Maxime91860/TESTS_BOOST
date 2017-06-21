// Serializable Members

class bus_stop
{
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & latitude;
        ar & longitude;
    }
    gps_position latitude;
    gps_position longitude;
protected:
    bus_stop(const gps_position & lat_, const gps_position & long_) :
    latitude(lat_), longitude(long_)
    {}
public:
    bus_stop(){}
    // See item # 14 in Effective C++ by Scott Meyers.
    // re non-virtual destructors in base classes.
    virtual ~bus_stop(){}
};