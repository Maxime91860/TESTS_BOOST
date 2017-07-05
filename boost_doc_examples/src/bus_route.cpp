// Pointers


class bus_route
{
    friend class boost::serialization::access;
    bus_stop * stops[10];
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        int i;
        for(i = 0; i < 10; ++i)
            ar & stops[i];
    }
public:
    bus_route(){}
};