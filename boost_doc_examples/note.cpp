// Exemple Open Classrooms

#include <fstream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

class Note
{
private:
        friend class boost::serialization::access;

        template<class Archive>
        void serialize(Archive & ar, const unsigned int version) {
                ar & numerateur;
                ar & denominateur;
        }

        int numerateur;
        int denominateur;
public:
        Note() {};
        Note(int n, int d) :
                        numerateur(n), denominateur(d) {}
};


int main()
{
        std::ofstream ofs("fichierDeSerialisation");

        // Vous avez vu comme je travaille bien ? :)
        const Note maNoteDePhysisque(20,20);


        {
                boost::archive::text_oarchive oa(ofs);
                oa << maNoteDePhysisque;
        }

        /** Quelque temps plus tard… ***/

        Note monAncienneNote;

        {
                std::ifstream ifs("fichierDeSerialisation");
                boost::archive::text_iarchive ia(ifs);

                ia >> monAncienneNote;
        }

        return 0;
}
