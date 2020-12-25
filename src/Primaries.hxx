#ifndef PRIMARIES_HXX
#define PRIMARIES_HXX

#include <CubeEvent.hxx>

extern Cube::Event* event;

/**
 * @brief Primaries class
 * @details Primaries class consists of std::vector<int> particle; \n
 * muon, anti muon, electron, gamma, pion, neutron, proton, others \n
 * The particles are primary particles
 */
class Primaries
{
    public:
        /**
         * default constructor
         */
        Primaries()
        {
            Initailize();
        }

        /**
         * @brief Get number of primary muon
         */
        const int GetNumberOfMuon() const;

        /**
         * @brief Get number of primary anti muon
         */
        const int GetNumberOfAntiMuon() const;

        /**
         * @brief Get number of primary electron, it contains position
         */
        const int GetNumberOfElectron() const;

        /**
         * @brief Get number of primary gamma
         */
        const int GetNumberOfGamma() const;

        /**
         * @brief Get number of primary pion, it contains pi+-, pi0
         */
        const int GetNumberOfPion() const;

        /**
         * @brief Get number of primary nuetron
         */
        const int GetNumberOfNeutron() const;

        /**
         * @brief Get number of primary proton
         */
        const int GetNumberOfProton() const;

        /**
         * @brief Get number of primary the others
         */
        const int GetNumberOfOther() const;

        /**
         * @brief Get number of all primary particles
         */
        const int GetNumberOfPrimary() const;

    private:
        /**
         * @brief assign trajectoryId to the vectors
         */
        void Initailize();

        /**
         * @brief std::vector<int> of primary muon, int trajectoryId
         */
        std::vector<int> mMuons;

        /**
         * @brief std::vector<int> of primary muon, int trajectoryId
         */
        std::vector<int> mAntiMuons;

        /**
         * @brief std::vector<int> of primary electron (e+-), int trajectoryId
         */
        std::vector<int> mElectrons;

        /**
         * @brief std::vector<int> of primary gamma, int trajectoryId
         */
        std::vector<int> mGammas;

        /**
         * @brief std::vector<int> of primary pions (pi+-, pi0), int trajectoryId
         */
        std::vector<int> mPions;

        /**
         * @brief std::vector<int> of primary neutron, int trajectoryId
         */
        std::vector<int> mNeutrons;

        /**
         * @brief std::vector<int> of primary proton, int trajectoryId
         */
        std::vector<int> mProtons;

        /**
         * @brief std::vector<int> of primary the others, int trajectoryId
         */
        std::vector<int> mOthers;

        /**
         * @brief std::vector<int> of all primary particles, int trajectoryId
         */
        std::vector<int> mPrimId;
};

#endif
