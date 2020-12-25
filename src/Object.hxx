#ifndef TRACK_HXX
#define TRACK_HXX

#include <CubeReconTrack.hxx>
#include <CubeReconNode.hxx>
#include <CubeReconCluster.hxx>

/**
 * @brief Object class for Event.
 * @details Object class is either Cube::Handle<Cube::ReconTrack> 
 * or Cube::Handle<Cube::ReconCluster> + pdg, parent pdg and parentId information \n
 * It has both of the two reconstructed object as a member but only one of them 
 * is activated when the Object class is initialized. \n
 */
class Object
{
    public:
        /**
         * default constructor
         */
        Object();

        /**
         * @brief track constructor
         * @param Cube::Handle<Cube::ReconTrack>& sourceTrack
         * @details copy sourceTrack as a member track. \n
         * Set position of this Object as a front position of sourceTrack.
         */
        Object(const Cube::Handle<Cube::ReconTrack>& sourceTrack)
        {
            mTrack = sourceTrack;
            IsTrack = true;
            mPosition = sourceTrack->GetPosition();
            Cube::Handle<Cube::TrackState> backState = sourceTrack->GetNodes().back()->GetState();
            mBackPosition = backState->GetPosition();
            mEdep = sourceTrack->GetEDeposit();
        }

        /**
         * @brief cluster constructor
         * @param Cube::Handle<Cube::ReconCluster>& sourceCluster
         * @details copy sourceTrack as a member cluster. \n
         * Set position of this Object as a front position of sourceCluster.
         */
        Object(const Cube::Handle<Cube::ReconCluster>& sourceCluster)
        {
            mCluster = sourceCluster;
            IsCluster = true;
            mPosition = sourceCluster->GetPosition();
            mBackPosition = mPosition;
            mEdep = sourceCluster->GetEDeposit();
        }

        /**
         * @brief assignment operator
         */
        Object& operator=(const Object& rhs);

        /**
         * @brief Set pdg code of this Object
         */
        void SetPdg(int inPdg);

        /**
         * @brief Get pdg code of this Object
         */
        const int GetPdg() const;

        /**
         * @brief Set parent pdg code of this Object
         */
        void SetParentPdg(int inParentPdg);

        /**
         * @brief Get parent pdg code of this Object
         */
        const int GetParentPdg() const;

        /**
         * @brief Set parentId of this Object
         */
        void SetParentId(int inParentId);

        /**
         * @brief Get parentId of this Object
         */
        const int GetParentId() const;

        /**
         * @brief Get energy deposit of this Object
         */
        const double GetEDeposit() const;

        /**
         * @brief Get position of this Object
         * @returns TLorentzVector mPosition
         */
        const TLorentzVector& GetPosition() const;

        /**
         * @brief Get back position of this Object
         * @returns TLorentzVector mBackPosition
         */
        const TLorentzVector& GetBackPosition() const;

        /**
         * @brief Show information of this Object
         * @details (x, y, z, t), pdg, parentId, parentpdg \n
         * When the Object is Track, also show back state position
         */
        const void Show() const;

        /**
         * @brief flag for track
         */
        bool IsTrack = false;

        /**
         * @brief flag for cluster
         */
        bool IsCluster = false;

        /**
         * @brief flag for first object
         */
        bool IsFirstObject = false;

        /**
         * @brief flag for vertex
         */
        bool IsVertex = false;

        /**
         * @brief flag for last object
         */
        bool IsLast = false;

        /**
         * @brief return member track
         * @returns Cube::Handle<Cube::ReconTrack> mTrack
         */
        const Cube::Handle<Cube::ReconTrack> GetTrack() const;

        /**
         * @brief return member track
         * @returns Cube::Handle<Cube::ReconCluster> mCluster
         */
        const Cube::Handle<Cube::ReconCluster> GetCluster() const;

    private:
        /**
         * @brief member track
         */
        Cube::Handle<Cube::ReconTrack> mTrack;

        /**
         * @brief member cluster
         */
        Cube::Handle<Cube::ReconCluster> mCluster;

        /**
         * @brief position of this Object
         */
        TLorentzVector mPosition;

        /**
         * @brief back position of this Object
         */
        TLorentzVector mBackPosition;

        /**
         * @brief pdg code of this Object
         */
        int mPdg;

        /**
         * @brief parent pdg code of this Object
         */
        int mParentPdg;

        /**
         * @brief parent id of this Object
         */
        int mParentId;

        /**
         * @brief energy deposit
         */
        double mEdep;

        /**
         * @brief for operator=
         */
        friend void swap(Object& first, Object& second) noexcept;
};

#endif
