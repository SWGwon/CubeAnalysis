#ifndef EVENT_HXX
#define EVENT_HXX

#include <ToolMainTrajectory.hxx>

#include "Primaries.hxx"
#include "Object.hxx"

/**
 * @brief Event class
 * @details Event class consists of std::vector<Object> and Object first object
 * in time, Object interaction vertex. \n
 * It also has Cube::Handle<Cube::ReconObjectContainer> and 
 * Cube::Event::G4TrajectoryContainer.
 */
class Event
{
    public:
        /**
         * @brief default constructor
         * @details Get Cube::Handle<Cube::ReconObjectContainer> and
         * Cube::Event::G4TrajectoryContainer from Cube::Event \n
         * make std::vector<Object> and sort it in time order. \n
         * Set first Object in time and interaction vertex.
         */
        Event()
        {
            SetTrajectories();
            SetReconObjects();
            Initialize();
            SortObjectByTime();
            SetVertex();
            SetFirstObject();
            SetNumberOfVertexAssociated();
        }

        /**
         * @brief Get Primaries of this event
         */
        const Primaries& GetPrimaryParticles() const;

        /**
         * @brief Get std::vector<Object> of this event
         */
        const std::vector<Object>& GetObjects() const;

        /**
         * @brief Show information of this event
         * @details individual Object, vertex, first Object in time
         */
        const void Show() const;

        /**
         * @brief Get first Object in time
         */
        const Object& GetFirstObject() const;

        /**
         * @brief Get interaction vertex of this event
         */
        const Object& GetVertex() const;

        const int GetNumberOfVertexAssociated() const;

    private:
        /**
         * @brief assign Object,
         */
        void Initialize();

        /**
         * @brief sort Object in time order
         */
        void SortObjectByTime();

        /**
         * @brief vector if Object in this event
         */
        std::vector<Object> mObjects;

        /**
         * @brief Primaries of this event
         */
        Primaries mPrimaryParticles;

        /**
         * @brief Get Cube::Event::G4TrajectoryContainer from Cube::Event
         */
        void SetTrajectories();

        /**
         * @brief Trajectories of this event \n
         * access to trajectory: mTtrajectories[trajectory id]
         */
        Cube::Event::G4TrajectoryContainer mTtrajectories;

        /**
         * @brief Get Cube::Handle<Cube::ReconObjectContainer> from Cube::Event
         */
        void SetReconObjects();
        /**
         * @brief Objects of this event
         */
        Cube::Handle<Cube::ReconObjectContainer> mReconObjects;

        /**
         * @brief Set first Object in time
         */
        void SetFirstObject();

        /**
         * @brief first Object in time
         */
        Object mFirstObject;

        /**
         * @brief Set interaction vertex of this event
         */
        void SetVertex();

        /**
         * @brief interaction vertex of this event
         */
        Object mVertex;

        void SetNumberOfVertexAssociated();
        int mNumberOfVertexAssociated;
};

#endif
