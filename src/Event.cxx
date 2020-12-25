#include "Event.hxx"

void Event::Initialize()
{
    for (Cube::ReconObjectContainer::iterator o = mReconObjects->begin();
            o != mReconObjects->end(); 
            ++o) 
    {
        Cube::Handle<Cube::ReconTrack> track = *o;
        if (track) 
        {
            int mainTraj = Cube::Tool::MainTrajectory(*event,*track);
            if (mainTraj == -1)
            {
                continue;
            }
            Cube::Handle<Cube::G4Trajectory> traj = mTtrajectories.at(mainTraj);
            if (traj) 
            {
                std::unique_ptr<Object> tempTrack = std::make_unique<Object> (track);
                tempTrack->SetPdg(traj->GetPDGCode());
                tempTrack->SetParentId(traj->GetParentId());
                if (traj->GetParentId() != -1)
                {
                    tempTrack->SetParentPdg(this->mTtrajectories.at(traj->GetParentId())->GetPDGCode());
                }
                else
                {
                    tempTrack->SetParentPdg(0);
                }
                this->mObjects.push_back(*tempTrack.get());
            }
            else
            {
                continue;
            }
        }
        Cube::Handle<Cube::ReconCluster> cluster = *o;
        if (cluster) 
        {
            int mainTraj = Cube::Tool::MainTrajectory(*event,*cluster);
            if (mainTraj == -1)
            {
                continue;
            }
            Cube::Handle<Cube::G4Trajectory> traj = mTtrajectories[mainTraj];
            if (traj) 
            {
                std::unique_ptr<Object> tempCluster = std::make_unique<Object> (cluster);
                tempCluster->SetPdg(traj->GetPDGCode());
                tempCluster->SetParentId(traj->GetParentId());
                if (traj->GetParentId() != -1)
                {
                    tempCluster->SetParentPdg(mTtrajectories[traj->GetParentId()]->GetPDGCode());
                }
                else
                {
                    tempCluster->SetParentPdg(0);
                }
                this->mObjects.push_back(*tempCluster);
            }
            else
            {
                continue;
            }
        }
    }
}

const Primaries& Event::GetPrimaryParticles() const
{
    return this->mPrimaryParticles;
}

void Event::SetTrajectories()
{
    this->mTtrajectories = event->G4Trajectories;
}

void Event::SetReconObjects()
{
    this->mReconObjects = event->GetObjectContainer();
    if (!mReconObjects)
    {
        throw std::runtime_error("event->GetObjectContainer() failed");
    }
}

const std::vector<Object>& Event::GetObjects() const
{
    return this->mObjects;
}

bool tSort(Object& object1, Object& object2)
{
    double time1 = object1.GetPosition().T();
    double time2 = object2.GetPosition().T();

    return time1 < time2;
}

void Event::SortObjectByTime()
{
    std::sort(mObjects.begin(), mObjects.end(), tSort);
    if (mObjects.size() != 0)
    {
        mObjects.at(mObjects.size()-1).IsLast = true;
    }
}

const void Event::Show() const
{
    std::cout << "|-Object" << std::endl;
    for (auto o : mObjects)
    {
        std::cout << "| |-";
        o.Show();
    }
    std::cout << "|-vertex" << std::endl;
    std::cout << "| |-";
    mVertex.Show();
    std::cout << "|-first object" << std::endl;
    std::cout << "  |-";
    if (mFirstObject.GetPosition().X() == 0)
    {
        std::cout << "no first object candidate" << std::endl;
    }
    else
    {
        mFirstObject.Show();
    }
    std::cout.width(102);
    std::cout.fill('=');
    std::cout << "=" << std::endl;
}

void Event::SetFirstObject()
{
    bool isThereFirstObject = false;
    for (auto o : this->GetObjects())
    {
        if (std::abs(o.GetPdg()) != 13)
        {
            this->mFirstObject = o;
            this->mFirstObject.IsFirstObject = true;
            isThereFirstObject = true;
            break;
        }
    }
    if (!isThereFirstObject)
    {
    }
}

void Event::SetVertex()
{
    for (auto o : this->GetObjects())
    {
        if (o.GetPdg() == -13 && o.GetParentId() == -1)
        {
            mVertex = o;
            mVertex.IsVertex = true;
            break;
        }
    }
}

const Object& Event::GetFirstObject() const
{
    return this->mFirstObject;
}

const Object& Event::GetVertex() const
{
    return this->mVertex;
}
