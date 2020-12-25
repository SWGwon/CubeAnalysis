#include <CubeReconNode.hxx>

#include "Object.hxx"

Object::Object()
{
    IsTrack = false;
    IsCluster = false;
    mPdg = -1;
    mParentPdg = -1;
    mParentId = -1;
}

Object& Object::operator=(const Object& rhs)
{
    if (this == &rhs)
    {
        return *this;
    }

    Object tempObject(rhs);
    swap(*this, tempObject);
    return *this;
}

void swap(Object& first, Object& second) noexcept
{
    std::swap(first.IsTrack, second.IsTrack);
    std::swap(first.IsCluster, second.IsCluster);
    std::swap(first.mPdg, second.mPdg);
    std::swap(first.mParentPdg, second.mParentPdg);
    std::swap(first.mParentId, second.mParentId);
    std::swap(first.mPosition, second.mPosition);
    std::swap(first.mTrack, second.mTrack);
    std::swap(first.mCluster, second.mCluster);
}

void Object::SetPdg(int inPdg)
{
    this->mPdg = inPdg;
}

const int Object::GetPdg() const
{
    return this->mPdg;
}

void Object::SetParentPdg(int inParentPdg)
{
    this->mParentPdg = inParentPdg;
}

const int Object::GetParentPdg() const
{
    return this->mParentPdg;
}

void Object::SetParentId(int inParentId)
{
    this->mParentId = inParentId;
}

const int Object::GetParentId() const
{
    return this->mParentId;
}

const void Object::Show() const
{
    if (this->IsCluster)
    {
        std::cout << "cluster:";
        std::cout << " (";
        std::cout.width(9);
        std::cout.fill(' ');
        std::cout << this->GetPosition().X();
        std::cout << ", ";
        std::cout.width(9);
        std::cout.fill(' ');
        std::cout << this->GetPosition().Y();
        std::cout << ", ";
        std::cout.width(9);
        std::cout.fill(' ');
        std::cout << this->GetPosition().Z();
        std::cout << ", ";
        std::cout.width(9);
        std::cout.fill(' ');
        std::cout << this->GetPosition().T();
        std::cout << ")";
        std::cout << ", pdg: ";
        std::cout.width(5);
        std::cout.fill(' ');
        std::cout << this->GetPdg();
        std::cout << ", ParentId: ";
        std::cout.width(3);
        std::cout.fill(' ');
        std::cout << this->GetParentId();
        std::cout << ", ParentPdg: ";
        std::cout.width(5);
        std::cout.fill(' ');
        std::cout << this->GetParentPdg() << std::endl;
    }
    if (this->IsTrack)
    {
        std::cout << "track" << std::endl;
        if (!this->IsVertex && !this->IsFirstObject && !this->IsLast)
        {
            std::cout << "| | |-";
        }
        else if (this->IsVertex || this->IsLast)
        {
            std::cout << "|   |-";
        }
        else
        {
            std::cout << "    |-";
        }
        std::cout << "front:";
        std::cout << " (";
        std::cout.width(9);
        std::cout.fill(' ');
        std::cout << this->GetPosition().X();
        std::cout << ", ";
        std::cout.width(9);
        std::cout.fill(' ');
        std::cout << this->GetPosition().Y();
        std::cout << ", ";
        std::cout.width(9);
        std::cout.fill(' ');
        std::cout << this->GetPosition().Z();
        std::cout << ", ";
        std::cout.width(9);
        std::cout.fill(' ');
        std::cout << this->GetPosition().T();
        std::cout << ")";
        std::cout << ", pdg: ";
        std::cout.width(5);
        std::cout.fill(' ');
        std::cout << this->GetPdg();
        std::cout << ", ParentId: ";
        std::cout.width(3);
        std::cout.fill(' ');
        std::cout << this->GetParentId();
        std::cout << ", ParentPdg: ";
        std::cout.width(5);
        std::cout.fill(' ');
        std::cout << this->GetParentPdg() << std::endl;
        if (!this->IsVertex && !this->IsFirstObject && !this->IsLast)
        {
            std::cout << "| | |-";
        }
        else if (this->IsVertex || this->IsLast)
        {
            std::cout << "|   |-";
        }
        else
        {
            std::cout << "    |-";
        }
        std::cout << "back: ";
        std::cout << " (";
        std::cout.width(9);
        std::cout.fill(' ');
        Cube::Handle<Cube::TrackState> backState = this->GetTrack()->GetNodes().back()->GetState();
        TLorentzVector pos_back = backState->GetPosition();
        std::cout << pos_back.X();
        std::cout << ", ";
        std::cout.width(9);
        std::cout.fill(' ');
        std::cout << pos_back.Y();
        std::cout << ", ";
        std::cout.width(9);
        std::cout.fill(' ');
        std::cout << pos_back.Z();
        std::cout << ", ";
        std::cout.width(9);
        std::cout.fill(' ');
        std::cout << pos_back.T();
        std::cout << ")";
        std::cout << std::endl;
    }
}

const TLorentzVector& Object::GetPosition() const
{
    return this->mPosition;
}

const Cube::Handle<Cube::ReconTrack> Object::GetTrack() const
{
    if (!this->IsTrack)
    {
        throw std::runtime_error("not a track, Object->GetTrack() failed");
    }
    return this->mTrack;
}

const Cube::Handle<Cube::ReconCluster> Object::GetCluster() const
{
    if (!this->IsCluster)
    {
        throw std::runtime_error("not a cluster, Object->GetCluster() failed");
    }
    return this->mCluster;
}
