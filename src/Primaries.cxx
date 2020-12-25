#include <ToolPrimaryId.hxx>

#include "Primaries.hxx"

Cube::Event *event;

void Primaries::Initailize()
{
    Cube::Event::G4TrajectoryContainer& trajectories = event->G4Trajectories;
    mPrimId = Cube::Tool::AllPrimaries(*event);
    for (std::vector<int>::iterator tr = mPrimId.begin();
            tr != mPrimId.end(); ++tr) 
    {
        Cube::Handle<Cube::G4Trajectory> traj = trajectories[*tr];
        if (!traj) 
        {
            std::cout << "This cannot happen!";
            throw std::runtime_error("Invalid primary id");
        }
        switch (traj->GetPDGCode()) {
            case 11: case -11: mElectrons.push_back(*tr); break;
            case 13:  mMuons.push_back(*tr); break;
            case -13: mAntiMuons.push_back(*tr); break;
            case 22: mGammas.push_back(*tr); break;
            case 111: case 211: case -211: mPions.push_back(*tr); break;
            case 2112: mNeutrons.push_back(*tr); break;
            case 2212: mProtons.push_back(*tr); break;
            default: mOthers.push_back(*tr); break;
        }
    }
}

const int Primaries::GetNumberOfMuon() const
{
    return this->mMuons.size();
}

const int Primaries::GetNumberOfAntiMuon() const
{
    return this->mAntiMuons.size();
}

const int Primaries::GetNumberOfElectron() const
{
    return this->mElectrons.size();
}

const int Primaries::GetNumberOfGamma() const
{
    return this->mGammas.size();
}

const int Primaries::GetNumberOfPion() const
{
    return this->mPions.size();
}

const int Primaries::GetNumberOfNeutron() const
{
    return this->mNeutrons.size();
}

const int Primaries::GetNumberOfProton() const
{
    return this->mProtons.size();
}

const int Primaries::GetNumberOfOther() const
{
    return this->mOthers.size();
}

const int Primaries::GetNumberOfPrimary() const
{
    return this->mPrimId.size();
}

