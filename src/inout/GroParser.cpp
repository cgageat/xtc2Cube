//
// Created by gageat on 16/11/16.
//

#include "GroParser.h"

using namespace std;

GroParser::GroParser(std::string fileName):
    m_fileName(fileName)
{
    parseFile();
}

void GroParser::parseFile()
{
    ifstream groInputFile (m_fileName);
    string line;

    if (groInputFile.is_open())
    {
        getline (groInputFile,line);// Comment line

        getline (groInputFile,line);// Number of Atoms
        boost::algorithm::trim(line);
        const int nAtoms = stoi(line);

        for(int i=0; i<nAtoms; i++)
        {
            getline (groInputFile,line);

            Atom atom;

            atom.globalIndex = i;

            atom.resIndex = stoi(line.substr(0,5));
            atom.resName = line.substr(5,5);
            boost::algorithm::trim(atom.resName);

            atom.atomName = line.substr(10,5);
            boost::algorithm::trim(atom.atomName);


            atom.atomIndex = stoi(line.substr(15,5));

            atom.coordinate = { stof(line.substr(20,8)),
                                stof(line.substr(28,8)),
                                stof(line.substr(36,8)) };

            if(line.size()>44) {
                atom.xVelocity = stof(line.substr(44,8));
                atom.yVelocity = stof(line.substr(52,8));
                atom.zVelocity = stof(line.substr(60,8));
            }

            if(atom.resName.compare("SOL") == 0) {
                m_solvent.push_back(atom);
            } else {
                m_solutes.push_back(atom);
            }

        }

        getline (groInputFile,line);// Number of Atoms
        boost::algorithm::trim(line);

        vector<string> splitVec;
        boost::algorithm::split( splitVec, line, boost::algorithm::is_any_of(" "), boost::algorithm::token_compress_on );

        assert(splitVec.size() == 3 && "ERROR: box size line bad formatted");

        m_box.x = stof(splitVec[0]);
        m_box.y = stof(splitVec[1]);
        m_box.z = stof(splitVec[2]);


        groInputFile.close();
    }
}


vector<int> GroParser:: getIndexForSolvent(std::string atomName)
{
    vector<int> indexVectors = vector<int>();

    for(auto atom: m_solvent) {
        if(atom.atomName.compare(atomName) == 0) {
            indexVectors.push_back(atom.globalIndex);
        }
    }
    return indexVectors;
}


