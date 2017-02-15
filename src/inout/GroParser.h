//
// Created by gageat on 16/11/16.
//

#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>

#include <boost/algorithm/string.hpp>

struct Box {
    float x;
    float y;
    float z;
};

struct Coordinate {
    float x;
    float y;
    float z;
};

class GroParser {

public:
    GroParser(std::string fileName);

    std::vector<int> getIndexForSolvent(std::string atomName);

    Box box() const { return m_box; };

    int nSolutes() const { return m_solutes.size(); };

    Coordinate getCoordinateForSolute(int const index) const {
        assert(index<m_solutes.size());
        return { m_solutes[index].coordinate };
    }

    std::string getTypeForSolute(int const index) const {
        assert(index<m_solutes.size());
        return m_solutes[index].atomName.substr(0,1);
    }

private:
    void parseFile();

    struct Atom {
        int resIndex;
        std::string resName;
        std::string atomName;
        int atomIndex;
        Coordinate coordinate;
        float xVelocity;
        float yVelocity;
        float zVelocity;
        int globalIndex;
    };

    std::vector<Atom> m_solutes;
    std::vector<Atom> m_solvent;
    std::string m_fileName;

    Box m_box;
};

