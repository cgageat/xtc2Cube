#include <cstdio>
#include <iostream>
#include <string>
#include <cmath>

#include "../xdrfile-1.1.4/include/xdrfile.h"
#include "../xdrfile-1.1.4/include/xdrfile_xtc.h"

#include "inout/GroParser.h"
#include <boost/program_options.hpp>

using namespace std;
namespace options = boost::program_options;

int main (int argc, char* argv[])
{

    /////////////////////////////////////////////
    ////////////////   OPTIONS   ////////////////
    /////////////////////////////////////////////

    options::options_description desc (string (argv[0]).append(" options"));
    desc.add_options()
            ("h", "Display this message")
            ("f", options::value<string>()->default_value("MD_Every10ps.xtc"), "Traj file (.xtc)")
            ("s", options::value<string>()->default_value("MD.gro"), "Structure file (.gro)")
            ("d", options::value<float>()->default_value(0.10f), "Delta grid in nm")
            ;

    options::variables_map args;
    options::store (options::command_line_parser (argc, argv).options (desc)
                            .style (options::command_line_style::default_style |
                                    options::command_line_style::allow_long_disguise)
                            .run (), args);
    options::notify (args);

    if (args.count ("h")) {
        cout << desc << endl;
        return 0; }




    const float deltaGrid {args["d"].as<float>()};

    typedef map<float, map<float, map<float, float>>> grid_type;


	// Parse gro file
	GroParser groParser(args["s"].as<string>());
    vector<int> indexs = groParser.getIndexForSolvent("OW");


    //Read trajectory
	string fullFileName = args["f"].as<string>();
	
    XDRFILE *xd;
    int lastFrame;
    float *prec = new float();
    int nAtoms {0};
    matrix box;
	int *step = new int();
	float *time = new float();


	char *fileName = (char*)fullFileName.c_str();

    int response = read_xtc_natoms(fileName, &nAtoms);

    rvec *positions = (float(*)[3])malloc(3*nAtoms*sizeof(rvec));



    cout << "initialisation ..." << endl;

    xd = xdrfile_open(fileName,"r");

    float minX { numeric_limits<float>::infinity() };
    float maxX { -numeric_limits<float>::infinity() };
    float minY { numeric_limits<float>::infinity() };
    float maxY { -numeric_limits<float>::infinity() };
    float minZ { numeric_limits<float>::infinity() };
    float maxZ { -numeric_limits<float>::infinity() };

    int nFrames {0};
	while (response!=exdrENDOFFILE) {
		response = read_xtc(xd, nAtoms, step, time, box, positions, prec);
		if(response == exdrOK) {
            nFrames++;
            if(nFrames%100==0) {
                cout << "Frame en cours: " << nFrames <<"\r" << flush;
            }
            for (auto index: indexs) {
                if (positions[index][0] < minX) {
                    minX = positions[index][0];
                }
                if (positions[index][0] > maxX) {
                    maxX = positions[index][0];
                }
                if (positions[index][1] < minY) {
                    minY = positions[index][1];
                }
                if (positions[index][1] > maxY) {
                    maxY = positions[index][1];
                }
                if (positions[index][2] < minZ) {
                    minZ = positions[index][2];
                }
                if (positions[index][2] > maxZ) {
                    maxZ = positions[index][2];
                }
            }
		}
	}

    xdrfile_close(xd);

    cout << nFrames << " frames lues!" << endl;


    cout << "Préparation de la grille ..." << endl;
    grid_type grid;
    const int nXGrid { static_cast<int>( ceil((maxX-minX)/deltaGrid) ) };
    const int nYGrid { static_cast<int>( ceil((maxY-minY)/deltaGrid) ) };
    const int nZGrid { static_cast<int>( ceil((maxZ-minZ)/deltaGrid) ) };



    float xGrid {0.0f};
    for(int nx=0; nx<nXGrid; nx++) {
        grid[nx] = map<float, map<float, float>>();

        for(int ny=0; ny<nYGrid; ny++) {
            grid[nx][ny] = map<float, float>();

            for(int nz=0; nz<nZGrid; nz++) {
                grid[nx][ny][nz] = 0.0f;
            }
        }
    }




    cout << "Lecture des positions ..." << endl;

    xd = xdrfile_open(fileName,"r");

    response=exdrOK;
    int nFramesTmp {0};
    while (response!=exdrENDOFFILE) {
        response = read_xtc(xd, nAtoms, step, time, box, positions, prec);
        if(response == exdrOK) {
            nFramesTmp++;
            if(nFramesTmp%100==0) {
                cout << nFramesTmp << "/" << nFrames <<"\r" << flush;
            }
            for (auto index: indexs) {
                const float x = positions[index][0];
                const float y = positions[index][1];
                const float z = positions[index][2];

                grid[floor((x-minX)/deltaGrid+0.5)]
                    [floor((y-minY)/deltaGrid+0.5)]
                    [floor((z-minZ)/deltaGrid+0.5)] += 1.0f;
            }
        }
    }



    //normalization des densités
    cout << "Normalisation des densités ..." << endl;
    const float voxelVolume {deltaGrid*deltaGrid*deltaGrid};
    constexpr float rhoZero{0.0000333f };
    const float normalizationFactor {1.0f/(static_cast<float>(nFrames)*voxelVolume*rhoZero*1000000.0f)};

    for(int nx=0; nx<nXGrid; nx++) {
        for(int ny=0; ny<nZGrid; ny++) {
            for(int nz=0; nz<nZGrid; nz++) {
                grid[nx][ny][nz] *= normalizationFactor;
            }
        }
    }

    xdrfile_close(xd);




    cout << "Ecriture du fichier cube ..." << endl;
    ofstream densityCubeFile;
    densityCubeFile.open("density.cube");

    map<string, int> nameToZ { {"H",1},
                               {"C",6},
                               {"N",7},
                               {"O",8}
                             };

    constexpr float nmToBohr {18.8971616463};

    if (densityCubeFile)
    {
        densityCubeFile << " CPMD CUBE FILE." << endl;
        densityCubeFile << " OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z" << endl;
        densityCubeFile << groParser.nSolutes() << " " << minX << " " << minY << " " << minZ << endl;

        densityCubeFile << nXGrid << " " << deltaGrid*nmToBohr << " 0.0 0.0" << endl;
        densityCubeFile << nYGrid << " 0.0 " << deltaGrid*nmToBohr << " 0.0" << endl;
        densityCubeFile << nZGrid << " 0.0 0.0 " << deltaGrid*nmToBohr << endl;


        for(int index=0; index<groParser.nSolutes(); index++)
        {
            const Coordinate coordinate = groParser.getCoordinateForSolute(index);
            const string type = groParser.getTypeForSolute(index);
            densityCubeFile << nameToZ[groParser.getTypeForSolute(index)] << " 0.0 " << (coordinate.x)*nmToBohr << " " << (coordinate.y)*nmToBohr << " " << (coordinate.z)*nmToBohr << endl;
        }

        for(int nx=0; nx<nXGrid; nx++) {
            for(int ny=0; ny<nYGrid; ny++) {
                for(int nz=0; nz<nZGrid; nz++) {
                    densityCubeFile << grid[nx][ny][nz]<< endl;
                }
            }
        }

        densityCubeFile.close();
    }
}









