/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    POD

Description
    Compute the SVD 

SourceFiles
    POD.C

\*---------------------------------------------------------------------------*/

#ifndef UTILITIES_H
#define UTILITIES_H

#include "fvCFD.H"
#include <Eigen/Dense>

namespace Foam
{
class utilities
{
    public:
        utilities() = default;                         // Constructor
        utilities(const utilities&) = delete;          // no copy constructor
        void operator=(const utilities&) = delete;     // no assignment
        virtual ~utilities() = default;                // default destructor


        /****************************************/
        /*********** MEMBER FUNCTIONS ***********/
        /****************************************/
        
        void debugMsg(const string& msg, const Switch dbg)
        {
            if(dbg)
            {
                std::cout << "**____DEBUG____** --> " << msg << nl
                           << "*****************" << nl;
            }
        }

        template<typename... Args>
        void debugMsg(const Switch dbg, Args&&... args)
        {
            if (dbg)
            {
                std::ostringstream os;
                (os << ... << args);   // fold expression (C++17)

                std::cout << "**____DEBUG____** --> " << os.str() << nl
                        << "*****************" << nl;
            }
        }


        void plotVector(const Eigen::VectorXd& vec, word xLabel, word yLabel, word title)
        {
            FILE* gp = popen("gnuplot -persistent", "w");
            
            if (!gp)
            {
                Foam::Warning << "Cannot open gnuplot pipe\n";
                return;
            }
            
            fprintf(gp,
                    // Terminal and window setup
                    "set terminal qt size 800,600 enhanced font 'Arial,12' title 'Vector Data Visualization'\n"
                    
                    // Margins and layout
                    "set bmargin 7\n"
                    "set lmargin 15\n"
                    "set rmargin 5\n"
                    "set tmargin 5\n"
                    
                    // Labels with better positioning
                    "set xlabel '%s' font 'Arial,14' offset 0,-1\n"
                    //"set ylabel 'Sigma' font 'Arial,14' offset -2,0\n"
                    "set ylabel '%s' font 'Arial,14'\n"
                    "set title '%s' font 'Arial,16' offset 0,1\n"
                    
                    // Grid styling
                    "set grid xtics ytics\n"
                    "set grid linewidth 1 linetype 0 linecolor rgb '#cccccc'\n"
                    
                    // Axis styling
                    "set border linewidth 2\n"
                    "set xtics font 'Arial,11'\n"
                    "set ytics font 'Arial,11'\n"
                    
                    // Set y axis to logarithmic scale
                    "set logscale y\n"

                    // Auto-scaling with some padding
                    "set autoscale\n"
                    //"set yrange [0:10000]\n"
                    //"set offsets 0.05, 0.05, 0.05, 0.05\n"  // Add 5% padding around data
                    
                    // Color and style
                    "set style line 1 linecolor rgb '#2E86AB' pointtype 7 pointsize 1.2 linewidth 2\n"
                    
                    // Plot command - choose based on data characteristics
                    "plot '-' with points pt 7 ps 1 title 'data'\n",
                    xLabel.c_str(),
                    yLabel.c_str(),
                    title.c_str()
                );
                // Stream the data: “index value” per line
                for (label i = 0; i < vec.size(); ++i)
                {
                    fprintf(gp, "%d %g\n", i, vec[i]);
                }

            // Signal end of data
            fprintf(gp, "e\n");
            fflush(gp);
            pclose(gp);
        }

        //helper function that creates a list of timefolders without the constant folder
        instantList getTimeDirs(const Time& runTime, bool excludeZero = false, bool debug = false)
        {
            instantList allTimes = runTime.times();
            instantList timeDirs;

            for (const instant& t : allTimes)
            {
                if (t.name() != "constant" && (!excludeZero || t.value() > 0))
                {
                    timeDirs.append(t);
                }
            }
            debugMsg("this is the first timeDirectory without constant and 0: " + string(timeDirs[0].name()) + " and this is the dimension: " + name(timeDirs.size()), debug);
            
            return timeDirs;
        }

        //computes dimension when giving first entry, last entry and separation of a wanted vector
        label computeDim(scalar firstTime, scalar lastTime, scalar separation, bool debug = false)
        {
            if (separation <= SMALL)
            {
                FatalErrorInFunction << "Separation must be > 0" << exit(FatalError);
            }

            if (lastTime < firstTime)
            {
                FatalErrorInFunction << "lastTime must be >= firstTime" << exit(FatalError);
            }

            label dim = std::floor((lastTime - firstTime)/separation) + 1;

            debugMsg("Computed dimension: " + name(dim) + " (from " + name(firstTime) + " to "   + name(lastTime) + " with separation " + name(separation) + ")", debug);

            return dim;
        }

        instantList createCustomTime(const Time& runTime, scalar firstTime, scalar lastTime, scalar separation, bool debug = false)
        {
            label dim = computeDim(firstTime, lastTime, separation, debug);
            instantList newTimes(dim);

            //fill the list
            for (label i = 0; i < dim; ++i)
            {
                scalar t = firstTime + i * separation;

                // Construct an instant with time value and name string
                newTimes[i] = instant(t, name(t));
            }
            debugMsg("created " + name(newTimes.size()) + " instants from " + name(firstTime) + " to " + name(lastTime) + " with separation " + name(separation), debug);
            
            return newTimes;
        }

    }; 
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
