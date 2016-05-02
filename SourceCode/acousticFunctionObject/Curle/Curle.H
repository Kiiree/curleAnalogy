/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::Curle

Description
    Implementation of the Curle analogy for incompressible flow as
    described in:
    \verbatim
        Johan Larsson, Lars Davidson, Magnus Olsson, Lars-Erik Eriksson (2004).
        Aeroacoustic Investigation of an Open Cavity at Low Mach Number.
        AIAA Journal, Vol. 42, No. 12.
    \endverbatim

    Example of function object specification:
    \verbatim
    Curle1
    {
        type                    Curle;
        functionObjectLibs      ("libAcousticAnalogy.so");
        log                     true;

        patches                 ( cylinder1 cylinder2 );
        cellZone                cellZoneName;
        rhoRef                  1.225;
        cRef                    343;

        observers
        {
            microphone-A
            {
                position    	(-0.4760595 1.58962725 0);
            }
        }
    }
    \endverbatim

    \heading Function object usage
    \table
        Property     | Description                                       | Required    | Default value
        type         | type name: Curle                                  | yes         |
        log          | write acoustic data to standard output            | no          | no
        patches      | patches included in the acoustic calculation      | yes         |
        cellZone     | name of cell zone included in acoustic calculaion | no          |
        rhoRef       | reference density                                 | yes         |
        cRef         | reference speed of sound                          | yes         |
        pName        | pressure field name                               | no          | p
        UName        | velocity field name                               | no          | U
    \endtable


\*---------------------------------------------------------------------------*/

#ifndef Curle_H
#define Curle_H

#include "functionObjectFile.H"
#include "volFieldsFwd.H"
#include "Switch.H"
#include "soundObserver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class objectRegistry;
class dictionary;
class polyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class Curle Declaration
\*---------------------------------------------------------------------------*/

class Curle
:
    public functionObjectFile
{
protected:

    // Protected data

        //- Name of this set of Curle,
        //  Also used as the name of the probes directory.
        word name_;

        const objectRegistry& obr_;

        //- On/off switch
        bool active_;

        //- Initialised flag
        bool initialised_;
        
        //- Switch for compressible
        bool compressible_;

        //- Switch to send output to Info as well as to file
        Switch log_;
        
        //- Patches to sample
        labelHashSet patches_;
        
        //- Name of cellZone
        word cellZoneName_;
        
        //- ID of cellZone
        label cellZoneID_;
        
        //- Name of pressure field (optional)
        word pName_;
        
        //- Name of density field (optional)
        word rhoName_;
        
        //- Name of velocity field (optional)
        word UName_;

        //- Reference density needed for incompressible cases
        scalar rhoRef_;
        
        //- Reference speed of sound
        scalar cRef_;

        //- List of observers
        List<SoundObserver> observers_;
        
    // Protected Member Functions

        //- Create file names for Curle
        wordList createFileNames(const dictionary& dict) const;

        //- Output file header information
        virtual void writeFileHeader(const label i);

        //- Helper function to write acoustic data
        void writeCurle();

        //- Initialise the fields
        void initialise();

        //- Disallow default bitwise copy construct
        Curle(const Curle&);

        //- Disallow default bitwise assignment
        void operator=(const Curle&);

        //- Return pressure
        tmp<volScalarField> p() const;

        //- Return pressure time derivative
        tmp<volScalarField> dpdt() const;

        //- Return Lighthill tensor
        tmp<volTensorField> Tij() const;

        //- Return first time derivative of Lighthill tensor 
        tmp<volTensorField> dTijdt() const;

        //- Return second time derivative of Lighthill tensor 
        tmp<volTensorField> d2Tijdt2() const;


public:

    //- Runtime type information
    TypeName("Curle");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        Curle
        (
            const word& name,
            const objectRegistry&,
            const dictionary&,
            const bool loadFromFiles = false,
            const bool readFields = true
        );


    //- Destructor
    virtual ~Curle();


    // Member Functions

        //- Return name of the set of Curle
        virtual const word& name() const
        {
            return name_;
        }

        //- Read the Curle data
        virtual void read(const dictionary&);

        //- Execute, currently does nothing
        virtual void execute();

        //- Execute at the final time-loop, currently does nothing
        virtual void end();

        //- Called when time was set at the end of the Time::operator++
        virtual void timeSet();

        //- Write the acoustic data
        virtual void write();

        //- Calculate the acoustic data
        virtual void calculate();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh&)
        {}

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh&)
        {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
