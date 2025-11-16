/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::functionObjects::MachNoRealH2

Group
    grpFieldFunctionObjects

Description
    Computes the Mach number as a \c volScalarField.

    Operands:
    \table
      Operand          | Type           | Location
      input            | volVectorField | $FOAM_CASE/\<time\>/\<inpField\>
      output file      | -              | -
      output field     | volScalarField | $FOAM_CASE/\<time\>/\<outField\>
    \endtable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    MachNumberRealH2
    {
        // Mandatory entries (unmodifiable)
        type            MachNoRealH2;
        libs            ("libfieldmyFunctionObjects.so");

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Req'd | Dflt
      type         | Type name: MachNoRealH2                  | word |  yes  | -
      libs         | Library name: libfieldmyFunctionObjects | word |  yes  | -
    \endtable

    The inherited entries are elaborated in:
     - \link functionObject.H \endlink
     - \link fieldExpression.H \endlink


SourceFiles
    MachNoRealH2.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_MachNoRealH2_H
#define functionObjects_MachNoRealH2_H

#include "fieldExpression.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class MachNoRealH2 Declaration
\*---------------------------------------------------------------------------*/

class MachNoRealH2
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the Mach number field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("MachNoRealH2");


    // Constructors

        //- Construct for given objectRegistry and dictionary.
        //  Allow the possibility to load fields from files
        MachNoRealH2
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- No copy construct
        MachNoRealH2(const MachNoRealH2&) = delete;

        //- No copy assignment
        void operator=(const MachNoRealH2&) = delete;


    //- Destructor
    virtual ~MachNoRealH2() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
