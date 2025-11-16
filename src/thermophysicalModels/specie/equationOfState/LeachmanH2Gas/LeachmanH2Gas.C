/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "LeachmanH2Gas.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::LeachmanH2Gas<Specie>::LeachmanH2Gas
(
    const dictionary& dict
)
:
    Specie(dict),
    Tc_(dict.subDict("equationOfState").get<scalar>("Tc")),
    Vc_(dict.subDict("equationOfState").get<scalar>("Vc")),
    Pc_(dict.subDict("equationOfState").get<scalar>("Pc"))
{

}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::LeachmanH2Gas<Specie>::write(Ostream& os) const
{
    Specie::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const LeachmanH2Gas<Specie>& pg
)
{
    pg.write(os);
    return os;
}


// ************************************************************************* //
