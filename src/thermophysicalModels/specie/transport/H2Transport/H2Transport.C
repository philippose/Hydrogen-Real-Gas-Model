/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "H2Transport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::H2Transport<Thermo>::H2Transport(const dictionary& dict)
:
    Thermo(dict),
    Tc_(dict.subDict("equationOfState").get<scalar>("Tc")),
    Vc_(dict.subDict("equationOfState").get<scalar>("Vc"))
{}

template<class Thermo>
Foam::H2Transport<Thermo>::H2Transport(const Thermo& t, const dictionary& dict)
:
    Thermo(t)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::H2Transport<Thermo>::write(Ostream& os) const
{
    os.beginBlock(this->specie::name());

    Thermo::write(os);

    // Entries in dictionary format
    {
        os.beginBlock("transport");
        os.endBlock();
    }

    os.endBlock();
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const H2Transport<Thermo>& st
)
{
    st.write(os);
    return os;
}

// ************************************************************************* //
