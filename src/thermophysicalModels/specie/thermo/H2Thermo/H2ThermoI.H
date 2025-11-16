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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::H2Thermo<EquationOfState>::H2Thermo
(
    const EquationOfState& st
)
:
    EquationOfState(st)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::H2Thermo<EquationOfState>::H2Thermo
(
    const word& name,
    const H2Thermo& ct
)
:
    EquationOfState(name, ct)
{}


template<class EquationOfState>
inline Foam::autoPtr<Foam::H2Thermo<EquationOfState>>
Foam::H2Thermo<EquationOfState>::clone() const
{
    return autoPtr<H2Thermo<EquationOfState>>::New(*this);
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::H2Thermo<EquationOfState>>
Foam::H2Thermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<H2Thermo<EquationOfState>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::H2Thermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::H2Thermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    return EquationOfState::Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::H2Thermo<EquationOfState>::Ha
(
    const scalar p, const scalar T
) const
{
    return Hs(p, T) + Hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::H2Thermo<EquationOfState>::Hs
(
    const scalar p, const scalar T
) const
{
    return EquationOfState::H(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::H2Thermo<EquationOfState>::Hc() const
{
    return 0.0;
}


template<class EquationOfState>
inline Foam::scalar Foam::H2Thermo<EquationOfState>::S
(
    const scalar p, const scalar T
) const
{
    return EquationOfState::S(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::H2Thermo<EquationOfState>::Gstd
(
    const scalar T
) const
{
    return EquationOfState::H(Pstd, T) - T*EquationOfState::S(Pstd, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::H2Thermo<EquationOfState>::dCpdT
(
    const scalar p, const scalar T
) const
{
    return EquationOfState::dCpdT(p, T);;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::H2Thermo<EquationOfState>::operator+=
(
    const H2Thermo<EquationOfState>& ct
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(ct);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        scalar Y2 = ct.Y()/this->Y();
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::H2Thermo<EquationOfState> Foam::operator+
(
    const H2Thermo<EquationOfState>& ct1,
    const H2Thermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
      + static_cast<const EquationOfState&>(ct2)
    );

    if (mag(eofs.Y()) < SMALL)
    {
        return H2Thermo<EquationOfState>
        (
            eofs
        );
    }
    else
    {
        return H2Thermo<EquationOfState>
        (
            eofs
        );
    }
}


template<class EquationOfState>
inline Foam::H2Thermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const H2Thermo<EquationOfState>& ct
)
{
    return H2Thermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ct)
    );
}


template<class EquationOfState>
inline Foam::H2Thermo<EquationOfState> Foam::operator==
(
    const H2Thermo<EquationOfState>& ct1,
    const H2Thermo<EquationOfState>& ct2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ct1)
     == static_cast<const EquationOfState&>(ct2)
    );

    return H2Thermo<EquationOfState>
    (
        eofs
    );
}


// ************************************************************************* //
