/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::H2Transport

Group
    grpSpecieTransport

Description

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s]

SourceFiles
    H2TransportI.H
    H2Transport.C

\*---------------------------------------------------------------------------*/

#ifndef H2Transport_H
#define H2Transport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class Thermo> class H2Transport;

template<class Thermo>
inline H2Transport<Thermo> operator+
(
    const H2Transport<Thermo>&,
    const H2Transport<Thermo>&
);

template<class Thermo>
inline H2Transport<Thermo> operator*
(
    const scalar,
    const H2Transport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const H2Transport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class H2Transport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class H2Transport
:
    public Thermo
{
public:


private:

    // Private Data

        //- Critical temperature [K]
        scalar Tc_;

        //- Critical volume [m^3/kmol]
        scalar Vc_;


public:

    // Constructors

        //- Construct from components
        inline H2Transport
        (
            const Thermo& t,
            const scalar& Tc,
            const scalar& Vc
        );

        //- Construct as named copy
        inline H2Transport(const word&, const H2Transport&);

        //- Construct from dictionary
        explicit H2Transport(const dictionary& dict);

        //- Construct from base thermo and dictionary
        H2Transport(const Thermo& t,const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<H2Transport> clone() const;

        // Selector from dictionary
        inline static autoPtr<H2Transport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "H2Transport<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const H2Transport&);

        inline void operator*=(const scalar);


    // Friend Operators

        friend H2Transport operator+ <Thermo>
        (
            const H2Transport&,
            const H2Transport&
        );

        friend H2Transport operator* <Thermo>
        (
            const scalar,
            const H2Transport&
        );


    // IOstream Operators

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const H2Transport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "H2TransportI.H"

#ifdef NoRepository
    #include "H2Transport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
