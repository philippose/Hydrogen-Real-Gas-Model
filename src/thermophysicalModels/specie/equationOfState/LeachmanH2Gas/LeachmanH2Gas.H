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

Class
    Foam::LeachmanH2Gas

Group
    grpSpecieEquationOfState

Description
    LeachmanH2Gas gas equation of state.

SourceFiles
    LeachmanH2GasI.H
    LeachmanH2Gas.C

\*---------------------------------------------------------------------------*/

#ifndef LeachmanH2Gas_H
#define LeachmanH2Gas_H

#include "autoPtr.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class LeachmanH2Gas;

template<class Specie>
inline LeachmanH2Gas<Specie> operator+
(
    const LeachmanH2Gas<Specie>&,
    const LeachmanH2Gas<Specie>&
);

template<class Specie>
inline LeachmanH2Gas<Specie> operator*
(
    const scalar,
    const LeachmanH2Gas<Specie>&
);

template<class Specie>
inline LeachmanH2Gas<Specie> operator==
(
    const LeachmanH2Gas<Specie>&,
    const LeachmanH2Gas<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const LeachmanH2Gas<Specie>&
);



/*---------------------------------------------------------------------------*\
                           Class LeachmanH2Gas Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class LeachmanH2Gas
:
    public Specie
{
public:

    // Public Data

private:

    // Private data

        //- Critical temperature [K]
        scalar Tc_;

        //- Critical volume [m^3/kmol]
        scalar Vc_;

        //- Critical pressure [Pa]
        scalar Pc_;

public:

    // Constructors

        //- Construct from components
        inline LeachmanH2Gas
        (
            const Specie& sp,
            const scalar& Tc,
            const scalar& Vc,
            const scalar& Pc
        );

        //- Construct from dictionary
        LeachmanH2Gas(const dictionary& dict);

        //- Construct as named copy
        inline LeachmanH2Gas(const word& name, const LeachmanH2Gas&);

        //- Construct and return a clone
        inline autoPtr<LeachmanH2Gas> clone() const;

        // Selector from dictionary
        inline static autoPtr<LeachmanH2Gas> New
        (
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "LeachmanH2Gas<" + word(Specie::typeName_()) + '>';
        }

        // Fundamental properties


            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = false;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K)]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return internal energy departure [J/kg]
            inline scalar E(const scalar p, const scalar T) const;

            //- Return Cv departure [J/(kg K)]
            inline scalar Cv(scalar p, scalar T) const;

            //- Return entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kg K)]
            inline scalar CpMCv(scalar p, scalar T) const;

            //- Temperature derivative of heat capacity at constant pressure [J/(kg K^2)]
            inline scalar dCpdT(const scalar p, const scalar T) const;

            //- Speed of sound [m/s]
            inline scalar w(const scalar p, const scalar T) const;

    // IO

            //- Write to Ostream
            void write(Ostream& os) const;

    // Member operators

        inline void operator+=(const LeachmanH2Gas&);
        inline void operator*=(const scalar);


    // Friend operators

        friend LeachmanH2Gas operator+ <Specie>
        (
            const LeachmanH2Gas&,
            const LeachmanH2Gas&
        );

        friend LeachmanH2Gas operator* <Specie>
        (
            const scalar s,
            const LeachmanH2Gas&
        );

        friend LeachmanH2Gas operator== <Specie>
        (
            const LeachmanH2Gas&,
            const LeachmanH2Gas&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const LeachmanH2Gas&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LeachmanH2GasI.H"

#ifdef NoRepository
    #include "LeachmanH2Gas.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
