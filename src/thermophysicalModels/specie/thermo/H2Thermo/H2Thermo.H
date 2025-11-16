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
    Foam::H2Thermo

Group
    grpSpecieThermo

Description


SourceFiles
    H2ThermoI.H
    H2Thermo.C

\*---------------------------------------------------------------------------*/

#ifndef H2Thermo_H
#define H2Thermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class H2Thermo;

template<class EquationOfState>
inline H2Thermo<EquationOfState> operator+
(
    const H2Thermo<EquationOfState>&,
    const H2Thermo<EquationOfState>&
);

template<class EquationOfState>
inline H2Thermo<EquationOfState> operator*
(
    const scalar,
    const H2Thermo<EquationOfState>&
);

template<class EquationOfState>
inline H2Thermo<EquationOfState> operator==
(
    const H2Thermo<EquationOfState>&,
    const H2Thermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const H2Thermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class H2Thermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class H2Thermo
:
    public EquationOfState
{

    // Private Member Functions

        //- Construct from components
        inline H2Thermo
        (
            const EquationOfState& st
        );

public:

    // Constructors

        //- Construct from dictionary
        H2Thermo(const dictionary& dict);

        //- Construct as named copy
        inline H2Thermo(const word&, const H2Thermo&);

        //- Construct and return a clone
        inline autoPtr<H2Thermo> clone() const;

        //- Selector from dictionary
        inline static autoPtr<H2Thermo> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "H2Thermo<" + EquationOfState::typeName() + '>';
        }

        //- Limit temperature to be within the range
        inline scalar limit(const scalar T) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/(kg K)]
            inline scalar Cp(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kg]
            inline scalar Ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar Hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            //- Entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            //- Gibbs free energy of the mixture in the standard state [J/kg]
            inline scalar Gstd(const scalar T) const;

            #include "HtoEthermo.H"


        // Derivative term used for Jacobian


            //- Temperature derivative of heat capacity at constant pressure
            inline scalar dCpdT(const scalar p, const scalar T) const;

        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const H2Thermo&);


    // Friend operators

        friend H2Thermo operator+ <EquationOfState>
        (
            const H2Thermo&,
            const H2Thermo&
        );

        friend H2Thermo operator* <EquationOfState>
        (
            const scalar,
            const H2Thermo&
        );

        friend H2Thermo operator== <EquationOfState>
        (
            const H2Thermo&,
            const H2Thermo&
        );


    // IOstream Operators

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const H2Thermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "H2ThermoI.H"

#ifdef NoRepository
    #include "H2Thermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
