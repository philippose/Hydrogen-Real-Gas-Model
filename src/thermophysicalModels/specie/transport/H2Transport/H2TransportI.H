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

#include "specie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::H2Transport<Thermo>::H2Transport
(
    const Thermo& t,
    const scalar& Tc,
    const scalar& Vc
)
:
    Thermo(t),
    Tc_(Tc),
    Vc_(Vc)
{}

template<class Thermo>
inline Foam::H2Transport<Thermo>::H2Transport
(
    const word& name,
    const H2Transport& st
)
:
    Thermo(name, st),
    Tc_(st.Tc_),
    Vc_(st.Vc_)
{}

template<class Thermo>
inline Foam::autoPtr<Foam::H2Transport<Thermo>>
Foam::H2Transport<Thermo>::clone() const
{
    return autoPtr<H2Transport<Thermo>>::New(*this);
}

template<class Thermo>
inline Foam::autoPtr<Foam::H2Transport<Thermo>>
Foam::H2Transport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<H2Transport<Thermo>>::New(dict);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::H2Transport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    const scalar rho = this->rho(p,T);
    const scalar M = RR/this->R()/1000.0;
    const scalar Tau = Tc_/T;
    const scalar Delta = rho/M/1000.0*Vc_;

    scalar A = 7.74547668e-01 + 1.29238538e+00/Tau + -6.84212063e-02/sqr(Tau) + 3.48747415e-03/sqr(Tau)/Tau + -7.67422945e-05/sqr(Tau)/sqr(Tau); // 0
    scalar B = Delta*(-1.46086756e+04*sqr(Tau)*sqr(Tau)*sqr(Tau) + 1.66933564e+04*sqr(Tau)*sqr(Tau)*Tau + -7.65723295e+03*sqr(Tau)*sqr(Tau) + 1.56405852e+03*sqr(Tau)*Tau + -1.67026246e+00*sqr(Tau) + -7.66928573e+01*Tau + 2.09739133e+01 + -2.83918911e+00/Tau + 2.21891496e-01/sqr(Tau) + -9.69992933e-03/sqr(Tau)/Tau + 1.81478140e-04/sqr(Tau)/sqr(Tau)); // 1
    scalar C = sqr(Delta)*(8.28934333e-01 + 2.65215204e-02/Tau + 1.36337564e-03/sqr(Tau)); // 2
    scalar D = sqr(Delta)*sqr(Delta)*(1.58943616e-01*Tau + -6.42377129e-03 + 3.02064008e-04/Tau); // 4
    scalar E = sqr(Delta)*sqr(Delta)*sqr(Delta)*(5.97459453e-02*sqr(Tau) + -4.25572314e-02*Tau + 4.36489790e-03); // 6
    scalar F = sqr(Delta)*sqr(Delta)*sqr(Delta)*sqr(Delta)*(4.56887436e-02*sqr(Tau)*Tau + -3.15245364e-02*sqr(Tau) + 9.20429561e-03*Tau); // 8 

    return (A+B+C+D+E+F)*1.0e-06;
}

template<class Thermo>
inline Foam::scalar Foam::H2Transport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    const scalar rho = this->rho(p,T);
    const scalar M = RR/this->R()/1000.0;
    const scalar Tr = T/Tc_;
    const scalar d = rho/M/1000.0*Vc_;
    
    // Zero-density limit
    scalar m_ = -3.40976e-1 + 4.58820e0*Tr + -1.45080e0*sqr(Tr) + 3.26394e-1*sqr(Tr)*Tr + 3.16939e-3*sqr(Tr)*sqr(Tr) + 1.90592e-4*sqr(Tr)*sqr(Tr)*Tr + -1.13900e-6*sqr(Tr)*sqr(Tr)*sqr(Tr);
    scalar n_ = 1.38497e2 + -2.21878e1*Tr + 4.57151e0*sqr(Tr) + 1.00000e0*sqr(Tr)*Tr;
    scalar kt0 = m_/n_;
    
    // Excess part   
    scalar kte = (3.63081e-2+1.83370e-3*Tr)*d + (-2.07629e-2+ -8.86716e-3*Tr)*sqr(d) + (3.14810e-2+1.58260e-2*Tr)*sqr(d)*d + (-1.43097e-2+ -1.06283e-2*Tr)*sqr(d)*sqr(d) + (1.74980e-3+2.80673e-3*Tr)*sqr(d)*sqr(d)*d;

    return kt0 + kte;
}

template<class Thermo>
inline Foam::scalar Foam::H2Transport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{

    return kappa(p, T)/this->Cp(p, T);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::H2Transport<Thermo>::operator+=
(
    const H2Transport<Thermo>& st
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(st);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        scalar Y2 = st.Y()/this->Y();

        Tc_ = Y1*Tc_ + Y2*st.Tc_;
        Vc_ = Y1*Vc_ + Y2*st.Vc_;
    }
}

template<class Thermo>
inline void Foam::H2Transport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::H2Transport<Thermo> Foam::operator+
(
    const H2Transport<Thermo>& st1,
    const H2Transport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );

    if (mag(t.Y()) < SMALL)
    {
        return H2Transport<Thermo>
        (
            t,
            0,
            st1.Tc_,
            st1.Vc_
        );
    }
    else
    {
        scalar Y1 = st1.Y()/t.Y();
        scalar Y2 = st2.Y()/t.Y();

        const scalar Tc = Y1*st1.Tc_ + Y2*st2.Tc_;
        const scalar Vc = Y1*st1.Vc_ + Y2*st2.Vc_;

        return H2Transport<Thermo>
        (
            t,
            Y1*st1.Tc_ + Y2*st2.Tc_,
            Y1*st1.Vc_ + Y2*st2.Vc_
        );
    }
}

template<class Thermo>
inline Foam::H2Transport<Thermo> Foam::operator*
(
    const scalar s,
    const H2Transport<Thermo>& st
)
{
    return H2Transport<Thermo>
    (
        s*static_cast<const Thermo&>(st),
        st.Tc_,
        st.Vc_
    );
}

// ************************************************************************* //