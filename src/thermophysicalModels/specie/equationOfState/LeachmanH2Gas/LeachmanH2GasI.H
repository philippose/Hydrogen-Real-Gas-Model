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
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::LeachmanH2Gas<Specie>::LeachmanH2Gas
(
    const Specie& sp,
    const scalar& Tc,
    const scalar& Vc,
    const scalar& Pc
)
:
    Specie(sp),
    Tc_(Tc),
    Vc_(Vc),
    Pc_(Pc)
{

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::LeachmanH2Gas<Specie>::LeachmanH2Gas
(
    const word& name,
    const LeachmanH2Gas& pg
)
:
    Specie(name, pg),
    Tc_(pg.Tc_),
    Vc_(pg.Vc_),
    Pc_(pg.Pc_)
{

}


template<class Specie>
inline Foam::autoPtr<Foam::LeachmanH2Gas <Specie>>
Foam::LeachmanH2Gas<Specie>::clone() const
{
    return autoPtr<LeachmanH2Gas<Specie>>::New(*this);
}


template<class Specie>
inline Foam::autoPtr<Foam::LeachmanH2Gas<Specie>>
Foam::LeachmanH2Gas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<LeachmanH2Gas<Specie>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::LeachmanH2Gas<Specie>::rho(scalar p, scalar T) const
{
    const scalar M = RR/this->R()/1000.0;
    const scalar Pi = p/Pc_;
    const scalar Tau = Tc_/T;

    scalar A = (3.29664/Pi + -0.000296503/sqrt(Pi) + -0.00271293 + -0.000307942*sqrt(Pi) + 0.000108365*Pi + -5.78755e-07*sqr(Pi))/Tau;
    scalar B = -0.000299925/Pi + 0.00480231/sqrt(Pi) + 0.352393 + 0.0107914*sqrt(Pi) + -0.00378648*Pi + 1.87303e-05*sqr(Pi);
    scalar C = (0.0145274/Pi + -0.200007/sqrt(Pi) + -0.522984 + -0.252758*sqrt(Pi) + 0.0507639*Pi + -0.000222019*sqr(Pi))*Tau;
    scalar D = (-0.239496/Pi + 3.13049/sqrt(Pi) + -5.52193 + 2.39076*sqrt(Pi) + -0.262458*Pi + 0.000805104*sqr(Pi))*sqr(Tau);
    return 1/(A+B+C+D)*1000*M/Vc_;
}

template<class Specie>
inline Foam::scalar Foam::LeachmanH2Gas<Specie>::H(scalar p, scalar T) const
{
    const scalar M = RR/this->R()/1000.0;
    const scalar Delta = Vc_*this->rho(p,T)/1000.0/M;
    const scalar Tau = Tc_/T;

    scalar A = (303.06192241/Tau + 2036.48916359/sqrt(Tau) + -8773.8794081 + 14130.92327428*sqrt(Tau) + -7559.95968243*Tau);
    scalar B = (63.6698282/Tau + -698.57300897/sqrt(Tau) + 2816.01869592 + -4920.31164185*sqrt(Tau) + 3153.77861237*Tau)*sqrt(Delta);
    scalar C = (-135.34402051/Tau + 1976.27750815/sqrt(Tau) + -7979.24744228 + 13042.77442704*sqrt(Tau) + -8357.12295415*Tau)*Delta;
    scalar D = (199.07836741/Tau + -2199.14125078/sqrt(Tau) + 9306.93434938 + -16824.36618175*sqrt(Tau) + 11169.70297796*Tau)*sqr(Delta);
    scalar E = (-124.08489154/Tau + 1445.30765548/sqrt(Tau) + -6065.70448805 + 11069.72708383*sqrt(Tau) + -7467.61831306*Tau)*sqr(Delta)*Delta;
    scalar F = (19.31080942/Tau + -245.91556346/sqrt(Tau) + 1107.87606047 + -2095.74171864*sqrt(Tau) + 1463.99140087*Tau)*sqr(Delta)*sqr(Delta);

    return 1000.0*(A+B+C+D+E+F); 
}

template<class Specie>
inline Foam::scalar Foam::LeachmanH2Gas<Specie>::Cp(scalar p, scalar T) const
{
    const scalar M = RR/this->R()/1000.0;
    const scalar Delta = Vc_*this->rho(p,T)/1000.0/M;
    const scalar Tau = Tc_/T;

    scalar A = (2.51784942e+00/Tau + -3.02199542e+01/sqrt(Tau) + 1.45586430e+02 + -2.37269737e+02*sqrt(Tau) + 1.42263678e+02*Tau);
    scalar B = (-1.09265065e+00/Tau + 1.31451913e+01/sqrt(Tau) + -5.90151640e+01 + 1.17076401e+02*sqrt(Tau) + -8.66693462e+01*Tau)*sqrt(Delta);
    scalar C = (2.77612445e+00/Tau + -3.36660433e+01/sqrt(Tau) + 1.52107624e+02 + -3.05173134e+02*sqrt(Tau) + 2.41113270e+02*Tau)*Delta;
    scalar D = (-2.52109955e+00/Tau + 3.06412477e+01/sqrt(Tau) + -1.38073553e+02 + 2.72938596e+02*sqrt(Tau) + -2.04863672e+02*Tau)*sqr(Delta);
    scalar E = (1.08594055e+00/Tau + -1.34328537e+01/sqrt(Tau) + 6.10807324e+01 + -1.19535356e+02*sqrt(Tau) + 8.50322317e+01*Tau)*sqr(Delta)*Delta;
    scalar F = (-1.41852995e-01/Tau + 1.87578589e+00/sqrt(Tau) + -8.99491325e+00 + 1.81843744e+01*sqrt(Tau) + -1.29059115e+01*Tau)*sqr(Delta)*sqr(Delta);

    return 1000.0*(A+B+C+D+E+F);
}

template<class Specie>
inline Foam::scalar Foam::LeachmanH2Gas<Specie>::E(scalar p, scalar T) const
{
    const scalar M = RR/this->R()/1000.0;
    const scalar Delta = Vc_*this->rho(p,T)/1000.0/M;
    const scalar Tau = Tc_/T;

    scalar A = (166.32634318/Tau + 2036.94035344/sqrt(Tau) + -8776.93781107 + 14140.60135103*sqrt(Tau) + -7570.7776267*Tau);
    scalar B = (64.03991777/Tau + -707.7740136/sqrt(Tau) + 2898.32759136 + -5213.93247713*sqrt(Tau) + 3491.41807805*Tau)*sqrt(Delta);
    scalar C = (-143.10545613/Tau + 1703.2612291/sqrt(Tau) + -7202.27250717 + 12616.04727441*sqrt(Tau) + -8484.05993196*Tau)*Delta;
    scalar D = (190.09186338/Tau + -2106.18428532/sqrt(Tau) + 8712.86905844 + -15822.51142415*sqrt(Tau) + 10708.84424757*Tau)*sqr(Delta);
    scalar E = (-118.74396454/Tau + 1345.9416812/sqrt(Tau) + -5593.07509049 + 10180.84983044*sqrt(Tau) + -6902.68557014*Tau)*sqr(Delta)*Delta;
    scalar F = (20.72132864/Tau + -245.43544098/sqrt(Tau) + 1054.68414088 + -1960.79495233*sqrt(Tau) + 1351.97794081*Tau)*sqr(Delta)*sqr(Delta);

    return 1000.0*(A+B+C+D+E+F);
}

template<class Specie>
inline Foam::scalar Foam::LeachmanH2Gas<Specie>::Cv(scalar p, scalar T) const
{
    const scalar M = RR/this->R()/1000.0;
    const scalar Delta = Vc_*this->rho(p,T)/1000.0/M;
    const scalar Tau = Tc_/T;

    scalar A = (2.48433355e+00/Tau + -2.98138623e+01/sqrt(Tau) + 1.39625614e+02 + -2.33601612e+02*sqrt(Tau) + 1.39534379e+02*Tau);
    scalar B = (-1.05737915e-01/Tau + 1.19805013e+00/sqrt(Tau) + -5.06874692e+00 + 9.50235450e+00*sqrt(Tau) + -6.67922407e+00*Tau)*sqrt(Delta);
    scalar C = (3.08025379e-01/Tau + -3.60605176e+00/sqrt(Tau) + 1.61548529e+01 + -3.10082627e+01*sqrt(Tau) + 2.29383702e+01*Tau)*Delta;
    scalar D = (-3.95008404e-01/Tau + 4.52838640e+00/sqrt(Tau) + -1.95632585e+01 + 3.82857951e+01*sqrt(Tau) + -2.85627895e+01*Tau)*sqr(Delta);
    scalar E = (2.58875704e-01/Tau + -2.91702883e+00/sqrt(Tau) + 1.22902914e+01 + -2.30972488e+01*sqrt(Tau) + 1.66438796e+01*Tau)*sqr(Delta)*Delta;
    scalar F = (-5.84333431e-02/Tau + 6.54612018e-01/sqrt(Tau) + -2.73146861e+00 + 5.05992916e+00*sqrt(Tau) + -3.56414701e+00*Tau)*sqr(Delta)*sqr(Delta);

    return 1000.0*(A+B+C+D+E+F);
}

template<class Specie>
inline Foam::scalar Foam::LeachmanH2Gas<Specie>::S(scalar p, scalar T) const
{
    const scalar M = RR/this->R()/1000.0;
    const scalar Delta = Vc_*this->rho(p,T)/1000.0/M;
    const scalar Tau = Tc_/T;

    scalar A = (-4.63085946e-07/Tau + 4.38983303e-06/sqrt(Tau) + -1.04136783e-05)/sqr(Delta);
    scalar B = (3.17938414e-05/Tau + -4.03692882e-04/sqrt(Tau) + 1.18029707e-03)/sqrt(Delta)/Delta;
    scalar C = (-3.67985180e-04/Tau + 1.23038756e-02/sqrt(Tau) + -5.22620140e-02)/Delta;
    scalar D = (-1.07799343e-02/Tau + -1.48187013e-01/sqrt(Tau) + 1.35461475e+00)/sqrt(Delta);
    scalar E = (-5.67119427e-01/Tau + 1.22574647e+01/sqrt(Tau) + 1.31791433e+01);
    scalar F = (-1.54770334e+00/Tau + -1.22244271e-01/sqrt(Tau) + -2.47217215e+01)*sqrt(Delta);
    scalar G = (4.21846210e+00/Tau + -5.09724931e+00/sqrt(Tau) + 1.90368593e+01)*Delta;
    scalar H = (-5.07759701e+00/Tau + 1.12461342e+01/sqrt(Tau) + -1.63655648e+01)*sqrt(Delta)*Delta;
    scalar I = (2.11854559e+00/Tau + -5.94802690e+00/sqrt(Tau) + 6.21796592e+00)*sqr(Delta);

    return 1000.0*(A+B+C+D+E+F+G+H+I);
}

template<class Specie>
inline Foam::scalar Foam::LeachmanH2Gas<Specie>::psi(scalar p, scalar T) const
{
    const scalar Z = this->Z(p, T);
    return 1.0/(Z*this->R()*T);
}

template<class Specie>
inline Foam::scalar Foam::LeachmanH2Gas<Specie>::Z(scalar p, scalar T) const
{
    return p/(this->rho(p,T)*this->R()*T);
}

template<class Specie>
inline Foam::scalar Foam::LeachmanH2Gas<Specie>::CpMCv(scalar p, scalar T) const
{
    return this->Cp(p,T) - this->Cv(p,T);
}

template<class Specie>
inline Foam::scalar Foam::LeachmanH2Gas<Specie>::dCpdT(scalar p, scalar T) const
{
    const scalar M = RR/this->R()/1000.0;
    const scalar Delta = Vc_*this->rho(p,T)/1000.0/M;
    const scalar Tau = Tc_/T;

    scalar A = (-3.19188449e-02/Tau + 4.04656403e-01/sqrt(Tau) + -1.86551235e+00 + 3.63554981e+00*sqrt(Tau) + -2.42520661e+00*Tau);
    scalar B = (4.08230761e-02/Tau + -4.80087240e-01/sqrt(Tau) + 2.10461470e+00 + -4.07813749e+00*sqrt(Tau) + 2.94990439e+00*Tau)*sqrt(Delta);
    scalar C = (-1.06154926e-01/Tau + 1.26574989e+00/sqrt(Tau) + -5.64339853e+00 + 1.11896145e+01*sqrt(Tau) + -8.40244236e+00*Tau)*Delta;
    scalar D = (8.69219842e-02/Tau + -1.03697673e+00/sqrt(Tau) + 4.58874814e+00 + -8.93738980e+00*sqrt(Tau) + 6.49448415e+00*Tau)*sqr(Delta);
    scalar E = (-3.04006648e-02/Tau + 3.74628506e-01/sqrt(Tau) + -1.68018962e+00 + 3.24593259e+00*sqrt(Tau) + -2.27179212e+00*Tau)*sqr(Delta)*Delta;
    scalar F = (1.28075966e-03/Tau + -2.48025207e-02/sqrt(Tau) + 1.39261388e-01 + -3.02580462e-01*sqrt(Tau) + 2.20447333e-01*Tau)*sqr(Delta)*sqr(Delta);

    return 1000.0*(A+B+C+D+E+F);
}

template<class Specie>
inline Foam::scalar Foam::LeachmanH2Gas<Specie>::w(scalar p, scalar T) const
{
    const scalar M = RR/this->R()/1000.0;
    const scalar Delta = Vc_*this->rho(p,T)/1000.0/M;
    const scalar Tau = Tc_/T;

    scalar A = (-2.71181100e+01/Tau + 7.76265357e+02/sqrt(Tau) + -1.52238557e+03 + 2.81681226e+03*sqrt(Tau) + -1.67759981e+03*Tau);
    scalar B = (1.32644598e+00/Tau + -1.72586937e+01/sqrt(Tau) + 4.33246438e+01 + 1.09521136e+02*sqrt(Tau) + -3.58363795e+02*Tau)*sqrt(Delta);
    scalar C = (-2.62136725e+01/Tau + 3.60252864e+02/sqrt(Tau) + -7.38311941e+02 + 6.20798522e+02*sqrt(Tau) + 8.08971987e+01*Tau)*Delta;
    scalar D = (1.24474593e+01/Tau + -2.20683480e+02/sqrt(Tau) + 1.20546492e+03 + -1.84072792e+03*sqrt(Tau) + 6.56727769e+02*Tau)*sqr(Delta);
    scalar E = (5.65468710e+00/Tau + 1.45713696e+01/sqrt(Tau) + -3.55818262e+02 + 9.45803639e+02*sqrt(Tau) + -4.29181971e+02*Tau)*sqr(Delta)*Delta;
    scalar F = (-6.81222200e+00/Tau + 5.41649906e+01/sqrt(Tau) + -1.26559998e+02 + 9.53260650e+01*sqrt(Tau) + -5.12543488e+01*Tau)*sqr(Delta)*sqr(Delta);

    return (A+B+C+D+E+F);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::LeachmanH2Gas<Specie>::operator+=
(
    const LeachmanH2Gas<Specie>& pg
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(pg);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = pg.Y()/this->Y();

        Tc_ = Y1*Tc_ + Y2*pg.Tc_;
        Vc_ = Y1*Vc_ + Y2*pg.Vc_;
        Pc_ = Y1*Pc_ + Y2*pg.Pc_;
    }
}


template<class Specie>
inline void Foam::LeachmanH2Gas<Specie>::operator*=(const scalar s)
{
     Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Specie>
Foam::LeachmanH2Gas<Specie> Foam::operator+
(
    const LeachmanH2Gas<Specie>& pg1,
    const LeachmanH2Gas<Specie>& pg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pg1)
      + static_cast<const Specie&>(pg2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return LeachmanH2Gas<Specie>
        (
            sp,
            pg1.Tc_,
            pg1.Vc_,
            pg1.Pc_
        );
    }
    else
    {
        const scalar Y1 = pg1.Y()/sp.Y();
        const scalar Y2 = pg2.Y()/sp.Y();

        const scalar Tc = Y1*pg1.Tc_ + Y2*pg2.Tc_;
        const scalar Vc = Y1*pg1.Vc_ + Y2*pg2.Vc_;
        const scalar Pc = Y1*pg1.Pc_ + Y2*pg2.Pc_;


        return LeachmanH2Gas<Specie>
        (
            sp,
            Tc,
            Vc,
            Pc
        );
    }
}

template<class Specie>
Foam::LeachmanH2Gas<Specie> Foam::operator*
(
    const scalar s,
    const LeachmanH2Gas<Specie>& pg
)
{
    return LeachmanH2Gas<Specie>
    (
        s*static_cast<const Specie&>(pg),
        pg.Tc_,
        pg.Vc_,
        pg.Pc_
    );
}

template<class Specie>
Foam::LeachmanH2Gas<Specie> Foam::operator==
(
    const LeachmanH2Gas<Specie>& pg1,
    const LeachmanH2Gas<Specie>& pg2
)
{
    Specie sp
    (
        static_cast<const Specie&>(pg1)
     == static_cast<const Specie&>(pg2)
    );

    const scalar Y1 = pg1.Y()/sp.Y();
    const scalar Y2 = pg2.Y()/sp.Y();

    const scalar Tc = Y2*pg2.Tc_ - Y1*pg1.Tc_;
    const scalar Vc = Y2*pg2.Vc_ - Y1*pg1.Vc_;
    const scalar Pc = Y2*pg2.Pc_ - Y1*pg1.Pc_;

    return LeachmanH2Gas<Specie>
    (
        sp,
        Tc,
        Vc,
        Pc
    );
}


// ************************************************************************* //
