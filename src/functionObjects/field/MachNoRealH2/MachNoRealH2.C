/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "MachNoRealH2.H"
#include "fluidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(MachNoRealH2, 0);
    addToRunTimeSelectionTable(functionObject, MachNoRealH2, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::MachNoRealH2::calc()
{
    if
    (
        foundObject<volVectorField>(fieldName_)
     && foundObject<fluidThermo>(fluidThermo::dictName)
    )
    {
        const fluidThermo& thermo =
            lookupObject<fluidThermo>(fluidThermo::dictName);

        const volVectorField& U = lookupObject<volVectorField>(fieldName_);

	const dimensionedScalar Vc_("Vc_",dimensionSet(0, 3, 0, 0, -1, 0 ,0), 0.064483);
	const dimensionedScalar Tc_("Tc_",dimensionSet(0, 0, 0, 1, 0, 0 ,0), 33.145);
    	const volScalarField Delta = Vc_*thermo.rho()/thermo.W();
    	const volScalarField Tau = Tc_/thermo.T();

    	const volScalarField A = (-2.71181100e+01/Tau + 7.76265357e+02/sqrt(Tau) + -1.52238557e+03 + 2.81681226e+03*sqrt(Tau) + -1.67759981e+03*Tau);
    	const volScalarField B = (1.32644598e+00/Tau + -1.72586937e+01/sqrt(Tau) + 4.33246438e+01 + 1.09521136e+02*sqrt(Tau) + -3.58363795e+02*Tau)*sqrt(Delta);
    	const volScalarField C = (-2.62136725e+01/Tau + 3.60252864e+02/sqrt(Tau) + -7.38311941e+02 + 6.20798522e+02*sqrt(Tau) + 8.08971987e+01*Tau)*Delta;
    	const volScalarField D = (1.24474593e+01/Tau + -2.20683480e+02/sqrt(Tau) + 1.20546492e+03 + -1.84072792e+03*sqrt(Tau) + 6.56727769e+02*Tau)*sqr(Delta);
    	const volScalarField E = (5.65468710e+00/Tau + 1.45713696e+01/sqrt(Tau) + -3.55818262e+02 + 9.45803639e+02*sqrt(Tau) + -4.29181971e+02*Tau)*sqr(Delta)*Delta;
    	const volScalarField F = (-6.81222200e+00/Tau + 5.41649906e+01/sqrt(Tau) + -1.26559998e+02 + 9.53260650e+01*sqrt(Tau) + -5.12543488e+01*Tau)*sqr(Delta)*sqr(Delta);

    	const volScalarField w = A+B+C+D+E+F;
	const dimensionedScalar meterPerSecond("meterPerSecond",dimensionSet(0, 1, -1, 0, 0, 0 ,0), 1.0);
	const volScalarField w_dim = w*meterPerSecond;

        return store
        (
            resultName_,
	    mag(U)/w_dim
        );
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::MachNoRealH2::MachNoRealH2
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, "U")
{
    setResultName("MaRealH2", "U");
}


// ************************************************************************* //
