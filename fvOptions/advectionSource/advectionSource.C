/*---------------------------------------------------------------------------*  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "advectionSource.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "fvMatrix.H"


#include "fvCFD.H"
#include "word.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

defineTypeNameAndDebug(advectionSource, 0);
addRemovableToRunTimeSelectionTable
(
    option,
    advectionSource,
    dictionary
);




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

advectionSource::advectionSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

advectionSource::~advectionSource()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void advectionSource::correct
(
    GeometricField<scalar, fvPatchField, volMesh>& fld
)
{
}


void advectionSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
 
            Info<< "\nAdding advection to solid region: "
                << mesh().name() << endl;

	    volScalarField rho1 =  mesh().lookupObject<volScalarField>("thermo:rho");
	    volScalarField h =  mesh().lookupObject<volScalarField>("h");
	    volScalarField betav_solid =  mesh().lookupObject<volScalarField>("betavSolid");
	    volScalarField alpha_solid =  mesh().lookupObject<volScalarField>("thermo:alpha");
            if(mesh().foundObject<volVectorField>("U"))
            {
                const volVectorField& U
                (
                    mesh().lookupObject<volVectorField>("U")
                );
                const volScalarField& h = eqn.psi();
                surfaceScalarField phi = fvc::flux(U);
                eqn -= fvm::div(phi,h);
            }
            else
            {
                volVectorField U
                (
                    IOobject
                    (
                        "U",
                        mesh().time().constant(),
                        mesh(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh()
                );

                const volScalarField& h = eqn.psi();
                surfaceScalarField phi = fvc::flux(rho1*U);
                eqn -= fvm::div(phi,h);

            }
		if(mesh().time().write())
		{
			rho1.write();
			alpha_solid.write();
          		h.write();
		}
}


void advectionSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
            Info<< "\nAdding advection to solid region: "
                << mesh().name() << endl;

	    volScalarField rho1 =  mesh().lookupObject<volScalarField>("thermo:rho");
	    volScalarField h =  mesh().lookupObject<volScalarField>("h");
	    volScalarField betav_solid =  mesh().lookupObject<volScalarField>("betavSolid");
	    volScalarField alpha_solid =  mesh().lookupObject<volScalarField>("thermo:alpha");
            if(mesh().foundObject<volVectorField>("U"))
            {
                const volVectorField& U
                (
                    mesh().lookupObject<volVectorField>("U")
                );
                const volScalarField& h = eqn.psi();
                surfaceScalarField phi = fvc::flux(U);
                eqn -= fvm::div(phi,h);
            }
            else
            {
                volVectorField U
                (
                    IOobject
                    (
                        "U",
                        mesh().time().constant(),
                        mesh(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh()
                );

                const volScalarField& h = eqn.psi();
                surfaceScalarField phi = fvc::flux(rho1*U);
                eqn -= fvm::div(phi,h);

            }
		if(mesh().time().write())
		{
			rho1.write();
			alpha_solid.write();
          		h.write();
		}
}


void advectionSource::constrain
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

} // End namespace fv
// ************************************************************************* //

