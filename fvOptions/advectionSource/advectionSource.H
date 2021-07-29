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

Description
    Adds an advection term to the governing equation using a pre-defined velocity field.
    To be used for the energy equation in chtMultiRegionFoam

Author
    Federico Municchi. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef advectionSource_H
#define advectionSource_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fv
{

/*---------------------------------------------------------------------------*                         A templated CodedFvOption
\*---------------------------------------------------------------------------*/

class advectionSource
:
    public cellSetOption
{
public:

    //- Runtime type information
    TypeName("advectionSource");


    // Constructors

        //- Construct from patch and internal field
        advectionSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~advectionSource();


    // Member Functions

            //- Correct field
            virtual void correct
            (
                GeometricField<scalar, fvPatchField, volMesh>&
            );

            //- Explicit and implicit matrix contributions
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

            //- Explicit and implicit matrix contributions for compressible
            //  equations
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

            //- Set value
            virtual void constrain
            (
                fvMatrix<scalar>& eqn,
                const label fieldi
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

