/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Class
    Foam::processedMeatSolidTransport

Description
    Constant properties Transport package.
    Templated into a given thermodynamics package (needed for thermal
    conductivity).

SourceFiles
    processedMeatSolidTransportI.H
    processedMeatSolidTransport.C

\*---------------------------------------------------------------------------*/

#ifndef processedMeatSolidTransport_H
#define processedMeatSolidTransport_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Thermo> class processedMeatSolidTransport;

template<class Thermo>
inline processedMeatSolidTransport<Thermo> operator*
(
    const scalar,
    const processedMeatSolidTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const processedMeatSolidTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                   Class processedMeatSolidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class processedMeatSolidTransport
:
    public Thermo
{
    // Private data

        //- Constant isotropic thermal conductivity
        scalar kappa_;

        //- Water content
        scalar Yw_;

        //- Freezing temperature
        scalar Tf_;


    // Private Member Functions

        //- Construct from components
        inline processedMeatSolidTransport(const Thermo& t, const scalar kappa);


public:

    // Constructors

        //- Construct as named copy
        inline processedMeatSolidTransport
        (
            const word&,
            const processedMeatSolidTransport&
        );

        //- Construct from dictionary
        processedMeatSolidTransport(const dictionary& dict);

        // Selector from dictionary
        inline static autoPtr<processedMeatSolidTransport> New
        (
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "processedMeat<" + Thermo::typeName() + '>';
        }

        //- Is the thermal conductivity isotropic
        static const bool isotropic = true;

        //- Isotropic thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Un-isotropic thermal conductivity [W/mK]
        inline vector Kappa(const scalar p, const scalar T) const;

         //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;


        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline void operator=(const processedMeatSolidTransport&);
        inline void operator+=(const processedMeatSolidTransport&);


    // Friend operators

        friend processedMeatSolidTransport operator* <Thermo>
        (
            const scalar,
            const processedMeatSolidTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const processedMeatSolidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "processedMeatSolidTransportI.H"

#ifdef NoRepository
    #include "processedMeatSolidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
