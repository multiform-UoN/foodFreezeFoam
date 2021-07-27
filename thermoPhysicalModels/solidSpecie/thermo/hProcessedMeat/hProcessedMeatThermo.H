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
    Foam::hProcessedMeatThermo

Description
    Constant properties thermodynamics package
    templated into the EquationOfState.

SourceFiles
    hProcessedMeatThermoI.H
    hProcessedMeatThermo.C

\*---------------------------------------------------------------------------*/

#ifndef hProcessedMeatThermo_H
#define hProcessedMeatThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class hProcessedMeatThermo;

template<class EquationOfState>
inline hProcessedMeatThermo<EquationOfState> operator+
(
    const hProcessedMeatThermo<EquationOfState>&,
    const hProcessedMeatThermo<EquationOfState>&
);

template<class EquationOfState>
inline hProcessedMeatThermo<EquationOfState> operator*
(
    const scalar,
    const hProcessedMeatThermo<EquationOfState>&
);

template<class EquationOfState>
inline hProcessedMeatThermo<EquationOfState> operator==
(
    const hProcessedMeatThermo<EquationOfState>&,
    const hProcessedMeatThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const hProcessedMeatThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class hProcessedMeatThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hProcessedMeatThermo
:
    public EquationOfState
{
    // Private data

        scalar Cp_;
        scalar Hf_;


    // Private Member Functions

        //- Construct from components
        inline hProcessedMeatThermo
        (
            const EquationOfState& st,
            const scalar cp,
            const scalar hf
        );


public:

    // Constructors

        //- Construct from dictionary
        hProcessedMeatThermo(const dictionary& dict);

        //- Construct as named copy
        inline hProcessedMeatThermo(const word&, const hProcessedMeatThermo&);

        //- Construct and return a clone
        inline autoPtr<hProcessedMeatThermo> clone() const;

        //- Selector from dictionary
        inline static autoPtr<hProcessedMeatThermo> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "hProcessedMeat<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
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


        // Derivative term used for Jacobian

            //- Derivative of Gibbs free energy w.r.t. temperature
            inline scalar dGdT(const scalar p, const scalar T) const;

            //- Temperature derivative of heat capacity at constant pressure
            inline scalar dCpdT(const scalar p, const scalar T) const;


        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const hProcessedMeatThermo&);


    // Friend operators

        friend hProcessedMeatThermo operator+ <EquationOfState>
        (
            const hProcessedMeatThermo&,
            const hProcessedMeatThermo&
        );

        friend hProcessedMeatThermo operator* <EquationOfState>
        (
            const scalar,
            const hProcessedMeatThermo&
        );

        friend hProcessedMeatThermo operator== <EquationOfState>
        (
            const hProcessedMeatThermo&,
            const hProcessedMeatThermo&
        );


    // IOstream Operators

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const hProcessedMeatThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hProcessedMeatThermoI.H"

#ifdef NoRepository
    #include "hProcessedMeatThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
