/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "rhoProcessedMeat.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoProcessedMeat<Specie>::rhoProcessedMeat
(
    const Specie& sp,
    const scalar Yw,
    const scalar Tf
)
:
    Specie(sp),
    Yw_(Yw),
    Tf_(Tf)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoProcessedMeat<Specie>::rhoProcessedMeat
(
    const word& name,
    const rhoProcessedMeat<Specie>& ico
)
:
    Specie(name, ico),
    Yw_(ico.Yw_),
    Tf_(ico.Tf_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::rhoProcessedMeat<Specie>>
Foam::rhoProcessedMeat<Specie>::clone() const
{
    return autoPtr<rhoProcessedMeat<Specie>>(new rhoProcessedMeat<Specie>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::rhoProcessedMeat<Specie>::rho(scalar p, scalar T) const
{

    scalar rho0(1053.);

    if(T>Tf_)
    {
        return rho0;
    }
    else
    {

        scalar TC(T-273.15);

        return (
            rho0
            /
            (
                scalar(0.98221)
              + scalar(0.1131)*Yw_
              + scalar(0.25746)*(scalar(1.)-Yw_)/TC
            )

        );
    }
}


template<class Specie>
inline Foam::scalar Foam::rhoProcessedMeat<Specie>::H(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoProcessedMeat<Specie>::Cp(scalar p, scalar T) const
{
    if(T>Tf_)
    {
        return (
            scalar(1448.)*(scalar(1.)-Yw_)
          + scalar(4187.)*Yw_
        );
    }
    else
    {
        scalar TC(T-273.15);

        return (
            scalar(3874.)
          - scalar(2534.)*Yw_
          + scalar(902893.)*(scalar(1.)-Yw_)/(TC*TC)

        );
    }
}


template<class Specie>
inline Foam::scalar Foam::rhoProcessedMeat<Specie>::S(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoProcessedMeat<Specie>::psi(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoProcessedMeat<Specie>::Z(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoProcessedMeat<Specie>::CpMCv(scalar p, scalar T) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::rhoProcessedMeat<Specie>::operator+=(const rhoProcessedMeat<Specie>& ico)
{
    scalar Y1 = this->Y();
    Specie::operator+=(ico);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        const scalar Y2 = ico.Y()/this->Y();

        Yw_ = Y1*Yw_ + Y2*ico.Yw_;
        Tf_ = Y1*Tf_ + Y2*ico.Tf_;
    }
}


template<class Specie>
inline void Foam::rhoProcessedMeat<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoProcessedMeat<Specie> Foam::operator+
(
    const rhoProcessedMeat<Specie>& ico1,
    const rhoProcessedMeat<Specie>& ico2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ico1)
      + static_cast<const Specie&>(ico2)
    );

    if (mag(sp.Y()) < small)
    {
        return rhoProcessedMeat<Specie>
        (
            sp,
            ico1.Yw_,
            ico1.Tf_
        );
    }
    else
    {
        const scalar Y1 = ico1.Y()/sp.Y();
        const scalar Y2 = ico2.Y()/sp.Y();

        return rhoProcessedMeat<Specie>
        (
            sp,
            Y1*ico1.Yw_ + Y2*ico2.Yw_,
            Y1*ico1.Tf_ + Y2*ico2.Tf_
        );
    }
}


template<class Specie>
inline Foam::rhoProcessedMeat<Specie> Foam::operator*
(
    const scalar s,
    const rhoProcessedMeat<Specie>& ico
)
{
    return rhoProcessedMeat<Specie>
        (s*static_cast<const Specie&>(ico), ico.Yw_);
}


template<class Specie>
inline Foam::rhoProcessedMeat<Specie> Foam::operator==
(
    const rhoProcessedMeat<Specie>& ico1,
    const rhoProcessedMeat<Specie>& ico2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ico1)
     == static_cast<const Specie&>(ico2)
    );

    const scalar Y1 = ico1.Y()/sp.Y();
    const scalar Y2 = ico2.Y()/sp.Y();

    return rhoProcessedMeat<Specie>
    (
        sp,
        Y2*ico2.Yw_ - Y1*ico1.Yw_,
        Y1*ico1.Tf_ + Y2*ico2.Tf_
    );
}


// ************************************************************************* //
