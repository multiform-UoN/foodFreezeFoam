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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class thermo>
inline Foam::processedMeatSolidTransport<thermo>::processedMeatSolidTransport
(
    const thermo& t,
    const scalar kappa
)
:
    thermo(t),
    kappa_(kappa)
{}


template<class thermo>
inline Foam::processedMeatSolidTransport<thermo>::processedMeatSolidTransport
(
    const word& name,
    const processedMeatSolidTransport& ct
)
:
    thermo(name, ct),
    kappa_(ct.kappa_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::processedMeatSolidTransport<Thermo>>
Foam::processedMeatSolidTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<processedMeatSolidTransport<Thermo>>
    (
        new processedMeatSolidTransport<Thermo>(dict)
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class thermo>
inline Foam::scalar Foam::processedMeatSolidTransport<thermo>::
kappa(const scalar p, const scalar T) const
{

    scalar TC(T-273.15);

    if(T>Tf_)
    {
        return
            scalar(0.0866)
          + scalar(0.501)*Yw_
          + scalar(5.052e-04)*Yw_*TC;
    }
    else
    {
        return
            scalar(0.378)
          + scalar(1.376)*Yw_
          + scalar(0.930)/TC;

    }
}

template<class thermo>
inline Foam::vector Foam::processedMeatSolidTransport<thermo>::
Kappa(const scalar p, const scalar T) const
{
    return vector(kappa_, kappa_, kappa_);
}


template<class thermo>
inline Foam::scalar Foam::processedMeatSolidTransport<thermo>::
mu(const scalar p, const scalar T) const
{
    NotImplemented;
    return scalar(0);
}


template<class thermo>
inline Foam::scalar Foam::processedMeatSolidTransport<thermo>::
alphah(const scalar p, const scalar T) const
{
    return kappa_/this->Cp(p, T);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class thermo>
inline void Foam::processedMeatSolidTransport<thermo>::operator=
(
    const processedMeatSolidTransport<thermo>& ct
)
{
    thermo::operator=(ct);
    kappa_ = ct.kappa_;
}


template<class thermo>
inline void Foam::processedMeatSolidTransport<thermo>::operator+=
(
    const processedMeatSolidTransport<thermo>& ct
)
{
    scalar Y1 = this->Y();
    thermo::operator+=(ct);

    Y1 /= this->Y();
    scalar Y2 = ct.Y()/this->Y();

    kappa_ = Y1*kappa_ + Y2*ct.kappa_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class thermo>
inline Foam::processedMeatSolidTransport<thermo> Foam::operator*
(
    const scalar s,
    const processedMeatSolidTransport<thermo>& ct
)
{
    return processedMeatSolidTransport<thermo>
    (
        s*static_cast<const thermo&>(ct),
        ct.kappa_
    );
}


// ************************************************************************* //
