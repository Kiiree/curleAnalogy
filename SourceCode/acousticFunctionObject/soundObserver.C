/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "soundObserver.H"

Foam::SoundObserver::SoundObserver(word name, vector pos)
:
name_(name),
position_(pos)
{
    // Initial fluctuating pressure to zero
    pPrime_ = 0.0;
}


Foam::SoundObserver::SoundObserver(const SoundObserver& so)
:
name_(so.name_),
position_(so.position_),
pPrime_(so.pPrime_)
{
    // Copy initial fluctuating pressure
    pPrime_ = so.pPrime();
}


Foam::SoundObserver::SoundObserver()
:
name_(Foam::word::null),
position_(vector::zero),
pPrime_(0.0)
{}


void Foam::SoundObserver::pPrime(scalar pPrime)
{
    pPrime_ = pPrime;
}

