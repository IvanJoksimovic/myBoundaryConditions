/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "fixedMapFvPatchField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fixedMapFvPatchField<Type>::fixedMapFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    sourcePatch_(),
    offset_(),
    globalPatchInd_()
{}


template<class Type>
Foam::fixedMapFvPatchField<Type>::fixedMapFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    sourcePatch_(dict.get<word>("sourcePatch")),
    offset_(dict.get<vector>("offset")),
    globalPatchInd_(1)
{


    this->patchType() = dict.getOrDefault<word>("patchType", word::null);

    Info << "CREATING A PATCH MAP ON THE PATCH: " << this->patch().patch().name()<< endl;

    const label sourcePatchID = this->patch().patch().boundaryMesh().findPatchID(sourcePatch_);

    const polyPatch& sourcePatch =  this->patch().patch().boundaryMesh()[sourcePatchID];

    const vectorField& mpfC = sourcePatch.faceCentres();

    const vectorField& thisPatchC =this->patch().patch().faceCentres();
    Info << "SIZE OF THE PATCH: " << thisPatchC.size() << endl;
    Info << "SIZE OF THE PATCH: " << mpfC.size() << endl;

    if(Pstream::parRun())
    {
    	// Send field to master processor
    	if(!Pstream::master)
    	{

    	}
    	else
    	{

    	}


    }
    else
    {
    	// Search in serial



    	const vectorField mpfCTransformed(mpfC - offset_);



    	globalPatchInd_.resize(thisPatchC.size());

    	forAll(thisPatchC,i)
    	{
    		scalar minDist = 1e20;
    		label minInd = -1;

    		forAll(mpfCTransformed,j)
    		{
    			scalar dist = mag(thisPatchC[i] - mpfCTransformed[j]);

    			if(dist<minDist)
    			{
    				minDist = dist;
    				minInd = j;
    			}
    		}

    		globalPatchInd_[i] = minInd;
    	}



    	forAll(thisPatchC,i)
    	{
    		Info << "thisPatchC = " << thisPatchC[i] <<"  , sourcePatchC = " << mpfC[globalPatchInd_[i]] << endl;
    	}

    }

}


template<class Type>
Foam::fixedMapFvPatchField<Type>::fixedMapFvPatchField
(
    const fixedMapFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    sourcePatch_(ptf.sourcePatch_),
    offset_(ptf.offset_),
    globalPatchInd_(1)
{}


template<class Type>
Foam::fixedMapFvPatchField<Type>::fixedMapFvPatchField
(
    const fixedMapFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    sourcePatch_(ptf.sourcePatch_),
    offset_(ptf.offset_),
    globalPatchInd_(1)
{}


template<class Type>
Foam::fixedMapFvPatchField<Type>::fixedMapFvPatchField
(
    const fixedMapFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    sourcePatch_(ptf.sourcePatch_),
    offset_(ptf.offset_),
    globalPatchInd_(1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fixedMapFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const GeometricField<Type, fvPatchField, volMesh>& f
    (
        dynamic_cast<const GeometricField<Type, fvPatchField, volMesh>&>
        (
            this->internalField()
        )
    );

    const fvPatch& p = this->patch();
    const label sourcePatchID = p.patch().boundaryMesh().findPatchID(sourcePatch_);

    if (sourcePatchID < 0)
    {
        FatalErrorInFunction
            << "Unable to find patch " << sourcePatch_
            << abort(FatalError);
    }

    //const fvPatch& sourcePatch = p.boundaryMesh()[sourcePatchID];

    const fvPatchField<Type>& sourcePatchField = f.boundaryField()[sourcePatchID];

    Info << "Average of sourcePatchField is " <<  gAverage(sourcePatchField ) << endl;

    auto tnewValues = tmp<Field<Type>>::New();
    auto& newValues = tnewValues.ref();

    newValues.setSize(sourcePatchField.size());

    forAll(newValues,i)
    {
      newValues[i] = sourcePatchField[globalPatchInd_[i]];
    }

    Info << "Average of newValues is " <<  gAverage(newValues) << endl;


    this->operator==(tnewValues);

    fixedValueFvPatchField<Type>::updateCoeffs();

    Info << "GOTOVO!!!!" << endl;

}


template<class Type>
void Foam::fixedMapFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeEntry("sourcePatch", sourcePatch_);
    os.writeEntry("offset", offset_);
    this->writeEntry("value", os);
}


// ************************************************************************* //
