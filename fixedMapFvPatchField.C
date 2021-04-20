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
    mapFrom_(),
    fieldName_(),
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
    mapFrom_(dict.get<word>("mapFrom")),
    fieldName_(dict.get<word>("fieldName")),
    offset_(dict.get<vector>("offset")),
    globalPatchInd_(1)
{
    this->patchType() = dict.getOrDefault<word>("patchType", word::null);

    const label mapFromPatchID = this->patch().patch().boundaryMesh().findPatchID(mapFrom_);

    const polyPatch& mapFromPatch =  this->patch().patch().boundaryMesh()[mapFromPatchID]; 

    const vectorField& mpfC = mapFromPatch.faceCentres();

    const vectorField& thisPatchC =this->patch().patch().faceCentres();

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

    	Info << "Creating a patchMap" << endl;

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

    	/*

    	forAll(thisPatchC,i)
    	{
    		Info << "thisPatchC = " << thisPatchC[i] <<"  , mapFromPatchC = " << mpfC[globalPatchInd_[i]] << endl;
    	}
    	*/
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
    mapFrom_(ptf.mapFrom_),
    fieldName_(ptf.fieldName_),
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
    mapFrom_(ptf.mapFrom_),
    fieldName_(ptf.fieldName_),
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
    mapFrom_(ptf.mapFrom_),
    fieldName_(ptf.fieldName_),
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

    //const polyPatch& thisPatch = this->patch().patch(); //hshs(); //.boundaryMesh().patches(); //.mesh().boundary().patches();

    wordList mapFroms = this->patch().patch().boundaryMesh().names();




    if(Pstream::parRun())
    {
    	//std::cout<<"\n UPDATING on processor:" << Pstream::myProcNo() << std::endl;
    	forAll(mapFroms,i)
    	{
    		//std::cout << "\n On processor " << Pstream::myProcNo() <<" , mapFrom[ " << i <<"]: " << mapFroms[i] << "\n"<<endl;
    		//std::cout << "mapFrom_ = " << mapFrom_ << endl;
    	}


    	
    	//PtrList<boundaryPatch> patches = this->patch().boundaryMesh().patches();
    }
    
    Info << "UPDATING MAPPED FIELD!!!!" << endl;

    const fvPatch& pt = this->patch();

    const fvPatchField<Type> & mapField = pt.lookupPatchField<GeometricField<Type, fvPatchField, volMesh>, Type>(fieldName_);

    Field<Type> newValues(this->patch().patch().size());

    forAll(newValues,i)
    {
    	newValues[i] = mapField[globalPatchInd_[i]];
    }

    this->operator==(newValues);
     
    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::fixedMapFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeEntry("mapFrom", mapFrom_);
    os.writeEntry("fieldName", fieldName_);
    os.writeEntry("offset", offset_);
    this->writeEntry("value", os);
}


// ************************************************************************* //
