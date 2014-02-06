/*---------------------------------------------------------------------------*\
    Copyright            : (C) 2013 Symscape
    Website              : www.symscape.com
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

Application
    mapDualFields

Description
    Maps volume and boundary face fields from one mesh to another, reading and
    interpolating all fields present in the time directory of the source case.
    Parallel cases need to be reconstructed.
    Based on mapFields.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "labelListIOList.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


class MapDualFields
{
public:
  MapDualFields
  (
   fvMesh & sourceMesh,
   fvMesh & targetMesh,
   const labelListList & cellMap,
   const labelListList & faceMap
   );

  void perform();

private:
  template<class Type>
  void mapVolFields(const IOobjectList& objects);

  template<class Type>
  tmp<GeometricField<Type, fvPatchField, volMesh> > interpolate
  (
   const GeometricField<Type, fvPatchField, volMesh>& sourceVf
   );

  template<class Type>
  void mapCellField
  (
   const Field<Type>& sourceF,
   Field<Type>& targetF
   );

  template<class Type>
  void mapFaceField
  (
   const fvPatchField<Type> & sourceF,
   const label sourceStart,
   fvPatchField<Type> & targetF,
   const label targetStart
   );

private:
  fvMesh & sourceMesh_;
  fvMesh & targetMesh_;
  const labelListList & cellMap_;
  const labelListList & faceMap_;
}; // class MapDualFields


MapDualFields::MapDualFields
(
 fvMesh & sourceMesh,
 fvMesh & targetMesh,
 const labelListList & cellMap,
 const labelListList & faceMap
 )
  : sourceMesh_(sourceMesh)
  , targetMesh_(targetMesh)
  , cellMap_(cellMap)
  , faceMap_(faceMap)
{}


template<class Type>
void
MapDualFields::mapCellField
(
 const Field<Type>& sourceF,
 Field<Type>& targetF
 )
{
  forAll(targetF, cellI)
    {      
      const labelList & sourceCells = cellMap_[cellI];
      targetF[cellI] = pTraits<Type>::zero;

      forAll(sourceCells, sourceCellI)
	{
	  targetF[cellI] += sourceF[sourceCells[sourceCellI]];
	}

      targetF[cellI] /= scalar(sourceCells.size());
    }
}


template<class Type>
void
MapDualFields::mapFaceField
(
 const fvPatchField<Type> & sourceF,
 const label sourceStart,
 fvPatchField<Type> & targetF,
 const label targetStart
 )
{
  forAll(targetF, faceI)
    {      
      const labelList & sourceFaces = faceMap_[faceI + targetStart];
      targetF[faceI] = pTraits<Type>::zero;

      forAll(sourceFaces, sourceFaceI)
	{
	  targetF[faceI] += sourceF[sourceFaces[sourceFaceI] - sourceStart];
	}

      targetF[faceI] /= scalar(sourceFaces.size());
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > 
MapDualFields::interpolate
(
 const GeometricField<Type, fvPatchField, volMesh>& sourceVf
 )
{
  Field<Type> internalField(targetMesh_.nCells());
  mapCellField(sourceVf, internalField);

  PtrList<fvPatchField<Type> > patchFields
    (
     sourceVf.boundaryField().size()
     );

  forAll(sourceVf.boundaryField(), patchI)
    {
      const label sourceStart = sourceMesh_.boundaryMesh()[patchI].start();
      const label targetStart = targetMesh_.boundaryMesh()[patchI].start();

      const fvPatchField<Type> & sourceP = sourceVf.boundaryField()[patchI];
      tmp<fvPatchField<Type> > targetP = fvPatchField<Type>::New
      (
       sourceP.type(),
       targetMesh_.boundary()[patchI],
       DimensionedField<Type, volMesh>::null()
       );

      mapFaceField(sourceP, sourceStart, targetP(), targetStart);

      patchFields.set
        (
	 patchI,
	 targetP
	 );
    }

  tmp<GeometricField<Type, fvPatchField, volMesh> > ttoF
    (
     new GeometricField<Type, fvPatchField, volMesh>
     (
      IOobject
      (
       "interpolated(" + sourceVf.name() + ')',
       targetMesh_.time().timeName(),
       targetMesh_,
       IOobject::NO_READ,
       IOobject::NO_WRITE
       ),
      targetMesh_,
      sourceVf.dimensions(),
      internalField,
      patchFields
      )
     );

  return ttoF;
}


template<class Type>
void  
MapDualFields::mapVolFields(const IOobjectList& objects)
{
  word fieldClassName
    (
     GeometricField<Type, fvPatchField, volMesh>::typeName
     );

  IOobjectList fields = objects.lookupClass(fieldClassName);

  forAllIter(IOobjectList, fields, fieldIter)
    {
      Info<< "    interpolating " << fieldIter()->name()
	  << endl;

      // Read field
      GeometricField<Type, fvPatchField, volMesh> sourceField
        (
	 *fieldIter(),
	 sourceMesh_
	 );

      IOobject targetFieldIOobject
        (
	 fieldIter()->name(),
	 targetMesh_.time().timeName(),
	 targetMesh_,
	 IOobject::NO_READ,
	 IOobject::AUTO_WRITE
	 );

      GeometricField<Type, fvPatchField, volMesh> targetField
	(
	 targetFieldIOobject,
	 interpolate(sourceField)
	 );

      targetField.write();
    }
}


void 
MapDualFields::perform()
{
  {
    IOobjectList objects(sourceMesh_, sourceMesh_.time().timeName());

    mapVolFields<scalar>(objects);
    mapVolFields<vector>(objects);
    mapVolFields<sphericalTensor>(objects);
    mapVolFields<symmTensor>(objects);
    mapVolFields<tensor>(objects);
  }
}


static
void
renumberElements
(
 const string & name,
 labelListList & elementMap,
 fvMesh & sourceMesh,
 const label elementCount
 )
{
  labelList oldToNewMap;

  {
    labelIOList newToOldMap
      (
       IOobject
       (
	name,
	sourceMesh.facesInstance(),
	polyMesh::meshSubDir,
	sourceMesh,
	IOobject::READ_IF_PRESENT,
	IOobject::NO_WRITE,
	false
	)
       );

    if (elementCount != newToOldMap.size()) return;

    Info<< "Renumber elements according to: " << name
	<< endl;

    oldToNewMap.setSize(newToOldMap.size());

    // Invert to give old -> new
    forAll(newToOldMap, elementI)
      {
	oldToNewMap[newToOldMap[elementI]] = elementI;
      }
  }

  forAll(elementMap, elementI)
    {
      labelList & sourceElements = elementMap[elementI];

      forAll(sourceElements, sourceElementI)
	{
	  sourceElements[sourceElementI] = oldToNewMap[sourceElements[sourceElementI]];
	}
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
  argList::addNote
    (
     "map volume and boundary faces fields from dual mesh to original mesh"
     );
  argList::noParallel();
  argList::validArgs.append("sourceCase");

  argList::addOption
    (
     "sourceTime",
     "scalar",
     "specify the source time"
     );

  argList args(argc, argv);

  if (!args.check())
    {
      FatalError.exit();
    }

  fileName rootDirTarget(args.rootPath());
  fileName caseDirTarget(args.globalCaseName());

  const fileName casePath = args[1];
  const fileName rootDirSource = casePath.path();
  const fileName caseDirSource = casePath.name();

  Info<< "Source: " << rootDirSource << " " << caseDirSource << endl;
  word sourceRegion = fvMesh::defaultRegion;

  Info<< "Target: " << rootDirTarget << " " << caseDirTarget << endl;
  word targetRegion = fvMesh::defaultRegion;

#include "createTimes.H" 
#include "setTimeIndex.H"

  Info<< "Create meshes\n" << endl;

  fvMesh sourceMesh
    (
     IOobject
     (
      sourceRegion,
      runTimeSource.timeName(),
      runTimeSource
      )
     );

  fvMesh targetMesh
    (
     IOobject
     (
      targetRegion,
      runTimeTarget.timeName(),
      runTimeTarget
      )
     );

  if (sourceMesh.boundary().size() != targetMesh.boundary().size())
    {
      FatalErrorIn
        (
	 "mapDualCells"
	 ) << "Incompatible meshes: different number of boundaries"
	   << exit(FatalError);
    }
  
  Info<< "Source mesh size: " << sourceMesh.nCells() << tab
      << "Target mesh size: " << targetMesh.nCells() << nl << endl;

 labelListCompactIOList cellMap
   (
    IOobject
    (
     "cellDualMap",
     sourceMesh.facesInstance(),
     polyMesh::meshSubDir,
     sourceMesh,
     IOobject::MUST_READ,
     IOobject::NO_WRITE,
     false
     )
    );
  
  if (cellMap.size() != targetMesh.nCells()) 
    {
      FatalErrorIn
        (
	 "mapDualFields"
	 ) << "cellDualMap does not match target mesh"
	   << exit(FatalError);     
    }

  renumberElements("cellMap", cellMap, sourceMesh, sourceMesh.nCells());

  labelListCompactIOList faceMap
    (
     IOobject
     (
      "faceDualMap",
      sourceMesh.facesInstance(),
      polyMesh::meshSubDir,
      sourceMesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE,
      false
      )
     );
  
  if (faceMap.size() != targetMesh.nFaces()) 
    {
      FatalErrorIn
        (
	 "mapDualFields"
	 ) << "faceDualMap does not match target mesh"
	   << exit(FatalError);     
    }

  renumberElements("faceMap", faceMap, sourceMesh, sourceMesh.nFaces());

  Info<< nl
      << "Creating and mapping fields for time "
      << sourceMesh.time().timeName() << nl << endl;

  MapDualFields mapDualFields(sourceMesh, 
			      targetMesh, 
			      cellMap, 
			      faceMap);
  mapDualFields.perform();

  Info<< "\nEnd\n" << endl;

  return 0;
}


// ************************************************************************* //
