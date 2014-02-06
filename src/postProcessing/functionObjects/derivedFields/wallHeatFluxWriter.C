/*---------------------------------------------------------------------------*\
    Copyright            : (C) 2012 Symscape
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

\*---------------------------------------------------------------------------*/

#include "wallHeatFluxWriter.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "wallFvPatch.H"

#include "incompressible/RAS/RASModel/RASModel.H"
#include "compressible/RAS/RASModel/RASModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
  defineTypeNameAndDebug(wallHeatFluxWriter, 0);
  const word wallHeatFluxWriter::FIELD_NAME("wallHeatFlux");
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::wallHeatFluxWriter::htc(const objectRegistry& obr,
			      bool & foundField) const
{
  if (obr.foundObject<incompressible::RASModel>("RASProperties"))
    {
      const incompressible::RASModel& ras
        = obr.lookupObject<incompressible::RASModel>("RASProperties");
      dimensionedScalar Pr(ras.transport().lookup("Pr"));
      dimensionedScalar Cp0(ras.transport().lookup("Cp"));

      return Cp0 * (ras.nu()/Pr + obr.lookupObject<volScalarField>("alphat"));
    }
  else if (obr.foundObject<compressible::RASModel>("RASProperties"))
    {
      const compressible::RASModel& ras
        = obr.lookupObject<compressible::RASModel>("RASProperties");

      return ras.thermo().Cp() * ras.alphaEff();
    }
  else
    {
      foundField = false;
      return volScalarField::null();
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallHeatFluxWriter::wallHeatFluxWriter(const dictionary & dict)
  : derivedFieldWriter(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallHeatFluxWriter::write(const objectRegistry& obr)
{
  const fvMesh& mesh = refCast<const fvMesh>(obr);

  bool foundField = true;
  word hName;
  tmp<volScalarField> thtc = htc(obr, foundField);


  volScalarField wallHeatFlux
    (
     IOobject
     (
      FIELD_NAME,
      obr.time().timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE
      ),
     mesh,
     dimensionedScalar
     (
      FIELD_NAME,
      dimensionSet(1, 0, -3, 0, 0, 0, 0), // W.m-2
      0.
      )
     );


  if (foundField) {
    const volScalarField & T = obr.lookupObject<volScalarField>("T");

    surfaceScalarField heatFlux
      (
       fvc::interpolate(thtc())*fvc::snGrad(T)
       );

    const surfaceScalarField::GeometricBoundaryField& patchHeatFlux =
      heatFlux.boundaryField();

    const fvPatchList& patches = mesh.boundary();
	    
    volScalarField::GeometricBoundaryField & wHeatFluxb =
      wallHeatFlux.boundaryField();
    
    forAll(patches, patchi)
      {
        const fvPatch& currPatch = patches[patchi];
        
        if (isA<wallFvPatch>(currPatch))
          {	    
	    scalarField & whfp = wHeatFluxb[patchi];
	    whfp = patchHeatFlux[patchi];
	    writeStatistics(whfp, currPatch.name(), FIELD_NAME);
	  }
      }
  }

  wallHeatFlux.write();
 }


// ************************************************************************* //
