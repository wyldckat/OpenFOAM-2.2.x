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

Class
    Foam::wallHeatFluxWriter

Description
    Write wallHeatFlux field

SourceFiles
    wallShearStress.C

\*---------------------------------------------------------------------------*/

#ifndef wallHeatFluxWriter_H
#define wallHeatFluxWriter_H

#include "derivedFieldWriter.H"
#include "volFieldsFwd.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  class dictionary;


  /*---------------------------------------------------------------------------*\
    Class wallHeatFluxWriter Declaration
  \*---------------------------------------------------------------------------*/

  class wallHeatFluxWriter : public derivedFieldWriter
  {
    // Private data

    // Private methods

    //- Return heat transfer coefficient
    tmp<volScalarField> htc(const objectRegistry& obr, 
			    bool & foundField) const;
  public:
    //- Runtime type information
    TypeName("wallHeatFluxWriter");

    //- Field file name.
    static const word FIELD_NAME;

    // Constructors

    wallHeatFluxWriter(const dictionary& dict);

    // Member Functions

    //- Write derived field.
    virtual void write(const objectRegistry& obr);
  };

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
