    autoPtr<radiation::radiationModel> radiation
    (
        radiation::radiationModel::New(T)
    );

    dimensionedScalar rhoCpRef
    (
        "rhoCpRef",
        dimDensity*dimEnergy/dimMass/dimTemperature,
        1.0
    );

    if (radiation->radiation())
    {
        IOdictionary transportProperties
        (
            IOobject
            (
                "transportProperties",
                runTime.constant(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false  // do not register!
            )
        );

        dimensionedScalar rhoRef(transportProperties.lookup("rhoRef"));
	// RJS 10/9/13 Prefer Cp over original CpRef for compatiblity with Cp in thermophysicalProperties
        dimensionedScalar CpRef(transportProperties.lookup("Cp"));

        rhoCpRef = rhoRef*CpRef;
    }
