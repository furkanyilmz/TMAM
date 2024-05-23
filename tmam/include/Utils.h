#pragma once

enum LinearUnwrappingType
{
	UNIFORM_WEIGHTS,
	HARMONIC_WEIGHTS,
	MEAN_VALUE_WEIGHTS,
	RANDOM_WEIGHTS
};

enum MappingType
{
	DISC,
	SQUARE,
	CONVEX_HULL,
	ARBITRARY
};
static bool isWireframeMode = false;