export module ConstantsAndUtils;
import std;
using namespace std;

export constexpr double sparseChance = 0.3;
export constexpr double denseChance = 0.8;
export constexpr double normalizationRatio = 0.999;
export constexpr double armijoRatio = 1e-12;
export constexpr double constError = 1e-2;
export double roundDouble(double value) {
	return round(value * 100.0) / 100.0;
}
