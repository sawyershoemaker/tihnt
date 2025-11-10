#pragma once

#include <cstdint>
#include <vector>

#include "board.hpp"

namespace solve {

enum class Mark : uint8_t {
	None = 0,
	Safe = 1,
	Mine = 2,
	Guess = 3,
	Chord = 4,              // chord center that still needs supporting flags placed
	FlagForChord = 5,       // cell that should be flagged to prepare a chord
	ChordReady = 6,         // chord center with all required flags already placed
	FlagForChordReady = 7   // cell already flagged and contributing to a chord
};

struct Overlay {
	std::vector<Mark> marks;
	bool hasGuaranteedSafe = false;
	std::vector<double> mineProbability;
};

Overlay compute_overlay(const game::Board& board, int totalMines=-1, bool enableChords=true, int threads=0);

}



