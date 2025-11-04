#pragma once

#include <cstdint>
#include <vector>

#include "board.hpp"

namespace solve {

enum class Mark : uint8_t { None = 0, Safe = 1, Mine = 2, Guess = 3, Chord = 4, FlagForChord = 5 };

struct Overlay {
	std::vector<Mark> marks;
	bool hasGuaranteedSafe = false;
	std::vector<double> mineProbability;
};

Overlay compute_overlay(const game::Board& board, int totalMines=-1, bool enableChords=true, int threads=0);

}



