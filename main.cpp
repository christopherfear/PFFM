#include "modeI.h"

int main() {
    // Everything is handled by the objects
    ModeISimulation sim;
    
    try {
        sim.Run();
    } catch (const std::exception& e) {
        std::cerr << "CRITICAL ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}