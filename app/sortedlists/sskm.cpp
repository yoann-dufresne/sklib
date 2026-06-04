#include "cli.hpp"
#include "commands.hpp"

int main(int argc, char* argv[]) {
    auto const parsed {parse_cli(argc, argv)};

    if (parsed.construct) {
        return run_construct(*parsed.construct);
    }
    if (parsed.query) {
        return run_query(*parsed.query);
    }
    if (parsed.experiment) {
        return run_experiment(*parsed.experiment);
    }

    return 0;
}
