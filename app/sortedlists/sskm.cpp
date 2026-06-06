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
    if (parsed.setop) {
        return run_setop(*parsed.setop);
    }

    return 0;
}
