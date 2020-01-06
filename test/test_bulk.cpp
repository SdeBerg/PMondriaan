#include <pmondriaan.hpp>

int main() {
    bulk::thread::environment env;

    env.spawn(env.available_processors(), [](bulk::world& world) {
        int s = world.rank();
        int p = world.active_processors();

        // hello world
        for (int t = 0; t < p; ++t) {
            if (s == t)
                std::cout << "Hello, world " << s << "/" << p << "\n";
            world.sync();
        }
    });

    return 0;
}
