#include <pmondriaan.hpp>

int main() {
    bulk::thread::environment env;

    env.spawn(env.available_processors(), [](bulk::world& world) {
        int s = world.rank();
        // int p = world.active_processors();

        auto coar = bulk::coarray<int>(world, 1);
        coar[0] = s + 1;
        auto f = coar(0)[0].get();
        world.sync();
        world.log("s %d: %d", s, f.value());
        world.sync();
        world.log("s %d: %d", s, f.value());
    });

    return 0;
}
