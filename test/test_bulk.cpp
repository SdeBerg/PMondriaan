#include <pmondriaan.hpp>

int main() {
    bulk::thread::environment env;

    env.spawn(env.available_processors(), [](bulk::world& world) {
        // int s = world.rank();
        // int p = world.active_processors();

        auto q = bulk::queue<long>(world);
        q(0).send(2);
        world.sync();
        q(0).send(5);
        world.sync();
        for (auto item : q) {
            world.log("Received %d", item);
        }
    });

    return 0;
}
