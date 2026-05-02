#include <cerrno>
#include <chrono>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <string_view>

#include <sys/resource.h>
#include <sys/time.h>
#include <sys/wait.h>
#include <unistd.h>

void usage(const char* program_name, std::ostream& out) {
    out << "usage: " << program_name << " [--help | -h] COMMAND ...\n";
}

long long int micros(const timeval& tv) {
    return tv.tv_sec * 1'000'000LL + tv.tv_usec;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        usage(argv[0], std::cerr);
        return 1;
    }
    const std::string_view first = argv[1];
    if (first == "-h" || first == "--help") {
        usage(argv[0], std::cout);
        return 0;
    }

    const auto before = std::chrono::steady_clock::now();

    const pid_t child = fork();
    if (child == pid_t(-1)) {
        const int error = errno;
        std::cerr << "Unable to fork process: " << std::strerror(error) << '\n';
        return 1;
    }
    if (child == 0) {
        // We're the child process.
        const int rc = execv(argv[1], &argv[1]);
        if (rc == -1) {
            const int error = errno;
            std::cerr << "Unable to execute program: " << std::strerror(error) << '\n';
            return 1;
        }
    }

    // We're the parent process.
    pid_t rc;
    int child_status;
    do {
        const int flags = 0;
        rc = waitpid(child, &child_status, flags);
    } while (rc == pid_t(-1) && errno == EINTR);

    const auto after = std::chrono::steady_clock::now();

    if (rc == pid_t(-1)) {
        const int error = errno;
        std::cerr << "Unable to wait for child process: " << std::strerror(error) << '\n';
        return 1;
    }

    // `child_status` isn't necessarily the exit status: we need to upack it first.
    if (WIFEXITED(child_status)) {
        child_status = WEXITSTATUS(child_status);
    } else if (WIFSIGNALED(child_status)) {
        child_status = -WTERMSIG(child_status);
    }

    struct rusage child_usage;
    if (-1 == getrusage(RUSAGE_CHILDREN, &child_usage)) {
        const int error = errno;
        std::cerr << "Unable to get resource usage info of child process: " << std::strerror(error) << '\n';
        return 1;
    }
    
    // TODO: `std::quoted` is not adequate for JSON, but it will work for arguments that don't
    // contain backslashes.
    std::cout << "{\"status\": " << child_status << ", \"command\": [" << std::quoted(argv[1]);
    for (int i = 2; i < argc; ++i) {
        std::cout << ", " << std::quoted(argv[i]);
    }
    std::cout << "], \"cpu_user_micros\": " << micros(child_usage.ru_utime);
    std::cout << ", \"cpu_system_micros\": " << micros(child_usage.ru_stime);
    std::cout << ", \"wall_micros\": " << std::chrono::duration_cast<std::chrono::microseconds>(after - before).count();
    std::cout << ", \"max_resident_size_kb\": " << child_usage.ru_maxrss;
    std::cout << "}\n";
}
