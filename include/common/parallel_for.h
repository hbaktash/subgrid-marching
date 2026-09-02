#pragma once

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

// Minimal data-parallel loop over an index range. std::thread only -- no OpenMP
// (AppleClang ships none by default) and no TBB, because the work here is coarse
// enough that a real scheduler buys nothing.

// Resolve a thread count: <= 0 means "all cores".
inline int resolve_thread_count(int requested) {
    if (requested > 0) return requested;
    const unsigned hc = std::thread::hardware_concurrency();
    return hc > 0 ? static_cast<int>(hc) : 1;
}

// Run fn(i) for every i in [begin, end).
//
// Chunks are handed out dynamically from an atomic cursor rather than split into
// fixed per-thread ranges: per-tet cost varies by orders of magnitude (an empty
// tet versus one with several boundary curves) and the non-empty tets cluster
// near the surface, so a static split would leave most threads idle.
//
// An exception thrown by fn is captured and rethrown on the calling thread; the
// first one wins and the remaining chunks are abandoned. Letting it escape a
// worker would call std::terminate.
template <typename Fn>
void parallel_for(size_t begin, size_t end, int num_threads, Fn&& fn, size_t grain = 64) {
    if (end <= begin) return;
    if (grain == 0) grain = 1;

    const size_t total = end - begin;
    int nt = resolve_thread_count(num_threads);
    nt = static_cast<int>(std::min<size_t>(static_cast<size_t>(nt), (total + grain - 1) / grain));
    if (nt <= 1) {
        for (size_t i = begin; i < end; ++i) fn(i);
        return;
    }

    std::atomic<size_t> cursor{begin};
    std::atomic<bool> failed{false};
    std::exception_ptr error;
    std::mutex error_mutex;

    auto worker = [&]() {
        for (;;) {
            if (failed.load(std::memory_order_relaxed)) return;
            const size_t lo = cursor.fetch_add(grain, std::memory_order_relaxed);
            if (lo >= end) return;
            const size_t hi = std::min(lo + grain, end);
            try {
                for (size_t i = lo; i < hi; ++i) fn(i);
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!error) error = std::current_exception();
                failed.store(true, std::memory_order_relaxed);
                return;
            }
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(static_cast<size_t>(nt) - 1);
    for (int t = 0; t < nt - 1; ++t) pool.emplace_back(worker);
    worker();  // the calling thread takes a share too
    for (auto& th : pool) th.join();

    if (error) std::rethrow_exception(error);
}
