if (TARGET CGAL::CGAL)
    return()
endif()

include(CPM)

# Pin the CGAL release so the exact (EPECK) query handler is reproducible and
# independent of whatever CGAL happens to be installed system-wide, matching how
# the other CPM deps are frozen. CGAL's core is header-only, so we only download
# and point find_package at the unpacked release (no install needed).
#
# NOTE: EPECK still relies on GMP, MPFR, and Boost, which CGAL's own config
# locates from the system. Those remain required system packages and are not
# vendored here (CGAL's exact kernels depend on them by design).
set(SUBGRID_CGAL_VERSION "6.0.1" CACHE STRING "Pinned CGAL release used by --cgal")

CPMAddPackage(
    NAME CGAL
    URL "https://github.com/CGAL/cgal/releases/download/v${SUBGRID_CGAL_VERSION}/CGAL-${SUBGRID_CGAL_VERSION}.tar.xz"
    DOWNLOAD_ONLY YES
)

if (CGAL_SOURCE_DIR)
    # The release tarball ships a ready CGALConfig.cmake at its root; use it
    # directly instead of any system-installed CGAL.
    find_package(CGAL REQUIRED PATHS "${CGAL_SOURCE_DIR}" NO_DEFAULT_PATH)
endif()
