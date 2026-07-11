if (TARGET fcpw::fcpw)
    return()
endif()

include(CPM)

# Use FCPW's Eigen fallback instead of Enoki on this toolchain.
set(FCPW_USE_ENOKI OFF CACHE BOOL "Disable Enoki backend for FCPW" FORCE)

# Enoki (required by fcpw)
CPMAddPackage(
    NAME enoki
    GITHUB_REPOSITORY mitsuba-renderer/enoki
    GIT_TAG 141cf4bd18eee674841a7c3e3c28f3db44adc6fa
    DOWNLOAD_ONLY YES
)
if (enoki_SOURCE_DIR)
    # Teach Enoki's ARM check about Apple Silicon's "arm64" (else it probes -march=native, which Apple clang rejects).
    file(READ "${enoki_SOURCE_DIR}/CMakeLists.txt" _enoki_cml)
    string(REPLACE "MATCHES \"aarch64\"" "MATCHES \"aarch64|arm64\"" _enoki_cml "${_enoki_cml}")
    file(WRITE "${enoki_SOURCE_DIR}/CMakeLists.txt" "${_enoki_cml}")
    add_subdirectory("${enoki_SOURCE_DIR}" "${enoki_BINARY_DIR}")
endif()

# Create a consistent alias if CPM provided an 'enoki' target
if (TARGET enoki)
    set_target_properties(enoki PROPERTIES FOLDER third_party)
    add_library(enoki::enoki ALIAS enoki)
endif()

# fcpw itself
CPMAddPackage(
    NAME fcpw
    GITHUB_REPOSITORY rohan-sawhney/fcpw
    GIT_TAG 61814ff0c7b69d61dac3b9725cda1541b1b3ec4f
)

# If CPM created a target called 'fcpw', expose an alias and attach Enoki if available.
if (TARGET fcpw)
    set_target_properties(fcpw PROPERTIES FOLDER third_party)
    add_library(fcpw::fcpw ALIAS fcpw)

    # If Enoki was provided, make sure fcpw links to it so include dirs propagate.
    if (TARGET enoki::enoki OR TARGET enoki)
        if (TARGET enoki::enoki)
            target_link_libraries(fcpw PUBLIC enoki::enoki)
        else()
            target_link_libraries(fcpw PUBLIC enoki)
            add_library(enoki::enoki ALIAS enoki) # create alias for consistency
        endif()
    endif()

# If CPM created 'fcpw::fcpw' target directly, ensure it links to Enoki as well.
elseif (TARGET fcpw::fcpw)
    set_target_properties(fcpw::fcpw PROPERTIES FOLDER third_party)
    if (TARGET enoki::enoki OR TARGET enoki)
        if (TARGET enoki::enoki)
            target_link_libraries(fcpw::fcpw PUBLIC enoki::enoki)
        else()
            target_link_libraries(fcpw::fcpw PUBLIC enoki)
            add_library(enoki::enoki ALIAS enoki)
        endif()
    endif()
else()
    message(WARNING "CPM did not create a target named 'fcpw' or 'fcpw::fcpw'. Inspect CPM results.")
endif()