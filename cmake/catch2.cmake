if (TARGET Catch2::Catch2WithMain)
    return()
endif()

include(CPM)

CPMAddPackage(
    NAME Catch2
    GITHUB_REPOSITORY catchorg/Catch2
    GIT_TAG v3.7.1
)
