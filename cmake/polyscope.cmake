if (TARGET polyscope::polyscope)
    return()
endif()

include(CPM)
CPMAddPackage(
    NAME polyscope
    GITHUB_REPOSITORY nmwsharp/polyscope
    GIT_TAG 59da72df6517cab8379865899bdffdbc96171301
)

set_target_properties(polyscope PROPERTIES FOLDER third_party)
add_library(polyscope::polyscope ALIAS polyscope)
