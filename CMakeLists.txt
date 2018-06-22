add_definitions(-std=c++11 -pthread)
set(CMAKE_CPP_FLAGS "-DNoMySQL -DNoXerces -Wno-deprecated-declarations -m64 -DNoQT")
set(CMAKE_CXX_FLAGS "-g -O2 -fPIC -fno-strict-aliasing -Wall -fno-strict-aliasing")
set(${CMAKE_LINK_LIBRARY_FLAG} -pthread)

add_subdirectory(math_h)
add_subdirectory(FitGen)
add_subdirectory(WASA-libs)
include_directories(${MATH_H_INC})
include_directories(${FITGEN_INC})
include_directories(${WASA_LIBS_INC})
link_directories(${ROOT_LIBRARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/WASA-RootSorterWrap/DataFiles DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/WASA-RootSorterWrap)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/scripts/config DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/scripts/WMC DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/scripts/RunPreselection.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/scripts/RunPreselectedAnalysis.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/scripts/RunAnalysis-mc.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/scripts/RunAnalysis-Reconstruction.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/scripts/RunWMC.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/.sorterrc DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/distributions DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/WASA-RootSorterWrap)
if(COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)
if(NOT $ENV{PLUTOSYS} STREQUAL "")
    link_directories($ENV{PLUTOSYS})
    include_directories(
        $ENV{PLUTOSYS}/src
        $ENV{PLUTOSYS}/plugins/brems
        $ENV{PLUTOSYS}/plugins/dalitz_mod
        $ENV{PLUTOSYS}/plugins/elementary
        $ENV{PLUTOSYS}/plugins/eta_decays
        $ENV{PLUTOSYS}/plugins/fairroot
        $ENV{PLUTOSYS}/plugins/hades
        $ENV{PLUTOSYS}/plugins/nucleus_fermi
        $ENV{PLUTOSYS}/plugins/pdg_unigen
        $ENV{PLUTOSYS}/plugins/pion_beam
        $ENV{PLUTOSYS}/plugins/scatter_mod
        $ENV{PLUTOSYS}/plugins/strangeness
        $ENV{PLUTOSYS}/plugins/w_decays
    )
    file(GLOB SIMULATIONS Simulation/sim-*.cpp)
    foreach(example ${SIMULATIONS})
        get_filename_component(ex_app ${example} NAME_WE)
        add_executable(${ex_app} ${example} Simulation/runpluto.cpp)
        target_link_libraries(${ex_app} Pluto ${ROOT_LIBRARIES} ${libs} experiment kinematics)
    endforeach()
    file(GLOB SIMULATIONS_ADV Simulation-advanced/sim-*.cpp)
    foreach(example ${SIMULATIONS_ADV})
	get_filename_component(ex_app ${example} NAME_WE)
	add_executable(${ex_app} ${example} Simulation-advanced/runsim.cpp Simulation-advanced/bound.cpp)
	target_link_libraries(${ex_app} Pluto ${ROOT_LIBRARIES} ${libs} experiment kinematics)
    endforeach()
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/Simulation-advanced/pp DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()
if(NOT $ENV{ROOTSORTERSYS} STREQUAL "")
    include_directories($ENV{ROOTSORTERSYS}/core/include)
    link_directories($ENV{ROOTSORTERSYS}/core/lib)
    link_directories($ENV{CERN_ROOT}/lib)
    set(WASA_LIBRARIES 
	WasaAnaRaw WasaRecMDC WasaAnaCalib WasaCalib 
	WasaParameter WasaEms WasaREvent WasaMain WasaAlig WasaCommon WasaBanks
	WasaRecFD WasaRecSE WasaRecPS WasaRecCD WasaRecFPC vt
	Core Cint RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread 
	m dl Thread XMLParser RSCore RSXml RSCommon RSParameter RSEms RSREvent z
	WasaRecAll MLP TreePlayer 
	kfit105 packlib gfortran
    )


    file(GLOB PRESELECTION_LIB WASA-RootSorterWrap/*.cpp)
    root_generate_dictionary(mainDict
	WASA-RootSorterWrap/analysismodule.hh
	LINKDEF WASA-RootSorterWrap/LinkDef.hh 
    )
    add_library(sorterwrap SHARED ${PRESELECTION_LIB} mainDict)
    target_link_libraries(sorterwrap
	${WASA_LIBRARIES} ${ROOT_LIBRARIES}
	FitGen experiment kinematics reconstructionfit
    )
    file(GLOB RECONSTRUCTION Reconstruction/*.cpp)
    add_library(reconstruction SHARED ${RECONSTRUCTION})
    target_link_libraries(reconstruction 
	sorterwrap ${WASA_LIBRARIES} ${ROOT_LIBRARIES}
    )
    file(GLOB REACTIONS Reactions/*.cpp)
    add_library(reactions SHARED ${REACTIONS})
    target_link_libraries(reactions reconstruction
	sorterwrap ${WASA_LIBRARIES} ${ROOT_LIBRARIES}
    )
    add_executable(rawanalysis rawanalysis.cpp)
    target_link_libraries(rawanalysis 
	sorterwrap reactions reconstruction
	${WASA_LIBRARIES} ${ROOT_LIBRARIES}
	FitGen experiment kinematics reconstructionfit
    )
endif()
