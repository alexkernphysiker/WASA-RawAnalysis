cmake_minimum_required(VERSION 3.0)
project(WASAatCOSY_etamesic_3He)

set(CMAKE_CPP_FLAGS "-DNoMySQL -DNoXerces -pthread -std=c++11 -Wno-deprecated-declarations -m64 -DNoQT")
set(CMAKE_CXX_FLAGS "-g -O2 -fPIC -fno-strict-aliasing -Wall -fno-strict-aliasing -std=c++11")
set(${CMAKE_LINK_LIBRARY_FLAG} -pthread)

add_subdirectory(WASA-libs)
include_directories(${WASA_LIBS_INC})
link_directories(${ROOT_LIBRARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/WASA-RootSorterWrap/DataFiles DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/WASA-RootSorterWrap)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/WASA-libs/RunPluto.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/WASA-libs)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/scripts/config DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/scripts/WMC DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/scripts/RunAnalysis-data.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/scripts/RunAnalysis-mc.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/scripts/RunAnalysis-Reconstruction.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/scripts/RunClear.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/scripts/RunWMC.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/rootlogon.C DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/.sorterrc DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

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

    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/WASA-RootSorterWrap)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/Reactions)

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
    file(GLOB REACTIONS Reactions/*.cpp)
    add_library(reactions SHARED ${REACTIONS})
    target_link_libraries(reactions 
	sorterwrap ${WASA_LIBRARIES} ${ROOT_LIBRARIES}
    )
    add_executable(rawanalysis rawanalysis.cpp)
    target_link_libraries(rawanalysis 
	sorterwrap reactions
	${WASA_LIBRARIES} ${ROOT_LIBRARIES}
	FitGen experiment kinematics reconstructionfit
    )
endif()