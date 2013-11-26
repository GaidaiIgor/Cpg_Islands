TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    template.cpp \
    main.cpp

QMAKE_CXXFLAGS += -msse3

OTHER_FILES += \
    input_sequence.txt \
    initial_probabilities.txt \
    transition_probabilities.txt \
    emission_probabilities.txt

CONFIG += console precompile_header

# Use Precompiled headers (PCH)
PRECOMPILED_HEADER  = precompiled.h

HEADERS += \
    sse_operator_traits.hpp \
    operator_traits.hpp \
    hmm.hpp \
    hmm_vector.hpp \
    hmm_table.hpp \
    hmm_matrix.hpp \
    float_traits.hpp \
    allocator_traits.hpp \
    precompiled.h
