cdef extern from "wrapper.cpp":
    pass

cdef extern from "wrapper.h" namespace "solvers":
    cdef cppclass Quick_solver:
        Quick_solver() except +
        Quick_solver(float, float, float, float, float, float, int) except +
        void alpha_evo_only(float, float, float, float, float, float, int, int)