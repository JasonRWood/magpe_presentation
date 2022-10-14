# distutils: language = c++
from solvers cimport Quick_solver

cdef class PySolver:
    
    cdef Quick_solver cpp_solver
    
    def __cinit__(self):
        self.cpp_solver = Quick_solver()
        return
    
    def alpha_ad_dyn(self, float beta_max, float alpha_max, float b, float q, float d, float gamma, int seed, int alpha_init):
        
        self.cpp_solver.alpha_evo_only(beta_max, alpha_max, b, q, d, gamma, seed, alpha_init)
        return
