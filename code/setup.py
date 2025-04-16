from setuptools import setup
from Cython.Build import cythonize
import numpy as np
from setuptools.extension import Extension

# Define extensions with OpenMP flags
extensions = [
    Extension(
        name="interpolation",
        sources=["interpolation.pyx"],
        extra_compile_args=["-fopenmp"],
        extra_link_args=["-fopenmp"],
        include_dirs=[np.get_include()],
    ),
    Extension(
        name="data_importer",
        sources=["data_importer.pyx"],
        extra_compile_args=["-fopenmp"],
        extra_link_args=["-fopenmp"],
        include_dirs=[np.get_include()],
    ),
    Extension(
        name="Fisher_calc",
        sources=["Fisher_calc.pyx"],
        extra_compile_args=["-fopenmp"],
        extra_link_args=["-fopenmp"],
        include_dirs=[np.get_include()],
    ),
    # Extension(
    #     name="Fisher_calc2",
    #     sources=["Fisher_calc2.pyx"],
    #     extra_compile_args=["-fopenmp"],
    #     extra_link_args=["-fopenmp"],
    #     include_dirs=[np.get_include()],
    # ),
    # Extension(
    #     name="Fisher_calc3",
    #     sources=["Fisher_calc3.pyx"],
    #     extra_compile_args=["-fopenmp"],
    #     extra_link_args=["-fopenmp"],
    #     include_dirs=[np.get_include()],
    # ),
    Extension(
    name="Fisher_calc_python_imp",
        sources=["Fisher_calc_python_imp.pyx"],
        extra_compile_args=["-fopenmp"],
        extra_link_args=["-fopenmp"],
        include_dirs=[np.get_include()],
    ),
    Extension(
    name="cmb_noise_fast",
        sources=["cmb_noise_fast.pyx"],
        include_dirs=[np.get_include()],
    )
]

# Setup
setup(
    ext_modules=cythonize(
        extensions,
        compiler_directives={
            "nonecheck": False,
            "boundscheck": False,
            "wraparound": False,
            "initializedcheck": False,
            "cdivision": True,
        },
        language_level=3,
    ),
)
