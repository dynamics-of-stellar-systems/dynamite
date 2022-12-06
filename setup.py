import os
import setuptools
import numpy.distutils.core
import distutils.command.install_lib

# this loads the version number from the dynamite/version.py module
ver_file = open("dynamite/_version.py")
version = ver_file.readlines()[-1].split()[-1].strip("\"'")
ver_file.close()

# load the readme as long description
with open("README.md", "r") as fh:
    long_description = fh.read()

# load the package requirements from requirements.txt
with open("requirements.txt", "r") as fp:
    required = fp.read().splitlines()

ext = numpy.distutils.core.Extension(name='dynamite.pyfort_GaussHerm',
                 sources=['legacy_fortran/sub/gausherm.f'],
                 f2py_options=['--quiet'],
                )

legacy_fortran = [
    "../legacy_fortran/modelgen",
    "../legacy_fortran/orbitstart",
    "../legacy_fortran/orblib",
    "../legacy_fortran/orblib_new_mirror",
    "../legacy_fortran/partgen",
    "../legacy_fortran/triaxmass",
    "../legacy_fortran/triaxmassbin",
    "../legacy_fortran/triaxnnls_CRcut",
    "../legacy_fortran/triaxnnls_noCRcut",
]

class my_install_lib(distutils.command.install_lib.install_lib):
  def run(self):
    distutils.command.install_lib.install_lib.run(self)
    for fn in self.get_outputs():
      if any([fn.endswith(f) for f in legacy_fortran]):
        # copied from distutils source - make the binaries executable
        mode = ((os.stat(fn).st_mode) | 0o555) & 0o7777
        distutils.log.info("changing mode of %s to %o", fn, mode)
        os.chmod(fn, mode)

setuptools.setup(
    name="dynamite",
    version=version,
    author="Prashin Jethwa, Sabine Thater, Thomas Maindl",
    author_email="prashin.jethwa@univie.ac.at",
    description="dynamics, age and metallicity indicators tracing evolution",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://dynamics.univie.ac.at/dynamite_docs/index.html",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    project_urls={
        "Source": "https://github.com/dynamics-of-stellar-systems/dynamite/",
        "Documentation": "https://dynamics.univie.ac.at/dynamite_docs/index.html",
    },
    python_requires=">=3.7",
    # use the already parsed requirements from requirements.txt
    install_requires=required,
    # extra requirements for testing
    extras_require={
        "testing": [
            "pytest",
            "coverage",
        ]
    }
)

numpy.distutils.core.setup(
    cmdclass={'install_lib':my_install_lib},
    name="dynamite",
    version=version,
    author="Prashin Jethwa, Sabine Thater, Thomas Maindl",
    author_email="prashin.jethwa@univie.ac.at",
    description="dynamics, age and metallicity indicators tracing evolution",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://dynamics.univie.ac.at/dynamite_docs/index.html",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    project_urls={
        "Source": "https://github.com/dynamics-of-stellar-systems/dynamite/",
        "Documentation": "https://dynamics.univie.ac.at/dynamite_docs/index.html",
    },
    python_requires=">=3.7",
    # extra requirements for testing
    extras_require={
        "testing": [
            "pytest",
            "coverage",
        ]
    },
    package_data={
        "dynamite": legacy_fortran
    },
    ext_modules=[ext]
)
