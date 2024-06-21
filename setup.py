import os
import setuptools

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

legacy_fortran = [
    "../legacy_fortran/orbitstart",
    "../legacy_fortran/orbitstart_bar",
    "../legacy_fortran/orblib",
    "../legacy_fortran/orblib_new_mirror",
    "../legacy_fortran/orblib_bar",
    "../legacy_fortran/partgen",
    "../legacy_fortran/triaxmass",
    "../legacy_fortran/triaxmass_bar",
    "../legacy_fortran/triaxmassbin",
    "../legacy_fortran/triaxmassbin_bar"
]
additional_ex = ["legacy_fortran/modelgen",
                 "legacy_fortran/triaxnnls_CRcut",
                 "legacy_fortran/triaxnnls_noCRcut",
                 "legacy_fortran/triaxnnls_bar"]
legacy_fortran.extend([f'../{e}' for e in additional_ex if os.path.isfile(e)])

setuptools.setup(
    name="dynamite",
    version=version,
    author="Prashin Jethwa, Sabine Thater, Thomas Maindl, Alice Zocchi, Edward Lilley",
    author_email="prashin.jethwa@univie.ac.at",
    description="dynamics, age and metallicity indicators tracing evolution",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://dynamics.univie.ac.at/dynamite_docs/index.html",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
    ],
    project_urls={
        "Source": "https://github.com/dynamics-of-stellar-systems/dynamite/",
        "Documentation": "https://dynamics.univie.ac.at/dynamite_docs/index.html",
    },
    python_requires=">=3.8",
    # use the already parsed requirements from requirements.txt
    install_requires=required,
    package_data={
        "dynamite": legacy_fortran
    },
    # extra requirements for testing
    extras_require={
        "cvxopt":
            "cvxopt>=1.2.6",
        "testing": [
            "pytest",
            "coverage",
        ]
    }
)
