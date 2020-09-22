import dynamite_src as dyn
from dynamite_src.parameter_space import Parameter

# 1) define components of the physical system

# ... 1a) stars

stellar_mge = dyn.mges.MGE_from_file(
    mge_filename='mge_file_of_stars.txt'
    )

# replace GH with BSplines...
# stellar_kinematics = dyn.data.GaussHermite(
#     filename='gh_file_of_stars.txt'
#     )
stellar_kinematics = dyn.data.BSplines(
    filename='gh_file_of_stars.txt'
    )

stellar_ages = dyn.data.PopulationMap(
    filename='age_map.txt'
    )

stellar_z = dyn.data.PopulationMap(
    filename='metallicity_map.txt'
    )

stellar_parameters = [
    Parameter('ml', lo=1, hi=10, step=1, fixed=False),
    Parameter('p', lo=0.2, hi=0.6, step=0.01, fixed=False),
    Parameter('q', lo=0.2, hi=0.6, step=0.01, fixed=False),
]

stars = dyn.physical_system.VisibleComponent(
    contributes_to_potential=True,
    symmetry='triaxial',
    mge=stellar_mge,
    kinematic_data=[stellar_kinematics],
    population_data=[stellar_ages, stellar_z],
    parameters=stellar_parameters
    )

# ... 1b) black hole

black_hole = dyn.physical_system.Plummer(
    Parameter('M', lo=1e-5, hi=1e-4, step=1e-5, fixed=False),
    Parameter('a', value=1e-6, fixed=True),
    )

# ... 1c) dark matter halo

dark_halo = dyn.physical_system.NFW(
    Parameter('M', lo=1, hi=10, step=1, fixed=False),
    Parameter('c', lo=5, hi=15, step=1, fixed=False),
    )

# ... 1d) globular clusters

gc_mge = dyn.mges.MGE_from_file(
    mge_filename='mge_file_of_gcs.txt'
    )

gc_kinematics = dyn.data.DiscreteLOS(
    filename='vLOS_file_of_gcs.txt'
    )

gc_z = dyn.data.DiscretePopulation(
    filename='gc_metallicities.txt'
    )

gcs = dyn.physical_system.VisibleComponent(
    contributes_to_potential=False,
    symmetry='spherical',
    mge=gc_mge,
    kinematic_data=[gc_kinematics],
    population_data=[gc_z],
    parameters=[] # spherical tracer population --> no parameters?
    )

# 2) define components of the physical system

system = dyn.physical_system.System(
    stars,
    black_hole,
    dark_halo,
    gcs
    )

system.n_pot
system.n_kin
system.n_pop
system.n_par

# 3) intialise parameter space

par_space = dyn.parameter_space.ParameterSpace(system)
par_space.n_par
par_space.n_par_fixed
par_space.n_par_free

# 4) generate parameters

pargen0 = dyn.parameter_space.GridWalk(par_space)
parset_list0 = pargen0.generate(current_models=None)

# 5) set weight_solver... but now use posterior solver

# weight_solver = dyn.weight_solvers.NNLS(
#     weight_solver_args={'regularisation':0.1}
#     )
# weight_solver.set_kinematics(system)

weight_solver = dyn.weight_solvers.LatentPosterior(
    weight_solver_args={'latent_dimension':10}
)
weight_solver.set_kinematics(system)

# 5) run models

schw_models = dyn.model.SchwarzschildModelLoop(
    system=system,
    parset_list=parset_list0,
    weight_solver=weight_solver
)
schw_models.models

# 5) generate more parameters... but now using a GPE

# parset_list1 = pargen0.generate(current_models=schw_models)
pargen1 = dyn.parameter_space.GaussianProcessEmulator(par_space)
parset_list1 = pargen1.generate(current_models=schw_models)

# 5) run more models

schw_models1 = dyn.model.SchwarzschildModelLoop(
    system=system,
    parset_list=parset_list1,
    weight_solver=weight_solver
)
schw_models.models


# end
