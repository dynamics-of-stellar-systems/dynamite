import dynamite_src as dyn
import dynamite_src.mges as mges
import dynamite_src.data as data
from dynamite_src.parameter_space import Parameter, StellarParameters, ConfigurationReaderYaml

# 1) define components of the physical system

# ... 1a) stars

stellar_mge = mges.MGE_from_file(
    mge_filename='NGC6278_rbandmge.txt'
    )

stellar_kinematics = data.GaussHermite(
    filename='gh_file_of_stars.txt'
    )

config = ConfigurationReaderYaml('config_example.yaml')
stellar_parameters = StellarParameters(name='q',**config.params['model_components']['stars']['parameters']['q'])
stellar_parameters.add(name='p',**config.params['model_components']['stars']['parameters']['p'])
stellar_parameters.add(name='u',**config.params['model_components']['stars']['parameters']['u'])

if (stellar_parameters.validate() is False):
    raise SystemExit('Not all stellar parameters q, p, and u specified')

# stellar_parameters = [
#     Parameter('ml', lo=1, hi=10, step=1, fixed=False),
#     Parameter('p', lo=0.2, hi=0.6, step=0.01, fixed=False),
#     Parameter('q', lo=0.2, hi=0.6, step=0.01, fixed=False),
# ]
raise SystemExit
stars = dyn.physical_system.VisibleComponent(
    contributes_to_potential=True,
    symmetry='triaxial',
    mge=stellar_mge,
    kinematic_data=[stellar_kinematics],
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

# 2) define components of the physical system

system = dyn.physical_system.System(
    stars,
    black_hole,
    dark_halo)

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

pargen0 = dyn.parameter_space.GridSearch(par_space)
parset_list0 = pargen0.generate(current_models=None)

# 5) set weight_solver

weight_solver = dyn.weight_solvers.NNLS(
    weight_solver_args={'regularisation':0.1}
    )
weight_solver.set_kinematics(system)

# 5) run models

schw_models = dyn.schwarzschild.SchwarzschildModelLoop(
    system=system,
    parset_list=parset_list0,
    weight_solver=weight_solver
)
schw_models.models

# 5) generate more parameters

parset_list1 = pargen0.generate(current_models=schw_models)

# 5) run more models

schw_models1 = dyn.schwarzschild.SchwarzschildModelLoop(
    system=system,
    parset_list=parset_list1,
    weight_solver=weight_solver
)
schw_models.models

# end
