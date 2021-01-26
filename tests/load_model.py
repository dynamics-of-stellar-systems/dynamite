# coding: utf-8

import sys
sys.path.insert(0, '../')

import dynamite as dyn
print('Using DYNAMITE version:', dyn.__version__)
print('Located at:', dyn.__path__)

# read configuration
fname = 'test_reimplement_nnls.yaml'
c = dyn.config_reader.Configuration(fname, silent=True)
parset0 = c.all_models.get_parset_from_row(2)
mod0 = dyn.model.LegacySchwarzschildModel(system=c.system,
                                          settings=c.settings,
                                          parspace=c.parspace,
                                          parset=parset0)
orblib0 = mod0.get_orblib()
ws0 = mod0.get_weights(orblib0)
orbmat_no_error = ws0.read_nnls_orbmat_noerror()
kinematics = c.system.cmp_list[2].kinematic_data[0]
orb_gh0 = kinematics.transform_orblib_to_observables(orblib0)
tmp = kinematics.get_observed_values_and_uncertainties()
obs_gh, obs_gh_error = tmp
orbmat_with_error, rhs, solution = ws0.read_nnls_orbmat_rhs_and_solution()

weight_solver = dyn.weight_solvers.PrashsCoolNewWeightSolver(
    system=c.system,
    settings=c.settings.weight_solver_settings,
    directory_noml=mod0.directory_noml)

orbmat_prash, rhs_prash = weight_solver.construct_nnls_matrix_and_rhs(orblib0)

solution_prash = weight_solver.solve(orblib0)

def get_3d_density_for_orb(i):
    x = orbmat_no_error[i,1:361]
    y = np.ravel(orblib0.density_3D[i])
    x = x/np.sum(x)
    fig, ax = plt.subplots(2, 1, sharex=True)
    ax[0].plot(x)
    ax[0].plot(y)
    ax[1].plot(x/y)

def get_2d_masses_for_orb(i):
    x = orbmat_no_error[i,361:361+152]
    y = np.sum(orblib0.losvd_histograms.y, 1)[i]
    x = x/np.sum(x)
    y = y/np.sum(y)
    print(np.sum(x), np.sum(y))
    fig, ax = plt.subplots(2, 1, sharex=True)
    ax[0].plot(x)
    ax[0].plot(y)
    ax[1].plot(x/y)

def get_2d_masses_for_orb(i):
    x = orbmat_no_error[i,361:361+152]
    y = np.sum(orblib0.losvd_histograms.y, 1)[i]
    x = x/np.sum(x)
    y = y/np.sum(y)
    print(np.sum(x), np.sum(y))
    fig, ax = plt.subplots(2, 1, sharex=True)
    ax[0].plot(x)
    ax[0].plot(y)
    ax[1].plot(x/y)

def plot_gh_for_orbit(i):
    x = orbmat_no_error[i,361+152:]
    y = orb_gh0[i,:,:]
    y = np.ravel(y.T)
    fig, ax = plt.subplots(2, 1, sharex=True)
    ax[0].plot(x)
    ax[0].plot(y, ls='--')
    ax[1].plot(x-y)
    ax[1].set_ylim(-1,1)
    ax[0].set_ylim(-0.001,0.001)
    fig.tight_layout()

def compare_columns_of_orbmat_noerr(i):
    x = orbmat_no_error[i,:]
    y = orbmat_noerr_prash[i,:]
    fig, ax = plt.subplots(2, 1, sharex=True)
    ax[0].plot(x)
    ax[0].plot(y, ls='--')
    ax[1].plot((y-x)/x)
    for ax0 in ax:
        ax0.axvline(360, ls=':', color='k')
        for i in range(4):
            ax0.axvline(360+152*i, ls=':', color='k')
    fig.tight_layout()

def compare_orbmat(no_error=True):
    is_close = np.isclose(orbmat_with_error, orbmat_prash.T)
    fig, ax = plt.subplots(1, 1)
    col = ax.imshow(is_close,
                    interpolation='none',
                    aspect='auto',
                    cmap=plt.cm.Greys_r)
    ax.set_ylabel('orbit')
    ax.set_xlabel('constraint')
    fig.colorbar(col)
    fig.tight_layout()

def compare_rhs():
    if np.allclose(rhs, rhs_prash):
        print('EVERYTHING WORKS PRASHIN YOU CAN REST NOW')
    else:
        print("you're not there yet but still rest anyway it's late")
