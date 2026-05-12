#!/usr/bin/env python3
import sys
import itertools
import numpy as np
from copy import deepcopy as copy
from os.path import basename, isfile
from pyFIT3D.common.stats import pdl_stats, _STATS_POS
from pyFIT3D.common.constants import __FWHM_to_sigma__, _MODELS_ELINE_PAR
from pyFIT3D.common.io import ReadArguments, sel_waves, trim_waves, get_data_from_fits
from pyFIT3D.common.io import get_wave_from_header, print_time, read_masks_file, remove_isfile
from pyFIT3D.common.gas_tools import output_emission_lines_spectra, output_emission_lines_parameters, fit_elines_main
from pyFIT3D.common.gas_tools import ConfigEmissionModel, create_emission_lines_parameters, append_emission_lines_parameters, load_spec

np.set_printoptions(precision=4, suppress=True, linewidth=200)



def calc_chi_sq(f_obs, f_mod, ef_obs, ddof=0):
    """
    Calculates the Chi Square of a fitted model.

    Parameters
    ----------
    f_obs : array_like
        Observed spectrum
    f_mod : array_like
        Modeled spectrum
    ef_obs : array_like
        Error of observed spectrum.

    Returns
    -------
    float
        The Chi Square of the fit.
    int
        The number of observations.
    """
    mask = ef_obs != 0
    N_obs = mask.sum()
    chi = np.divide(
        f_obs - f_mod, ef_obs,
        where=mask,
        out=np.zeros_like(f_obs)
    )
    chi_sq = np.sum(chi**2)
    chi_sq_red = ((chi_sq / (N_obs - ddof)) if N_obs - ddof > 0 else chi_sq)
    return chi_sq_red, N_obs

class ReadArgumentsLocal(ReadArguments):
    """
    Argument parser for script that fits the emission-lines of a FITS cube.

    To add an argument to the script:
        Argument name in `__mandatory__` or `__optional__` list.
        Argument conversion function in `__mandatory_conv_func___` or `__optional_conv_func___` list.
        Argument default value (if not mandatory) in `__def_optional__`
    """
    # class static configuration:
    # arguments names and conversion string to number functions
    __script_name__ = basename(sys.argv[0])
    __mandatory__ = ['spec_file', 'config_file', 'mask_file', 'w_min', 'w_max', 'out_file', 'spatial_guide']
    __optional__ = ['central_coordinate', 'n_MC', 'n_loops', 'plot', 'scale_ini', 'prefix', 'memo', 'vel_map', 'sigma_map']
    __arg_names__ = __mandatory__ + __optional__
    __N_tot_args__ = len(__arg_names__)
    # default values of optional arguments with __optional__ as keys
    __def_optional__ = {
        'central_coordinate': '0,0', 'n_MC': '50', 'n_loops': '15', 'plot': '0', 'scale_ini': '0.15', 'memo': '0',
        # 'out_file': 'out.fit_spectra',
        # 'out_mod_res_final': __script_name__.replace('.py', '.out')
    }

    # parse functions
    __conv_func_mandatory__ = {'spec_file': str, 'config_file': str, 'mask_file': str, 'out_file': str, 'spatial_guide': str}
    __conv_func_optional__ = {'prefix': str, 'vel_map': str, 'sigma_map': str, 'central_coordinate': str}
    __conv_func__ = __conv_func_mandatory__.copy()
    __conv_func__.update(__conv_func_optional__)

    # usage message
    __usage_msg__ = 'usage: {} SPECTRUM_FILE[,ERROR_FILE] CONFIG MASK_LIST W_MIN W_MAX OUTFILE SPATIAL_GUIDE'.format(__script_name__)
    __usage_msg__ += ' [CENTRAL_COORDINATE={:s}]'.format(__def_optional__['central_coordinate'])
    __usage_msg__ += ' [N_MC={:d}]'.format(eval(__def_optional__['n_MC']))
    __usage_msg__ += ' [N_LOOPS={:d}]]'.format(eval(__def_optional__['n_loops']))
    __usage_msg__ += ' [PLOT] [SCALE_INI={:.2f}]'.format(eval(__def_optional__['scale_ini']))
    __usage_msg__ += ' [PREFIX] [MEMO=0/1]'
    __usage_msg__ += ' [input_vel_map.fits,input_vel_mask_now,FIXED=0,1]'
    __usage_msg__ += ' [input_disp_map.fits,FIXED=0,1 | FWHM=2.354*sigma]'
    # __usage_msg__ += ' [OUTPUT_FILENAME] [OUTPUT_MOD_RES_FILENAME]'

    def __init__(self, args_list=None, verbose=False):
        ReadArguments.__init__(self, args_list, verbose=verbose)
        self.redefine_max = 0
        self._parse_filenames()
        if self.central_coordinate is not None:
            _tmp = self.central_coordinate.split(',')
            # Integer as mandatory
            self.central_coordinate = [np.int(i) for i in _tmp]

    def _parse_filenames(self):
        self.error_file = None
        spec_file_args = self.spec_file.split(',')
        if len(spec_file_args) == 2:
            self.spec_file = spec_file_args[0]
            self.error_file = spec_file_args[1]
        self.vel_fixed = None
        self.guided = False
        if self.vel_map is not None:
            self.guided = True
            vel_map = self.vel_map.split(',')
            self.vel_map_file = vel_map[0]
            self.vel_mask_file = vel_map[1]
            self.vel_fixed = int(vel_map[2])
        self.sigma_fixed = None
        self.guided_sigma = False
        if self.sigma_map is not None:
            self.guided_sigma = True
            sigma_map = self.sigma_map.split(',')
            if isfile(sigma_map[0]):
                self.sigma_map_file = sigma_map[0]
            else:
                self.sigma_val = eval(sigma_map[0])*__FWHM_to_sigma__
            if len(sigma_map) > 1:
                self.sigma_fixed = int(sigma_map[1])

def parse_spatial_guide(spatial_guide, central_coordinate=None, nx=0, ny=0):
    spx_list = np.asarray(list(itertools.product(range(nx), range(ny))))
    if spatial_guide == 'from_center':
        x0, y0 = nx//2, ny//2
        if central_coordinate is not None:
            x0, y0 = central_coordinate
        ix, iy = np.indices((nx, ny))
        _r = np.sqrt((x0 - ix)**2 + (y0 - iy)**2)
        _r_list = np.asarray([_r[ixy] for ixy in itertools.product(range(nx), range(ny))])
        spx_list = spx_list[np.argsort(_r_list)]
    else:
        if spatial_guide is not None:
            if isfile(spatial_guide):
                spx_list = np.loadtxt(spatial_guide, delimiter=',')
            else:
                print(f'{basename}: {spatial_guide}: spatial guide file not found.')
    spx_list = [(ixy[0], ixy[1]) for ixy in spx_list]
    return spx_list

def kin_cube_elines_main(
    wavelength,
    cube_flux,
    config,
    out_file=None,
    cube_eflux=None,
    run_mode='RND',
    vel_map=None,
    vel_mask_map=None,
    vel_fixed=2,
    sigma_map=None,
    sigma_fixed=2,
    memo=False,
    n_MC=50,
    n_loops=15,
    scale_ini=0.15,
    redefine_max=True,
    max_factor=2,
    redefine_min=True,
    min_factor=0.012,
    plot=0,
    oversize_chi=False,
    spatial_guide=None,
    ):
    """
    The main function of the kin_cube_elines script. It will run the `fit_elines_main`
    over all observed spectra `cube_flux`, with shape ``(NW, NY, NX)``.
    All spectra has to be sampled over the same `wavelength` with shape ``(NW)``.

    Parameters
    ----------
    wavelength : array
        The observed wavelengths.

    cube_flux : array like
        The observed flux spectra with shape ``(NW, NY, NX)``.

    config : :class:`ConfigEmissionModel`
        The :class:`ConfigEmissionModel` instance that will configure the fit.

    out_file : str, optional
        The filename of the file where the output of :class:`EmissionLines` result will be recorded.

    cube_eflux : array like, optional
        The error in `flux` with shape ``(NW, NY, NX)``.

    run_mode : 'RND', 'LM', 'both', optional
        It will configure which :class:`EmissionLines` it will run.
        'RND': The RND algorithm of fit.
        'LM': Levemberg-Marquadt algorithm of fit.
        'both': It will run first the 'RND' method and inputs the `final_fit_params` of 'RND' at the 'LM' method.
        Defaults to RND method.

    vel_map : array like, optional
        Sets up the guide velocity map. (see: `guide_vel`).

    vel_mask_map : array like, optional
        Sets up the mask map of the guide velocity. (see: `guide_vel`).

    vel_fixed : int or boolean, optional
        Sets up the to_fit parameter of the guide velocity. (see: `guide_vel`).

    sigma_map : array like, optional
        Sets up the guide sigma map. (see: `guide_sigma`).

    sigma_fixed : int or boolean
        Sets up the to_fit parameter of the guide sigma map. (see: `guide_sigma`).

    memo : boolean, optional
        While fitting the cube, memorizes the last result and inputs to the next spectrum
        in the loop.

    n_MC : int, optional
        Number of Monte-Carlo iterations in the fit. Defaults to 50.

    n_loops : int, optional
        Number of loops of the fit. Defaults to 15.

    scale_ini : float, optional
        The scale_ini parameter of fit. Defaults to 0.15

    redefine_max : int or boolean, optional
        The redefine_max parameter of the fit. Defaults to True.

    max_factor : float, optional
        Set up the max flux factor when redefine_max is on. Defaults to 2.

    redefine_min : boolean, optional
        If redefine_max is on, it will also redefine the min flux. Defaults
        to True.

    min_factor : float, optional
        Set up the min flux factor when redefine_max is on. Defaults to 0.012.

    plot : int, optional
        If 1 it will plot fit. If 2 it will plot to a file. None it will be set as 0, i.e.
        do not plot.

    oversize_chi : bool, optional
        if True, will triplify the first chi of each MC round in RND mode. Defaults to False.

    spatial_guide : list of tuples, optional
        A list of pairs of coordinates (X, Y) guiding the spatial order of the fit.

    Returns
    -------
    tuple
        Constructed as::

            list of array like
                Observed spectra, model spectra and residual spectra, all with same shape as `cube_flux`.

            dict
                Output maps of fitted parameters
                (see: `create_emission_lines_parameters`, `append_emission_lines_parameters`
                and `output_emission_lines_parameters`).

    """
    nw, ny, nx = cube_flux.shape
    out_file = basename(sys.argv[0]).replace('.py', '.out') if out_file is None else out_file

    # create models output
    output_el_models = create_emission_lines_parameters(config, shape=(ny, nx))
    # create spectra output
    org__wyx = np.zeros_like(cube_flux)
    mod__wyx = np.zeros_like(cube_flux)
    res__wyx = np.zeros_like(cube_flux)

    #for ixy in itertools.product(range(nx), range(ny)):
    first_spaxel = spatial_guide[0]
    for ixy in spatial_guide:
        ix, iy = ixy
        print(f'# ID {ix}/{nx - 1},{iy}/{ny - 1}')
        current_spaxel = ixy
        current_guide_vel = None if vel_map is None else vel_map[iy, ix]
        current_guide_vel_mask = None if vel_mask_map is None else vel_mask_map[iy, ix]
        current_guide_sigma = None if sigma_map is None else sigma_map[iy, ix]
        flux = cube_flux[:, iy, ix]
        # setbadtoval(0)
        flux[~np.isfinite(flux)] = 0
        st_f = pdl_stats(flux)
        st_f_median = np.abs(st_f[_STATS_POS['median']])
        st_f_sigma = np.abs(st_f[_STATS_POS['pRMS']])
        if cube_eflux is not None:
            eflux = cube_eflux[:, iy, ix]
            eflux[~np.isfinite(eflux)] = 0
        else:
            eflux = np.ones_like(flux)*0.1
            if (st_f_median != 0) and (st_f_sigma != 0):
                eflux *= (st_f_median + st_f_sigma**2)
            # setvaltobad(0)
            eflux[eflux == 0] = np.nan
            if st_f_median != 0:
                eflux[np.isnan(eflux)] = st_f_median
            else:
                eflux[np.isnan(eflux)] = 0.001
        # enter if all flux is zero or if is in the first spaxel
        # if not ((ix or iy) and flux.any()):
        cf = copy(config)
        if (ixy[0] == first_spaxel[0]) and (ixy[1] == first_spaxel[1]):
            print('first spaxel')
        else:
            if flux.any():
                last_chi_sq, _ = calc_chi_sq(EL.flux, EL.model, EL.sigma_flux,
                                             EL.config.n_models + 1)
                # Memorize last result to input the next spectrum
                # if memo and (last_chi_sq < 2*EL.config.chi_goal):

                if memo and (np.abs(last_chi_sq - EL.config.chi_goal) < 0.5):
                    print('memo config')

                    # guess = last_spaxel_values*(1 +/- 0.15)
                    pars = np.asarray(EL.final_fit_params_mean)
                    pars_0 = copy(pars)
                    pars_1 = copy(pars)
                    sel_model = (np.asarray(cf.model_types) == 'eline')
                    delta = pars*0.15
                    delta[:, _MODELS_ELINE_PAR['v0']] *= 0.5
                    pars_0[sel_model] = pars[sel_model] - delta[sel_model]
                    pars_1[sel_model] = pars[sel_model] + delta[sel_model]
                    cf.update_ranges(min_values=pars_0, max_values=pars_1)

                    cf.guess = cf._set_linkings(EL.final_fit_params_mean)
                    cf._correct_guess()

                    # cf.update_config(EL.final_fit_params_mean, update_ranges=True, frac_range=0.15)
            #     else:
            #         cf = copy(config)
            # else:
            #     cf = copy(config)

        flux_max = np.abs(st_f[_STATS_POS['max']] - st_f[_STATS_POS['median']])
        EL = fit_elines_main(wavelength, flux, eflux, cf,
                             vel_guide=current_guide_vel, sigma_guide=current_guide_sigma,
                             vel_fixed=vel_fixed, sigma_fixed=sigma_fixed,
                             vel_mask=current_guide_vel_mask,
                             n_MC=n_MC, plot=plot, n_loops=n_loops,
                             scale_ini=scale_ini, fine_search=True,
                             redefine_max=redefine_max, redefine_min=redefine_min,
                             max_factor=max_factor, min_factor=min_factor, flux_max=flux_max,
                             run_mode=run_mode, randomize_flux=True, oversize_chi=oversize_chi)
        EL.output(filename_output=out_file, append_output=True, spec_id=ixy)
        org__wyx[:, iy, ix] = EL.flux
        mod__wyx[:, iy, ix] = EL.model
        res__wyx[:, iy, ix] = EL.flux - EL.model
        append_emission_lines_parameters(EL, output_el_models, ixy)
        # ixy, output_el_models)
    output_el_spectra = np.array([org__wyx, mod__wyx, res__wyx])
    return output_el_spectra, output_el_models

def kin_cube_elines(
    spec_file,
    config_file,
    w_min,
    w_max,
    out_file=None,
    mask_file=None,
    error_file=None,
    prefix=None,
    run_mode='RND',
    memo=False,
    vel_map_file=None,
    vel_mask_file=None,
    vel_fixed=2,
    sigma_map_file=None,
    sigma_fixed=2,
    sigma_val=None,
    n_MC=50,
    n_loops=15,
    scale_ini=0.15,
    redefine_max=True,
    max_factor=2,
    redefine_min=False,
    min_factor=0.012,
    plot=0,
    seed=None,
    spatial_guide=None,
    central_coordinate=None,
    ):
    """
    The kin_cube_elines script (i.e. a well prepared wrap of `kin_cube_elines_main`).

    It will run the `run_mode` fit the the emission models defined by the `config_file`
    on the cube of spectra in `spec_file` masked by `mask_file`, inside the define wavelength
    range (`w_min`, `w_max`).

    If `run_mode` is::

        'RND': mimics fit_elines_rnd script.
        'LM': mimics fit_elines_LM script.
        'both': It will run the 'RND' first and inputs the results of the fit in the 'LM' fit.


    Parameters
    ----------
    spec_file : str
        The filename of the FITS file containing the wavelengths and observed flux spectra.

    config_file : str
        The filename of the :class:`ConfigEmissionModel` input file.

    w_min : int
        The minimun (bluer) wavelength which defines the fitted wavelength range.

    w_max : int
        The maximum (redder) wavelength which defines the fitted wavelength range.

    out_file : str, optional
        The name of the output results of the fit. Defaults to
        ``basename(sys.argv[0]).replace('py', 'out')``.

    mask_file : str, optional
        If mask_file is a valid file it will mask the wavelengths.

    error_file : str
        The filename of the FITS file containing the errors in observed flux spectra.

    prefix : str, optional
        It will define the prefix used in `output_emission_lines_parameters` and
        `output_emission_lines_spectra`. Defaults to basename of `out_file`.

    run_mode : str {'RND', 'LM', 'both'}, optional
        It will configure which :class:`EmissionLines` it will run.
        'RND': The RND algorithm of fit.
        'LM': Levemberg-Marquadt algorithm of fit.
        'both': It will run first the 'RND' method and inputs the `final_fit_params` of
        'RND' at the 'LM' method.
        Defaults to the RND method.

    memo : boolean, optional
        While fitting the cube, memorizes the last result and inputs to the next spectrum
        in the loop.

    vel_map_file : str, optional
        The filename of the guide velocity map. (see: `guide_vel`).

    vel_mask_file : array like, optional
        The filename of the mask map of the guide velocity. (see: `guide_vel`).

    vel_fixed : int or boolean, optional
        Sets up the to_fit parameter of the guide velocity. (see: `guide_vel`).

    sigma_map_file : array like, optional
        The filename of the guide sigma map. (see: `guide_sigma`).

    sigma_fixed : int or boolean
        Sets up the to_fit parameter of the guide sigma map. (see: `guide_sigma`).

    sigma_val : float, optional
        Rewrites the sigma_map using a single value for the guide sigma map. (see: `guide_sigma`)

    n_MC : int, optional
        Number of Monte-Carlo iterations in the fit. Defaults to 50.

    n_loops : int, optional
        Number of loops of the fit. Defaults to 15.

    scale_ini : float, optional
        The scale_ini parameter of fit. Defaults to 0.15

    redefine_max : int or boolean, optional
        The redefine_max parameter of the fit. Defaults to True.

    max_factor : float, optional
        Set up the max flux factor when redefine_max is on. Defaults to 2.

    redefine_min : boolean, optional
        If redefine_max is on, it will also redefine the min flux. Defaults
        to False.

    min_factor : float, optional
        Set up the min flux factor when redefine_max is on. Defaults to 0.012.

    plot : int, optional
        If 1 it will plot fit. If 2 it will plot to a file. None it will be set as 0, i.e.
        do not plot.

    seed : int, optional
        It will define the input seed. Defaults to the ``int(time.time())``.

    spatial_guide : str, optional
        Configures the spatial order of the fit of cube spectra.
        If `spatial_guide` is `from_center` the analysis will start from `central_coordinate`
        and will follow the spaxels ordered by the distance from the center (NOTE: the option `central_coordinate` should be also input). Otherwise the function will understand
        `spatial_guide` as a txt file containing a list of spaxels coordinates (X,Y), e.g.::

            0,0
            1,1
            4,2
            7,3
            ...
            25,23
            25,24

    central_coordinate : tuple, optional
        The central spaxel coordinate (X, Y).

    See also
    --------
    :class:`ConfigEmissionModel`, `fit_elines_main`, `kin_cube_elines_main`, `output_emission_lines_parameters`, `output_emission_lines_spectra`

    """
    seed = print_time() if seed is None else print_time(time_ini=seed)
    np.random.seed(seed)
    time_ini_run = print_time(print_seed=False, get_time_only=True)
    memo = False if memo is None else bool(memo)
    vel_fixed = 2 if vel_fixed is None else vel_fixed
    sigma_fixed = 2 if sigma_fixed is None else sigma_fixed
    out_file = basename(sys.argv[0]).replace('.py', '.out') if out_file is None else out_file
    if prefix is None:
        _s = out_file.split('.')
        prefix = '.'.join(_s[0:-1]) if len(_s) > 1 else f'{out_file}'
    remove_isfile(out_file)
    _runs = ['RND', 'LM', 'both']
    if (run_mode is None) or (run_mode not in _runs):
        run_mode = 'RND'
    if vel_map_file is None or not isfile(vel_map_file):
        vel_fixed = 2
    if sigma_map_file is None or not isfile(sigma_map_file):
        sigma_fixed = 2
    # read masks file
    masks, n_masks = read_masks_file(mask_file)
    print(f'{n_masks} regions to mask')
    error_extension = 0
    if error_file is not None:
        _tmp = error_file.split('[')
        if len(_tmp) > 1:
            try:
                error_file = _tmp[0]
                error_extension = int(_tmp[-1].split(']')[0])
            except:
                print(f'kin_cube_elines: error extension not found, using {error_extension}')
    wave__w, flux__wyx, header, eflux__wyx = load_spec(filename=spec_file,
                                                       error_filename=error_file,
                                                       error_extension=error_extension)

    nW, ny, nx = flux__wyx.shape
    spx_list = parse_spatial_guide(spatial_guide, central_coordinate=central_coordinate, nx=nx, ny=ny)
    sel_wl_range = trim_waves(wave__w, [w_min, w_max])
    sel_wl_goods = sel_waves(masks, wave__w)
    sel_wl = sel_wl_range & sel_wl_goods
    org_cf = ConfigEmissionModel(config_file)

    # EL: 2020-09-04: bugfix - If there is no avaible spectral information for the fit
    # This 0 could be changed to another value for a better fit.
    if sel_wl.astype('int').sum() < 2:
        print('[kin_cube_elines]: n_wavelength < 2: No avaible spectra to perform the configured analysis.')
        output_el_models = create_emission_lines_parameters(org_cf, shape=(ny, nx))
    else:
        wave_msk__w = wave__w[sel_wl]
        flux_msk__wyx = flux__wyx[sel_wl]
        if eflux__wyx is not None:
            eflux_msk__wyx = eflux__wyx[sel_wl]
        else:
            eflux_msk__wyx = None
        nw, ny, nx = flux_msk__wyx.shape
        # read config FOR TESTING PURPOSES -------------------------------------------------------------
        vel_map__yx = None
        if (vel_map_file is not None) and (isfile(vel_map_file)):
            vel_map__yx = get_data_from_fits(vel_map_file)
        vel_mask_map__yx = None
        if (vel_mask_file is not None) and (isfile(vel_mask_file)):
            vel_mask_map__yx = get_data_from_fits(vel_mask_file)
        sigma_map__yx = None
        if (sigma_map_file is not None) and (isfile(sigma_map_file)):
            sigma_map__yx = get_data_from_fits(sigma_map_file)*__FWHM_to_sigma__
        if sigma_val is not None:
            sigma_map__yx = np.ones((ny, nx), dtype='float')*sigma_val
        r = kin_cube_elines_main(
            wave_msk__w, flux_msk__wyx, org_cf, out_file,
            cube_eflux=eflux_msk__wyx,
            vel_map=vel_map__yx, vel_mask_map=vel_mask_map__yx, vel_fixed=vel_fixed,
            sigma_map=sigma_map__yx, sigma_fixed=sigma_fixed,
            memo=memo, n_MC=n_MC, n_loops=n_loops, redefine_max=redefine_max,
            scale_ini=scale_ini, run_mode=run_mode,
            plot=plot, spatial_guide=spx_list,
        )
        output_el_spectra, output_el_models = r
        output_emission_lines_spectra(wave_msk__w, output_el_spectra, header, prefix)
    output_emission_lines_parameters(prefix, org_cf, output_el_models)
    time_end = print_time(print_seed=False)
    time_total = time_end - time_ini_run
    print(f'# SECONDS = {time_total}')

if __name__ == "__main__":
    # parse arguments
    pa = ReadArgumentsLocal()
    kin_cube_elines(spec_file=pa.spec_file, config_file=pa.config_file,
                    w_min=pa.w_min, w_max=pa.w_max, out_file=pa.out_file,
                    mask_file=pa.mask_file, error_file=pa.error_file,
                    n_MC=pa.n_MC, n_loops=pa.n_loops, scale_ini=pa.scale_ini,
                    prefix=pa.prefix, memo=pa.memo, plot=pa.plot,
                    vel_map_file=pa.vel_map_file, vel_mask_file=pa.vel_mask_file,
                    sigma_map_file=pa.sigma_map_file, sigma_val=pa.sigma_val,
                    sigma_fixed=pa.sigma_fixed, vel_fixed=pa.vel_fixed,
                    spatial_guide=pa.spatial_guide, central_coordinate=pa.central_coordinate,
                    run_mode='both', seed=None)
