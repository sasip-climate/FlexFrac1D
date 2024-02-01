#!/usr/bin/env python3

"""Basic value tests on the attributes and methods of a Floe instance"""

import hashlib
import itertools
import numpy as np
import polars as pl
import polars.testing as pltesting
import pytest
import scipy.sparse as sparse

from flexfrac1d.wave import Wave
from flexfrac1d.ice import Floe


SRC_TARGET = "gen_end_to_end/floe"
BIG_NUM = 2**64
DUM_WAVE = Wave(1, 100)


def get_floe(comb, wave=None):
    thickness, left_edge, length, disp_rel, dx = comb
    kwargs = dict()
    if disp_rel is not None:
        kwargs["dispersion"] = disp_rel
    if dx is not None:
        kwargs["dx"] = dx

    floe = Floe(thickness, left_edge, length, **kwargs)
    if wave is not None:
        floe.setWPars(wave)
    return floe


def sub_dict(dct, keys):
    return {k: dct[k] for k in keys}


def get_digest(comb):
    h = hashlib.shake_256()
    for _t in comb:
        h.update(str(_t).encode("utf-8"))
    return h.hexdigest(16)


def gen_method_arguments():
    method_arguments = {
        "kw": (0, 0.6, 1),  # 0: default parameter
        "iFracs": (20, (10, 15, 21)),
        "wave": (Wave(1, 50), Wave(0.6, 60)),
        "t": (465, 987, 3654),
        "EType": ("Disp", "Flex"),
        "verbose": (False,),
        "recompute": (False,),
        "istart": (5, 9),
        "iend": (13, 19),
        "multiFrac": (False, True),
        "maxFracs": (1, 2, 3),
        "x_fracs": (46, 49.3, (60, 70.1)),
    }

    return method_arguments


method_arguments = gen_method_arguments()
# FindE_min and FindE_minVerbose both call computeEnergyIfFrac.
# In the latter, there seems to be a call to the floe.a0 before it is
# initialised. These three methods thus systematically fail, and might be
# broken at least for isolated Floe objects.
# Idem for computeEnergySubFloe, called by FindE_min.
methods = {
    "calc_alpha": sub_dict(method_arguments, ("kw",)),
    "calc_curv": dict(),
    "calc_du": dict(),
    "calc_Eel": sub_dict(method_arguments, ("EType",)),
    "calc_strain": dict(),
    "calc_w": dict(),
    "fracture": sub_dict(method_arguments, ("x_fracs",)),
    "mslf_int": dict(),
}

constructor_arguments = {
    "thickness": (1, 1.1),
    "left_edge": (-10, 40, 45.2),
    "length": (100, 120.1),
    "dispersion": (None, "Open", "ML", "El", "ElML"),
    "dx": (None, 0.1, 1),
}
constructor_combinations = list(itertools.product(*constructor_arguments.values()))

old_to_new_map = {
    "h": "thickness",
    "x0": "left_edge",
    "L": "length",
    "dx": "dx",
    "hw": "draft",
    "ha": "freeboard",
    "I": "quad_moment",
    "k": "frac_toughness",
    "E": "youngs_modulus",
    "v": "poissons_ratio",
    "DispType": "dispersion",
}
new_to_old_map = {v: k for k, v in old_to_new_map.items()}


def test_attributes():
    attributes = [
        "h",
        "x0",
        "L",
        "dx",
        "hw",
        "ha",
        "I",
        "k",
        "E",
        "v",
        "DispType",
        "xF",
        "A",
    ]
    n_combinations = len(constructor_combinations)
    df_dict = dict(zip(constructor_arguments.keys(), zip(*constructor_combinations)))
    problematic_attributes = ("DispType", "xF", "A")
    df_dict["DispType"] = n_combinations * [None]
    df_dict |= {
        att: np.full(n_combinations, np.nan)
        for att in attributes
        if att not in problematic_attributes
    }

    for i, comb in enumerate(constructor_combinations):
        _h = get_digest(comb)
        floe = get_floe(comb)

        for att in attributes:
            if att not in problematic_attributes or att == "DispType":
                df_dict[att][i] = getattr(floe, old_to_new_map[att])
            else:
                if att == "xF":
                    truth = np.loadtxt(f"{SRC_TARGET}/att_{att}_{_h}.csv")
                    assert np.allclose(truth, floe.xF)
                else:
                    is_sp = floe.xF.size > 250
                    if is_sp:
                        arr = floe.disp_matrix_sparse
                        truth = sparse.load_npz(f"{SRC_TARGET}/att_{att}_{_h}.npz")
                        assert (
                            np.all(truth.indptr == arr.indptr)
                            and np.all(truth.indices == arr.indices)
                            and np.allclose(truth.data, arr.data)
                        )
                    else:
                        arr = floe.disp_matrix
                        truth = np.loadtxt(f"{SRC_TARGET}/att_{att}_{_h}.csv")
                        assert np.allclose(truth, arr)
    df = pl.from_dict(df_dict)

    df_src = pl.read_parquet(f"{SRC_TARGET}/attributes_reference.parquet")
    pltesting.assert_frame_equal(df_src, df)


@pytest.mark.parametrize("method, args", methods.items())
def test_methods(method, args):
    def init_eel(floe, wvf, _args):
        floe.calc_w(wvf)
        _kwargs = dict(zip(("wvf", "EType"), (wvf, _args[-1])))
        floe.calc_Eel(**_kwargs)
        return floe

    n_combinations = len(constructor_combinations)
    argument_combinations = list(itertools.product(*args.values()))
    n_met_args = len(argument_combinations)
    df_keys = sum(
        map(lambda _d: tuple(_d.keys()), (constructor_arguments, args)), tuple()
    ) + (method,)
    df_dict = {k: [None] * n_combinations * len(argument_combinations) for k in df_keys}

    # iteration on different instances
    for i, cstr_args in enumerate(constructor_combinations):
        floe = get_floe(cstr_args, DUM_WAVE)
        wvf = 0.1 * np.sin(0.4 * floe.xF + 0.36)  # arbitrary "wave realisation"
        floe.calc_w(wvf)

        # iteration on the combinations of method arguments
        for j, _args in enumerate(argument_combinations):
            idx = i * n_met_args + j
            # populate the dataframe keys: constructor and method parameters
            for k, v in zip(constructor_arguments, cstr_args):
                df_dict[k][idx] = v
            for k, v in zip(args, _args):
                if k == "x_fracs":
                    # cast to list necessary for Polars not to scream
                    df_dict[k][idx] = list(np.atleast_1d(v))
                else:
                    df_dict[k][idx] = v

            # populate the dataframe values: method output
            try:
                df_dict[method][idx] = getattr(floe, method)(*_args)
            except AttributeError as e:
                if method in ("calc_curv", "calc_du", "calc_strain"):
                    floe.calc_w(wvf)
                    df_dict[method][idx] = getattr(floe, method)()
                    if method == "calc_strain":
                        # call to getattr first so the attribute is created
                        df_dict[method][idx] = floe.strain
                else:
                    raise e
            except TypeError as e:
                if method == "calc_w":
                    floe.calc_w(wvf)
                    df_dict[method][idx] = floe.w
                elif method == "mslf_int":
                    df_dict[method][idx] = getattr(floe, method)(wvf)
                elif method == "calc_Eel":
                    floe = init_eel(floe, wvf, _args)
                    df_dict[method][idx] = floe.Eel
                elif method == "FindE_minVerbose":
                    floe = init_eel(floe, wvf, _args)
                    df_dict[method][idx] = getattr(floe, method)(*_args[:-1])[2]
                else:
                    raise e
            if method == "fracture":
                df_dict[method][idx] = [
                    [floe.thickness, floe.left_edge, floe.length]
                    for floe in df_dict[method][idx]
                ]

    df = pl.from_dict(df_dict)
    df_src = pl.read_parquet(f"{SRC_TARGET}/" f"met_{method}_reference.parquet")
    pltesting.assert_frame_equal(df_src, df)
