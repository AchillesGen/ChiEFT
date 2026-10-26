"""The pyNUISANCE calls: published data tables and binned predictions.

Needs pyNUISANCE (NUISANCE3 with the NUISANCE2 record), i.e. the container of
docker/Dockerfile.
"""
from dataclasses import dataclass

import numpy as np

PB_TO_CM2 = 1e-36
_record = None


def analysis(name):
    global _record
    import pyNUISANCE as pn
    if _record is None:
        _record = pn.RecordFactory().make_record({"type": "nuisance2"})
    return _record.analysis(name)


def data_table(name):
    """(values, covariance, bin edges or None) of a published measurement."""
    import pyNUISANCE as pn
    a = analysis(name)
    binned = a.get_data()[0]
    values = np.asarray(binned.values, float).reshape(-1)
    errors = np.asarray(binned.errors, float).reshape(-1)
    cov = np.asarray(a.get_covariance_matrix(), float)
    if cov.shape != (values.size, values.size):
        cov = np.diag(errors ** 2)
    else:
        # NUISANCE gives the covariance in the sample's published units while the
        # values are absolute; recover the factor from the per-bin errors.
        sd = np.sqrt(np.clip(np.diag(cov), 0.0, None))
        ok = (sd > 0) & (errors > 0)
        if ok.any():
            cov = cov * float(np.median(errors[ok] / sd[ok])) ** 2
    try:
        edges = np.asarray(pn.Binning.get_bin_edges1D(binned.binning.bins), float)
    except Exception:
        edges = None
    return values, cov, edges


@dataclass
class Sample:
    """Selected events of one measurement: bin, weight (in the published
    cross-section units, bin width divided) and the extra columns."""
    bins: np.ndarray
    weights: np.ndarray
    nbins: int
    extra: dict


def histogram(path, names, divisor, add_columns, extra_names, block_size=50_000):
    """Bin every measurement in `names` in one pass over the event file `path`.

    divisor: per-atom cross section / divisor = the published convention.
    add_columns(fg) registers the `extra_names` columns, returned per event."""
    import pyNUISANCE as pn
    fg = pn.EventFrameGen(pn.EventSource(path), block_size)
    specs = []
    for i, name in enumerate(names):
        a = analysis(name)
        projs = a.get_projections()
        fg.add_int_column(f"sel{i}", a.get_selection().op)
        for j, p in enumerate(projs):
            fg.add_double_column(f"proj{i}_{j}", p.op)
        binned = a.get_data()[0]
        specs.append(dict(name=name, binning=binned.binning,
                          nbins=int(np.asarray(binned.values).size),
                          projs=[f"proj{i}_{j}" for j in range(len(projs))],
                          bins=[], weights=[], extra={n: [] for n in extra_names}))
    add_columns(fg)

    fatx = None
    block = fg.first()
    while block is not None:
        table = np.asarray(block.table, float)
        if table.size == 0:
            break
        cols = {n: i for i, n in enumerate(block.column_names)}
        fatx = float(table[-1, cols["fatx_per_sumw.pb_per_target.estimate"]])
        for i, s in enumerate(specs):
            p = s["projs"]
            idx = np.asarray(s["binning"].find_bins(block, p[0] if len(p) == 1 else p),
                             dtype=np.int64).reshape(-1)
            ok = (table[:, cols[f"sel{i}"]] != 0) & (idx >= 0) & (idx < s["nbins"])
            s["bins"].append(idx[ok])
            s["weights"].append(table[ok, cols["weight.cv"]])
            for n in extra_names:
                s["extra"][n].append(table[ok, cols[n]])
        block = fg.next()

    out = {}
    for s in specs:
        widths = np.asarray(list(s["binning"].bin_sizes()), float).reshape(-1)[:s["nbins"]]
        scale = np.divide(fatx * PB_TO_CM2 / divisor, widths, out=np.zeros(s["nbins"]),
                          where=widths != 0)
        bins = np.concatenate(s["bins"]).astype(int)
        out[s["name"]] = Sample(bins, np.concatenate(s["weights"]) * scale[bins], s["nbins"],
                                {n: np.concatenate(v) for n, v in s["extra"].items()})
    return out
