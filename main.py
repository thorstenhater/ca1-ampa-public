#!/usr/bin/env python3

import arbor as A
from arbor import units as U
import subprocess as sp
from pathlib import Path
from time import perf_counter as pc
import polars as pl
from argparse import ArgumentParser
import numpy as np

here = Path(__file__).parent

args = ArgumentParser(
    "ca1-ampa",
    description="Simulate Hay 2011 model with AMPA synapses",
)

args.add_argument("-T", "--final-time", default=1000, help="Stop time [ms]", type=float)
args.add_argument(
    "-t", "--time-step", default=0.0025, help="Time step [ms]", type=float
)
args.add_argument(
    "-c",
    "--correlated-noise",
    default=False,
    help="Correlated noise on all synapses.",
    action="store_true",
)

args.add_argument(
    "-g",
    "--max-conductance",
    default=0.0035,
    help="Peak conductance.",
    type=float,
)

args.add_argument(
    "-w",
    "--weight",
    default=1.0,
    help="Weight.",
    type=float,
)


args.add_argument(
    "-s",
    "--synapse-sites",
    default="terminals",
    choices=["sphere", "random5"],
    help="""Where to place synapses""",
)


opts = args.parse_args()
g = opts.max_conductance
correlate = opts.correlated_noise
syns = opts.synapse_sites
weight = opts.weight

dn = here / f"results-syns={syns}-weight={weight}-correlate={correlate}-g={g}"
dn.mkdir(exist_ok=True)

n_syns = 10

ends = []
if syns == "sphere":
    ids = [532, 302, 301, 57, 3595, 529, 707, 3651, 775, 86]
    ends = ids[:n_syns]
elif syns == "random5":
    ids = [3668, 3201,  541, 3066, 2958,  675, 2315,  997, 2254,  401]
    ends = ids[:n_syns]
else:
    raise RuntimeError(f"Unknown synapse site choice {syns}.")

def make_cell(g, ends):
    cid = "L5PC"
    mrf = "morphology_L5PC"
    nml = A.neuroml(f"{here}/mrf/{mrf}.nml").morphology(mrf, allow_spherical_root=True)
    lbl = A.label_dict()
    lbl.append(nml.metadata.groups())
    lbl["all"] = "(all)"
    syn = A.synapse("local_exp2syn", g=g)
    dec = A.load_component(f"{here}/acc/{cid}.acc").component
    dec.place('(location 0 0.5)', A.threshold_detector(-10*U.mV), 'det')
    for end in ends:
        dec.place(f"(on-components 0.5 (segment {end}))", syn, f"syn-{end}")

    return A.cable_cell(nml.morphology, dec, lbl, A.cv_policy_fixed_per_branch(5))


cell = make_cell(g, ends)

def compile(fn, cat):
    fn = fn.resolve()
    cat = cat.resolve()
    recompile = False
    if fn.exists():
        for src in cat.glob("*.mod"):
            src = Path(src).resolve()
            if src.stat().st_mtime > fn.stat().st_mtime:
                recompile = True
                break
    else:
        recompile = True
    if recompile:
        sp.run(f"arbor-build-catalogue local {cat}", shell=True, check=True)
    return A.load_catalogue(fn)


class recipe(A.recipe):
    def __init__(self):
        A.recipe.__init__(self)
        self.props = A.neuron_cable_properties()
        cat = compile(here / "local-catalogue.so", here / "cat")
        self.props.catalogue.extend(cat, "local_")

    def num_cells(self):
        return 1

    def cell_kind(self, _):
        return A.cell_kind.cable

    def cell_description(self, gid):
        return cell

    def probes(self, _):
        return [A.cable_probe_membrane_voltage("(location 0 0.5)", tag="Um")]

    def global_properties(self, _):
        return self.props

    def event_generators(self, _):
        if correlate:
            return [
                A.event_generator(
                    f"syn-{end}", weight, A.poisson_schedule(0.01 * U.kHz, seed=42)
                )
                for end in ends
            ]
        else:
            return [
                A.event_generator(
                    f"syn-{end}", weight, A.poisson_schedule(0.01 * U.kHz, seed=end)
                )
                for end in ends
            ]

ctx = A.context()
if A.config()["profiling"]:
    A.profiler_initialize(ctx)

mdl = recipe()
ddc = A.partition_load_balance(mdl, ctx)
sim = A.simulation(mdl, ctx, ddc)
sim.record(A.spike_recording.all)
scd = A.regular_schedule(0.1 * U.ms)
hdl = sim.sample((0, "Um"), scd)

print("Running simulation for 1s...")
t0 = pc()
sim.run(1 * U.s, 0.0025 * U.ms)
t1 = pc()
print(f"Simulation done, took: {t1 - t0:.3f}s")
if A.config()["profiling"]:
    print(A.profiler_summary(1.0))

res = {}
for d, m in sim.samples(hdl):
    res["time"] = d[:, 0]
    res["Um"] = d[:, 1]
pl.DataFrame(res).write_parquet(dn / "um.parquet")

spks = sim.spikes()
pl.DataFrame({"time": spks["time"], "gid": spks["source"]["gid"]}).write_parquet(
    dn / "spikes.parquet"
)

inp = []
if correlate:
    inp += A.poisson_schedule(0.01 * U.kHz, seed=42).events(0*U.ms, 1*U.s)
else:
    for end in ends:
        inp += A.poisson_schedule(0.01 * U.kHz, seed=42).events(0*U.ms, 1*U.s)

pl.DataFrame({"time": np.array(inp)}).write_parquet(
    dn / "input.parquet"
)
