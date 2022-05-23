"""
Generates shell scripts for `idpconfgen build` job submission on SLURM
managed HPCs. Optimized for the ComputeCanada Graham cluster.

Handle generation of job scripts/input for other workload managers
other than SLURM is possible upon request.

Generate as many job-scripts as required by the user. Especially if
multiple nodes are required for scalable parallelization.

Generates "all*" and "cancel*" shell scripts to submit or cancel all jobs.

USAGE:
    $ idpconfgen sethpc \\
        --destination <SAVE_DIRECTORY> \\
        --account <ACCOUNT> \\
        --job-name <NAME> \\
        --nodes <#> \\
        --ntasks-per-node <#> \\
        --mem <#> \\
        --time-per-node <d-hh:mm:ss> \\
        --mail-user <@> \\
        <IDPConfGen Build>
"""
import argparse
import os
import re
from copy import deepcopy
from pathlib import Path

from idpconfgen import log
from idpconfgen.cli_build import ap as build_ap
from idpconfgen.libs import libcli
from idpconfgen.libs.libio import make_folder_or_cwd
from idpconfgen.logger import S, T, init_files


LOGFILESNAME = 'idpconfgen_sethpc'
_name = 'sethpc'
_help = 'Set up job scripts for SLURM managed HPC.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)


ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    '-des',
    '--destination',
    help='Directory to save job files. Optional.',
    type=str,
    default=None,
    )

ap.add_argument(
    '--account',
    help='Account associated with cluster. Required.',
    required=True,
    type=str,
    )

ap.add_argument(
    '--job-name',
    help='Prefix for the job name(s). Required.',
    required=True,
    type=str,
    )

ap.add_argument(
    '--nodes',
    help='Number of nodes/systems required. Default to 1.',
    default=1,
    type=int,
    )

ap.add_argument(
    '--ntasks-per-node',
    help=(
        'Number of physical CPUs per node required. '
        'Defaults to 32 cores.'
        ),
    default=32,
    type=int,
    )

ap.add_argument(
    '--mem',
    help=(
        'Number of physical CPUs/node required. '
        'Defaults to 32 GB.'
        ),
    default='32G',
    type=str,
    )

ap.add_argument(
    '--time-per-node',
    help=(
        'Number of time per node required '
        'in the format d-hh:mm:ss. Required.'
        ),
    required=True,
    type=str,
    )

ap.add_argument(
    '--mail-user',
    help='Optional E-mail address for job notifications.',
    type=str,
    )

# The following grouping of arguments is inspired from @joaomcteixeira's
# suggestion in PR#202
build_ap_group = ap.add_argument_group(
    "IDPConfGen Build",
    "Commands for `idpconfgen build`",
    )

for action in build_ap._actions[1:]:
    args = deepcopy(action.__dict__)
    positional= args.pop("option_strings")
    args.pop("container")

    if isinstance(action, argparse._StoreTrueAction):
        args["action"] = "store_true"
        args.pop("nargs")
        args.pop("const")
        args.pop("type")
        args.pop("choices")
        args.pop("metavar")

    if isinstance(action, argparse._StoreFalseAction):
        args["action"] = "store_false"
        args.pop("nargs")
        args.pop("const")
        args.pop("type")
        args.pop("choices")
        args.pop("metavar")

    build_ap_group.add_argument(*positional, **args)


def main(
        destination,
        account,
        job_name,
        time_per_node,
        mail_user,
        nodes=1,
        ntasks_per_node=32,
        mem='32G',
        **kwargs,
        ):
    """
    Perform main logic of the the script.

    Parameters
    ----------
    destination : string or path, optional
        Path to saving the shell scripts to run the jobs.

    account : string, required
        Account to run the jobs off of.

    job_name : string, required
        Name of the set of jobs to run.

    time_per_node : string, required
        Time for each job to run in parallel.

    mail_user : string, optional
        E-mail address of the user to recieve job notifications.

    nodes : int, optional
        Number of nodes required for job.
        Defaults to 1.

    ntasks_per_node : int, optional
        Number of cores/workers to use.
        Defaults to 32.

    mem : string, optional.
        Amount of RAM to use per node.
        Defaults to '32G'.
    """
    init_files(log, LOGFILESNAME)
    destination = make_folder_or_cwd(destination)
    log.info(T('Writing #SBATCH headers'))

    # simple checking for common mistakes
    if not re.fullmatch(r"\d*[-]\d*[:]\d*[:]\d*", time_per_node):
        emsg = (
            f'WARNING: the default format for `--time={time_per_node}` '
            'is d-hh:mm:ss.'
            )
        log.info(S(emsg))

    if not mem[-1].isalpha():
        log.info(S(f'WARNING: the default for `--mem={mem}` is MB, not GB.'))

    # writes SBATCH headers
    __header = (
        "#!/bin/bash",
        f"#SBATCH --account={account}",
        f"#SBATCH --job-name={job_name}",
        "#SBATCH --nodes=1",
        f"#SBATCH --ntasks-per-node={ntasks_per_node}",
        f"#SBATCH --mem={mem}",
        f"#SBATCH --time={time_per_node}",
        )
    _header = os.linesep.join(__header) + os.linesep


    if mail_user:
        _header += (
            f"#SBATCH --mail-user={mail_user}{os.linesep}"
            f"#SBATCH --mail-type=ALL{os.linesep}"
            )

    log.info(S('done'))

    log.info(T('Writing job file contents'))
    seeds = [kwargs['random_seed']]
    remaining = 0

    if nodes > 1:
        seeds = [kwargs['random_seed']+i for i in range(nodes)]
        confs_per_node = int(kwargs['nconfs'] / nodes)
        remaining = kwargs['nconfs'] % nodes
        log.info(S(f'Producing {nodes} job files with random seeds {seeds}...'))
        log.info(S(
            f"Splitting {kwargs['nconfs']} conformers to {confs_per_node} "
            f"confs. per job with {remaining} additional confs. "
            "for the last job."
            ))
        kwargs['nconfs'] = confs_per_node

    # assumes you've followed the Graham installaion steps for idpconfgen
    _tmp = (
        "module load scipy-stack dssp boost",
        "cd",
        "source idpconfgen/bin/activate",
        'start_time="$(date -u +%s)"',
        "idpconfgen build \ ",
        "\t",
        )

    _header += os.linesep.join(_tmp)

    # trimming build arguments
    kwargs.pop('func', None)
    kwargs['ncores'] = ntasks_per_node
    filtered = {k: v for k, v in kwargs.items() if v is not None}
    kwargs.clear()
    kwargs.update(filtered)

    # organizing output information
    output=[]
    job_names=[]
    of = kwargs['output_folder']
    for s in seeds:
        if s == seeds[-1]:
            kwargs['nconfs'] += remaining
        kwargs['random_seed'] = s
        kwargs['output_folder'] += f'_rs{s}'
        job_names.append(f'{job_name}_rs{s}.sh')
        _output = _header

        for arg in kwargs:
            if isinstance(kwargs[arg], bool):
                if kwargs[arg] == True:
                    _output += f"--{arg} \ {os.linesep}\t"
            else:
                _output += f"--{arg} {kwargs[arg]} \ {os.linesep}\t"
        _output = _output[:-4]
        _output += os.linesep + os.linesep
        _output += (
            f'end_time="$(date -u +%s)"{os.linesep}'
            f'elapsed="$(($end_time-$start_time))"{os.linesep}'
            'echo "Total of $elapsed seconds elapsed to process: '
            f'{job_name} rs {s}."{os.linesep}'
            )
        output.append(_output)
        kwargs['output_folder'] = of
    log.info(S('done'))

    # saving shell scripts
    log.info(T('Writing job files to working directory'))
    i=0
    for job in output:
        with open(Path(destination, job_names[i]), mode="w") as wfile:
            wfile.write(job)
        i+=1

    if nodes > 1:
        with open(Path(destination, f'run_{job_name}_all.sh'), mode="w") as wfile:
            wfile.write(f"#!/bin/bash{os.linesep}")
            for names in job_names:
                wfile.write(f"sbatch {names}{os.linesep}")

        with open(Path(destination, f'cancel_{job_name}_all.sh'), mode="w") as wfile:
            wfile.write(f"#!/bin/bash{os.linesep}")
            for names in job_names:
                wfile.write(f"scancel -n {names}{os.linesep}")
        log.info(S('done'))
        log.info(S(f'Tip: to queue up {nodes} jobs, run: bash /{destination}/run_{job_name}_all.sh'))
    else:
        log.info(S('done'))
        log.info(S(f'Tip: to queue up the job, run: sbatch /{destination}/{job_name}.sh'))

    return
