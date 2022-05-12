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
        <IDPConfGen Build> \\
"""
import argparse, re

from copy import deepcopy

from idpconfgen.libs.libio import make_folder_or_cwd
from idpconfgen.cli_build import ap as build_ap
from idpconfgen.libs import libcli
from idpconfgen import log
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
    help=('Directory to save job files. Optional.'
          ),
    type=str,
    default=None,
)

ap.add_argument(
    '--account',
    help=('Account associated with cluster. Required.'
        ),
    required=True,
    type=str,
    )

ap.add_argument(
    '--job-name',
    help=('Prefix for the job name(s). Required.'
          ),
    required=True,
    type=str,
    )

ap.add_argument(
    '--nodes',
    help=('Number of nodes/systems required. '
          'Default is 1.'
          ),
    default=1,
    type=int,
    )

ap.add_argument(
    '--ntasks-per-node',
    help=(
        'Number of physical CPUs per node required. '
        'Default is 32 cores.'
        ),
    default=32,
    type=int,
    )

ap.add_argument(
    '--mem',
    help=(
        'Number of physical CPUs/node required. '
        'The default is 32 GB.'
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
    help=(
        'Optional E-mail address for job notifications.'
        ),
    type=str,
    )

# The following grouping of arguments is inspired from @joaomcteixeira's suggestion in PR#202
build_ap_group = ap.add_argument_group("IDPConfGen Build", "Commands for `idpconfgen build`")
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
    init_files(log, LOGFILESNAME)
    destination = make_folder_or_cwd(destination)
    log.info(T('Writing #SBATCH headers'))
    
    if not re.search("\d*[-]\d*[:]\d*[:]\d*", time_per_node): 
        log.info(S(f'WARNING: the default format for `--time={time_per_node}` is d-hh:mm:ss.'))
    if not mem[-1].isalpha(): 
        log.info(S(f'WARNING: the default for `--mem={mem}` is MB, not GB.'))
    
    _header = ("#!/bin/bash\n"
               f"#SBATCH --account={account}\n"
               f"#SBATCH --job-name={job_name}\n"
               "#SBATCH --nodes=1\n"
               f"#SBATCH --ntasks-per-node={ntasks_per_node}\n"
               f"#SBATCH --mem={mem}\n"
               f"#SBATCH --time={time_per_node}\n"
               )
    if mail_user:
        _header += (f"#SBATCH --mail-user={mail_user}\n"
                    "#SBATCH --mail-type=ALL\n"
                    )
    log.info(S('done'))
    
    log.info(T('Writing job file contents'))
    seeds = [kwargs['random_seed']]
    if nodes > 1:
        seeds = [kwargs['random_seed']+i for i in range(nodes)]
        log.info(S(f'Producing {nodes} job files with random seeds {seeds}...'))
        
    # assumes you've followed the Graham installaion steps for idpconfgen
    _header += ("module load scipy-stack dssp boost\ncd\n"
                "source idpconfgen/bin/activate\n"
                "idpconfgen build \ \n\t"                
                )
    
    # trimming build arguments
    kwargs.pop('func', None)
    kwargs['ncores'] = ntasks_per_node
    filtered = {k: v for k, v in kwargs.items() if v is not None}
    kwargs.clear()
    kwargs.update(filtered)
    
    output=[]
    job_names=[]
    of = kwargs['output_folder']
    for s in seeds:
        kwargs['random_seed'] = s
        kwargs['output_folder'] += f'_rs{s}'
        job_names.append(f'{job_name}_rs{s}.sh')
        _output = _header
        
        for arg in kwargs:
            if isinstance(kwargs[arg], bool):
                if kwargs[arg] == True: _output += f"--{arg} \ \n\t"
            else:
                _output += f"--{arg} {kwargs[arg]} \ \n\t"

        output.append(_output)
        kwargs['output_folder'] = of
    log.info(S('done'))
    
    log.info(T('Writing job files to working directory'))
    i=0
    for job in output:
        with open(f"{destination}/{job_names[i]}", mode="w") as wfile:
            wfile.write(job)
        i+=1
    
    if nodes > 1:
        with open(f'{destination}/run_{job_name}_all.sh', mode="w") as wfile:
            wfile.write("#!/bin/bash\n")
            for names in job_names:
                wfile.write(f"sbatch {names}\n")
        
        with open(f'{destination}/cancel_{job_name}_all.sh', mode="w") as wfile:
            wfile.write("#!/bin/bash\n")
            for names in job_names:
                wfile.write(f"scancel -n {names}\n")
        log.info(S('done'))
        log.info(S(f'Tip: to queue up {nodes} jobs, run: bash /{destination}/run_{job_name}_all.sh'))
    else:
        log.info(S('done'))
        log.info(S(f'Tip: to queue up the job, run: sbatch /{destination}/{job_name}.sh'))

    return
