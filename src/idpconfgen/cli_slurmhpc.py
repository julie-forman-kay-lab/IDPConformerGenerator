"""
Generates shell scripts for `idpconfgen build` job submission on SLURM 
managed HPCs. Optimized for the ComputeCanada Graham cluster.

Any requests in the future, this subclient will be renamed to "sethpc"
to handle generation of job scripts/input for other management software other than SLURM 

Generates as many job-scripts as required by the user. Especially if
multiple nodes are required for scalable parallelization.

Generates "all*" and "cancel*" shell scripts to submit or cancel all jobs.

USAGE:
    $ idpconfgen slurmhpc \
        --account <ACCOUNT> \
        --job-name <NAME> \
        --nodes <#> \
        --ntasks-per-node <#> \
        --mem <#> \
        --time-per-node <d-hh:mm:ss> \
        --mail-user <@> \
        -db <DATABASE> \
        -seq <SEQUENCE> \
        -subs <SUBSTITUTIONS> \
        -dr <SEC_STR_FILTERS> \
        -nc <NUM_CONFORMERS> \
        -rs <BASE_RANDOM_SEED> \
        -xp <XMER_PROBS> \
        -lj <BACKBONE/SIDECHAIN/METHODS/LOG> \
        -scms <SIDECHAIN_OPTIONS> \
        -misc <MISCELLANEOUS> \
        -of <OUTPUT_FOLDER> \
"""
import argparse

from idpconfgen import log
from idpconfgen.libs import libcli
from idpconfgen.core import help_docs
from idpconfgen.core.definitions import seq_shortcuts, edssmat50_dict
from idpconfgen.logger import S, T, init_files

LOGFILESNAME = 'idpconfgen_slurmhpc'
_name = 'slurmhpc'
_help = 'Set-up job scripts for SLURM managed HPC.'

_prog, _des, _usage = libcli.parse_doc_params(__doc__)


ap = libcli.CustomParser(
    prog=_prog,
    description=libcli.detailed.format(_des),
    usage=_usage,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    '--account',
    help=('Account associated with cluster. Required'
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

ap.add_argument(
    '-db',
    '--database',
    help='Path to the IDPConfGen database.',
    required=True,
    )

ap.add_argument(
    '-seq',
    '--input-seq',
    help='The conformer residue sequence. String or FASTA file.',
    required=True,
    nargs='?',
    )

ap.add_argument(
    '-subs',
    '--residue-substitutions',
    help=help_docs.residue_substitutions_cli_help,
    default=None,
    )

ap.add_argument(
    '-dr',
    '--dfilters',
    help=('Secondary structure filters. '
          'E.g. -dr dloop-off csss=/path/ '
          'Defaults to loops only.'
          ),
    nargs='*',
    default=None,
    )

ap.add_argument(
    '-nc',
    '--nconfs',
    help='Total number of conformers to build.',
    default=1,
    type=int,
    )

libcli.add_argument_random_seed(ap)

ap.add_argument(
    '-xp',
    '--xmer-probs',
    help='Path to the xmer-probs file.',
    )

ap.add_argument(
    '-lj',
    '--energy-params',
    help=('Parameters for energy considerations. '
          "E.g. -lj '{'etbb':100, 'etss':1000, 'et':'pairs'}'"
          ),
    default=None,
    action=libcli.ReadDictionary,
    )

ap.add_argument(
    '-scms',
    '--sidechain-methods',
    help=('List of method(s) for side-chain addition. '
          'E.g. -scm dsd, -scm faspr, -scm mcsce mcsce-n_trials=128 mcsce-mode=simple '
          'Defaults to loops only.'
          ),
    nargs='*',
    default=None,
    )

ap.add_argument(
    '-misc',
    '--miscellaneous',
    help=('Parameters for force-fields and bgeo-path. '
          "E.g. -misc '{'ff':'{Amberff14SB}', 'bgeo_path':'/path'}'"
          ),
    default=None,
    action=libcli.ReadDictionary,
    )

libcli.add_argument_header_folder(ap)

def main(
    account,
    job_name,
    time_per_node,
    mail_user,
    database,
    input_seq,
    output_folder,
    residue_substitutions,
    xmer_probs,
    dfilters,
    energy_params,
    sidechain_methods,
    miscellaneous,
    nodes=1,
    ntasks_per_node=32,
    mem='32G',
    nconfs=1,
    random_seed=0,
    **kwargs,
    ):
    init_files(log, LOGFILESNAME)
    log.info(T('Initializing SBATCH job file(s)'))
    log.info(S('writing #SBATCH headers...'))
    
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
        
    seeds = [random_seed]
    if nodes > 1:
        seeds = [random_seed+i for i in range(nodes)]
        log.info(S(f'Producing {nodes} job files with random seeds {seeds}...'))
        
    # assumes you've followed the Graham installaion steps for idpconfgen
    _header += ("module load scipy-stack dssp boost\ncd"
                "source idpconfgen/bin/activate\n"                
                )
    output=[]
    for s in seeds:
        _output = _header
        _output += ("idpconfge build \ \n\t"
                    f"-db {database} \ \n\t"
                    f"-seq {input_seq} \ \n\t"
                    f"-nc {nconfs} \ \n\t"
                    f"-n {ntasks_per_node} \ \n\t"
                    f"-rs {s} \ \n\t"
                    f"-of {output_folder} \ \n\t"
                    )
        if residue_substitutions: _output += f"-subs {residue_substitutions} \ \n\t"
        if xmer_probs: _output += f"-xp {xmer_probs} \ \n\t"
        #TODO: dfilters, energy_params, sidechain_methods, miscellaneous
        output.append(_output)

    
    
    log.info(S('done')) 
    
    log.info(S(f'Tip: to queue up {nodes} jobs, simply: `sbatch run_{job_name}_all.sh'))
    log.info(S('done'))
    return