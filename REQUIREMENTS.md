# Requirements for Conformer Generator

_under construction_

## User Interface

User interface should be that of a shell command-line interface. Graphical interfaces can be considered at a latter stage.

### Usage

The following describes examples of usage:

```
./conformer_generator \
    -n NUMBER_OF_CONFORMERS_TO_GENERATE \
    -m MINIMUM_PROTEIN_CHUNK_SIZE \
    -x MAXIMUM_ALLOWED_MISMATCH \
    -l MAXIMUM_%_OF_ANGLES_FROM_LOOP_REGIONS \ ?? still considering
    -h MAXIMUM_%_OF_ANGLES_FROM_HELICES \ ??
    -e MAXIMUM_%_OF_ANGLES_FROM_SHEETS \ ??
    -o FOLDER_WHERE_TO_DEPOSIT_THE_STRUCTURES, defaults to CWD
    

```

