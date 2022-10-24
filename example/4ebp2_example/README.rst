Building Folded Regions - Beads on a String Model (FLDR/S)
==========================================================

.. start-description

The following example will walk you through building a folded domain within the
internal coordinate space upon which IDPConformerGenerator natively operates within.

For this exercise, please enter the :code:`example/4ebp2_example` folder where you
will find the FASTA sequence: ``4ebp2.fasta``, and a tarball containing snapshots of the
folded domain extracted from PDB ID 2MX4: ``2MX4_4ebp2.tar``.

Steps from now will assume you're in the working director of ``example/4ebp2_example``
and have already compiled your preferred reusable IDPConformerGenerator database. For
instructions on the database, please visit the previous exercise "A Real Case Scenario".

Create your folded domain database using the following commands::

    idpconfgen sscalc 2MX4_4ebp2.tar -n
    idpconfgen torsions sscalc_splitted.tar -sc sscalc.json -o torsions.json -n
    idpconfgen bgeodb sscalc_splitted.tar -sc torsions.json -o 4ebp2_fld.json -n

In this case, the folded domain resides between residue positions 18-62. We will
specify this using the ``-fldr`` flag. To build 20 conformers of 4EBP2 with the folded
domain, execute the following build function::

    idpconfgen build \
        -db <PATH TO DATABASE.JSON> \
        -seq 4ebp2.fasta \
        -etbb 100 \
        -etss 250 \
        -nc 20 \
        -fldr 18-62 \
        -flds 4ebp2_fld.json \
        -of ./4ebp2_fldrs_L+_faspr \
        -n

Please note, this method currently is great for folded regions that do not
have significant contacts with other folded regions (i.e. beads of folded regions
on a string of disordered protein). It is also capable of building multi-folded region
domains. Simply have all the structures containing your folded region within the
reusable secondary database prepared in step 1 (``4ebp2_fld.json``) and delimit
the ranges by comma. For example if we have 2 folded regions between residues 1-10 and 40-60::

    idpconfgen build \
        ... \
        -fldr 1-10,40-60 \
        -flds custom_database.json \
 
