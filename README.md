# raman_vib_analysis

Raman_vib_analysis is a script used for analysing the Raman vibrational modes of a given calculation and estimation of contribution of each atom or a group to the specific vibrational mode. The usage of the program is as follows:
'''
usage: vib_analysis filename [options]

Extracting vibrational data from ORCA, GAUSSIAN and TURBOMOLE output files and
analyzing it. Program already creates output files so there's no need to pipe
the output to a file.

positional arguments:
  filename              Name of the file to be analyzed. It can be either
                        ORCA, GAUSSIAN or TURBOMOLE (aoforce) output file.
                        When analyzing TURBOMOLE output file, the control file
                        must be specified with the '-c' argument.

optional arguments:
  -h, --help            show this help message and exit
  -c CONTROL, --control CONTROL
                        Name of the control file in TM vibrational analysis.
  -r RANGE RANGE, --range RANGE RANGE
                        Frequency range to be included in the output. Default
                        is 1000-2000 cm^-1
  -s SCALING, --scaling SCALING
                        Scaling factor for the frequencies. Default is 1.00
  -x [INDICES [INDICES ...]], --indices [INDICES [INDICES ...]]
                        List of atoms to be included in the functional group
                        contribution analysis. When not specified, the program
                        will try to guess the functional groups based on the
                        atom types and proximity. The supported groups are CO,
                        NH2, and NH.
  -int, --intensities   Calculates Raman intensities based on the Raman
                        activities and save them to a file.
  -xf INDICES_FILE, --indices-file INDICES_FILE
                        File containing the list of atoms to be included in
                        the functional group contribution analysis. Each group
                        should be in a separate line. The atoms in the group
                        should be separated by spaces.
  -p, --print           Prints the contributions of functional groups to the
                        console.
  -P, --print-all       Prints the contributions of the functional groups and
                        individual atoms to the console.
  -T T                  Temperature in Kelvin (default: 298.15)
  -WL WL                Excitation wavelength in cm^-1 (default: 18797 cm^-1 /
                        532 nm)
  -sr SORT, --sort SORT
                        Sorts the output frequencies by the contributions of
                        inputed group (given as a index of a group).

Enjoy the program! :)
'''

Wheb analysing the Turbomole Raman calculation, aoforce file should be specified as the positional argument "filename" and the control file has to be given with the use of "-c" optional argument.
