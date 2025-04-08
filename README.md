# raman_vib_analysis

Raman_vib_analysis is a script used for analysing the Raman vibrational modes of a given calculation and estimation of contribution of each atom or a group to the specific vibrational mode. The usage of the program is as follows:
```
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
```

Wheb analysing the Turbomole Raman calculation, aoforce file should be specified as the positional argument "filename" and the control file has to be given with the use of "-c" optional argument.

## Examples

### Water

If we use a the script on the water.out file (orca Raman calculation given in the test subdirectory) with the following command:
```
python3 vib_analysis.py water.out -int -p -r 1000 4000
```
we get the following output:
```
	*** GROUP CONTRIBUTIONS ***

     idx       v     Int

                         

       7  1558.9  1691.7

       8  3755.7 10000.0

       9  3878.6 3553.59
```
The script is not able to find any functional groups because it only looks for C=O, NH2 and NH groups.

We can use the "-x" argument to specify the exact indices we are interested in in the following manner:
```
python3 vib_analysis.py water.out -int -p -r 1000 4000 -x 1 2
```
Output:
```
	*** GROUP CONTRIBUTIONS ***

     idx       v     Int        OH

                               1,2

       7  1558.9  1691.7      53.8

       8  3755.7 10000.0      51.3

       9  3878.6 3553.59      54.5
```


### Phenol

For analysis of more complicated structures we can use a list file to specify the exact group indices we are interested in. An example of a list file is as follows:
```
12 13 #this is an OH group
3 12 #this is the C-O bond
1 2 3 4 5 6 "ring" #these are the carbon atoms of the benzene ring
```

To use the list in the analysis we just have to specify its name with "-xf" argument:
```
python3 vib_analysis.py orca.out -r 1200 4000 -s 0.95 -p -int -xf list
```
and we get the following output:
```
	*** GROUP CONTRIBUTIONS ***

     idx       v     Int          OH          CO        ring

                                12,13        3,12     1,...,6

      27  1246.0 2577.75        18.1        63.7        63.2

      28  1287.8   38.17         8.1        10.3        80.5

      29  1323.3   44.74        14.1         8.1        21.3

      30  1459.8   53.43         5.2        13.3        61.7

      31  1488.1  108.91         2.8        17.1        63.5

      32  1617.1  1724.8         1.2        20.2        91.3

      33  1626.9 1828.62         1.9        19.2        90.8

      34  3078.4 1808.47         0.0         0.0         8.3

      35  3086.2 3180.19         0.0         0.0         8.2

      36  3094.1 2816.26         0.0         0.0         8.6

      37  3101.4  869.41         0.0         0.0         8.9

      38  3111.2 10000.0         0.0         0.0         9.2

      39  3712.6 1739.92       100.0         5.9         0.0
```

We can also sort the vibrational modes by a contribution of a certain group. For example with we would want to sort the frequencies by the contribution of the C-O bond, we have to add "-sr 2" (because the C-O bond is the second printed group). The output changes to:
```
	*** GROUP CONTRIBUTIONS ***

     idx       v     Int          OH          CO        ring

                                12,13        3,12     1,...,6

      27  1246.0 2577.75        18.1        63.7        63.2

      32  1617.1  1724.8         1.2        20.2        91.3

      33  1626.9 1828.62         1.9        19.2        90.8

      31  1488.1  108.91         2.8        17.1        63.5

      30  1459.8   53.43         5.2        13.3        61.7

      28  1287.8   38.17         8.1        10.3        80.5

      29  1323.3   44.74        14.1         8.1        21.3

      39  3712.6 1739.92       100.0         5.9         0.0

      37  3101.4  869.41         0.0         0.0         8.9

      34  3078.4 1808.47         0.0         0.0         8.3

      36  3094.1 2816.26         0.0         0.0         8.6

      38  3111.2 10000.0         0.0         0.0         9.2

      35  3086.2 3180.19         0.0         0.0         8.2
```
