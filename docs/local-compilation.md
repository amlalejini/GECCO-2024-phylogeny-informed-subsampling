# Compile and run locally

You will need a C++ compiler that supports at least C++17.
We used g++13 for all local compilations.

First, clone `GECCO-2024-phylogeny-informed-subsampling` repository that contains the code needed to run our experiment software: <https://github.com/amlalejini/GECCO-2024-phylogeny-informed-subsampling>.

Once cloned, `cd` into your local `GECCO-2024-phylogeny-informed-subsampling` repository directory.
Then, initialize and update all of the git submodules:

```
git submodule update --init --recursive
```

Once the submodules are updated, you should be able to compile either the diagnostics or program synthesis experiment code.
To specify which experiment you would like to compile, adjust the `PROJECT` variable in the `Makefile`.

- `PROJECT := prog_synth` for compiling the program synthesis code
- `PROJECT := diagnostics` for compiling the selection scheme diagnostics code

To compile in debug mode, run `make debug` from repository directory, and to compile in release mode (all optimizations turned on), run `make native`.

If you get the following error:

```
third-party/Empirical/include/emp/matching/../../../third-party/robin-hood-hashing/src/include/robin_hood.h:54:14: fatal error: sys/auxv.h: No such file or directory
   54 | #    include <sys/auxv.h> // for getauxval
```

you can fix it by doing the following:

```
cd third-party/Empirical/third-party/robin-hood-hashing
git checkout master
```

Once you have an executable, you can generate a configuration file by running:

- `./prog_synth --gen prog_synth.cfg` or
- `./diagnostics --gen diagnostics.cfg`