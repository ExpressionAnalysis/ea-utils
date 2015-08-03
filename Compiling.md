# Introduction #

Compiling & installing.


# Details #

You should be able to run

```
> make
```
and
```
> make test
```

... to install.

Some caveats are:

The sparse hash library may need updating.

The new version of varcall requires the GNU scientific library to be installed in order to compile.

On UBUNTU :

apt-get install libgsl0-dev

On CENTOS/REDHAT :

rpm -i gsl-devel

On WINDOWS:

Use MinGW, and use the [Windows port of GSL](http://gnuwin32.sourceforge.net/packages/gsl.htm)