# Unus

A Perl package for phylogenomic analyses

## Requirements

### Perl modules

*   [BioPerl](http://www.bioperl.org), including (at least) bioperl-core and bioperl-run.
    Be sure you have Muscle or ClustalW (depending on your preference) and BLAST installed.
    Muscle or ClustalW (or both) are required to be in your PATH.  BLAST can be located
    somewhere else, but you should specify the `blastbins` option.

*   `Getopt::Long`

*   `File::Basename`

*   `Pod::Usage`

*   `Log::Log4perl`

*   `Switch`

*   `Parallel::ForkManager`

## Installation

After installing all the requirements, just locate the Unus libraries in a convenient place.
You can specify the location later.

## Run

### Configuration files

A real example is available at `conf/XsPaperBSR_tblastx.conf`.  Please note that many configuration files
at `conf/` are required by the example (or any practical file).  Therefore, it is usually a good idea to
link them from the working directory (wherever you are going to call unus from):

    ln -s /path/to/unus/conf/* .

Please see the example for the syntaxis.  It is basically a key-value pairs file, supporting `#` comments.
Please note that all values are assumed to be text, even booleans.  Therefore, you can not specify a negative
boolean (for example `tblastx=0` **is illegal** and can cause Unus to call tblastx).  All the booleans are
false by default, this means that you can just delete (or comment) boolean options to make them false.

### Calling Unus

#### Option 1: perl oneliner

Ideally, you should use a configuration file.  All the parameters can be passed as command-line arguments,
but this can result in endless commands.  If you already built your configuration file, you can just run:

    perl -MUnus::Unus -e 'Unus::Unus->new->run' -- -conf configuration_file.conf

Note that you can also provide a custom location for the Unus library:

    perl -I$HOME/my.libs/unus2 -MUnus::Unus -e 'Unus::Unus->new->run' -- -conf configuration_file.conf

You can save this in your ~/.bashrc file as an alias, for example adding the following line:

    alias unus2="perl -I$HOME/my.libs/unus2 -MUnus::Unus -e 'Unus::Unus->new->run' --"

And then (after re-loging or running `source ~/.bashrc`) run:

    unus2 -conf configuration_file.conf

#### Option 2: perl script

As an alternative, you can build a small Perl script like this:

```perl
    BEGIN { push @INC, '~/my.libs/unus2 }
    use Unus::Unus;
    my $unus = Unus::Unus->new();
    $unus->run;
```

And save it somewhere in your `$PATH`.  Mine is located at `/usr/local/bin/unus2`, but you can also save it, for
example, at `~/bin/unus2`.  Do not forget to give it execution permits:

    chmod +x /usr/local/bin/unus2 # Or wherever you saved it.

And then, again, just:

   unus2 -conf configuration_file.conf


