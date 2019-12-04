## Release following the email by Kurt Hornik
This correct failures on Debian by transforming tests like 
`class(.) == *` to `is(., *)`.

## Test environments
* ubuntu 14.04 (travis-ci), release/devel
* Windows Server 2012 R2 x64 (appveyor), R 3.5.3

## R CMD check results
There were no ERROR, WARNING or NOTE.

## Reverse dependencies
There are no reverse dependencies.
