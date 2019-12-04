## Release following the email by Kurt Hornik
This release correct failures on Debian by transforming tests like 
`class(.) == *` into `is(., *)`.

## Test environments
* ubuntu 16.04 (travis-ci), oldrel/release/devel
* Windows Server 2012 R2 x64 (appveyor), R 3.6.1

## R CMD check results
There were no ERROR, WARNING or NOTE.

## Reverse dependencies
There are no reverse dependencies.
