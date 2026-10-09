# Removed placeholder

`Dockerfile` and `run_introme.sh` are gone. Introme is a native WDL port of
[CCICB/introme](https://github.com/CCICB/introme) `nextflow_fns` at
`d6da48bc84abda66f588470f79abcf1244acfb21` (`workflow/wdl/Introme.wdl`).

The old placeholder recorded source pin `523c052b1d87d41f7556a3bcda3a04cf663b7f0d`
and did not install Introme. `introme_docker` is still accepted and ignored so
older Terra JSON parses. Provenance tool is `introme`, version `v2@d6da48b`.
