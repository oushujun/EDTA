print ([[

Set the environment variable `GT_MEM_BOOKKEEPING=on` to enable memory
bookkeeping (e.g., like this: `env GT_MEM_BOOKKEEPING=on gt`).

Set the environment variable `GT_ENV_OPTIONS=-spacepeak` to show a spacepeak
after program run.
Set the environment variable `GT_ENV_OPTIONS=-showtime` to show processing times
for some program parts if implemented.

Set the environment variable `GT_SEED` to an integer value to supply a seed for
the random number generator. Can be overridden by the `-seed` option.

Combinations are possible. Running the `gt` binary with `GT_ENV_OPTIONS=-help`
shows all possible "environment options".]])
